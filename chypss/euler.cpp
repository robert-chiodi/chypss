//                         MFEM Example 18 - Parallel Version
//
// Compile with: make ex18
//
// Sample runs:
//
//       mpirun -np 4 ex18p -p 1 -rs 2 -rp 1 -o 1 -s 3
//       mpirun -np 4 ex18p -p 1 -rs 1 -rp 1 -o 3 -s 4
//       mpirun -np 4 ex18p -p 1 -rs 1 -rp 1 -o 5 -s 6
//       mpirun -np 4 ex18p -p 2 -rs 1 -rp 1 -o 1 -s 3
//       mpirun -np 4 ex18p -p 2 -rs 1 -rp 1 -o 3 -s 3
//
// Description:  This example code solves the compressible Euler system of
//               equations, a model nonlinear hyperbolic PDE, with a
//               discontinuous Galerkin (DG) formulation.
//
//               Specifically, it solves for an exact solution of the equations
//               whereby a vortex is transported by a uniform flow. Since all
//               boundaries are periodic here, the method's accuracy can be
//               assessed by measuring the difference between the solution and
//               the initial condition at a later time when the vortex returns
//               to its initial location.
//
//               Note that as the order of the spatial discretization increases,
//               the timestep must become smaller. This example currently uses a
//               simple estimate derived by Cockburn and Shu for the 1D RKDG
//               method. An additional factor can be tuned by passing the --cfl
//               (or -c shorter) flag.
//
//               The example demonstrates user-defined bilinear and nonlinear
//               form integrators for systems of equations that are defined with
//               block vectors, and how these are used with an operator for
//               explicit time integrators. In this case the system also
//               involves an external approximate Riemann solver for the DG
//               interface flux. It also demonstrates how to use GLVis for
//               in-situ visualization of vector grid functions.
//
//               We recommend viewing examples 9, 14 and 17 before viewing this
//               example.

// Classes FE_Evolution, RiemannSolver, DomainIntegrator and FaceIntegrator
// shared between the serial and parallel version of the example.
#include "euler.hpp"

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "mfem.hpp"

// Choice for the problem setup. See InitialCondition in ex18.hpp.
int problem;

// Equation constant parameters.
const double specific_heat_ratio = 1.4;
const double gas_constant = 1.0;

// Maximum characteristic speed (updated by integrators)
double max_char_speed;

int main(int argc, char* argv[]) {
  // 1. Initialize MPI.
  MPI_Session mpi(argc, argv);

  // 2. Parse command-line options.
  problem = 1;
  const char* mesh_file = "../data/periodic-square.mesh";
  int ser_ref_levels = 0;
  int par_ref_levels = 1;
  int order = 3;
  int ode_solver_type = 4;
  double t_final = 2.0;
  double dt = -0.01;
  double cfl = 0.3;
  int vis_steps = 50;
  bool visit = false;

  int precision = 8;
  cout.precision(precision);

  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&problem, "-p", "--problem",
                 "Problem setup to use. See options in velocity_function().");
  args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                 "Number of times to refine the mesh uniformly before parallel"
                 " partitioning, -1 for auto.");
  args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                 "Number of times to refine the mesh uniformly after parallel"
                 " partitioning.");
  args.AddOption(&order, "-o", "--order",
                 "Order (degree) of the finite elements.");
  args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                 "ODE solver: 1 - Forward Euler,\n\t"
                 "            2 - RK2 SSP, 3 - RK3 SSP, 4 - RK4, 6 - RK6.");
  args.AddOption(&t_final, "-tf", "--t-final", "Final time; start time is 0.");
  args.AddOption(&dt, "-dt", "--time-step",
                 "Time step. Positive number skips CFL timestep calculation.");
  args.AddOption(&cfl, "-c", "--cfl-number",
                 "CFL number for timestep calculation.");
  args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                 "Visualize every n-th timestep.");
  args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                 "--no-visit-datafiles",
                 "Save data files for VisIt (visit.llnl.gov) visualization.");

  args.Parse();
  if (!args.Good()) {
    if (mpi.Root()) {
      args.PrintUsage(cout);
    }
    return 1;
  }

  if (visit) {
    if (mpi.Root()) {
      args.PrintOptions(cout);
    }
  }

  // 3. Read the mesh from the given mesh file. This example requires a 2D
  //    periodic mesh, such as ../data/periodic-square.mesh.
  Mesh mesh(mesh_file, 1, 1);
  const int dim = mesh.Dimension();

  MFEM_ASSERT(dim == 2,
              "Need a two-dimensional mesh for the problem definition");

  // 4. Define the ODE solver used for time integration. Several explicit
  //    Runge-Kutta methods are available.
  ODESolver* ode_solver = NULL;
  switch (ode_solver_type) {
    case 1:
      ode_solver = new ForwardEulerSolver;
      break;
    case 2:
      ode_solver = new RK2Solver(1.0);
      break;
    case 3:
      ode_solver = new RK3SSPSolver;
      break;
    case 4:
      ode_solver = new RK4Solver;
      break;
    case 6:
      ode_solver = new RK6Solver;
      break;
    default:
      if (mpi.Root()) {
        cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
      }
      return 3;
  }

  // 5. Refine the mesh in serial to increase the resolution. In this example
  //    we do 'ser_ref_levels' of uniform refinement, where 'ser_ref_levels' is
  //    a command-line parameter.
  for (int lev = 0; lev < ser_ref_levels; lev++) {
    mesh.UniformRefinement();
  }

  // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
  //    this mesh further in parallel to increase the resolution. Once the
  //    parallel mesh is defined, the serial mesh can be deleted.
  ParMesh pmesh(MPI_COMM_WORLD, mesh);
  mesh.Clear();
  for (int lev = 0; lev < par_ref_levels; lev++) {
    pmesh.UniformRefinement();
  }

  // 7. Define the discontinuous DG finite element space of the given
  //    polynomial order on the refined mesh.
  DG_FECollection fec(order, dim);
  // Finite element space for a scalar (thermodynamic quantity)
  ParFiniteElementSpace fes(&pmesh, &fec);
  // Finite element space for a mesh-dim vector quantity (momentum)
  ParFiniteElementSpace dfes(&pmesh, &fec, dim, Ordering::byNODES);
  // Finite element space for all variables together (total thermodynamic state)
  ParFiniteElementSpace vfes(&pmesh, &fec, num_equation, Ordering::byNODES);

  // This example depends on this ordering of the space.
  MFEM_ASSERT(fes.GetOrdering() == Ordering::byNODES, "");

  HYPRE_Int glob_size = vfes.GlobalTrueVSize();
  if (mpi.Root()) {
    cout << "Number of unknowns: " << glob_size << endl;
  }

  // 8. Define the initial conditions, save the corresponding mesh and grid
  //    functions to a file. This can be opened with GLVis with the -gc option.

  // The solution u has components {density, x-momentum, y-momentum, energy}.
  // These are stored contiguously in the BlockVector u_block.
  Array<int> offsets(num_equation + 1);
  for (int k = 0; k <= num_equation; k++) {
    offsets[k] = k * vfes.GetNDofs();
  }
  BlockVector u_block(offsets);

  // Momentum grid function on dfes for visualization.
  ParGridFunction mom(&dfes, u_block.GetData() + offsets[1]);

  // Initialize the state.
  VectorFunctionCoefficient u0(num_equation, InitialCondition);
  ParGridFunction sol(&vfes, u_block.GetData());
  sol.ProjectCoefficient(u0);

  // Output the initial solution.
  VisItDataCollection* visit_dc = nullptr;
  std::array<ParGridFunction, num_equation> uk_num;
  if (visit) {
    visit_dc = new mfem::VisItDataCollection(fes.GetComm(), "Euler", &pmesh);
    std::array<std::string, num_equation> names{
        {"density", "x-momentum", "y_momentum", "energy"}};
    for (int k = 0; k < num_equation; ++k) {
      // Setup mesh/variable for export
      uk_num[k].SetSpace(&fes);
      // Fun fact, this copy assignment only copies the data, not the space!
      // That is why we need to set the space above
      uk_num[k].SetFromTrueDofs(u_block.GetBlock(k));
      // Visit specific database for export
      visit_dc->RegisterField(names[k], &(uk_num[k]));
    }
    visit_dc->SetCycle(0);
    visit_dc->SetTime(0.0);
    visit_dc->Save();
  }

  // 9. Set up the nonlinear form corresponding to the DG discretization of the
  //    flux divergence, and assemble the corresponding mass matrix.
  MixedBilinearForm Aflux(&dfes, &fes);
  Aflux.AddDomainIntegrator(new DomainIntegrator(dim));
  Aflux.Assemble();

  ParNonlinearForm A(&vfes);
  RiemannSolver rsolver;
  A.AddInteriorFaceIntegrator(new FaceIntegrator(rsolver, dim));

  // 10. Define the time-dependent evolution operator describing the ODE
  //     right-hand side, and perform time-integration (looping over the time
  //     iterations, ti, with a time-step dt).
  FE_Evolution euler(vfes, A, Aflux.SpMat());

  // Determine the minimum element size.
  double hmin;
  if (cfl > 0) {
    double my_hmin = pmesh.GetElementSize(0, 1);
    for (int i = 1; i < pmesh.GetNE(); i++) {
      my_hmin = min(pmesh.GetElementSize(i, 1), my_hmin);
    }
    // Reduce to find the global minimum element size
    MPI_Allreduce(&my_hmin, &hmin, 1, MPI_DOUBLE, MPI_MIN, pmesh.GetComm());
  }

  // Start the timer.
  tic_toc.Clear();
  tic_toc.Start();

  double t = 0.0;
  euler.SetTime(t);
  ode_solver->Init(euler);

  if (cfl > 0) {
    // Find a safe dt, using a temporary vector. Calling Mult() computes the
    // maximum char speed at all quadrature points on all faces.
    max_char_speed = 0.;
    Vector z(sol.Size());
    A.Mult(sol, z);
    // Reduce to find the global maximum wave speed
    {
      double all_max_char_speed;
      MPI_Allreduce(&max_char_speed, &all_max_char_speed, 1, MPI_DOUBLE,
                    MPI_MAX, pmesh.GetComm());
      max_char_speed = all_max_char_speed;
    }
    dt = cfl * hmin / max_char_speed / (2 * order + 1);
  }

  // Integrate in time.
  bool done = false;
  for (int ti = 0; !done;) {
    double dt_real = min(dt, t_final - t);

    ode_solver->Step(sol, t, dt_real);
    if (cfl > 0) {
      // Reduce to find the global maximum wave speed
      {
        double all_max_char_speed;
        MPI_Allreduce(&max_char_speed, &all_max_char_speed, 1, MPI_DOUBLE,
                      MPI_MAX, pmesh.GetComm());
        max_char_speed = all_max_char_speed;
      }
      dt = cfl * hmin / max_char_speed / (2 * order + 1);
    }
    ti++;

    done = (t >= t_final - 1e-8 * dt);

    if (visit) {
      if (done || ti % vis_steps == 0) {
        if (mpi.Root()) {
          cout << "time step: " << ti << ", time: " << t << endl;
        }
        // Update ParGridFunction from u_block data.
        for (int k = 0; k < num_equation; ++k) {
          uk_num[k].SetFromTrueDofs(u_block.GetBlock(k));
        }
        visit_dc->SetCycle(ti);
        visit_dc->SetTime(t);
        visit_dc->Save();
      }
    }
  }

  tic_toc.Stop();
  if (mpi.Root()) {
    cout << " done, " << tic_toc.RealTime() << "s." << endl;
  }

  // 12. Compute the L2 solution error summed for all components.
  if (t_final == 2.0) {
    const double error = sol.ComputeLpError(2, u0);
    if (mpi.Root()) {
      cout << "Solution error: " << error << endl;
    }
  }

  // Free the used memory.
  delete ode_solver;

  return 0;
}
