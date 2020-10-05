

#include <mfem/mfem.hpp>

#include "chyps/input_parser.hpp"
#include "chyps/precice_adapter.hpp"

namespace chyps {

int main(int argc, char** argv) {
  mfem::MPI_Session mpi(argc, argv);
  PreciceAdapter("Heat_Solver", "TestXML", mpi.WorldRank(), mpi.WorldSize());

  InputParser input;

  return 0;
}

}  // namespace chyps

int main(int argc, char** argv) { return chyps::main(argc, argv); }
