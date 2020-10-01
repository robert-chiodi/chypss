

#include <mfem.hpp>

#include "chypss/precice_adapter.hpp"

int main(int argc, char** argv) {
  mfem::MPI_Session mpi(argc, argv);
  chyps::PreciceAdapter("Heat_Solver", "TestXML", mpi.WorldRank(),
                        mpi.WorldSize());
}
