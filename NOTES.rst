For large problems, HYPRE must be built with --enable-bigint, which makes it use long long int for offsets. 
Otherwise, MFEM will produce the error:

Verification failed: ((*offsets[i])[0] >= 0 && (*offsets[i])[1] >= 0) is false:
 --> overflow in offsets
 ... in function: void mfem::ParMesh::GenerateOffsets(int, int *, mfem::Array<int> **) const
 ... in file: /home1/03602/rmc298/mfem/mesh/pmesh.cpp:1661

BUILDING BOOST
===== USING BOOST 1.68.0 ====

./bootstrap.sh --with-libraries=log,thread,system,filesystem,program_options,test --prefix=/Users/rmc298/lib_builds/boost_built --with-toolset=darwin

Modify darwin entry in project-config.jam to  using darwin : 10.2 : g++-10 ; (in the if condition)

Modify tools/build/src/tools/darwin.jam  to remove -fcolasce flag

./b2 target-os=darwin architecture=x86 address-model=64 binary-format=mach-o link=static,shared variant=release --without-python define=_GLIBCXX_USE_CXX11_ABI=1 cxxflags="-fPIC -std=c++17" -j2 --prefix=/Users/rmc298/lib_builds/boost_built/ install

RUN BELOW SCRIPT ON INSTALLED LIB LOCATION. Rewrites their RPATH with an absolute RPATH
#!/bin/bash

# Modify the absolute dylib paths baked into the libraries
for i in *.dylib
do
FULLPATH=`pwd`/$i
install_name_tool -id $FULLPATH $i
echo -change $i $FULLPATH
done > changes
for i in *.dylib
do
install_name_tool `cat changes` $i
done
rm changes

Building METIS
make config prefix=/Users/rmc298/lib_builds/metis_built

Building PETSC
./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90  --with-fortran-bindings=0 COPTFLAGS='-O3 -march=native -mtune=native' CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native' --with-debugging=0 --prefix=/Users/rmc298/lib_builds/petsc_built
make PETSC_DIR=/Users/rmc298/lib_builds/petsc PETSC_ARCH=arch-darwin-c-opt all
make PETSC_DIR=/Users/rmc298/lib_builds/petsc PETSC_ARCH=arch-darwin-c-opt install
make PETSC_DIR=/Users/rmc298/lib_builds/petsc_built PETSC_ARCH="" check


Compiling PRECICE
cmake -DBUILD_SHARED_LIBS=ON -D PRECICE_PythonActions=OFF -D Boost_DIR=/Users/rmc298/lib_builds/boost_built/ -D Boost_INCLUDE_DIR=/Users/rmc298/lib_builds/boost_built/include/ -D PRECICE_PETScMapping=OFF -D CMAKE_INSTALL_PREFIX=/Users/rmc298/lib_builds/precice_built/  -D CMAKE_BUILD_TYPE=Release ..
make install


HYPRE
./configure --prefix=/Users/rmc298/lib_builds/hypre

MFEM -- DEBUG
cmake -DMFEM_USE_MPI=YES -D HYPRE_DIR=/Users/rmc298/lib_builds/hypre -DMETIS_DIR=/Users/rmc298/lib_builds/metis_built -DMFEM_DEBUG=YES -DMFEM_USE_METIS_5=YES -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=/Users/rmc298/lib_builds/mfem_built ..

SPDLOG
cmake -DCMAKE_INSTALL_PREFIX=/Users/rmc298/lib_builds/spdlog_built -DCMAKE_BUILD_TYPE=Release ..


Building CHyPS
cmake -D CMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON -D Boost_DIR=/Users/rmc298/lib_builds/boost_built/ -D Boost_INCLUDE_DIR=/Users/rmc298/lib_builds/boost_built/include -D Precice_DIR=/Users/rmc298/lib_builds/precice_built/ -D Petsc_DIR=/Users/rmc298/lib_builds/petsc_built/ -D Spdlog_DIR=/Users/rmc298/lib_builds/spdlog_built/ -D Hypre_DIR=/Users/rmc298/lib_builds/hypre -D Metis_DIR=/Users/rmc298/lib_builds/metis_built -D Mfem_DIR=/Users/rmc298/lib_builds/mfem_built -DCMAKE_CXX_FLAGS="-Wall -Werror" -D GCOV_PATH=/usr/local/bin/gcov-10 ..
