For large problems, HYPRE must be built with --enable-bigint, which makes it use long long int for offsets. 
Otherwise, MFEM will produce the error:

Verification failed: ((*offsets[i])[0] >= 0 && (*offsets[i])[1] >= 0) is false:
 --> overflow in offsets
 ... in function: void mfem::ParMesh::GenerateOffsets(int, int *, mfem::Array<int> **) const
 ... in file: /home1/03602/rmc298/mfem/mesh/pmesh.cpp:1661
