#cmake -S . -B ../bin-X/ -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_APPLICATIONS=ON -DUSE_CGAL_LOCAL=ON -DCGAL_DIR=/sfw/cgal/gcc14.3.0/lib64/cmake/CGAL -DUSE_HYPRE=ON -DUSE_PE=OFF -DENABLE_FBM_ACCELERATION=OFF -DCMAKE_INSTALL_PREFIX=../INSTALL-X

module purge 
module load cmake/3.28.3 gcc/latest-v13 openmpi/4.1.6 python/3.13.5 cgal/6.0.1 
cmake -S . -B ../bin-CGAL6 -DCMAKE_BUILD_TYPE=Release -DBUILD_APPLICATIONS=ON -DUSE_HYPRE=ON -DUSE_CGAL_LOCAL=ON -DCGAL_DIR=/sfw/cgal/gcc14.3.0/lib64/cmake/CGAL -DUSE_PE=OFF -DCMAKE_INSTALL_PREFIX=../INSTALL-CGAL6
cmake --build ../bin-CGAL6 --target q2p1_sse metis -j 8
cmake --build ../bin-CGAL6 --target install


