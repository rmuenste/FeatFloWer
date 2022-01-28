#!/usr/bin/env bash

CGALPATH=/home/user/omierka/nobackup1/FeatFlower_202008/bin-intel-202110/extern/libraries/cgal-install
INSTALLPATH=/home/user/omierka/nobackup1/FeatFlower_202008/install-intel-202201

cmake -DUSE_CGAL=True -DUSE_CGAL_LOCAL=True -DCGAL_INCLUDE_DIR=${CGALPATH}/include -DCGAL_LIBRARIES="${CGALPATH}/lib64/libCGAL_Core.a;${CGALPATH}/lib64/libCGAL.a" -DCMAKE_INSTALL_PREFIX=${INSTALLPATH} ../Feat_FloWer

exit 0