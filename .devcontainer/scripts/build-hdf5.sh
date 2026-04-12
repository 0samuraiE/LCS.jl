#!/bin/bash
set -e

# TSUBAME4.0 1.14.3
HDF5_VERSION=1.14.6
INSTALL_PREFIX=/usr/local/hdf5
MPI_PREFIX=/usr/local/openmpi

apt-get update && apt-get install -y --no-install-recommends \
  cmake \
  zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

cd /tmp
wget https://github.com/HDFGroup/hdf5/releases/download/hdf5_${HDF5_VERSION}/hdf5-${HDF5_VERSION}.tar.gz
tar xzf hdf5-${HDF5_VERSION}.tar.gz
mkdir -p hdf5-${HDF5_VERSION}/build
cd hdf5-${HDF5_VERSION}/build

cmake .. \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
  -DCMAKE_C_COMPILER=${MPI_PREFIX}/bin/mpicc \
  -DCMAKE_Fortran_COMPILER=${MPI_PREFIX}/bin/mpifort \
  -DMPI_C_COMPILER=${MPI_PREFIX}/bin/mpicc \
  -DMPI_Fortran_COMPILER=${MPI_PREFIX}/bin/mpifort \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON \
  -DHDF5_ENABLE_PARALLEL=ON \
  -DHDF5_ENABLE_Z_LIB_SUPPORT=ON \
  -DHDF5_ENABLE_SZIP_SUPPORT=OFF \
  -DHDF5_ENABLE_SUBFILING_VFD=OFF \
  -DHDF5_ENABLE_THREADSAFE=OFF \
  -DHDF5_BUILD_FORTRAN=ON \
  -DHDF5_BUILD_STATIC_LIBS=ON \
  -DHDF5_BUILD_HL_LIB=ON \
  -DHDF5_BUILD_TOOLS=ON \
  -DHDF5_BUILD_CPP_LIB=OFF \
  -DHDF5_BUILD_JAVA=OFF \

make -j$(nproc)
make install

${INSTALL_PREFIX}/bin/h5pcc -showconfig | grep -E "Parallel HDF5|Libraries"

cd /tmp
rm -rf hdf5-${HDF5_VERSION} hdf5-${HDF5_VERSION}.tar.gz
