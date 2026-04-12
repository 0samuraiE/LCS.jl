#!/bin/bash
set -e

OPENMPI_VERSION=5.0.2
INSTALL_PREFIX=/usr/local/openmpi
CUDA_DIR=/usr/local/cuda

apt-get update && apt-get install -y --no-install-recommends \
  build-essential \
  python3 \
  perl \
  gfortran \
  wget \
  && rm -rf /var/lib/apt/lists/*

cd /tmp
wget https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-${OPENMPI_VERSION}.tar.gz
tar xzf openmpi-${OPENMPI_VERSION}.tar.gz
cd openmpi-${OPENMPI_VERSION}

./configure \
  --prefix=${INSTALL_PREFIX} \
  --with-cuda=${CUDA_DIR} \
  --enable-static \
  --enable-prte-prefix-by-default \
  LIBS=-ldl

make -j$(nproc)
make install

${INSTALL_PREFIX}/bin/ompi_info | head -5
${INSTALL_PREFIX}/bin/ompi_info --all | grep -E "built_with_cuda|extensions"

cd /tmp
rm -rf openmpi-${OPENMPI_VERSION} openmpi-${OPENMPI_VERSION}.tar.gz
