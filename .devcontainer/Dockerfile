FROM nvidia/cuda:12.2.2-cudnn8-devel-ubuntu22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

COPY scripts/build-openmpi.sh /tmp/build-openmpi.sh
RUN ln -s /usr/local/cuda/lib64/stubs/libcuda.so \
          /usr/local/cuda/lib64/libcuda.so && \
    ln -s /usr/local/cuda/lib64/stubs/libcuda.so \
          /usr/local/cuda/lib64/libcuda.so.1
RUN bash /tmp/build-openmpi.sh

COPY scripts/build-hdf5.sh /tmp/build-hdf5.sh
RUN bash /tmp/build-hdf5.sh

RUN  rm /usr/local/cuda/lib64/libcuda.so /usr/local/cuda/lib64/libcuda.so.1

FROM nvidia/cuda:12.2.2-cudnn8-devel-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-numpy \
    libpython3.10 \
    libosmesa6 \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    ca-certificates \
    git \
    git-lfs \
    lsyncd \
    make \
    ssh \
    wget \
    vim \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/openmpi /usr/local/openmpi
COPY --from=builder /usr/local/hdf5 /usr/local/hdf5

ENV PATH="/usr/local/openmpi/bin:/usr/local/hdf5/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/openmpi/lib:/usr/local/hdf5/lib:${LD_LIBRARY_PATH:-}"

WORKDIR /workspace
