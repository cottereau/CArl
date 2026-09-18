FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

# Base tools
RUN apt update && apt install -y \
    build-essential \
    gfortran \
    cmake \
    vim \
    git \
    wget \
    python \
    python-dev \
    curl \
    libopenmpi-dev \
    openmpi-bin \
    autoconf \
    automake \
    libtool \
    autoconf-archive \
    pkg-config \
    libboost-all-dev \
    libcgal-dev \
    libeigen3-dev \
    libgmp-dev \
    libmpfr-dev \
    libvtk6-dev	   
    
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Install TetGen 1.5 (compatible with libMesh 1.1.0)
RUN cd /opt && \
    wget http://wias-berlin.de/software/tetgen/1.5/src/tetgen1.5.0.tar.gz && \
    tar -xzf tetgen1.5.0.tar.gz && \
    cd tetgen1.5.0 && \
    make    

RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py \
    && python get-pip.py

# ======================
# Install PETSc 3.6.2
# ======================

WORKDIR /opt
RUN wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.6.2.tar.gz \
    && tar -xzf petsc-3.6.2.tar.gz

WORKDIR /opt/petsc-3.6.2

RUN ./configure \
    --prefix=/opt/petsc-install \
    --with-cc=mpicc \
    --with-cxx=mpicxx \
    --with-fc=mpif90 \
    --with-debugging=0 \
    --download-fblaslapack=1 \
    && make \
    && make install

ENV PETSC_DIR=/opt/petsc-install
ENV PETSC_ARCH=""

#ENV PETSC_DIR=/opt/petsc-3.6.2
#ENV PETSC_ARCH=arch-linux-c-opt

# ======================
# Install libMesh 1.1.0
# ======================

WORKDIR /opt
RUN git clone https://github.com/libMesh/libmesh.git

WORKDIR /opt/libmesh
RUN git checkout v1.1.0

RUN ./bootstrap \
    && ./configure \
        --enable-mpi \
        --with-petsc=$PETSC_DIR \
        --enable-tetgen \
        --with-tetgen=/opt/tetgen1.5.0 \
        --enable-vtk \
        --disable-strict-lgpl \
        --disable-debugging \
    && make -j4 \
    && make install

ENV LIBMESH_DIR=/usr/local
RUN echo "/usr/local/lib" > /etc/ld.so.conf.d/libmesh.conf
RUN ldconfig

WORKDIR /workspace