# CArl

![Status](https://img.shields.io/badge/status-under%20development-orange)
![Platform](https://img.shields.io/badge/platform-Linux-blue)
![Docker](https://img.shields.io/badge/docker-supported-2496ED?logo=docker&logoColor=white)
![License](https://img.shields.io/badge/license-open%20source-green)
![C++](https://img.shields.io/badge/language-C%2B%2B17-blue?logo=cplusplus&logoColor=white)
![MPI](https://img.shields.io/badge/parallel-MPI-blueviolet)

> ⚠️ **This branch (`dyncoup`) is under active development and has not been fully tested. Use with caution.**

## PRESENTATION

This project is focused on the development of a software based on the [Arlequin multi-model coupling method](https://www.sciencedirect.com/science/article/pii/S0045782508003630). The main interest of this software is to allow, by its specific structure, the easy interfacing of different third-party softwares (developed and maintained outside of this project), and adapted to each of the models appearing in the coupling.

This software is mainly developed at the laboratory [LMPS (Laboratoire Mécanique de Paris-Saclay)](https://lmps.ens-paris-saclay.fr/en), of which MSSMat laboratory (École Centrale Paris - CNRS), that originally developed the [CArl](https://github.com/cottereau/CArl) software is now part of.

* 📧 contact : [Regis Cottereau](mailto:cottereau@lma.cnrs-mrs.fr)
* 👥 contributors (by order of first commit): R. Cottereau, C. Zaccardi, Y. Le Guennec, D. Neron, T. M. Schlittler, F. Gatti, C. Luo, R. Ruyssen

More detail on usage and examples can be found on the [related help web page](https://cottereau.github.io/CArl/).

---

## ⚙️ C++ / MPI IMPLEMENTATION

The C++ / MPI implementation of the CArl software can be found in the directory `Cpp`. It is capable of interfacing with external solvers based on the [PETSc](http://www.mcs.anl.gov/petsc/) toolkit (including the libMesh solvers, when compiled with PETSc support). Two solvers are available:

| Solver           | Description                                 | Status                                                              |
| ---------------- | ------------------------------------------- | ------------------------------------------------------------------- |
| **CArl-Static**  | Stationary solution of a coupled system     | ![status](https://img.shields.io/badge/-not%20tested-red)           |
| **CArl-Dynamic** | Time-dependent solution of a coupled system | ![status](https://img.shields.io/badge/-under%20development-orange) |

---

## 📦 REQUIREMENTS

The following third-party libraries are required. All are installed automatically by the provided Dockerfile.

| Library                                                                   | Version              | Status                                                      |
| ------------------------------------------------------------------------- | -------------------- | ----------------------------------------------------------- |
| [CMake](https://cmake.org)                                                | system               | ![status](https://img.shields.io/badge/-tested-brightgreen) |
| [Boost](http://www.boost.org)                                             | system (≥ 1.66)      | ![status](https://img.shields.io/badge/-tested-brightgreen) |
| [CGAL](http://www.cgal.org) (Core component)                              | system (≥ 5.0)       | ![status](https://img.shields.io/badge/-tested-brightgreen) |
| [OpenMPI](https://www.open-mpi.org)                                       | system               | ![status](https://img.shields.io/badge/-tested-brightgreen) |
| [HDF5](https://www.hdfgroup.org/solutions/hdf5/) — parallel/MPI           | system               | ![status](https://img.shields.io/badge/-tested-brightgreen) |
| [PETSc](http://www.mcs.anl.gov/petsc/) — latest release, with MPI + HDF5  | compiled from source | ![status](https://img.shields.io/badge/-tested-brightgreen) |
| [libMesh](https://libmesh.github.io) — master, with PETSc + HDF5 + TetGen | compiled from source | ![status](https://img.shields.io/badge/-tested-brightgreen) |

> All libraries must be built with the **same MPI implementation**.

---

## 🐳 DOCKER-BASED INSTALLATION (recommended)

![Docker](https://img.shields.io/badge/docker-recommended-2496ED?logo=docker&logoColor=white)
![Status](https://img.shields.io/badge/status-under%20development-orange)

The recommended way to build and run CArl is via the provided Dockerfile `carl_dyncoup.dockerfile`, included in the repository. The image is based on **Ubuntu 24.04** and builds the full dependency chain automatically.

### Quick start

Clone the repository and run the install script:

```bash
git clone --branch dyncoup https://github.com/cottereau/CArl.git
cd CArl
bash install_carl.sh
```

The install script wraps the following Docker commands:

```bash
# Build the image
docker buildx build . -t carl -f carl_dyncoup.dockerfile

# Run the container
docker run --name carl -p 8880:8000 carl
```

### Connect to the running container

```bash
docker exec -it carl /bin/bash
```

CArl executables are available on the `PATH` inside the container under `/opt/CArl/build/bin/`.

### Rebuild without cache

To force a full rebuild from scratch:

```bash
docker buildx build --no-cache . -t carl -f carl_dyncoup.dockerfile
```

To rebuild only the CArl step while keeping cached PETSc and libMesh layers, increment the `CACHE_BUST` argument:

```bash
docker buildx build --build-arg CACHE_BUST=2 . -t carl -f carl_dyncoup.dockerfile
```

---

## 🔧 WHAT THE DOCKERFILE DOES

![Status](https://img.shields.io/badge/status-under%20development-orange)

The Dockerfile proceeds in four stages:

### 1. System packages (Ubuntu 24.04)

Installs: `build-essential`, `cmake`, `git`, `gfortran`, `libopenmpi-dev`, `openmpi-bin`, `libgmp-dev`, `libmpfr-dev`, `libcgal-dev`, `libboost-all-dev`, `autoconf`, `autoconf-archive`, `automake`, `libtool`, `m4`, `pkg-config`, `libhdf5-openmpi-dev`, `hdf5-tools`.

Key environment paths set at this stage:

| Variable           | Value                                    |
| ------------------ | ---------------------------------------- |
| `MPI_DIR`          | `/usr/lib/x86_64-linux-gnu/openmpi`      |
| `HDF5_DIR`         | `/usr/lib/x86_64-linux-gnu/hdf5/openmpi` |
| `HDF5_INCLUDE_DIR` | `/usr/include/hdf5/openmpi`              |
| `BOOST_ROOT`       | `/usr`                                   |

### 2. PETSc — latest release branch

Cloned from `https://gitlab.com/petsc/petsc.git` into `/opt/petsc`, configured with:
- MPI wrappers (`mpicc`, `mpicxx`, `mpif90`)
- Downloaded BLAS/LAPACK, Metis, ParMetis
- System parallel HDF5
- Optimised build (`-O3`, no debugging)

### 3. libMesh — master branch

Cloned from `https://github.com/libMesh/libmesh.git` into `/opt/libmesh-src`, installed to `/opt/libmesh`, configured with:
- PETSc support (pointing to step 2)
- System HDF5 and Boost
- TetGen enabled (via bundled submodule)
- Optimised build (`opt` method only)

> Note: source and install directories are kept separate to avoid a known `make install` conflict.

### 4. CArl — `dyncoup` branch

Cloned from `https://github.com/cottereau/CArl.git` into `/opt/CArl`, built with CMake pointing to all of the above. Runtime paths set:

| Variable          | Value                                                                                     |
| ----------------- | ----------------------------------------------------------------------------------------- |
| `PATH`            | `/opt/CArl/build/bin` prepended                                                           |
| `LD_LIBRARY_PATH` | `/opt/libmesh/lib:/opt/petsc/arch-linux-c-opt/lib:/usr/lib/x86_64-linux-gnu/hdf5/openmpi` |

---

## 🩹 SOURCE CODE PATCHES

The `dyncoup` branch requires the following patches to compile against current versions of libMesh and PETSc. These are already applied in the repository:

| File                                          | Patch                                                                                                           | Reason                                                   |
| --------------------------------------------- | --------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------- |
| `Cpp/src/include/common_functions.h`          | `reduced_system_init()` rewritten to use `System::add_variables()`                                              | `DofMap::add_variable_group()` removed in modern libMesh |
| `Cpp/CMakeLists.txt`                          | Added `link_directories` for libMesh/PETSc; explicit linking of `mesh_opt`, `timpi_opt`, `petsc` to all targets | Modern linker requires explicit DSO listing              |
| `Cpp/src/include/mesh_intersection_methods.h` | Added `#include "libmesh/tetgen_mesh_interface.h"`                                                              | TetGen interface header not transitively included        |

---

## ✅ TESTED CONFIGURATION

| Component | Version                   |
| --------- | ------------------------- |
| OS        | Ubuntu 24.04              |
| Compiler  | GCC 13 + OpenMPI          |
| PETSc     | latest release branch     |
| libMesh   | master branch             |
| CGAL      | 5.6 (system)              |
| Boost     | 1.83.0 (system)           |
| HDF5      | parallel/openmpi (system) |

---

## 📚 REFERENCES

Original references for the Arlequin method:

1. H. Ben Dhia. Multiscale mechanical problems: the Arlequin method, _Comptes Rendus de l'Academie des Sciences - Series IIB 326_ (1998), pp. 899-904.
2. H. Ben Dhia, G. Rateau. The Arlequin method as a flexible engineering design tool, _Int. J. Numer. Meths. Engr._ 62 (2005), pp. 1442-1462.

A list of scientific papers that make use of the CArl software can be found in [references.rtf](references.rtf).