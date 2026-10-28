# Building CalculiX

CalculiX can be build using either make or CMake. 

Building with make requires manual adjustment of `src/Makefile` to set the correct paths to the required dependencies.

Building with CMake is currently tested on Linux only but CMake will automatically discover dependencies if they are installed to standard paths. In this case the CMake build works out of the box without any further configuration.

## 1. Prerequisites

Building CalculiX requires the following dependencies installed on your system:

- C and Fortran compilers with OpenMP support (e.g. GCC)
- a BLAS and LAPACK implementation (e.g. OpenBLAS)
- ARPACK/ARPACK-ng

Additionally, at least one of the following sparse linear solvers must be provided:

- SPOOLES (enabled as default)
- PARDISO (provided as part of Intel oneAPI)
- PaStiX (either upstream v6.4.0 or PaStiX4CalculiX, see below for more information)

If you build on Linux, you might be able to install the dependencies from your distributions package manager. For example, on Ubuntu 24.04 you can install all required dependencies for a build with SPOOLES via
```
sudo apt install gfortran libopenblas-dev libarpack2-dev libspooles-dev 
```
and for building with CMake you also need to install CMake via
```
sudo apt install cmake
```

## 2. Building with CMake

CalculiX ships with CMake presets to provide you with several preconfigured build specifications for GCC or Intel (LLVM) compilers and 32- or 64-bit integers.

__Note:__ 64-bit integer support is only needed for *very* large models where e.g. node indices exceed the limits of 32-bit integers. If you enable 64-bit support you need to also provide 64-bit enabled versions of the dependencies (BLAS, ARPACK, and the solvers). SPOOLES does not support 64-bit integers.

Using the `gnu-lp64` preset (GCC, 32-bit integers and SPOOLES solver) distributed with CalculiX, you may configure and then build by running the following commands at the project root:
```
cmake --preset gnu-lp64
cmake --build --preset gnu-lp64 -j
```
Note that the build directory is defined in the preset, see `CMakePresets.json` for more details. Therefore, the executable is found in this case at `build/gnu-lp64/CalculiX`.

Standard CMake workflows apply, for example a debug build next to the release build may be configured and built by
```
cmake -B build/gnu-lp64-debug --preset gnu-lp64 -DCMAKE_BUILD_TYPE=Debug
cmake --build build/gnu-lp64-debug -j
```

**Note:** If the dependencies are not installed to standard paths, CMake might fail discovering them automatically. See section 2.3 below on how to tell CMake where to look in such cases. 

### 2.1 Building with Intel oneAPI compilers

**Note:** If you build using CMake, you need to use Intel oneMKL version 2021.3 or later, since Intel shipped the CMake config module MKLConfig.cmake which we rely on starting from this oneAPI version.

Source your Intel oneAPI installation's `setvars.sh` script and pick one of the two Intel presets from `CMakePresets.json`, i.e. to build CalculiX with the MKL PARDISO solver:
```
source </path/to/intel/compiler/root>/setvars.sh
cmake --preset intel-llvm-lp64 -DPARDISO=ON
cmake --build --preset intel-llvm-lp64 -j
```
You need to substitute the correct path to `setvars.sh` for the Intel oneAPI installation you wish to use. Furthermore, note that the build path is `build/intel-llvm-lp64` as defined in the `intel-llvm-lp64` preset, see `CMakePresets.json` for more details.

### 2.2 Solvers

The CMake build supports the SPOOLES, PARDISO and PaStiX solvers by enabling their corresponding CMake options at CMake configure time. SPOOLES is set as default but you may build CalculiX with any combination of the supported solvers.

#### 2.2.1 SPOOLES

Configure the CalculiX build with SPOOLES with `cmake --preset gnu-lp64 -DSPOOLES=ON`.

#### 2.2.2 PARDISO

Building CalculiX with PARDISO requires the Intel oneAPI MKL, see above instructions for building with Intel oneAPI compilers.

#### 2.2.3 PaStiX

The CMake build supports building CalculiX with either one of currently two variants of PaStiX.

- PaStiX v6.4.0 as released by [Inria](https://gitlab.inria.fr/solverstack/pastix/-/releases/v6.4.0)
- PaStiX4CalculiX ([commit d6b493b](https://github.com/Dhondtguido/PaStiX4CalculiX/tree/d6b493bc431ecf860b3578a2ee6a25f392eb44be)) which has been forked from PaStiX v6.0.1

__Note:__ PaStiX4CalculiX contains some optimizations for reusing reorderings of the sparse matrix structure in case of repeated solves of an equation system where only the values but not the sparsity patterns change. While these optimizations have been [shown](https://mediatum.ub.tum.de/1575696) in the past to improve performance e.g. in some contact problems, they have not been merged in upstream PaStiX.

__Note:__ Inria's official PaStiX v6.4.0 release is supported to allow easier building of CalculiX by

- using the CMake config package provided by upstream PaStiX for automatic discovery 
- supporting PaStiX builds without CUDA
- supporting PaStiX builds with CUDA on current CUDA versions
- supporting PaStiX binary installations installed from package managers (e.g. `pacman` or `Spack`).

To build CalculiX with PaStiX4CalculiX, configure the build with `cmake --preset gnu-ilp64 -DPASTIX4CALCULIX=ON` before building.

To build CalculiX with Inrias official PaStiX v6.4.0, configure the build with `cmake --preset gnu-ilp64 -DPASTIX=ON` before building. CMake will automatically swap in the correct binding for this PaStiX release for the build from the `cmake/compat/pastix/v6.4.0` folder without touching `src/`.

__Note:__ PaStiX is only tested with 64-bit integer support, so use the ILP64-flavored presets `gnu-ilp64` or `intel-llvm-ilp64` and make sure the dependencies you provide to the build also support 64-bit integers.

### 2.3 Dependencies in non-standard paths

If you have some of the dependencies installed manually into non-standard paths, CMake might not be able to find them on its own. In such cases, you can specify the path(s) to their installation(s) directly in the call to CMake, for example:

- __OpenBLAS:__ Use `CMAKE_PREFIX_PATH`, for example with OpenBLAS installed to `/opt/openblas/` invoke CMake with `cmake --preset gnu-lp64 -DCMAKE_PREFIX_PATH=/opt/openblas/`.
- __ARPACK-ng:__ Use `-Darpackng_DIR`, for example with ARPACK-ng installed to `/opt/arpackng`, invoke CMake with `cmake --preset gnu-lp64 -Darpackng_DIR=/opt/arpackng/lib64/cmake/arpackng`.
- (Legacy) __ARPACK:__ Use `-DARPACK_ROOT` or `-DARPACK_LIBRARY_DIR` to point at your ARPACK installation's root or library directory.
- __SPOOLES:__ Use either `-DSPOOLES_ROOT` or `-DSPOOLES_LIBRARY_DIR` and `-DSPOOLES_INCLUDE_DIR` for your SPOOLES installation's root or library/include directories, respectively.
- __PARDISO__/__MKL__/__Intel oneAPI:__ No intervention needed - handled by Intel oneAPI's `setvars.sh` script, see section 2.1 above. 
- __PaStiX:__ Use `-DPASTIX_DIR` and `-DSCOTCH_DIR`, for example `cmake --preset gnu-ilp64 -DPASTIX_DIR=/opt/pastix/lib/cmake/pastix -DSCOTCH_DIR=/opt/scotch/lib64/cmake/scotch`.
- __PaStiX4CalculiX:__ Append the `pkgconfig` path of your PaStiX4CalculiX installation to `PKG_CONFIG_PATH` and use `-D<libname>_DIR` for the transitive dependencies, for example `PKG_CONFIG_PATH="/opt/pastix4calculix/lib/pkgconfig:${PKG_CONFIG_PATH}" cmake --preset gnu-ilp64 -DSCOTCH_DIR=/opt/scotch/lib64/cmake/scotch`.

### 2.4 Patching hooks

The CMake build exposes a set of patching hooks that allow you to modify or expand on the CalculiX sources without directly overwriting the files in `src/` or patching `CMakeLists.txt`.

They are plain CMake variables, so they can be set via `-D` at configure time (e.g. in a toolchain or custom preset).

| Variable                    | Description                                                                |
|-----------------------------|----------------------------------------------------------------------------|
| `CCX_SOURCE_EXCLUDES`       | List of `src/<basename>` entries to remove from the build (e.g. `utemp.f`) |
| `CCX_SOURCE_ADDITIONS`      | List of absolute paths of extra source files to compile                    |
| `CCX_EXTRA_COMPILE_DEFS`    | Additional compile definitions (e.g. for preprocessor macros )             |
| `CCX_EXTRA_COMPILE_OPTIONS` | Additional compile options (e.g. `-march=native`)                          |
| `CCX_EXTRA_INCLUDE`         | Additional include directories                                             |
| `CCX_EXTRA_LIBS`            | Additional link libraries                                                  |

For example, overriding an existing source file can be done using CMake:
```
cmake --preset gnu-lp64 \
  -DCCX_SOURCE_EXCLUDES=myfile.f \
  -DCCX_SOURCE_ADDITIONS=/path/to/myfile_alt.f
```
