# Installation

This section describes how libgran can be installed into a project.
If your project is a git repository, libgran can be added to it as a git submodule. To install libgran and
its dependency libtimestep in a directory named deps under your project root, execute the following commands:
```Bash
git submodule add https://github.com/egor-demidov/libgran deps/libgran
git submodule add https://github.com/egor-demidov/libtimestep deps/libtimestep
```
Since libgran and libtimestep are header-only libraries, they can be included in your project without linking against
an object. Simply, add the include directories of libgran and libtimestep in your **CMakeLists.txt** file:
```cmake
include_directories(deps/libgran/include)
include_directories(deps/libtimestep/include)
```
By default, libgran uses C++ 17 parallel algorithms for computation of binary interactions. In case you would like to
use OpenMP instead, `LIBGRAN_USE_OMP` compile definition needs to be added. In your **CMakeLists.txt**, add the following
line:
```cmake
add_compile_definitions(LIBGRAN_USE_OMP)
```
Note that for improved performance / parallelization capabilities, additional compiler flags might be required on your
system. Platform- and toolchain-specific instructions are provided below. Also, the optimal number of threads depends on
the size of the granular system. Using too many threads can be detrimental to performance.

### Linux & GNU C++ compiler

The following compiler flags are recommended for best performance:
```
-O3 -flto=auto -march=native
```

<tabs>
<tab title="Parallelization with C++ 17 algorithms">
If you opted for C++ 17 parallel algorithms, Intel TBB library needs to be installed and linked against.
On Ubuntu, TBB can be installed with:

```Bash
sudo apt install libtbb2-dev
```

Then, in your CMakeLists.txt file, add:

```cmake
find_package(TBB REQUIRED)

# Set up your targets...

target_link_libraries(<your target name> PRIVATE TBB::tbb)
```

</tab>
<tab title="Parallelization with OpenMP">
If you opted for OpenMP, the following flag needs to be added:

```
-fopenmp
```

Then, every time you run your simulation, the number of threads can be set prior to execution with:

```Bash
export OMP_NUM_THREADS=4
```

</tab>
</tabs>

## Windows & MSVC compiler

The following compiler flags are recommended for best performance:
```
/O2 /EHsc /GL /fp:except
```
It is important that you **build configuration is set to Release**. Otherwise, MSVC produces
very slow binaries.

<tabs>
<tab title="Parallelization with C++ 17 algorithms">
C++ 17 parallel algorithms work out of the box under MSVC without the need for additional compiler
flags or libraries.

</tab>
<tab title="Parallelization with OpenMP">
If you opted for OpenMP, the following flag needs to be added:

```
/openmp
```

Then, every time you run your simulation, the number of threads can be set prior to execution with:

```shell
set OMP_NUM_THREADS=4
```

</tab>
</tabs>

## macOS & LLVM/clang compiler

LLVM/clang compiler has to be used with OpenMP parallelization.
LLVM/clang and OpenMP need to be installed:
```shell
brew install llvm libomp
```
To use LLVM/clang instead of AppleClang with CMake, run:
```shell
export CC=/usr/local/opt/llvm/bin/clang
export CXX=/usr/local/opt/llvm/bin/clang++
cmake <path to CMakeLists.txt>
```
The following compiler flags are recommended for best performance:
```
-O3 -fopenmp -flto -march=native
```
The following linker flags are recommended for best performance:
```
-lomp -flto
```
Every time you run your simulation, the number of threads can be set prior to execution with:
```shell
export OMP_NUM_THREADS=4
```

## macOS & AppleClang compiler

The default AppleClang compiler bundled with XCode **is not supported**, because it does not support
C++ 17 parallel algorithms or OpenMP.

<seealso>
<category ref="related">
    <a href="Overview.md">Overview</a>
    <a href="Tutorials.md">Tutorials</a>
    <a href="Class-reference.md">Class reference</a>
</category>
<category ref="external">
    <a href="https://github.com/egor-demidov/libgran">libgran on GitHub</a>
    <a href="https://github.com/egor-demidov/libtimestep">libtimestep on GitHub</a>
</category>
</seealso>
