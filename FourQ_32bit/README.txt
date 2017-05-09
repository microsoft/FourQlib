
# FourQlib v3.0 (C Edition)
# 32-bit implementation
 
## Contents

The `FourQ_32bit` folder contains:

* Visual Studio/FourQ/: folder with Visual Studio 2015 solution and project files for compilation in Windows.
* Visual Studio/fp_tests/: folder with Visual Studio project files for testing field arithmetic functions in Windows.
* Visual Studio/ecc_tests/: folder with Visual Studio project files for testing ECC functions in Windows.
* Visual Studio/crypto_tests/: folder with Visual Studio project files for testing cryptographic functions, specifically key exchange and signatures, in Windows.
* makefile: Makefile for compilation using GNU GCC or clang compilers on Linux. 
* Main .c and .h files: library and header files. Public API for ECC scalar multiplication, key exchange and signatures is located in FourQ_api.h        
* generic/: folder with library files for 32-bit implementation.
* tests/: test files.
* README.md: this readme file.


## Supported platforms

This implementation is supported on 32-bit platforms such as x86 and ARM-based processors running Windows or Linux. We have tested the library with Microsoft Visual Studio 2015, GNU GCC v4.9 and clang v3.8. 

See instructions below to choose an implementation option and compile on one of the supported platforms.


## Complementary crypto functions

Random values are generated with `/dev/urandom` in the case of Linux, and with the function `BCryptGenRandom()` in the case of Windows.
  
The library includes an implementation of SHA-512 which is used by default by SchnorrQ signatures.

Users can experiment with different options by replacing functions in the [`random`](random/) and [`sha512`](sha512/) folders and 
applying the corresponding changes to the settings in [`FourQ_32bit\FourQ.h`](FourQ_32bit\FourQ.h). 


## Instructions for Windows

### Building the library with Visual Studio

Open the solution file ([`FourQ.sln`](FourQ_32bit\Visual Studio\FourQ\FourQ.sln)) in Visual Studio 2015, select the "Generic" configurations from the Solution Configurations menu (Win32 should appear as Solution Platform). 

By default, `USE_ENDO=true` is defined. To modify this configuration, go to the property window of the FourQ project, go to `Configuration Properties > C/C++ > Preprocessor`. Make any suitable changes, e.g., `USE_ENDO=true` or `false`. Repeat these steps for the `fp_tests`, `ecc_tests` and `crypto_tests` projects.

Finally, select "Build Solution" from the "Build" menu. 

### Running the tests

After building the solution, run `fp_tests.exe`, `ecc_tests.exe` and `crypto_tests.exe`.

### Using the library

After building the solution, add the `FourQ.lib` file to the set of References for a project, and add [`FourQ.h`](FourQ_32bit\FourQ.h) and [`FourQ_api.h`](FourQ_32bit\FourQ_api.h) to the list of header files of a project.

## Instructions for Linux

### Building the library and executing the tests with GNU GCC or clang

To compile on Linux using the GNU GCC compiler or the clang compiler, execute the following command from the command prompt:
 
```sh 
$ make ARCH=[x86/ARM] CC=[gcc/clang] USE_ENDO=[TRUE/FALSE] EXTENDED_SET=[TRUE/FALSE] CACHE_MEM=[TRUE/FALSE]
```

After compilation, run `fp_tests`, `ecc_tests` or `crypto_tests`.

By default GNU GCC is used, as well as endomorphisms and extended settings. Similarly, `CACHE_MEM=TRUE` is set by default indicating that the targeted platform contains a cache memory.

For example, to compile using clang with the efficient endomorphisms on an x86 machine, execute:

```sh
$ make ARCH=x86 CC=clang
```

As another example, to compile using GNU GCC with the efficient endomorphisms on an ARM machine, execute:

```sh
$ make ARCH=ARM
```

By default `EXTENDED_SET` is enabled, which sets the following compilation flags: `-fwrapv -fomit-frame-pointer -march=native`. To disable this, use `EXTENDED_SET=FALSE`. Users are encouraged to experiment with the different flag options.
