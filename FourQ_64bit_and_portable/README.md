# FourQlib v3.0 (C Edition): portable and 64-bit optimized implementation

## Contents

The `FourQ_64bit_and_portable` folder contains:

* [`FourQ_64bit_and_portable/Visual Studio/`](Visual%20Studio/): folder with Visual Studio 2015 solution and 
project files for compilation and testing in Windows.
* [`FourQ_64bit_and_portable/makefile`](makefile): Makefile for compilation using GNU GCC or clang compilers 
on Linux. 
* Main .c and .h files: library and header files. Public API for ECC scalar multiplication, key exchange and signatures is in 
[`FourQ_64bit_and_portable/FourQ_api.h`](FourQ_api.h).        
* [`FourQ_64bit_and_portable/AMD64/`](AMD64/): folder with library files for optimized x64 implementation.
* [`FourQ_64bit_and_portable/ARM64/`](ARM64/): folder with library files for optimized 64-bit ARM 
implementation.
* [`FourQ_64bit_and_portable/generic/`](generic/): folder with library files for portable implementation.
* [`FourQ_64bit_and_portable/tests/`](tests/): test files.
* [`FourQ_64bit_and_portable/README.md`](README.md): this readme file.

## Supported platforms

This implementation is supported in a wide range of platforms including x64, x86, 32-bit ARM and 64-bit ARM,
running Windows or Linux. We have tested the library with Microsoft Visual Studio 2015, GNU GCC v4.9 and 
clang v3.8. 

See instructions below to choose an implementation option and compile on one of the supported platforms. 

## Complementary crypto functions

Random values are generated with `/dev/urandom` in the case of Linux, and with the function `BCryptGenRandom()` in the case of Windows.

The library includes an implementation of SHA-512 which is used by default by SchnorrQ signatures.

Users can experiment with different options by replacing functions in the `random` and `sha512` folders and 
applying the corresponding changes to the settings in [`FourQ.h`](FourQ.h). 

## Implementation options

The following compilation options are available for the `FourQ_64bit_and_portable` implementation:

* A portable implementation (enabled by the "GENERIC" option).
* Optimized implementations for x64 and 64-bit ARM (ARMv8). Note that the rest of platforms are only supported
  by the generic implementation. 
* Use of AVX or AVX2 instructions enabled by defining `_AVX_` or `_AVX2_` (Windows) or by the "AVX" and "AVX2" 
  options (Linux).
* Optimized x64 assembly implementations in Linux.
* Use of fast endomorphisms enabled by the "USE_ENDO" option.

Follow the instructions below to configure these different options.

## Instructions for Windows

### Building the library with Visual Studio

Open the solution file ([`FourQ.sln`](Visual%20Studio/FourQ/FourQ.sln)) in Visual Studio 2015, select 
one of the available configurations from
the Solution Configurations menu ("Release" corresponding to the high-speed x64 implementation and "Generic" 
corresponding to the portable implementation) and select one of the Solution Platforms (x64 or Win32). Note 
that Win32 is only supported with the "Generic" solution configuration.

By default, `USE_ENDO=true` and (for x64) `_AVX_` is defined. To modify this configuration, go to the property 
window of the FourQ project, go to `Configuration Properties > C/C++ > Preprocessor`. Make any suitable changes, 
e.g., delete `_AVX_` if AVX instructions are not supported, replace `_AVX_` by `_AVX2_` if AVX2 instructions
are supported, or set `USE_ENDO=true` or `false`. Repeat these steps for the `fp_tests`, `ecc_tests` and `crypto_tests` projects.

Finally, select "Build Solution" from the "Build" menu. 

### Running the tests

After building the solution, run `fp_tests.exe`, `ecc_tests.exe` and `crypto_tests.exe`.

### Using the library

After building the solution, add the `FourQ.lib` file to the set of References for a project, and add 
[`FourQ.h`](FourQ.h) and [`FourQ_api.h`](FourQ_api.h) to the list of header files of a project.

## Instructions for Linux

### Building the library and executing the tests with GNU GCC or clang

To compile on Linux using the GNU GCC compiler or the clang compiler, execute the following command from the 
command prompt:

```sh
$ make ARCH=[x64/x86/ARM/ARM64] CC=[gcc/clang] ASM=[TRUE/FALSE] AVX=[TRUE/FALSE] AVX2=[TRUE/FALSE] 
     EXTENDED_SET=[TRUE/FALSE] USE_ENDO=[TRUE/FALSE] GENERIC=[TRUE/FALSE] SERIAL_PUSH=[TRUE/FALSE] 
```

After compilation, run `fp_tests`, `ecc_tests` or `crypto_tests`.

By default GNU GCC is used, as well as the endomorphisms and the extended settings.

In the case of x64, AVX2 instructions and the high-speed assembly implementation are enabled by default.
In the case of x86 and ARM, the portable ("GENERIC") implementation is used by default.

For example, to compile the optimized x64 implementation in assembly with GNU GCC using the efficient
endomorphisms on a machine with AVX2 support (e.g, Intel's Haswell or Broadwell), execute:

```sh
$ make ARCH=x64
```

For example, to compile the optimized ARM64 implementation with GNU GCC using the efficient endomorphisms, 
execute:

```sh
$ make ARCH=ARM64
```

As another example, to compile the portable implementation with clang using the efficient endomorphisms 
on an x86 machine, execute:

```sh
$ make ARCH=x86 CC=clang
```

`SERIAL_PUSH` can be enabled in some platforms (e.g., AMD without AVX2 support) to boost performance.

By default `EXTENDED_SET` is enabled, which sets the following compilation flags: `-fwrapv -fomit-frame-pointer 
-march=native`. To disable this, use `EXTENDED_SET=FALSE`.
Users are encouraged to experiment with the different flag options.

Whenever an unsupported configuration is applied, the following message will be displayed: `#error -- "Unsupported configuration". 
For example, the use of assembly or any of the AVX options is not supported when selecting the portable implementation 
(i.e., if `GENERIC=TRUE` or if `ARCH=[x86/ARM]`). 
