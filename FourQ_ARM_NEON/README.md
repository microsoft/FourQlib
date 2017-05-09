# FourQlib v3.0 (C Edition): 
# Optimized implementation for 32-bit ARM using NEON

## Contents

The `FourQ_ARM_NEON` folder contains:

* [`FourQ_ARM_NEON/makefile`](makefile): Makefile for compilation on ARM processors with NEON support using GNU GCC on Linux.
* Main .c and .h files: library and header files. Public API for ECC scalar multiplication, key exchange and signatures is in 
[`FourQ_ARM_NEON/FourQ_api.h`](FourQ_api.h).        
* [`FourQ_ARM_NEON/ARM/`](ARM/): folder with library files implementing low-level arithmetic for ARM.
* [`FourQ_ARM_NEON/tests/`](tests/): test files.
* [`FourQ_ARM_NEON/README.md`](README.md): this readme file.

## Supported platforms

This implementation is supported on 32-bit ARM platforms that contain the SIMD engine called NEON and run Linux. 
For example, platforms with a NEON engine include many cores with the ARMv7 architecture. The implementation 
has been optimized for ARM Cortex-A7, Cortex-A8, Cortex-A9 and Cortex-A15 based processors.

See instructions below to choose an implementation option and compile on one of the supported platforms.

## Complementary crypto functions

Random values are generated with `/dev/urandom`.
  
The library includes an implementation of SHA-512 which is used by default by SchnorrQ signatures.

Users can experiment with different options by replacing functions in the `random` and `sha512` folders and applying the 
corresponding changes to the settings in [`FourQ.h`](FourQ.h). 

## Instructions to build the library and execute the tests with GNU GCC or clang

To compile on Linux using the GNU GCC compiler or the clang compiler, execute the following command from the 
command prompt:

```sh 
$ make CC=[gcc/clang] USE_ENDO=[TRUE/FALSE] EXTENDED_SET=[TRUE/FALSE] INTERLEAVE=[TRUE/FALSE] MIX_ARM_NEON=[TRUE/FALSE]
```

After compilation, run `fp_tests`, `ecc_tests` or `crypto_tests`.

By default GNU GCC is used, as well as endomorphisms and extended settings. 

There are two special optimizations that can be exploited. `INTERLEAVE` improves performance on some platforms
by interleaving load/store instructions with other non-memory instructions. This optimization is recommended
for Cortex-A7, Cortex-A8 and Cortex-A9. `MIX_ARM_NEON` improves performance on some platforms by mixing ARM and
NEON instructions. This optimization is recommended for Cortex-A7, Cortex-A9 and Cortex-A15.
By default, `INTERLEAVE` is turned off and `MIX_ARM_NEON` is turned on.

For example, to compile using GNU GCC with the efficient endomorphisms on an ARM Cortex-A15 device, execute:

```sh 
$ make
```

As another example, to compile using clang with the efficient endomorphisms on an ARM Cortex-A8 device, execute:

```sh 
$ make CC=clang INTERLEAVE=TRUE MIX_ARM_NEON=FALSE
```

By default `EXTENDED_SET` is enabled, which sets the following compilation flags: `-fwrapv -fomit-frame-pointer 
-funroll-loops`. To disable this, use `EXTENDED_SET=FALSE`.
Users are encouraged to experiment with the different flag options.
