## FourQlib v3.0 (C Edition)

**FourQlib** implements essential elliptic curve and cryptographic functions based on FourQ, 
a high-security, high-performance elliptic curve that targets the 128-bit security level [1]. At the 
high level, **FourQlib** consists of a set of implementations targeting different platforms with different
levels of portability and performance. The cryptographic and elliptic curve API is common to all the 
implementations.   

The library was developed by [Microsoft Research](http://research.microsoft.com/) and is available under the MIT License. 

## Contents

Version 3.0 includes the following implementations:
 
* [`FourQ_32bit`](FourQ_32bit/): a portable implementation especially tailored for 32-bit platforms.
* [`FourQ_64bit_and_portable`](FourQ_64bit_and_portable/): a portable implementation for 32-bit and 64-bit platforms with optional
  optimizations for x64 and 64-bit ARMv8 platforms.
* [`FourQ_ARM`](FourQ_ARM/): an optimized implementation for 32-bit ARMv6, ARMv7 and ARMv7-M (Cortex-M4) platforms.
* [`FourQ_ARM_side_channel`](FourQ_ARM_side_channel/): an optimized implementation for 32-bit ARMv6, ARMv7 and ARMv7-M (Cortex-M4),
  including strong countermeasures against a wide variety of side-channel attacks.
* [`FourQ_ARM_NEON`](FourQ_ARM_NEON/): an optimized implementation for 32-bit ARM platforms with NEON support. 

The elliptic curve and crypto API can be found in `FourQ_api.h`, which is available per implementation. 

The [`FourQ-Magma`](FourQ-Magma/) folder includes easy-to-read scripts written in Magma. 

### Complementary cryptographic functions

Random values are generated with `/dev/urandom` in the case of Linux, and with the function `BCryptGenRandom()`
in the case of Windows. Check the [`random`](random/) folder for details.
  
The library includes an implementation of SHA-512 which is used by default by SchnorrQ signatures (see [`sha512`](sha512/)).

Users can provide their own PRNG and hash implementations by replacing the functions in the [`random`](random/) and [`sha512`](sha512/)folders, and applying the corresponding changes to the settings in `FourQ.h` (in a given implementation). 
Refer to [2] for the security requirements for the cryptographic hash function. 

## What's new in version 3.0

* New support for co-factor ECDH and SchnorrQ signatures. 
* New implementations for 32-bit processors, 32-bit ARMv6 and ARMv7 processors, and 32-bit ARM Cortex-M4 microcontroller. 
* New implementation for ARMv6/ARMv7/ARMv7-M with strong countermeasures against several side-channel attacks.
* New support for 64-bit ARMv8 processors.  

## Main features
   
* Support for co-factor Elliptic Curve Diffie-Hellman (ECDH) key exchange [3].
* Support for the SchnorrQ digital signature scheme [2]. 
* Support for 3 core elliptic curve operations: variable-base, fixed-base and double-scalar multiplications.
* Support for Windows using Microsoft Visual Studio and Linux using GNU GCC or clang.    
* Includes a basic implementation using portable C to enable support on a wide range of platforms including x64, x86 
  and ARM, Windows and Linux. 
* Includes optimized implementations for 64-bit ARMv8 and x64 platforms with optional, high-performance x64 assembly for Linux [1].
* Includes high-performance implementations for 32-bit ARM processors with NEON support [4], for 32-bit ARMv6 and 
  ARMv7 processors, and for 32-bit ARM Cortex-M4 microcontrollers [5].
* Includes side-channel secure implementations for 32-bit ARMv6/ARMv7 and for ARMv7-M (Cortex-M4) microcontrollers [5].
* Includes testing and benchmarking code for field arithmetic, elliptic curve and cryptographic functions. 
* All functions evaluating secret data have regular, constant-time execution, protecting against timing and cache attacks.
* Includes an option to disable the use of the fast endomorphisms.

## Quick start

### Building the library and executing the tests on Linux
    
One can quickly test a given implementation by executing from the corresponding folder and using a supported architecture:

```sh
$ make ARCH=[x64/x86/ARM/ARM64] 
```

GNU GCC is used by default. After compilation, run `fp_tests`, `ecc_tests` or `crypto_tests`. 

Below are the architectures supported by each implementation:

* [`FourQ_32bit`](FourQ_32bit/): x86 and ARM
* [`FourQ_64bit_and_portable`](FourQ_64bit_and_portable/): x64, x86, ARM and ARM64
* [`FourQ_ARM`](FourQ_ARM/), [`FourQ_ARM_side_channel`](FourQ_ARM_side_channel/) and [`FourQ_ARM_NEON`](FourQ_ARM_NEON/): ARM

For example, to compile the optimized x64 implementation using assembly with GNU GCC, using the efficient endomorphisms on a machine with AVX2 support (e.g, Intel's Haswell or Skylake), execute:

```sh
$ cd FourQ_64bit_and_portable
$ make ARCH=x64
```

Additional compilation options are available. Refer to the `README` files in a given implementation folder for complete details.

**NOTE:** the above instructions apply to all the "processor-class" implementations. For instructions on how to compile on an ARM Cortex-M (ARMv7-M) microcontroller, refer to the `README` files in [`FourQ_ARM_side_channel`](FourQ_ARM_side_channel/) or [`FourQ_ARM`](FourQ_ARM/).  
	  
### Building the library and executing the tests on Windows

`FourQ_32bit` and `FourQ_64bit_and_portable` include Visual Studio solutions for compilation on Windows. Refer to the corresponding `README` files for instructions.

## License

**FourQlib** is licensed under the MIT License; see [`License`](LICENSE) for details.

Files `stm32f4_wrapper.c` and `stm32f4_wrapper.h` in the [`FourQ_ARM`](FourQ_ARM/) and [`FourQ_ARM_side_channel`](FourQ_ARM_side_channel/) folders are by Joost Rijneveld and are released under the CC0 1.0 Universal license.

Files in the folder [`FourQ_ARM_side_channel/libopencm3`](FourQ_ARM_side_channel/libopencm3/) and [`FourQ_ARM/libopencm3`](FourQ_ARM/libopencm3/) are from the `libopencm3` project and are under the GNU LGPL v3.0 license.
 
The SHA-512 implementation is by D.J. Bernstein and is released to the public domain.

# References

[1]   Craig Costello and Patrick Longa, "FourQ: four-dimensional decompositions on a Q-curve over the Mersenne prime". Advances in Cryptology - ASIACRYPT 2015, 2015. 
The extended version is available [`here`](http://eprint.iacr.org/2015/565).

[2]   Craig Costello and Patrick Longa. "SchnorrQ: Schnorr signatures on FourQ". MSR Technical Report, 2016. 
Available [`here`](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/07/SchnorrQ.pdf).

[3]   Watson Ladd, Patrick Longa and Richard Barnes, "Curve4Q". Internet-Draft, draft-ladd-cfrg-4q-01, 2017.
Available [`here`](https://www.ietf.org/id/draft-ladd-cfrg-4q-01.txt).

[4]   Patrick Longa, "FourQNEON: faster elliptic curve scalar multiplications on ARM processors". Selected Areas in Cryptography (SAC 2016), 2016.
Preprint available [`here`](http://eprint.iacr.org/2016/645).

[5]   Zhe Liu, Patrick Longa, Geovandro Pereira, Oscar Reparaz and Hwajeong Seo, "FourQ on embedded devices with strong countermeasures against side-channel attacks".
Preprint available [`here`](http://eprint.iacr.org/2017/434).

# Contributing

This project has adopted the [Microsoft Open Source Code of Conduct](https://opensource.microsoft.com/codeofconduct/). For more information see the [Code of Conduct FAQ](https://opensource.microsoft.com/codeofconduct/faq/) or contact [opencode@microsoft.com](mailto:opencode@microsoft.com) with any additional questions or comments.
