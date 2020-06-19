# FourQlib v3.0 (C Edition): 
# Optimized implementation for 32-bit ARM and ARM Cortex-M4 with side-channel countermeasures

This is an **experimental** side-channel secure implementation. DO NOT USE AS IS IN PRODUCTION. See the "Security notes" below.  

This implementation includes scalar multiplication, ECDH and digital signature algorithms protected with a set of efficient
countermeasures that have been especially tailored for FourQ to minimize the risk of timing attacks, simple and differential 
side-channel analysis (SSCA/DSCA), correlation and collision attacks, including specialized attacks such the doubling attack, 
the refined power attack (RPA), zero-value point attacks (ZVP), same value attacks (SVA), exceptional procedure attacks, invalid 
point attacks, and small subgroup attacks. 

More details can be found in: 

**"FourQ on embedded devices with strong countermeasures against side-channel attacks"**, CHES 2017.                         
Zhe Liu, Patrick Longa, Geovandro Pereira, Oscar Reparaz, and Hwajeong Seo.                            
Preprint available [`here`](http://eprint.iacr.org/2017/434).

**SECURITY NOTES:** 
* No software implementation is able to guarantee 100% side-channel security. In some cases, certain powerful attacks 
such as template attacks can be carried out using a single target trace, making any randomization or masking technique useless.
Moreover, the issue gets more complicated for embedded devices that lack access to a good source of randomness. Since many SCA attacks
closely depend on the underlying hardware, it is recommended to include additional countermeasures at the software and hardware levels
depending on the targeted platform. Also, note that hardware countermeasures are usually required to properly deal with most
sophisticated invasive attacks.                                                          
* The hash function implementation in the `sha512` folder, which is used by SchnorrQ, is NOT protected against side-channel attacks such as DPA.
 
## Contents

The `FourQ_ARM_side_channel` folder contains:

* [`FourQ_ARM_side_channel/makefile`](makefile): Makefile for compilation on ARM processors (ARMv6 and ARMv7) using GNU GCC on Linux.
* [`FourQ_ARM_side_channel/makefile_Cortex-M4`](makefile_Cortex-M4): Makefile for compilation on ARM Cortex-M4 (STM32F4xx series) 
using GNU GCC on Linux.
* Main .c and .h files: library and header files. Public API for ECC scalar multiplication, key exchange and signatures is in 
[`FourQ_ARM_side_channel/FourQ_api.h`](FourQ_api.h).        
* [`FourQ_ARM_side_channel/ARM/`](ARM/): folder with library files implementing low-level arithmetic for ARM.
* [`FourQ_ARM_side_channel/libopencm3/`](libopencm3/): folder with firmware library files for ARM Cortex-M microcontrollers.
* [`FourQ_ARM_side_channel/random/`](random/): folder with pseudo-random generation function for ARM Cortex-M4.
* [`FourQ_ARM_side_channel/tests/`](tests/): test files for 32-bit ARM.
* [`FourQ_ARM_side_channel/tests_Cortex-M4/`](tests_Cortex-M4/): test files for ARM Cortex-M4.
* [`FourQ_ARM_side_channel/README.md`](README.md): this readme file.

`stm32f4_wrapper.c` and `stm32f4_wrapper.h` are by Joost Rijneveld and can be found 
[`here`](https://github.com/joostrijneveld/STM32-getting-started).

Files in the [`libopencm3`](libopencm3/) folder are from the [libopencm3 project](https://github.com/libopencm3/libopencm3).

## Supported platforms

This implementation is supported on ARM platforms and includes two variants: 

*  Implementation for ARM processors based on ARMv6 and ARMv7 architectures. This implementation was optimized
     for a first generation Raspberry Pi using a 700 MHz ARM1176JZF-S processor (ARMv6 architecture).
* Implementation for ARM Cortex-M4 processors based on the ARMv7-M architecture. This implementation was 
     developed and optimized on a STM32F4Discovery development board containing a Cortex-M4 STM32F407VG microcontroller (ARMv7-M architecture). It should be possible to extend the support to Cortex-M3 and Cortex-M7 based devices with small modifications.   

See instructions below to choose an implementation option and compile on one of the supported platforms.
   
## Complementary crypto functions

Random values are generated with `/dev/urandom` in the case of the 32-bit ARM implementation, and with the function
`random_int()` in the case of the ARM Cortex-M4 implementation.
  
The library includes an implementation of SHA-512 which is used by default by SchnorrQ signatures.

Users can experiment with different options by replacing functions in the `random`, `FourQ_ARM/random` and `sha512` folders 
and applying the corresponding changes to the settings in [`FourQ.h`](FourQ.h). 

## Instructions

### Building the library for ARMv6 or ARMv7

To compile on Linux using the GNU GCC compiler or the clang compiler, execute the following command from the command prompt:

```sh
$ make CC=[gcc/clang] EXTENDED_SET=[TRUE/FALSE]
```

After compilation, run `fp_tests`, `ecc_tests` or `crypto_tests`.

By default GNU GCC is used, as well as the extended settings. 
For example, to compile using GNU GCC, execute:

```sh
$ make
```

As another example, to compile using clang, execute:

```sh
$ make CC=clang
```

By default `EXTENDED_SET` is enabled, which sets the following compilation flags: `-fwrapv -fomit-frame-pointer
-funroll-loops`. To disable this, use `EXTENDED_SET=FALSE`.
Users are encouraged to experiment with the different flag options.

### Building the library for Cortex-M4 on the STM32F4DISCOVERY board

The following instructions have been tested on a Ubuntu 16.04 Linux machine.

First, install the ARM GNU GCC cross-compiler on the server machine:

```sh
$ sudo apt-get install gcc-arm-none-eabi libc6-dev-i386
```

Then, download, build and install [`stlink`](https://github.com/texane/stlink).

```sh
$ sudo apt-get install libusb-1.0-0-dev
$ git clone https://github.com/texane/stlink.git
$ cd stlink
$ make
$ cd build/Release/ && sudo make install
```

To compile the code, execute the following command from the FourQ_ARM_side_channel folder on the server machine:
 
```sh
$ make -f makefile_Cortex-M4
```

Power the STM32F4DISCOVERY board (with a USB to mini-USB cable) and connect it to the server machine via a 
USB-TTL converter as follows:

```sh
$ VDD -> VDD
GND -> GND 
TX  -> PA3 
RX  -> PA2 
```

Then, run from the server machine:

```sh
$ sudo ./tests_Cortex-M4/monitor.sh
```

From a different terminal window on the server machine, program the device with one of the following commands
from the `FourQ_ARM_side_channel` folder:

```sh
$ st-flash write tests_Cortex-M4/fp_tests.bin 0x8000000
$ st-flash write tests_Cortex-M4/ecc_tests.bin 0x8000000
$ st-flash write tests_Cortex-M4/crypto_tests.bin 0x8000000
```

The tests should begin to run on the first terminal window.

## Additional side-channel countermeasure

Some attacks try to target potential leakage when manipulating precomputed values during the scalar multiplication. 
To increase the resilience against this class of attacks, it is recommended to randomize the full table before extracting a point.

This countermeasure can be enabled in the implementation by uncommenting ``#define FULL_TABLE_RANDOMIZATION`` in [`FourQ.h`](FourQ.h). 

Note that this countermeasure is relatively expensive, so there is a security/performance trade-off to consider.
