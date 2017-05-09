
                                        FourQlib v3.0 (C Edition)
                                        =========================
					       Optimized implementation for 32-bit ARM using NEON
	                       ==================================================
 
1. CONTENTS:
   --------

The "FourQ_ARM_NEON" folder contains:

makefile                        - Makefile for compilation on ARM processors with NEON support using 
                                  GNU GCC or clang on Linux. 
*.c, *.h                        - Library and header files. Public API for ECC scalar multiplication, key
                                  exchange and signatures is located in FourQ_api.h        
ARM/                            - Folder with library files implementing low-level arithmetic for ARM. 
tests/                          - Test files.
README.txt                      - This readme file.


2. SUPPORTED PLATFORMS:
   -------------------

This implementation is supported on 32-bit ARM platforms that contain the SIMD engine called NEON and run Linux. 
For example, platforms with a NEON engine include many cores with the ARMv7 architecture. The implementation 
has been optimized for ARM Cortex-A7, Cortex-A8, Cortex-A9 and Cortex-A15 based processors.

See instructions below to choose an implementation option and compile on one of the supported platforms.


3. COMPLEMENTARY CRYPTO FUNCTIONS:
   ------------------------------

Random values are generated with /dev/urandom.
  
The library includes an implementation of SHA-512 which is used by default by SchnorrQ signatures.

Users can experiment with different options by replacing functions in the folders "random" and "sha512" and 
applying the corresponding changes to the settings in FourQ.h. 


4. INSTRUCTIONS TO BUILD THE LIBRARY AND EXECUTE THE TESTS WITH GNU GCC OR CLANG:
     ---------------------------------------------------------------------------

To compile on Linux using the GNU GCC compiler or the clang compiler, execute the following command from the 
command prompt:
 
make CC=[gcc/clang] USE_ENDO=[TRUE/FALSE] EXTENDED_SET=[TRUE/FALSE] INTERLEAVE=[TRUE/FALSE] 
     MIX_ARM_NEON=[TRUE/FALSE]

After compilation, run fp_tests, ecc_tests or crypto_tests.

By default GNU GCC is used, as well as endomorphisms and extended settings. 

There are two special optimizations that can be exploited. INTERLEAVE improves performance on some platforms
by interleaving load/store instructions with other non-memory instructions. This optimization is recommended
for Cortex-A7, Cortex-A8 and Cortex-A9. MIX_ARM_NEON improves performance on some platforms by mixing ARM and
NEON instructions. This optimization is recommended for Cortex-A7, Cortex-A9 and Cortex-A15.

By default, INTERLEAVE is turned off and MIX_ARM_NEON is turned on.

For example, to compile using GNU GCC with the efficient endomorphisms on an ARM Cortex-A15 device, execute:

make

As another example, to compile using clang with the efficient endomorphisms on an ARM Cortex-A8 device, 
execute:

make CC=clang INTERLEAVE=TRUE MIX_ARM_NEON=FALSE

By default EXTENDED_SET is enabled, which sets the following compilation flags: -fwrapv -fomit-frame-pointer 
-funroll-loops. To disable this, use EXTENDED_SET=FALSE.
Users are encouraged to experiment with the different flag options.