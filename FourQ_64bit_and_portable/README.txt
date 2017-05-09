
                                        FourQlib v3.0 (C Edition)
                                        =========================
								   64-bit and portable implementation
								   ==================================
 
1. CONTENTS:
   --------

The "FourQ_64bit_and_portable" folder contains:

Visual Studio/FourQ/            - Folder with Visual Studio 2015 solution and project files for compilation
                                  in Windows.
Visual Studio/fp_tests/         - Folder with Visual Studio project files for testing field arithmetic functions
                                  in Windows.
Visual Studio/ecc_tests/        - Folder with Visual Studio project files for testing ECC functions in Windows.
Visual Studio/crypto_tests/     - Folder with Visual Studio project files for testing cryptographic functions,
                                  specifically key exchange and signatures, in Windows.
makefile                        - Makefile for compilation using GNU GCC or clang compilers on Linux. 
*.c, *.h                        - Library and header files. Public API for ECC scalar multiplication, key
                                  exchange and signatures is located in FourQ_api.h.                 
AMD64/                          - Folder with library files for optimized x64 implementation.             
ARM64/                          - Folder with library files for optimized 64-bit ARM implementation.
generic/                        - Folder with library files for portable implementation.
tests/                          - Test files.
README.txt                      - This readme file.


2. SUPPORTED PLATFORMS:
   -------------------

This implementation is supported in a wide range of platforms including x64, x86, 32-bit ARM and 64-bit ARM,
running Windows or Linux OS. We have tested the library with Microsoft Visual Studio 2015, GNU GCC v4.9 and 
clang v3.8. See instructions below to choose an implementation option and compile on one of the supported 
platforms. 


3. COMPLEMENTARY CRYPTO FUNCTIONS:
   ------------------------------

Random values are generated with /dev/urandom in the case of Linux, and with the function BCryptGenRandom() 
in the case of Windows.
  
The library includes an implementation of SHA-512 which is used by default by SchnorrQ signatures.

Users can experiment with different options by replacing functions in the folders "random" and "sha512" and 
applying the corresponding changes to the settings in FourQ.h. 


4. IMPLEMENTATION OPTIONS:
   ----------------------

The following compilation options are available for the "FourQ_64bit_and_portable" implementation:

- A portable implementation (enabled by the "GENERIC" option).

- Optimized implementations for x64 and 64-bit ARM (ARMv8). Note that the rest of platforms are only supported
  by the generic implementation. 

- Use of AVX or AVX2 instructions enabled by defining _AVX_ or _AVX2_ (Windows) or by the "AVX" and "AVX2" 
  options (Linux).

- Optimized x64 assembly implementations in Linux.

- Use of fast endomorphisms enabled by the "USE_ENDO" option.

Follow the instructions in Subsection 3.1, "INSTRUCTIONS FOR WINDOWS OS" or Subsection 3.2, "INSTRUCTIONS FOR 
LINUX OS", to configure these different options.


4.1 INSTRUCTIONS FOR WINDOWS OS:
    ---------------------------

BUILDING THE LIBRARY WITH VISUAL STUDIO:
---------------------------------------

Open the solution file (FourQ.sln) in Visual Studio 2015, select one of the available configurations from
the Solution Configurations menu ("Release" corresponding to the high-speed x64 implementation and "Generic" 
corresponding to the portable implementation) and select one of the Solution Platforms (x64 or Win32). Note 
that Win32 is only supported with the "Generic" solution configuration.

By default, USE_ENDO=true and (for x64) _AVX_ is defined. To modify this configuration, go to the property 
window of the FourQ project, go to Configuration Properties > C/C++ > Preprocessor. Make any suitable chan-
ges, e.g., delete _AVX_ if AVX instructions are not supported, replace _AVX_ by _AVX2_ if AVX2 instructions
are supported, or set USE_ENDO=true or false. Repeat these steps for the fp_tests, ecc_tests and crypto_tests
projects.

Finally, select "Build Solution" from the "Build" menu. 

RUNNING THE TESTS:
-----------------

After building the solution, run fp_tests.exe, ecc_tests.exe and crypto_tests.exe.

USING THE LIBRARY:
-----------------

After building the solution, add the FourQ.lib file to the set of References for a project, and add FourQ.h 
and FourQ_api.h to the list of Header Files of a project.


4.2. INSTRUCTIONS FOR LINUX OS:
     -------------------------

BUILDING THE LIBRARY AND EXECUTING THE TESTS WITH GNU GCC OR CLANG:
------------------------------------------------------------------

To compile on Linux using the GNU GCC compiler or the clang compiler, execute the following command from the 
command prompt:
 
make ARCH=[x64/x86/ARM/ARM64] CC=[gcc/clang] ASM=[TRUE/FALSE] AVX=[TRUE/FALSE] AVX2=[TRUE/FALSE] 
     EXTENDED_SET=[TRUE/FALSE] USE_ENDO=[TRUE/FALSE] GENERIC=[TRUE/FALSE] SERIAL_PUSH=[TRUE/FALSE] 

After compilation, run fp_tests, ecc_tests or crypto_tests.

By default GNU GCC is used, as well as the endomorphisms and the extended settings.

In the case of x64, AVX2 instructions and the high-speed assembly implementation are enabled by default.
In the case of x86 and ARM, the portable ("GENERIC") implementation is used by default.

For example, to compile the optimized x64 implementation in assembly with GNU GCC using the efficient
endomorphisms on a machine with AVX2 support (e.g, Intel's Haswell or Broadwell), execute:

make ARCH=x64

For example, to compile the optimized ARM64 implementation with GNU GCC using the efficient endomorphisms, 
execute:

make ARCH=ARM64

As another example, to compile the portable implementation with clang using the efficient endomorphisms 
on an x86 machine, execute:

make ARCH=x86 CC=clang

SERIAL_PUSH can be enabled in some platforms (e.g., AMD without AVX2 support) to boost performance.

By default EXTENDED_SET is enabled, which sets the following compilation flags: -fwrapv -fomit-frame-pointer 
-march=native. To disable this, use EXTENDED_SET=FALSE.
Users are encouraged to experiment with the different flag options.

Whenever an unsupported configuration is applied, the following message will be displayed: #error -- "Unsu-
pported configuration". For example, the use of assembly or any of the AVX options is not supported when se-
lecting the portable implementation (i.e., if GENERIC=TRUE or if ARCH=[x86/ARM]). 