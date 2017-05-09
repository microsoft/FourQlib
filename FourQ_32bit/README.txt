
                                        FourQlib v3.0 (C Edition)
                                        =========================
									      32-bit implementation
									 ===============================
 
1. CONTENTS:
   --------

The "FourQ_32bit" folder contains:

Visual Studio/FourQ/            - Folder with Visual Studio 2015 solution and project files for compilation
                                  in Windows.
Visual Studio/fp_tests/         - Folder with Visual Studio project files for testing field arithmetic functions
                                  in Windows.
Visual Studio/ecc_tests/        - Folder with Visual Studio project files for testing ECC functions in Windows.
Visual Studio/crypto_tests/     - Folder with Visual Studio project files for testing cryptographic functions,
                                  specifically key exchange and signatures, in Windows.
makefile                        - Makefile for compilation using GNU GCC or clang compilers on Linux. 
*.c, *.h                        - Library and header files. Public API for ECC scalar multiplication, key
                                  exchange and signatures is located in FourQ_api.h        
generic/                        - Folder with library files for 32-bit implementation.
random/                         - Folder with pseudo-random generation functions.
sha512/                         - Folder with SHA-512 implementation.
tests/                          - Test files.
README.txt                      - This readme file.


2. SUPPORTED PLATFORMS:
   -------------------

This implementation is supported on 32-bit platforms such as x86 and ARM-based processors running Windows or 
Linux OS. We have tested the library with Microsoft Visual Studio 2015, GNU GCC v4.9 and clang v3.8. 
See instructions below to choose an implementation option and compile on one of the supported platforms.


3. COMPLEMENTARY CRYPTO FUNCTIONS:
   ------------------------------

Random values are generated with /dev/urandom in the case of Linux, and with the function BCryptGenRandom() 
in the case of Windows.
  
The library includes an implementation of SHA-512 which is used by default by SchnorrQ signatures.

Users can experiment with different options by replacing functions in the folders "random" and "sha512" and 
applying the corresponding changes to the settings in FourQ.h. 


4. INSTRUCTIONS FOR WINDOWS OS:
   ---------------------------

BUILDING THE LIBRARY WITH VISUAL STUDIO:
---------------------------------------

Open the solution file (FourQ.sln) in Visual Studio 2015, select the "Generic" configurations from the
Solution Configurations menu (Win32 should appear as Solution Platform). 

By default, USE_ENDO=true is defined. To modify this configuration, go to the property window of the FourQ 
project, go to Configuration Properties > C/C++ > Preprocessor. Make any suitable changes, e.g., USE_ENDO=true
or false. Repeat these steps for the fp_tests, ecc_tests and crypto_tests projects.

Finally, select "Build Solution" from the "Build" menu. 

RUNNING THE TESTS:
-----------------

After building the solution, run fp_tests.exe, ecc_tests.exe and crypto_tests.exe.

USING THE LIBRARY:
-----------------

After building the solution, add the FourQ.lib file to the set of References for a project, and add FourQ.h 
and FourQ_api.h to the list of Header Files of a project.


5. INSTRUCTIONS FOR LINUX OS:
   -------------------------

BUILDING THE LIBRARY AND EXECUTING THE TESTS WITH GNU GCC OR CLANG:
------------------------------------------------------------------

To compile on Linux using the GNU GCC compiler or the clang compiler, execute the following command from the 
command prompt:
 
make ARCH=[x86/ARM] CC=[gcc/clang] USE_ENDO=[TRUE/FALSE] EXTENDED_SET=[TRUE/FALSE] CACHE_MEM=[TRUE/FALSE]

After compilation, run fp_tests, ecc_tests or crypto_tests.

By default GNU GCC is used, as well as endomorphisms and extended settings. Similarly, CACHE_MEM=TRUE is set
by default indicating that the targeted platform contains a cache memory.

For example, to compile using clang with the efficient endomorphisms on an x86 machine, execute:

make ARCH=x86 CC=clang

As another example, to compile using GNU GCC with the efficient endomorphisms on an ARM machine, execute:

make ARCH=ARM

By default EXTENDED_SET is enabled, which sets the following compilation flags: -fwrapv -fomit-frame-pointer 
-march=native. To disable this, use EXTENDED_SET=FALSE.
Users are encouraged to experiment with the different flag options.