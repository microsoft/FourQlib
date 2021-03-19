####  Makefile for compilation using GNU GCC or clang  ####

OPT=-O3     # Optimization option by default

CC=gcc
ifeq "$(CC)" "gcc"
    COMPILER=gcc
else ifeq "$(CC)" "clang"
    COMPILER=clang
endif

ifeq "$(ARCH)" "x64"
    ARCHITECTURE=_AMD64_
ifeq "$(GENERIC)" "TRUE"
    USE_GENERIC=-D _GENERIC_
endif
ifeq "$(ASM)" "FALSE"
else ifeq "$(ASM)" "TRUE"
    USE_ASM=-D _ASM_
    ASM_var=yes
else
ifneq "$(GENERIC)" "TRUE"
    USE_ASM=-D _ASM_
    ASM_var=yes
endif
endif
ifeq "$(AVX)" "FALSE"
else ifeq "$(AVX)" "TRUE"
    USE_AVX=-D _AVX_
    SIMD=-mavx
else
ifneq "$(GENERIC)" "TRUE"
    USE_AVX=-D _AVX_
    SIMD=-mavx
endif  	
endif
ifeq "$(AVX2)" "FALSE"
else ifeq "$(AVX2)" "TRUE"
    USE_AVX2=-D _AVX2_
    SIMD=-mavx2
    AVX2_var=yes
else
ifneq "$(GENERIC)" "TRUE"
    USE_AVX2=-D _AVX2_
    SIMD=-mavx2
    AVX2_var=yes
endif  	
endif

else ifeq "$(ARCH)" "ARM64"
    ARCHITECTURE=_ARM64_
    ARM_SETTING=-lrt
ifeq "$(GENERIC)" "TRUE"
    USE_GENERIC=-D _GENERIC_
endif

else

USE_GENERIC=-D _GENERIC_
ifeq "$(GENERIC)" "FALSE"
    USE_GENERIC=
endif
ifeq "$(ASM)" "TRUE"
    USE_ASM=-D _ASM_
endif
ifeq "$(AVX)" "TRUE"
    USE_ASM=-D _ASM_
endif
ifeq "$(AVX2)" "TRUE"
    USE_ASM=-D _ASM_
endif
ifeq "$(ARCH)" "x86"
    ARCHITECTURE=_X86_
else ifeq "$(ARCH)" "ARM"
    ARCHITECTURE=_ARM_
    ARM_SETTING=-lrt
endif
endif

ADDITIONAL_SETTINGS=-fwrapv -fomit-frame-pointer -march=native
ifeq "$(EXTENDED_SET)" "FALSE"
    ADDITIONAL_SETTINGS=
endif

USE_ENDOMORPHISMS=-D USE_ENDO
ifeq "$(USE_ENDO)" "FALSE"
    USE_ENDOMORPHISMS=
endif

ifeq "$(SERIAL_PUSH)" "TRUE"
    USE_SERIAL_PUSH=-D PUSH_SET
endif

SHARED_LIB_TARGET=libFourQ.so
ifeq "$(SHARED_LIB)" "TRUE"
    DO_MAKE_SHARED_LIB=-fPIC
	SHARED_LIB_O=$(SHARED_LIB_TARGET)
endif

cc=$(COMPILER)
CFLAGS=-c $(OPT) $(ADDITIONAL_SETTINGS) $(SIMD) -D $(ARCHITECTURE) -D __LINUX__ $(USE_AVX) $(USE_AVX2) $(USE_ASM) $(USE_GENERIC) $(USE_ENDOMORPHISMS) $(USE_SERIAL_PUSH) $(DO_MAKE_SHARED_LIB)
LDFLAGS=
ifdef ASM_var
ifdef AVX2_var
    ASM_OBJECTS=fp2_1271_AVX2.o
else
    ASM_OBJECTS=fp2_1271.o
endif 
endif
OBJECTS=eccp2.o eccp2_no_endo.o eccp2_core.o $(ASM_OBJECTS) crypto_util.o schnorrq.o hash_to_curve.o kex.o sha512.o random.o 
OBJECTS_FP_TEST=fp_tests.o $(OBJECTS) test_extras.o 
OBJECTS_ECC_TEST=ecc_tests.o $(OBJECTS) test_extras.o 
OBJECTS_CRYPTO_TEST=crypto_tests.o $(OBJECTS) test_extras.o 
OBJECTS_ALL=$(OBJECTS) $(OBJECTS_FP_TEST) $(OBJECTS_ECC_TEST) $(OBJECTS_CRYPTO_TEST)

all: crypto_test ecc_test fp_test $(SHARED_LIB_O)

ifeq "$(SHARED_LIB)" "TRUE"
    $(SHARED_LIB_O): $(OBJECTS)
	    $(CC) -shared -o $(SHARED_LIB_O) $(OBJECTS)
endif

crypto_test: $(OBJECTS_CRYPTO_TEST)
	$(CC) -o crypto_test $(OBJECTS_CRYPTO_TEST) $(ARM_SETTING)

ecc_test: $(OBJECTS_ECC_TEST)
	$(CC) -o ecc_test $(OBJECTS_ECC_TEST) $(ARM_SETTING)

fp_test: $(OBJECTS_FP_TEST)
	$(CC) -o fp_test $(OBJECTS_FP_TEST) $(ARM_SETTING)

eccp2_core.o: eccp2_core.c AMD64/fp_x64.h
	$(CC) $(CFLAGS) eccp2_core.c

eccp2.o: eccp2.c
	$(CC) $(CFLAGS) eccp2.c

eccp2_no_endo.o: eccp2_no_endo.c
	$(CC) $(CFLAGS) eccp2_no_endo.c
    
ifdef ASM_var
ifdef AVX2_var
    AMD64/consts.s: AMD64/consts.c
	    $(CC) $(CFLAGS) -S -o $@ $<
	    sed '/.globl/d' -i $@
    fp2_1271_AVX2.o: AMD64/fp2_1271_AVX2.S AMD64/consts.s
	    $(CC) $(CFLAGS) -o $@ $<
else
    fp2_1271.o: AMD64/fp2_1271.S
	    $(CC) $(CFLAGS) AMD64/fp2_1271.S
endif
endif

schnorrq.o: schnorrq.c
	$(CC) $(CFLAGS) schnorrq.c

hash_to_curve.o: hash_to_curve.c
	$(CC) $(CFLAGS) hash_to_curve.c

kex.o: kex.c
	$(CC) $(CFLAGS) kex.c

crypto_util.o: crypto_util.c
	$(CC) $(CFLAGS) crypto_util.c

sha512.o: ../sha512/sha512.c
	$(CC) $(CFLAGS) ../sha512/sha512.c

random.o: ../random/random.c
	$(CC) $(CFLAGS) ../random/random.c

test_extras.o: tests/test_extras.c
	$(CC) $(CFLAGS) tests/test_extras.c

crypto_tests.o: tests/crypto_tests.c
	$(CC) $(CFLAGS) tests/crypto_tests.c

ecc_tests.o: tests/ecc_tests.c
	$(CC) $(CFLAGS) tests/ecc_tests.c

fp_tests.o: tests/fp_tests.c
	$(CC) $(CFLAGS) tests/fp_tests.c

.PHONY: clean

clean:
	rm -rf $(SHARED_LIB_TARGET) crypto_test ecc_test fp_test *.o AMD64/consts.s

