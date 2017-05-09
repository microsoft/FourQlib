####  Makefile for compilation using GNU GCC or clang on ARM-based processors with NEON support  ####

OPT=-O3     # Optimization option by default

CC=gcc
ifeq "$(CC)" "gcc"
    COMPILER=gcc
else ifeq "$(CC)" "clang"
    COMPILER=clang
endif

ADDITIONAL_SETTINGS=-fwrapv -fomit-frame-pointer -funroll-loops
ifeq "$(EXTENDED_SET)" "FALSE"
    ADDITIONAL_SETTINGS=
endif

USE_ENDOMORPHISMS=-D USE_ENDO
ifeq "$(USE_ENDO)" "FALSE"
    USE_ENDOMORPHISMS=
endif

INLINING_SETTINGS=
ifeq "$(CC)" "gcc"
    INLINING_SETTINGS=-finline-functions -finline-limit=300
endif

ifeq "$(INTERLEAVE)" "TRUE"
    USE_INTERLEAVING=-D _INTERLEAVE_
endif

USE_ARM_NEON_MIX=-D _MIX_ARM_NEON_
ifeq "$(MIX_ARM_NEON)" "FALSE"
    USE_ARM_NEON_MIX=
endif

CFLAGS=-c $(OPT) $(ADDITIONAL_SETTINGS) -D _ARM_ -D __LINUX__ -mfloat-abi=hard -mfpu=neon $(USE_ENDOMORPHISMS) $(USE_INTERLEAVING) $(USE_ARM_NEON_MIX) $(INLINING_SETTINGS)
LDFLAGS=
OBJECTS=eccp2.o fp2_1271_NEON.o eccp2_no_endo.o crypto_util.o schnorrq.o kex.o sha512.o random.o
OBJECTS_ECC_TEST=ecc_tests.o test_extras.o $(OBJECTS) 
OBJECTS_FP_TEST=$(OBJECTS) fp_tests.o test_extras.o 
OBJECTS_CRYPTO_TEST=crypto_tests.o $(OBJECTS) test_extras.o 
OBJECTS_ALL=$(OBJECTS) $(OBJECTS_ECC_TEST) $(OBJECTS_FP_TEST) $(OBJECTS_CRYPTO_TEST)

all: crypto_test ecc_test fp_test

crypto_test: $(OBJECTS_CRYPTO_TEST)
	$(CC) -o crypto_test $(OBJECTS_CRYPTO_TEST) -lrt

ecc_test: $(OBJECTS_ECC_TEST)
	$(CC) -o ecc_test $(OBJECTS_ECC_TEST) -lrt

fp_test: $(OBJECTS_FP_TEST)
	$(CC) -o fp_test $(OBJECTS_FP_TEST) -lrt

eccp2.o: eccp2.c
	$(CC) $(CFLAGS) eccp2.c

eccp2_no_endo.o: eccp2_no_endo.c
	$(CC) $(CFLAGS) eccp2_no_endo.c

fp2_1271_NEON.o: ARM/fp2_1271_NEON.c
	$(CC) $(CFLAGS) ARM/fp2_1271_NEON.c

schnorrq.o: schnorrq.c
	$(CC) $(CFLAGS) schnorrq.c

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
	rm -f crypto_test ecc_test fp_test fp2_1271_NEON.o $(OBJECTS_ALL)

