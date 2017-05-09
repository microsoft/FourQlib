####  Makefile for compilation using GNU GCC or clang on ARM-based processors without using NEON  ####

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
    INLINING_SETTINGS=-finline-functions -finline-limit=100
endif

MEM=-D _CACHE_MEM_
ifeq "$(CACHE_MEM)" "FALSE"
    MEM=-D _NO_CACHE_MEM_
endif

CFLAGS=-c $(OPT) $(ADDITIONAL_SETTINGS) -D _ARM_ -D __LINUX__ $(USE_ENDOMORPHISMS) $(MEM) $(INLINING_SETTINGS)
LDFLAGS=
OBJECTS=eccp2.o eccp2_no_endo.o fp2_1271_arm.o crypto_util.o schnorrq.o kex.o sha512.o random.o  
OBJECTS_FP_TEST=$(OBJECTS) fp_tests.o test_extras.o 
OBJECTS_ECC_TEST=ecc_tests.o test_extras.o $(OBJECTS) 
OBJECTS_CRYPTO_TEST=crypto_tests.o $(OBJECTS) test_extras.o 
OBJECTS_ALL=$(OBJECTS) $(OBJECTS_FP_TEST) $(OBJECTS_ECC_TEST) $(OBJECTS_CRYPTO_TEST)

all: crypto_test ecc_test fp_test

crypto_test: $(OBJECTS_CRYPTO_TEST)
	$(CC) -o crypto_test $(OBJECTS_CRYPTO_TEST) -lrt

ecc_test: $(OBJECTS_ECC_TEST)
	$(CC) -o ecc_test $(OBJECTS_ECC_TEST) -lrt

fp_test: $(OBJECTS_FP_TEST)
	$(CC) -o fp_test $(OBJECTS_FP_TEST) -lrt

eccp2.o: eccp2.c FourQ_internal.h ARM/fp_arm.h
	$(CC) $(CFLAGS) eccp2.c

eccp2_no_endo.o: eccp2_no_endo.c
	$(CC) $(CFLAGS) eccp2_no_endo.c
    
fp2_1271_arm.o: ARM/fp2_1271_arm.S
	$(CC) $(CFLAGS) ARM/fp2_1271_arm.S

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

fp_tests.o: tests/fp_tests.c FourQ_internal.h
	$(CC) $(CFLAGS) tests/fp_tests.c

.PHONY: clean

clean:
	rm -f crypto_test ecc_test fp_test fp2_1271_arm.o $(OBJECTS_ALL)

