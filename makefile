SHELL := /bin/bash
PLATFORM := $(shell uname -s)

CC = gcc
BLDFLAGS = -Wall -Wextra -mavx -std=c11
CFLAGS = -D__USE_POSIX -D__STDC_CONSTANT_MACROS -D__STDINT_MACROS -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE=1 -O3
INCLUDES = /usr/include
LIBS = -lm

signal_fn_gz = 10k.786_master_list_v011117a.starch.bed.gz
signal_fn = 10k.786_master_list_v011117a.starch.bed
bs_non_avx_fn = 10k.786_master_list_v011117a.starch.bed.non_avx.bs
bs_avx_fn = 10k.786_master_list_v011117a.starch.bed.avx.bs
binary = pearson-test

ifeq ($(PLATFORM),Darwin)
	CC = clang
	FLAGS += -Weverything
endif

all: extract test

extract:
	gunzip -c ${signal_fn_gz} > ${signal_fn}

test: build test-non-avx test-avx

build:
	@echo "build..."
	$(CC) -g $(BLDFLAGS) $(CFLAGS) -c pearson-test.c -o pearson-test.o
	$(CC) -g $(BLDFLAGS) $(CFLAGS) -I$(INCLUDES) pearson-test.o -o ${binary} $(LIBS)

test-non-avx: build
	@echo "test (non-AVX)..."
	time -p ./${binary} --input ${signal_fn} > ${bs_non_avx_fn}

test-avx: build
	@echo "test (AVX)..."
	time -p ./${binary} --avx --input ${signal_fn} > ${bs_avx_fn}

clean:
	@echo "clean..."
	rm -f ${signal_fn}
	rm -f ${bs_non_avx_fn}
	rm -f ${bs_avx_fn}
	rm -f pearson-test.o
	rm -f ${binary}
