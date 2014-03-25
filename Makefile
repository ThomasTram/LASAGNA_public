#Makefile for LASAGNA
#Thomas Tram, 25.08.2011
#tram@phys.au.dk

MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build
OUTPUT = $(MDIR)/output
.base:
	if ! [ -d $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base
	if ! [ -d $(OUTPUT) ]; then mkdir $(OUTPUT) ; fi;
	touch output/.base

vpath %.c source:tools:test:main
vpath %.o build
vpath .base build

use_superlu=no
ifeq ($(use_superlu),yes)
SuperLU_root = /media/sf_Shared/Neutrinos/SuperLU_MT_2.0
include $(SuperLU_root)/make.inc
HEADERSLU=$(SuperLU_root)/SRC
LIBSLU = $(SuperLU_root)/lib
endif

#CC = gcc
CC = icc

ifeq ($(use_superlu),yes)
INCLUDES = ../include -I$(HEADERSLU)
CCFLAG = $(CFLAGS) $(CDEFS) $(BLASDEF) -D _SUPERLU
LDFLAG = $(LOADOPTS) -D _SUPERLU
else
INCLUDES = ../include 
CCFLAG   = -O3 -axavx -msse3 -g $(DEFSLU)
LDFLAG   = -O3 -axavx -msse3 -g $(DEFSLU)
#CCFLAG   = -O4 -Wall -g --fast-math
#LDFLAG   = -O4 -Wall -g --fast-math
#CCFLAG   = -O1 -g
#LDFLAG   = -O1 -g
endif

%.o:  %.c .base
	cd $(WRKDIR);$(CC) $(CCFLAG) $(CDEFS) $(BLASDEF) -I$(INCLUDES) -c ../$< -o $*.o

ifeq ($(use_superlu),yes)
EVO_TOOLS  = multimatrix.o evolver_common.o sparse.o linalg_wrapper_dense_NR.o linalg_wrapper_sparse.o linalg_wrapper_SuperLU.o
else
EVO_TOOLS  = multimatrix.o evolver_common.o sparse.o linalg_wrapper_dense_NR.o linalg_wrapper_sparse.o
EXTRA_FILES = tools/linalg_wrapper_SuperLU.c include/linalg_wrapper_SuperLU.h
endif 
IO_TOOLS = mat_io.o parser.o
TOOLS = $(IO_TOOLS) $(EVO_TOOLS) newton.o evolver_ndf15.o  arrays.o evolver_rk45.o evolver_radau5.o  

TEST_WRAPPER_DENSE = test_wrapper_dense.o

TEST_WRAPPER_SPARSE = test_wrapper_sparse.o

LASAGNA = lasagna.o

LASAGNA_LYA = lasagna_lya.o

QKE_EQUATIONS = qke_equations.o

LYA_EQUATIONS = lya_equations.o

BACKGROUND = background.o

TEST_RKODE = test_rkode.o

TEST_MATIO = test_matio.o

TEST_PROFILE = test_profile.o

EXTRACT_MATRIX = extract_matrix.o

INPUT = input.o

LYA_INPUT = lya_input.o

C_TOOLS =  $(addprefix tools/, $(addsuffix .c,$(basename $(TOOLS))))
C_SOURCE = $(addprefix source/, $(addsuffix .c,$(basename $(QKE_EQUATIONS) $(LYA_EQUATIONS) $(BACKGROUND) $(INPUT) $(LYA_INPUT) )))
C_TEST = $(addprefix test/, $(addsuffix .c,$(basename $(TEST_MATIO) $(TEST_PROFILE) $(TEST_RKODE) )))
C_MAIN = $(addprefix main/, $(addsuffix .c,$(basename $(LASAGNA) $(LASAGNA_LYA) $(EXTRACT_MATRIX))))
C_ALL = $(C_MAIN) $(C_TOOLS) $(C_SOURCE) $(C_TEST)
H_ALL = $(addprefix include/, common.h $(addsuffix .h, $(basename $(notdir $(C_ALL)))))
MISC_FILES = make_loop_dir.sh main/prepare_job.c test/test_wrapper_sparse.c test/test_wrapper_dense.c load_and_plot.m lepton_number.m evolve_in_time.m dsdofHP_B.dat parameters.ini SuperLUpatch.tar.gz README.txt Makefile

all: lasagna lasagna_lya extract_matrix

ifeq ($(use_superlu),yes)
LINKSLU = $(LIBSLU)/$(SUPERLULIB) $(BLASLIB) $(MPLIB)
else
LINKSLU =
endif
lasagna: $(TOOLS) $(QKE_EQUATIONS) $(BACKGROUND) $(INPUT) $(LASAGNA)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) $(LINKSLU) -lm

lasagna_lya: $(TOOLS) $(QKE_EQUATIONS) $(LYA_EQUATIONS) $(BACKGROUND) $(INPUT) $(LYA_INPUT) $(LASAGNA_LYA)	
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) $(LINKSLU) -lm

extract_matrix: $(IO_TOOLS) $(EXTRACT_MATRIX)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_profile: $(TOOLS) $(QKE_EQUATIONS) $(BACKGROUND) $(TEST_PROFILE)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_rkode:$(TOOLS) $(TEST_RKODE)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_matio: $(TOOLS) $(TEST_MATIO)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_wrapper_dense: $(EVO_TOOLS)$ $(TEST_WRAPPER_DENSE) 
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_wrapper_sparse: $(EVO_TOOLS)$ $(TEST_WRAPPER_SPARSE)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) $(LIBSLU)/$(SUPERLULIB) $(BLASLIB) $(MPLIB) -lm

tar:
	tar czvf lasagna_1.0.tar.gz $(C_ALL) $(H_ALL) $(MISC_FILES) $(EXTRA_FILES)

clean: .base
	rm -rf $(WRKDIR);
