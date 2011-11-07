#Makefile for LASAGNA
#Thomas Tram, 25.08.2011

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

CC       = gcc

#CCFLAG   = -O1 -Wall -ggdb -pg -fopenmp
#LDFLAG   = -O1 -Wall -ggdb -pg -fopenmp
#CCFLAG   = -O0 -Wall -ggdb -fopenmp
#LDFLAG   = -O0 -Wall -ggdb -fopenmp
#CCFLAG   = -O0 -Wall -ggdb
#LDFLAG   = -O0 -Wall -ggdb
CCFLAG   = -O4 -Wall -fopenmp
LDFLAG   = -O4 -Wall -fopenmp

INCLUDES = ../include

%.o:  %.c .base
	cd $(WRKDIR);$(CC) $(CCFLAG) -I$(INCLUDES) -c ../$< -o $*.o

TOOLS = newton.o evolver_ndf15.o sparse.o arrays.o evolver_rk45.o mat_io.o

LASAGNA = lasagna.o

LASAGNA_LOOP = lasagna_loop.o

QKE_EQUATIONS = qke_equations.o

BACKGROUND = background.o

TEST_RKODE = test_rkode.o

TEST_MATIO = test_matio.o

TEST_PROFILE = test_profile.o

EXTRACT_MATRIX = extract_matrix.o

C_TOOLS =  $(addprefix tools/, $(addsuffix .c,$(basename $(TOOLS))))
C_SOURCE = $(addprefix source/, $(addsuffix .c,$(basename $(QKE_EQUATIONS) $(BACKGROUND))))
C_TEST = $(addprefix test/, $(addsuffix .c,$(basename $(TEST_MATIO) $(TEST_PROFILE) $(TEST_RKODE) )))
C_MAIN = $(addprefix main/, $(addsuffix .c,$(basename $(LASAGNA) $(LASAGNA_LOOP) $(EXTRACT_MATRIX))))
C_ALL = $(C_MAIN) $(C_TOOLS) $(C_SOURCE) $(C_TEST)
H_ALL = $(addprefix include/, common.h $(addsuffix .h, $(basename $(notdir $(C_ALL)))))
MISC_FILES = load_and_plot.m dsdofHP_B.dat Makefile

all: lasagna lasagna_loop extract_matrix

lasagna: $(TOOLS) $(QKE_EQUATIONS) $(BACKGROUND) $(LASAGNA)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

lasagna_loop: $(TOOLS) $(QKE_EQUATIONS) $(BACKGROUND) $(LASAGNA_LOOP)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_profile: $(TOOLS) $(QKE_EQUATIONS) $(BACKGROUND) $(TEST_PROFILE)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

extract_matrix: $(TOOLS) $(EXTRACT_MATRIX)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_rkode: $(TOOLS) $(TEST_RKODE)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_matio: $(TOOLS) $(TEST_MATIO)
	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

#test_2D_quadrature: $(TOOLS) $(TEST_2D_QUADRATURE)
#	$(CC) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

tar:
	tar czvf lasagna.tar.gz $(C_ALL) $(H_ALL) $(MISC_FILES)

clean: .base
	rm -rf $(WRKDIR);
