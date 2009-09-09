#
# test         Run unit tests
# clean        Clean directory
# test_ppml    Compile the unit tests
# dep          Generate dependencies
# 

TOP = ../../..

INC_LOCAL  = -I$(TOP)/include
LIB_LOCAL  = -L$(TOP)/lib -lcello

OBJS = \
	cello_hydro.o \
	Grid_SolveMHDEquations.o \
	PPML_Conservative.o \
	PPML_HLLD.o \
	PPML_MAIN.o \
	PPML_monot.o \
	PPML_potokx.o \
	PPML_potoky.o \
	PPML_potokz.o \
	PPML_Primitive.o \
	PPML_QDD6.o \
	PPML_sub.o \
	PPML_TimeStep.o \
	test_ppml.o


FFLAGS_LOCAL   = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -O3
CXXFLAGS_LOCAL = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -O3

FFLAGS   = $(FFLAGS_LOCAL)
CXXFLAGS = $(CXXFLAGS_LOCAL) $(INC_LOCAL)
LIB      = $(LIB_LOCAL) -lgfortran -lm -lhdf5


test_ppml: $(OBJS) test_ppml.o
	$(CXX) -o test_ppml test_ppml.o $(OBJS) $(LIB)

obj: $(OBJS)

test: test_ppml
	test_ppml

.PHONY: clean
clean:
	-@rm *.o test_ppml *~

include $(TOP)/src/MAKE/Make.include

include DEPEND
