#
# ( ) test:         Run unit tests
# ( ) test_timer:   Compile the unit tests
# ( ) libperf.a:
# ( ) clean:        Clean directory
# ( ) dep:          Update dependencies
# ( ) install:      Install lib, include, and bin files
# 

TOP = ../..
SRC = ..

TEST_OBJS = test_timer.o
OBJ       = 

#=======================================================================
INC_LOCAL  = -I$(SRC)/data
OBJ_EXTERN = $(SRC)/data/array.o

OBJS = $(OBJ) $(OBJ_EXTERN)
#=======================================================================

#=======================================================================
# TARGETS
#=======================================================================

#-----------------------------------------------------------------------

.PHONY: test
test: test_timer
	test_timer

#-----------------------------------------------------------------------

test_timer: $(TEST_OBJS) $(OBJS)
	$(CPP) -o test_timer $(TEST_OBJS) $(OBJS) $(INC)

time_matvec:  time_matvec.o $(OBJS)
	$(CPP) -o time_matvec  time_matvec.o $(OBJS) $(INC)


#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-rm *.o *~ test_timer time_matvec out.make.dep Makefile.bak lib*.a *.unit

#-----------------------------------------------------------------------

.PHONY: install
install: install-include test_timer $(OBJS)
	cp ./test_timer $(TOP)/bin
	ar rcu $(TOP)/lib/libcello.a $(OBJS)

#-----------------------------------------------------------------------

.PHONY: install-include

install-include:
	ln -sf `pwd`/timer.hpp  $(TOP)/include/perf_timer_.hpp

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/Make.include




