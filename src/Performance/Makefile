#
# ( ) test:         Run unit tests
# ( ) test_timer:   Compile the unit tests
# ( ) libperformance.a:
# ( ) clean:        Clean directory
# ( ) dep:          Update dependencies
# ( ) install:      Install lib, include, and bin files
# 

TOP = ../..
SRC = ..

PACKAGE   = performance
TEST_OBJS = test_timer.o
OBJ       = 

#=======================================================================
INC_LOCAL  = 
OBJ_EXTERN = 

OBJS = $(OBJ) $(OBJ_EXTERN)
#=======================================================================

#=======================================================================
# TARGETS
#=======================================================================

.PHONY: all
all: test time_matvec

#-----------------------------------------------------------------------

.PHONY: test
test: test_timer
	./test_timer

#-----------------------------------------------------------------------

test_timer: $(TEST_OBJS) $(OBJS)
	$(CXX) -o test_timer $(TEST_OBJS) $(OBJS) $(OPT) $(INC)

PERFORMANCE_FLAGS = 
time_matvec: time_matvec.cpp
	$(CXX) -o time_matvec.o time_matvec.cpp $(PERFORMANCE_FLAGS) -O $(INC) -c
	$(CXX) -o $@            time_matvec.o   $(PERFORMANCE_FLAGS) -O $(INC)
time_matvec-O0:  time_matvec.cpp
	$(CXX) -o time_matvec.o time_matvec.cpp $(PERFORMANCE_FLAGS) -O0 $(INC) -c
	$(CXX) -o $@            time_matvec.o   $(PERFORMANCE_FLAGS) -O0 $(INC)
time_matvec-O1:    time_matvec.cpp
	$(CXX) -o time_matvec.o time_matvec.cpp $(PERFORMANCE_FLAGS) -O1 $(INC) -c
	$(CXX) -o $@            time_matvec.o   $(PERFORMANCE_FLAGS) -O1 $(INC)
time_matvec-O2:    time_matvec.cpp
	$(CXX) -o time_matvec.o time_matvec.cpp $(PERFORMANCE_FLAGS) -O2 $(INC) -c
	$(CXX) -o $@            time_matvec.o   $(PERFORMANCE_FLAGS) -O2 $(INC)
time_matvec-O3:    time_matvec.cpp
	$(CXX) -o time_matvec.o time_matvec.cpp $(PERFORMANCE_FLAGS) -O3 $(INC) -c
	$(CXX) -o $@            time_matvec.o   $(PERFORMANCE_FLAGS) -O3 $(INC)
time_matvec-O4:    time_matvec.cpp
	$(CXX) -o time_matvec.o time_matvec.cpp $(PERFORMANCE_FLAGS) -O4 $(INC) -c
	$(CXX) -o $@            time_matvec.o   $(PERFORMANCE_FLAGS) -O4 $(INC)
time_matvec-O5:    time_matvec.cpp
	$(CXX) -o time_matvec.o time_matvec.cpp $(PERFORMANCE_FLAGS) -O5 $(INC) -c
	$(CXX) -o $@            time_matvec.o   $(PERFORMANCE_FLAGS) -O5 $(INC)


#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-@rm *.o *~ test_timer time_matvec time_matvec-O? out.make.dep Makefile.bak lib*.a *.unit *.png *.eps

#-----------------------------------------------------------------------

.PHONY: install
install: install-include test_timer $(OBJS)
	cp ./test_timer $(TOP)/bin
	ar rcu $(TOP)/lib/libcello.a $(OBJS)
	ar rcu $(TOP)/lib/libcello-$(PACKAGE).a $(OBJS)   

#-----------------------------------------------------------------------

INCLUDE =$(TOP)/include

.PHONY: install-include
install-include:
	ln -sf `pwd`/performance.hpp        $(INCLUDE)
	ln -sf `pwd`/performance_timer_.hpp $(INCLUDE)

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include




