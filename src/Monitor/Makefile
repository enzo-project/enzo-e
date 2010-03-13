
TOP     = ../..
PACKAGE = monitor
OBJS    = monitor.o pngwriter.o
TESTS   = test_monitor


include $(TOP)/src/MAKE/Make.include


CPPFLAGS_LOCAL = -DNO_FREETYPE

#=======================================================================
# TARGETS
#=======================================================================

test_monitor: test_monitor.o
	$(CXX) -o test_monitor test_monitor.o $(LIB) -lpng $(LIB_PNG)

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/monitor.hpp $(INCLUDE)
	ln -sf `pwd`/monitor_monitor.hpp $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:

