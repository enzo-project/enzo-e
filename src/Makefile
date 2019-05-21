#-----------------------------------------------------------------------
#
# @file    Makefile
# @brief   Top-level makefile
# @author  James Bordner
#
#-----------------------------------------------------------------------

#=======================================================================
# SCONS HOOKS
#=======================================================================

.PHONY: all clean

all:
	scons -u
clean:
	scons -u -c

#=======================================================================
# DOXYGEN
#=======================================================================

.PHONY: dox dox-html dox-latex dox-xml

dox: dox-html dox-latex dox-xml

dox-html:
	doxygen doxygen/Doxyfile-html
dox-xml:
	doxygen doxygen/Doxyfile-xml
dox-latex:
	doxygen doxygen/Doxyfile-latex
	@echo;echo "SEE Makefile FOR MANUAL UPDATES TO latex/refman.tex"; echo


#=======================================================================
# CCCC C and C++ Code Counter (output in .cccc/)
#=======================================================================

.PHONY: cccc
cccc:
	cccc */*.cpp
