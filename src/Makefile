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

.PHONY: doc doc-html doc-latex doc-xml

doc: doc-html doc-latex doc-xml

doc-html:
	doxygen doxygen/Doxyfile-html
doc-xml:
	doxygen doxygen/Doxyfile-xml
doc-latex:
	doxygen doxygen/Doxyfile-latex
	@echo;echo "SEE Makefile FOR MANUAL UPDATES TO latex/refman.tex"; echo


#=======================================================================
# CCCC C and C++ Code Counter (output in .cccc/)
#=======================================================================

.PHONY: cccc
cccc:
	cccc */*.cpp
