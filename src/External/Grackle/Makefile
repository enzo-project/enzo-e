#=======================================================================
#
# FILE:        Makefile
#
# SUMMARY:     Configurable Makefile for grackle
#              (originally written by James Border for Enzo)
#
# DESCRIPTION
#              See 'make help' for definitive description of targets
#
#              Makefile.config includes the following files:
# 
#              Make.config.settings   default configuration settings
#              Make.config.override   local user settings
#              Make.config.assemble   maps 'config' settings to 'flag' settings
#              Make.config.targets    configure targets
#              Make.mach.*            all machine-dependent settings
#              Make.config.objects    list of object files
#              DEPEND                 Make-generated dependencies
#
#              Make.mach.* should be the only file that one should
#              edit when porting grackle to a new machine.
#
#              Make.config.override should be the only file that
#              one should edit when defining local user settings.
#              Preferably, this is done implicitly through
#              the available make targets (e.g. "make precision-32").
#              Use "make help-config" for a list of all configuration
#              settings.  These make targets do error-checking; hand-editing 
#              Make.config.override is more error-prone.
#
#=======================================================================

# Use bash since sh on datastar does not recognize ">&" used in dep: target

SHELL    = /bin/bash

TOP_DIR  = ../..
DEFAULT_INSTALL_PREFIX = /usr/local
OUTPUT   = out.compile

GRACKLE_DIR = .
VERBOSE  = 0

SVN      = hg

#-----------------------------------------------------------------------
# Make.config.settings is used for setting default values to all compile-time 
# configuration settings.
#-----------------------------------------------------------------------

include $(GRACKLE_DIR)/Make.config.settings

#-----------------------------------------------------------------------
# Make.config.machine is used for setting which Make.mach.* file to use
#-----------------------------------------------------------------------

MAKE_CONFIG_MACHINE  = $(GRACKLE_DIR)/Make.config.machine
include $(GRACKLE_DIR)/Make.config.machine

#-----------------------------------------------------------------------
# Make.config.override is used for overriding the default settings in
# Make.config.settings.  This was made separate from the default settings 
# to enable easily interfacing grackle with a software testing environment 
# like lcatest.
#-----------------------------------------------------------------------

MAKE_CONFIG_OVERRIDE = $(GRACKLE_DIR)/Make.config.override

include $(MAKE_CONFIG_OVERRIDE)

#-----------------------------------------------------------------------
# Make.mach.<machine-name> defines all machine-dependent settings.
#
# mqk 03/22/2013: moved this to before Make.config.assemble, so the
# MACH_INSTALL_* variables are defined there.
#-----------------------------------------------------------------------

-include $(GRACKLE_DIR)/Make.mach.$(CONFIG_MACHINE)
-include $(HOME)/.grackle/Make.mach.$(CONFIG_MACHINE)

#-----------------------------------------------------------------------
# Make.config.assemble takes the settings in the Make.config.settings
# and Make.config.override, and generates the appropriate make variables
# required by this makefile.  E.g. $(CXX), $(CXXFLAGS), etc.
#-----------------------------------------------------------------------

include Make.config.assemble

#=======================================================================
# OBJECT FILES
#=======================================================================

include Make.config.objects

#-----------------------------------------------------------------------
# MAKE LIBGRACKLE BY DEFAULT
#-----------------------------------------------------------------------

# I believe this has to be done here, as it is target-specific.
lib: CFLAGS   += $(MACH_SHARED_FLAGS) 
lib: CXXFLAGS += $(MACH_SHARED_FLAGS) -DSHARED_LIBRARY
lib: FFLAGS   += $(MACH_SHARED_FLAGS)
lib: F90FLAGS += $(MACH_SHARED_FLAGS)
lib: LDFLAGS  += $(MACH_SHARED_FLAGS)
lib: libgrackle.$(MACH_SHARED_EXT)

libgrackle.$(MACH_SHARED_EXT): $(MODULES) autogen dep $(OBJS_LIB) MACHNOTES
	@rm -f $@
	@echo "Linking"
	-@$(LD) $(LDFLAGS) $(LDOUTPUT_FLAGS) $@ $(OBJS_LIB) $(LIBS) $(LIB_INSTALL_FLAGS) >& $(OUTPUT)
	@(if [ -e $@ ]; then \
	   echo "Success!"; \
	else \
	   echo "$(LD) $(LDFLAGS) $(LDOUTPUT_FLAGS) $*.$(MACH_SHARED_EXT) $*.o $(OBJS_LIB) $(LIBS) $(LIB_INSTALL_FLAGS)" >> temp1; \
	   cat temp1 $(OUTPUT) > temp2; \
	   rm -f temp1; \
	   mv -f temp2 $(OUTPUT); \
	   echo "Failed! See $(OUTPUT) for error messages"; \
	fi)

MACHNOTES: 
	@echo -e $(MACHINE_NOTES)

lib%.$(MACH_SHARED_EXT): $(MODULES) autogen dep $(OBJS_LIB)
	@rm -f $@
	@echo "Linking"
	-@$(LD) $(LDFLAGS) $(LDOUTPUT_FLAGS) $@ $(OBJS_LIB) $(LIBS) $(LIB_INSTALL_FLAGS) >& $(OUTPUT)
	@(if [ -e $@ ]; then \
	   echo "Success!"; \
	else \
	   echo "$(LD) $(LDFLAGS) $(LDOUTPUT_FLAGS) $*.$(MACH_SHARED_EXT) $*.o $(OBJS_LIB) $(LIBS) $(LIB_INSTALL_FLAGS)" >> temp1; \
	   cat temp1 $(OUTPUT) > temp2; \
	   rm -f temp1; \
	   mv -f temp2 $(OUTPUT); \
	   echo "Failed! See $(OUTPUT) for error messages"; \
	fi)

#-----------------------------------------------------------------------
# WRITE ALL COMPILER OUTPUT TO FILE
#-----------------------------------------------------------------------

.PHONY: verbose
verbose: VERBOSE = 1

#-----------------------------------------------------------------------
# Implicit rules
#-----------------------------------------------------------------------

.SUFFIXES: .C .F .o

# Inhibit removing any *.o files after compiling

.PRECIOUS: %.o %.C %.F

.F.o:
	@echo "Compiling $<"
	@rm -f $@
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  $(FC) -c -o $@ $(FFLAGS) $(DEFINES) $*.F >& $(OUTPUT) ; \
	  if [ ! -e $@ ]; then \
             echo; \
             echo "$(FC) -c -o $@ $(FFLAGS) $(DEFINES) $*.F"; \
             echo; \
             $(FC) -c -o $@ $(FFLAGS) $(DEFINES) $*.F; \
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(FC) -c -o $@ $(FFLAGS)  $(DEFINES) $*.f >> $(OUTPUT) 2>&1 ; \
	  if [ ! -e $@ ]; then \
	     echo "See $(OUTPUT) for error messages"; \
	     exit 1; \
	  fi ; \
	fi)

.C.o:
	@rm -f $@
	@echo "Compiling $<"
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  $(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C \
	    >& $(OUTPUT) ; \
	  if [ ! -e $@ ]; then \
             echo; \
             echo "$(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C"; \
             echo; \
             $(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C;\
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C \
	    >> $(OUTPUT) 2>&1 ; \
	  if [ ! -e $@ ]; then \
	     echo "See $(OUTPUT) for error messages"; \
	     exit 1; \
	  fi ; \
	fi)

#-----------------------------------------------------------------------
# Generate all make-generated source files
#-----------------------------------------------------------------------

.PHONY: autogen
autogen: auto_show_config.C auto_show_flags.C auto_show_version.C

# Force update of auto_show_config.C

.PHONY: auto_show_config.C
auto_show_config.C:
	-@$(MAKE) -s show-config  >& temp.show-config
	-@awk 'BEGIN {print "#include <stdio.h>\nvoid auto_show_config(FILE *fp) {"}; {print "   fprintf (fp,\""$$0"\\n\");"}; END {print "}"}' < temp.show-config > auto_show_config.C

# Force update of auto_show_flags.C

.PHONY: auto_show_flags.C
auto_show_flags.C:
	-@$(MAKE) -s show-flags  >& temp.show-flags
	-@awk 'BEGIN {print "#include <stdio.h>\nvoid auto_show_flags(FILE *fp) {"}; {print "   fprintf (fp,\""$$0"\\n\");"}; END {print "}"}' < temp.show-flags > auto_show_flags.C

# Force update of auto_show_version.C

.PHONY: auto_show_version.C
auto_show_version.C:
	-@$(MAKE) -s show-version  >& temp.show-version
	-@awk 'BEGIN {print "#include <stdio.h>\nvoid auto_show_version(FILE *fp) {"}; {print "   fprintf (fp,\""$$0"\\n\");"}; END {print "}"}' < temp.show-version > auto_show_version.C

#-----------------------------------------------------------------------
# Generate dependency file
#-----------------------------------------------------------------------

.PHONY: dep
dep:
	@echo "Updating DEPEND"
	-@(makedepend $(DEFINES) $(INCLUDES) -a -fDEPEND -o.o -m -- -- *.C) >> out.make.DEPEND 2>&1
	-@(makedepend $(DEFINES) $(INCLUDES) -a -fDEPEND -o.o -m -- -- *.F) >> out.make.DEPEND 2>&1
	-@(makedepend $(DEFINES) $(INCLUDES) -a -fDEPEND -o.o -m -- -- *.h) >> out.make.DEPEND 2>&1

include DEPEND

#-----------------------------------------------------------------------
# HELP TARGET
#-----------------------------------------------------------------------

help:
	@echo
	@echo "========================================================================"
	@echo "   Grackle Makefile Help"
	@echo "========================================================================"
	@echo
	@echo "   make                Compile and generate librackle"
	@echo "   make install        Copy the library somewhere"
	@echo "   make help           Display this help information"
	@echo "   make clean          Remove object files, executable, etc."
	@echo "   make dep            Create make dependencies in DEPEND file"
	@echo
	@echo "   make show-version   Display revision control system branch and revision"
	@echo "   make show-diff      Display local file modifications"
	@echo
	@echo "   make help-config    Display detailed help on configuration make targets"
	@echo "   make show-config    Display the configuration settings"
	@echo "   make show-flags     Display specific compilation flags"
	@echo "   make default        Reset the configuration to the default values"
	@echo

#-----------------------------------------------------------------------
# INSTALLATION TARGET
#-----------------------------------------------------------------------

install:
	@echo "Installing libgrackle.$(MACH_SHARED_EXT) to $(INSTALL_LIB_DIR)."
	@cp libgrackle.$(MACH_SHARED_EXT) $(INSTALL_LIB_DIR)/

	@echo "Installing grackle header files to $(INSTALL_INCLUDE_DIR)."
	@(if [ ! -d $(INSTALL_INCLUDE_DIR) ]; then \
		mkdir $(INSTALL_INCLUDE_DIR); \
	fi)
	@cp grackle.h grackle_macros.h grackle_types.h code_units.h chemistry_data.h $(INSTALL_INCLUDE_DIR)

#-----------------------------------------------------------------------

clean:
	-@rm -f *.$(MACH_SHARED_EXT) *.o DEPEND.bak *~ $(OUTPUT) *.exe auto_show*.C DEPEND
	-@touch DEPEND

#-----------------------------------------------------------------------
# Include configuration targets
#-----------------------------------------------------------------------

include $(GRACKLE_DIR)/Make.config.targets
# DO NOT DELETE
