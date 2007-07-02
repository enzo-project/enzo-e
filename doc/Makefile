ALL = project requirements design testing algorithms

.PHONY: all
all: $(ALL)

.PHONY: $(ALL)

project:
	cd project; ln -sf ../Makefile.1 makefile
	$(MAKE) -C project
requirements: 
	cd requirements; ln -sf ../Makefile.1 makefile
	$(MAKE) -C requirements
design: 
	cd design; ln -sf ../Makefile.1 makefile
	$(MAKE) -C design
testing: 
	cd testing; ln -sf ../Makefile.1 makefile
	$(MAKE) -C testing
algorithms: 
	cd algorithms; ln -sf ../Makefile.1 makefile
	$(MAKE) -C algorithms

#-----------------------------------------------------------------------

.PHONY: clean
clean: $(patsubst %,%-clean,$(ALL))
	-rm -f *~

%-clean:
	$(MAKE) -C $* clean 
