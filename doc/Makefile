ALL = project requirements design testing guide-devel guide-user algorithms

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
guide-devel: 
	cd guide-devel; ln -sf ../Makefile.1 makefile
	$(MAKE) -C guide-devel
guide-user: 
	cd guide-user; ln -sf ../Makefile.1 makefile
	$(MAKE) -C guide-user

#-----------------------------------------------------------------------

.PHONY: clean
clean: $(patsubst %,%-clean,$(ALL))
	-rm -f *~

%-clean:
	$(MAKE) -C $* clean

#-----------------------------------------------------------------------

.PHONY: distclean
distclean: $(patsubst %,%-distclean,$(ALL))
	-rm -f *~
%-distclean:
	$(MAKE) -C $* clean
	rm $*/include.tex
	rm $*/makefile
