ALL = project requirements design testing developer user algorithms

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
developer: 
	cd developer; ln -sf ../Makefile.1 makefile
	$(MAKE) -C developer
user: 
	cd user; ln -sf ../Makefile.1 makefile
	$(MAKE) -C user

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
