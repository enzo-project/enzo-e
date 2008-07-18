ALL = vision \
      performance \
      project \
      requirements \
      design \
      testing \
      manual-devel \
      manual-user \
      manual-reference

.PHONY: all
all: $(ALL)

.PHONY: $(ALL)

vision:
	cd vision; ln -sf ../Makefile.1 makefile
	$(MAKE) -C vision
performance:
	cd performance; ln -sf ../Makefile.1 makefile
	$(MAKE) -C performance
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
manual-devel: 
	cd manual-devel; ln -sf ../Makefile.1 makefile
	$(MAKE) -C manual-devel
manual-user: 
	cd manual-user; ln -sf ../Makefile.1 makefile
	$(MAKE) -C manual-user
manual-reference: 
	cd manual-reference; ln -sf ../Makefile.1 makefile
	$(MAKE) -C manual-reference

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
