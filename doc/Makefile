DOCS = \
      design \
      devel-uml \
      manual-devel \
      manual-reference \
      manual-user \
      notes \
      paper-tree \
      performance \
      project \
      requirements \
      talk-intro \
      testing \
      vision

DIRS = $(DOCS) uml
 

.PHONY: all
all: $(DOCS)

.PHONY: $(DOCS)

notes:
	cd notes; ln -sf ../Makefile.1 makefile
	$(MAKE) -C notes
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
devel-uml: 
	cd devel-uml; ln -sf ../Makefile.1 makefile
	$(MAKE) -C devel-uml
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
paper-tree: 
	cd paper-tree; ln -sf ../Makefile.1 makefile
	$(MAKE) -C paper-tree
talk-intro: 
	cd talk-intro; ln -sf ../Makefile.1 makefile
	$(MAKE) -C talk-intro

#-----------------------------------------------------------------------

.PHONY: clean
clean: $(patsubst %,%-clean,$(DIRS))
	-rm -f *~

%-clean:
	$(MAKE) -C $* clean

#-----------------------------------------------------------------------

.PHONY: distclean
distclean: $(patsubst %,%-distclean,$(DIRS))
	-@rm -f *~
%-distclean:
	$(MAKE) -C $* distclean

