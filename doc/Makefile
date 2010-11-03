DOCS = \
      depend \
      design \
      devel-uml \
      manual-devel \
      manual-reference \
      manual-user \
      paper-tree \
      performance \
      project \
      proposal-ascr \
      proposal-nsf \
      proposal-nsf-phy \
      report-nsf \
      requirements \
      talk-2009-09 \
      talk-2010-06 \
      testing \
      vision

DIRS = $(DOCS) uml

.PHONY: all
all: $(DOCS)

.PHONY: $(DOCS)

vision:
	cd vision; ln -sf ../Makefile.1 makefile
	$(MAKE) -C vision pdf
performance:
	cd performance; ln -sf ../Makefile.1 makefile
	$(MAKE) -C performance
project:
	cd project; ln -sf ../Makefile.1 makefile
	$(MAKE) -C project
proposal-nsf:
#       use own Makefile
	cd proposal-nsf; ln -sf ../papers.bib .
	$(MAKE) -C proposal-nsf
proposal-nsf-phy:
#       use own Makefile
	cd proposal-nsf-phy; ln -sf ../papers.bib .
	$(MAKE) -C proposal-nsf-phy
proposal-ascr:
	cd proposal-ascr; ln -sf ../Makefile.1 makefile
	cd proposal-ascr; ln -sf ../papers.bib .
	$(MAKE) -C proposal-ascr
report-nsf:
#       use own Makefile
	cd report-nsf; ln -sf ../papers.bib .
	$(MAKE) -C report-nsf
requirements: 
	cd requirements; ln -sf ../Makefile.1 makefile
	$(MAKE) -C requirements
depend: 
	$(MAKE) -C depend
design: 
	cd design; ln -sf ../Makefile.1 makefile
	$(MAKE) -C design pdf
devel-uml: 
	cd devel-uml; ln -sf ../Makefile.1 makefile
	$(MAKE) -C devel-uml pdf
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
talk-2009-09: 
	cd talk-2009-09; ln -sf ../Makefile.1 makefile
	$(MAKE) -C talk-2009-09 pdf2
talk-2010-06: 
	cd talk-2010-06; ln -sf ../Makefile.1 makefile
	$(MAKE) -C talk-2010-06 pdf2

#-----------------------------------------------------------------------

.PHONY: clean
clean: $(patsubst %,%-clean,$(DIRS))
	-rm `find . -type l`
	-rm -f *~

%-clean:
	$(MAKE) -C $* clean

#-----------------------------------------------------------------------

.PHONY: distclean
distclean: $(patsubst %,%-distclean,$(DIRS))
	-@rm -f *~
%-distclean:
	$(MAKE) -C $* distclean

