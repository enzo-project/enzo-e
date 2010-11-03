DOCS = \
      devel-depend \
      devel-uml \
      doc-design \
      doc-project \
      doc-requirements \
      doc-testing \
      doc-vision \
      manual-devel \
      manual-reference \
      manual-user \
      paper-tree \
      performance \
      report-nsf \
      talk-2009-09 \
      talk-2010-06

DIRS = $(DOCS) uml

.PHONY: all
all: $(DOCS)

.PHONY: $(DOCS)

doc-vision:
	cd doc-vision; ln -sf ../Makefile.1 makefile
	$(MAKE) -C doc-vision pdf
performance:
	cd performance; ln -sf ../Makefile.1 makefile
	$(MAKE) -C performance
doc-project:
	cd doc-project; ln -sf ../Makefile.1 makefile
	$(MAKE) -C doc-project
report-nsf:
#       use own Makefile
	cd report-nsf; ln -sf ../papers.bib .
	$(MAKE) -C report-nsf
doc-requirements: 
	cd doc-requirements; ln -sf ../Makefile.1 makefile
	$(MAKE) -C doc-requirements
devel-depend: 
	$(MAKE) -C devel-depend
doc-design: 
	cd doc-design; ln -sf ../Makefile.1 makefile
	$(MAKE) -C doc-design pdf
devel-uml: 
	cd devel-uml; ln -sf ../Makefile.1 makefile
	$(MAKE) -C devel-uml pdf
doc-testing: 
	cd doc-testing; ln -sf ../Makefile.1 makefile
	$(MAKE) -C doc-testing
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

