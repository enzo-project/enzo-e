DOCS = \
      devel-design \
      devel-project \
      devel-requirements \
      devel-testing \
      devel-vision \
      doc-depend \
      doc-performance \
      doc-uml \
      manual-devel \
      manual-reference \
      manual-user \
      paper-tree \
      report-nsf \
      talk-2009-09 \
      talk-2010-06 \
      talk-2010-12 \
      uml

DIRS = $(DOCS) uml

.PHONY: all
all: $(DOCS)

.PHONY: $(DOCS)

devel-vision:
	cd devel-vision; ln -sf ../Makefile.1 makefile
	$(MAKE) -C devel-vision
doc-performance:
	cd doc-performance; ln -sf ../Makefile.1 makefile
	$(MAKE) -C doc-performance
devel-project:
	cd devel-project; ln -sf ../Makefile.1 makefile
	$(MAKE) -C devel-project
report-nsf:
	cd report-nsf;  ln -sf ../Makefile.1 makefile
	cd report-nsf;  ln -sf ../papers.bib .
	$(MAKE) -C report-nsf
devel-requirements: 
	cd devel-requirements; ln -sf ../Makefile.1 makefile
	$(MAKE) -C devel-requirements
doc-depend: 
	$(MAKE) -C doc-depend
devel-design: 
	cd devel-design; ln -sf ../Makefile.1 makefile
	$(MAKE) -C devel-design
doc-uml: 
	cd doc-uml; ln -sf ../Makefile.1 makefile
	$(MAKE) -C doc-uml
devel-testing: 
	cd devel-testing; ln -sf ../Makefile.1 makefile
	$(MAKE) -C devel-testing
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
	$(MAKE) -C talk-2009-09
talk-2010-06: 
	cd talk-2010-06; ln -sf ../Makefile.1 makefile
	$(MAKE) -C talk-2010-06
talk-2010-12: 
	cd talk-2010-12; ln -sf ../Makefile.1 makefile
	$(MAKE) -C talk-2010-12

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
	-@rm `find . -type l`
	rm -rf html pdf eps png
%-distclean:
	$(MAKE) -C $* distclean

