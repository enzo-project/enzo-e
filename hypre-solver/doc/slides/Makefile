SHELL = /bin/sh

current:pdf
PAPER = slides

PS:dvi
	dvips -o $(PAPER).ps $(PAPER).dvi
pdf:dvi
	dvipdf $(PAPER).dvi $(PAPER).pdf
bib: mgnet.bib papers.bib
	bibtex $(PAPER)

dvi: 
	-latex $(PAPER).tex
mgnet.bib:
	-ln -s ~/references/mgnet.bib .
papers.bib:
	-ln -s ~/references/papers.bib .
clean:
	-rm -f *~ *.dvi *.ind *.idx *.ilg *.aux *.bbl *.blg *.log *.toc *.ps *.pdf *.out
tar:
	make clean; \
	set dir = $(PAPER); \
	cd ..; \
	tar cf $$dir.tar $$dir; \
	gzip $$dir.tar; \
	mv $$dir.tar.gz ~/Outgoing


