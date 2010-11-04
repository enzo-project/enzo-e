PAPER = report-nsf

.PHONY: all
all: $(PAPER).pdf 

.PHONY: bibliography
bibliography: $(PAPER).pdf
	bibtex $(PAPER)
	pdflatex $(PAPER)
#	pdftk $(PAPER).pdf cat 1-15   output $(PAPER)2.pdf
#	pdftk $(PAPER).pdf cat 16-end output bibliography.pdf
#A	mv $(PAPER)2.pdf $(PAPER).pdf

.PHONY: figures
figures:
	-@ln -s ../figures/* .

%.pdf: figures %.tex newcommand.tex
	pdflatex $*

.PHONY: clean
clean:
	-@rm -f *~ *.aux \
        $(PAPER).dvi \
        $(PAPER).log \
        $(PAPER).ps \
        $(PAPER).pdf \
        $(PAPER).html \
        $(PAPER).snm \
        $(PAPER).nav \
        $(PAPER).out \
        $(PAPER).blg \
        $(PAPER).bbl \
        $(PAPER).toc

.PHONY: distclean
distclean: clean

