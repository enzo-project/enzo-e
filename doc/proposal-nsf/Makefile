PAPER = proposal-nsf

.PHONY: all
all: figures bio-bordner.pdf $(PAPER).pdf summary.pdf bibliography

.PHONY: clean
clean:
	rm -f $(PAPER).pdf summary.pdf bio-bordner.pdf

.PHONY: bibliography
bibliography: $(PAPER).pdf
	bibtex $(PAPER)
	pdflatex $(PAPER)
	pdftk $(PAPER).pdf cat 1-15   output $(PAPER)2.pdf
	pdftk $(PAPER).pdf cat 16-end output bibliography.pdf
	mv $(PAPER)2.pdf $(PAPER).pdf

.PHONY: figures
figures:
	-@ln -s ../figures/* .

%.pdf: %.tex newcommand.tex
	pdflatex $*
