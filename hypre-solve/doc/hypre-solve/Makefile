all: hypre-solve.pdf  hypre-solve.ps

PAPER = hypre-solve
hypre-solve:  $(PAPER).ps 

FIG = fig/discret1.eps fig/discret2.eps fig/discret3.eps \
      fig/discret4.eps fig/discret5.eps fig/discret6.eps \
      fig/stencil0.eps fig/stencil1.eps fig/stencil2.eps \
      fig/discret-step-1.eps fig/stencil.eps


########################################################################
$(PAPER).ps:$(PAPER).dvi $(FIG)
	dvips -o $(PAPER).ps $(PAPER).dvi
$(PAPER).pdf:$(PAPER).dvi $(FIG)
	dvipdf $(PAPER).dvi
$(PAPER).dvi:$(PAPER).tex include.tex UML
	-latex $(PAPER)
	-latex $(PAPER)
UML:
#		make -C uml
html: $(PAPER).tex
	tth -e2 $(PAPER).tex > $(PAPER).html
$(PAPER)2.ps:$(PAPER).ps
	a2ps -M Letter -2 $(PAPER).ps -o $(PAPER)2.ps
pdf.$(PAPER):$(PAPER).dvi
	dvipdf $(PAPER).dvi $(PAPER).pdf
$(PAPER)2.pdf:$(PAPER)2.ps
	ps2pdf $(PAPER)2.ps $(PAPER)2.pdf

########################################################################

clean:
	-@rm -f *~ *.dvi *.ind *.idx *.ilg *.aux *.bbl *.blg *.log *.toc *.ps *.pdf *.glo `find . -name "*" -type l`
	$(MAKE) -C uml clean




