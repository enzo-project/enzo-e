all: hypre-solve.pdf  hypre-solve.ps

hypre-solve:  hypre-solve.ps 

FIG = fig/discret1.eps fig/discret2.eps fig/discret3.eps \
      fig/discret4.eps fig/discret5.eps fig/discret6.eps \
      fig/stencil0.eps fig/stencil1.eps fig/stencil2.eps \
      fig/discret-step-1.eps fig/stencil.eps


########################################################################
hypre-solve.ps:hypre-solve.dvi $(FIG)
	dvips -o hypre-solve.ps hypre-solve.dvi
hypre-solve.pdf:hypre-solve.dvi $(FIG)
	dvipdf hypre-solve.dvi
hypre-solve.dvi:hypre-solve.tex include.tex UML
	-latex hypre-solve
	-latex hypre-solve
UML:
#		make -C uml
hypre-solve2.ps:hypre-solve.ps
	a2ps -M Letter -2 hypre-solve.ps -o hypre-solve2.ps
pdf.hypre-solve:hypre-solve.dvi
	dvipdf hypre-solve.dvi hypre-solve.pdf
hypre-solve2.pdf:hypre-solve2.ps
	ps2pdf hypre-solve2.ps hypre-solve2.pdf

include.tex:
	ln -sf ~/latex/include-common.tex include.tex
########################################################################

clean:
	-@rm -f *~ *.dvi *.ind *.idx *.ilg *.aux *.bbl *.blg *.log *.toc *.ps *.pdf *.glo `find . -name "*" -type l`
	$(MAKE) -C uml clean




