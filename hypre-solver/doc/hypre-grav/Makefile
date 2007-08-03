all: hypre-grav.pdf

hypre-grav:  hypre-grav.ps 

FIG = fig/discret1.eps fig/discret2.eps fig/discret3.eps \
      fig/discret4.eps fig/discret5.eps fig/discret6.eps \
      fig/stencil0.eps fig/stencil1.eps fig/stencil2.eps \
      fig/discret-step-1.eps fig/stencil.eps


########################################################################
hypre-grav.ps:hypre-grav.dvi $(FIG)
	dvips -o hypre-grav.ps hypre-grav.dvi
hypre-grav.pdf:hypre-grav.dvi $(FIG)
	dvipdf hypre-grav.dvi
hypre-grav.dvi:hypre-grav.tex include.tex UML
	-latex hypre-grav
	-latex hypre-grav
UML:
	make -C uml
hypre-grav2.ps:hypre-grav.ps
	a2ps -M Letter -2 hypre-grav.ps -o hypre-grav2.ps
pdf.hypre-grav:hypre-grav.dvi
	dvipdf hypre-grav.dvi hypre-grav.pdf
hypre-grav2.pdf:hypre-grav2.ps
	ps2pdf hypre-grav2.ps hypre-grav2.pdf

include.tex:
	ln -sf ~/latex/include-common.tex include.tex
########################################################################

clean:
	-@rm -f *~ *.dvi *.ind *.idx *.ilg *.aux *.bbl *.blg *.log *.toc *.ps *.pdf *.glo `find . -name "*" -type l`
	$(MAKE) -C uml clean




