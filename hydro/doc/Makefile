all: hydro-scale.pdf  hydro-scale.ps

hydro-scale:  hydro-scale.ps 

FIG = 


########################################################################
hydro-scale.ps:hydro-scale.dvi $(FIG)
	dvips -o hydro-scale.ps hydro-scale.dvi
hydro-scale.pdf:hydro-scale.dvi $(FIG)
	dvipdf hydro-scale.dvi
hydro-scale.dvi:hydro-scale.tex include.tex
	-latex hydro-scale
	-latex hydro-scale
hydro-scale2.ps:hydro-scale.ps
	a2ps -M Letter -2 hydro-scale.ps -o hydro-scale2.ps
pdf.hydro-scale:hydro-scale.dvi
	dvipdf hydro-scale.dvi hydro-scale.pdf
hydro-scale2.pdf:hydro-scale2.ps
	ps2pdf hydro-scale2.ps hydro-scale2.pdf

include.tex:
	ln -sf ~/latex/include-common.tex include.tex
########################################################################

clean:
	-@rm -f *~ *.dvi *.ind *.idx *.ilg *.aux *.bbl *.blg *.log *.toc *.ps *.pdf *.glo `find . -name "*" -type l`




