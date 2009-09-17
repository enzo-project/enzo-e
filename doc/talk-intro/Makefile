cello-talk.pdf: cello-talk.tex
	pdflatex cello-talk.tex
	pdflatex cello-talk.tex
cello-talk.tex: cello-talk.txt
	perl outlinebeamer.pl < cello-talk.txt > cello-talk.tex

include makefile
