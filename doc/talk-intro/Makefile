PNG = dots.png dots-4-0.png dots-4-1.png dots-4-2.png dots-4-3.png dots-4-4.png dots-4-5.png dots-16-0.png dots-16-1.png dots-16-2.png dots-16-3.png dots-16-4.png dots-16-5.png cosmo2.png cosmo2-4-0.png cosmo2-4-1.png cosmo2-4-2.png cosmo2-4-3.png cosmo2-4-4.png cosmo2-4-5.png cosmo2-16-0.png cosmo2-16-1.png cosmo2-16-2.png cosmo2-16-3.png cosmo2-16-4.png cosmo2-16-5.png norman-10.png

cello-talk.pdf: cello-talk.tex $(PNG)
	pdflatex cello-talk.tex
cello-talk.tex: cello-talk.txt  $(PNG)
	perl outlinebeamer.pl < cello-talk.txt > cello-talk.tex

include makefile
