FIGS = array arrays amr task
PDF  = $(patsubst %,%.pdf, $(FIGS))
PS   = $(patsubst %,%.ps, $(FIGS))
OBJS = $(patsubst %,%.1, $(FIGS))
SRC  = $(patsubst %,%.mp,$(FIGS))


all: $(PDF)

.ps.pdf:
	epstopdf $*.ps
.1.ps:
	mv $*.1 $*.ps
.mp.1:
	mpost $*.mp

.SUFFIXES: .pdf .ps .mp .1

.PHONY: clean
clean:
	-@rm -rf *.log *~ 
distclean: clean
	-@rm -rf *.1  *.ps *.pdf *.png


