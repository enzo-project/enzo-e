FIGS = array arrays amr task

PDF  = $(patsubst %,%.pdf, $(FIGS))
PS   = $(patsubst %,%.ps,  $(FIGS))
OBJS = $(patsubst %,%.1,   $(FIGS))
SRC  = $(patsubst %,%.mp,  $(FIGS))
JPG  = $(patsubst %,%.jpg, $(FIGS))

all: $(PDF) $(JPG)




.ppm.jpg:
	ppmtojpeg $*.ppm > $*.jpg
.pdf.ppm:
	pdftoppm $*.pdf > $*.ppm
.ps.pdf:
	epstopdf $*.ps
.1.ps:
	mv $*.1 $*.ps
.mp.1:
	mpost $*.mp

.SUFFIXES: .pdf .ps .mp .1 .ppm .jpg

.PHONY: clean
clean:
	-@rm -rf *.log *~ 
distclean: clean
	-@rm -rf *.1  *.ps *.pdf *.jpg


