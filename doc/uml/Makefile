FIGS = array arrays amr core-components
PS   = $(patsubst %,%.ps, $(FIGS))
OBJS = $(patsubst %,%.1, $(FIGS))
SRC  = $(patsubst %,%.mp,$(FIGS))


all: $(PS)

.1.ps:
	mv $*.1 $*.ps
.mp.1:
	mpost $*.mp

.SUFFIXES: .ps .mp .1

.PHONY: clean
clean:
	-@rm -rf *.log *~
distclean: clean
	-@rm -rf *.1  *.ps


