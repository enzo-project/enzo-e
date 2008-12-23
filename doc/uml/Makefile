FIGS = array arrays amr core-components
OBJS = $(patsubst %,%.1, $(FIGS))
SRC  = $(patsubst %,%.mp,$(FIGS))


all: $(OBJS)
.mp.1:
	mpost $*.mp

.SUFFIXES: .mp .1

.PHONY: clean
clean:
	-@rm -rf *.log *~
distclean: clean
	-@rm -rf *.1 


