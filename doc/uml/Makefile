FIGS = array amr core-components
OBJS = $(patsubst %,%.1, $(FIGS))
SRC  = $(patsubst %,%.mp,$(FIGS))


all: $(OBJS)
.mp.1:
	mpost $*.mp

.SUFFIXES: .mp .1

.PHONY: clean
clean:
	rm -rf *.1 *.log *~

