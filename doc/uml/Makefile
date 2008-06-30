FIGS = array amr core-packages
OBJS = $(patsubst %,%.1, $(FIGS))
SRC  = $(patsubst %,%.mp,$(FIGS))


all: $(OBJS)
.mp.1:
	mpost $*.mp

.SUFFIXES: .mp .1
