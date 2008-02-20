DIR = array

all:   $(patsubst %,%-all,$(DIR))
clean: $(patsubst %,%-clean,$(DIR))
test:  $(patsubst %,%-test,$(DIR))

%-all:
	$(MAKE) -C $*
%-clean:
	-rm -f *~
	$(MAKE) -C $* clean
%-test:
	$(MAKE) -C $* test
doc:
	doxygen
