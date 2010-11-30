.PHONY: all src doc clean 

all: src

src:
	$(MAKE) -C src
doc:
	$(MAKE) -C doc

clean: clean-bin clean-src clean-include clean-lib
	$(MAKE) -C src     clean
	$(MAKE) -C bin     clean
	$(MAKE) -C include clean
	$(MAKE) -C lib     clean
	$(MAKE) -C doc     clean

