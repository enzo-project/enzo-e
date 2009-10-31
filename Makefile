.PHONY: all
all:
	make -C src
	make -C doc
	make -C src doc

.PHONY: clean

clean: clean-doc clean-src
	$(MAKE) -C doc     clean

.PHONY: clean-src

clean-src:
	$(MAKE) -C src     clean
	$(MAKE) -C include clean
	$(MAKE) -C lib     clean

.PHONY: clean-doc

clean-doc:
	$(MAKE) -C doc     clean

