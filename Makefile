.PHONY: all
all:
	$(MAKE) -C src install-include
	$(MAKE) -C src dep
	$(MAKE) -C src
	$(MAKE) -C doc
	$(MAKE) -C src doc

.PHONY: clean

clean: clean-src clean-include clean-lib

.PHONY: clean-src
clean-src:
	$(MAKE) -C src     clean

.PHONY: clean-include
clean-include:
	$(MAKE) -C include clean

.PHONY: clean-lib
clean-lib:
	$(MAKE) -C lib     clean

.PHONY: clean-doc
clean-doc:
	$(MAKE) -C doc     clean

