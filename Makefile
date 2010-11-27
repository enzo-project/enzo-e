.PHONY: all src doc

all: src
src:
	scons -u

doc:
	$(MAKE) -C doc
.PHONY: dot dot-png dot-eps


.PHONY: clean

clean: clean-bin clean-src clean-include clean-lib

.PHONY: clean-src
clean-src:
	$(MAKE) -C src     clean

.PHONY: clean-bin
clean-bin:
	$(MAKE) -C bin     clean

.PHONY: clean-include
clean-include:
	$(MAKE) -C include clean

.PHONY: clean-lib
clean-lib:
	$(MAKE) -C lib     clean

.PHONY: clean-doc
clean-doc:
	$(MAKE) -C doc     clean

