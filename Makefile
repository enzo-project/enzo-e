.PHONY: src clean 

src:
	$(MAKE) -C src

clean: clean-bin clean-src clean-include clean-lib
	$(MAKE) -C src     clean
	$(MAKE) -C bin     clean
	$(MAKE) -C lib     clean
	$(MAKE) -C include clean

