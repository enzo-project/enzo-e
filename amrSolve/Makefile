install:
	$(MAKE) -C docs install
	$(MAKE) -C include install
	$(MAKE) -C lib install
	$(MAKE) -C src install
clean:
	-rm -f *~
	$(MAKE) -C docs clean
	$(MAKE) -C include clean
	$(MAKE) -C lib clean
	$(MAKE) -C src clean

links:
	ln -sf docs/*.ps .
