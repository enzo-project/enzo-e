all:
	$(MAKE) -C doc
	$(MAKE) -C src
	$(MAKE) -C bin
clean:
	rm -f `find . -name index.html`
	$(MAKE) -C doc clean
	$(MAKE) -C src clean
	$(MAKE) -C bin clean


