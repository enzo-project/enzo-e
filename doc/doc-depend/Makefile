.PHONY: depend depend-eps depend-png
depend: depend-eps depend-png
depend-eps:
	dot -Teps depend-libraries.dot > depend-libraries.eps
	dot -Teps depend-includes.dot > depend-includes.eps
	$(MAKE) install-eps
depend-png:
	dot -Tpng depend-libraries.dot > depend-libraries.png
	dot -Tpng depend-includes.dot > depend-includes.png
	$(MAKE) install-png

#-----------------------------------------------------------------------
# INSTALL
#-----------------------------------------------------------------------

.PHONY: install-%

install-eps:
	if [ ! -e ../eps ]; then mkdir ../eps; fi
	cp *.eps ../eps
install-png:
	if [ ! -e ../png ]; then mkdir ../png; fi
	cp *.png ../png
#-----------------------------------------------------------------------
.PHONY: clean
clean:
	rm -rf *png *eps
#-----------------------------------------------------------------------
.PHONY: distclean
distclean: clean
#-----------------------------------------------------------------------
