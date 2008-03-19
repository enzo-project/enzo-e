.PHONY: all links clean


all: links include.tex
	make -C project
	make -C requirements
	make -C design
	make -C testing

include.tex:
	ln -s ~/latex/include-common.tex include.tex

links: project requirements design testing
	cd project;      rm -f Makefile;    ln -sf ../Makefile.1 Makefile
	cd requirements; rm -f Makefile;    ln -sf ../Makefile.1 Makefile
	cd design;       rm -f Makefile;    ln -sf ../Makefile.1 Makefile 
	cd testing;      rm -f Makefile;    ln -sf ../Makefile.1 Makefile
	cd project;      rm -f include.tex; ln -sf ../include.tex
	cd requirements; rm -f include.tex; ln -sf ../include.tex
	cd design;       rm -f include.tex; ln -sf ../include.tex
	cd testing;      rm -f include.tex; ln -sf ../include.tex
clean: links
	-@rm -f *~ *.ps
	make -C project clean
	make -C requirements clean
	make -C design clean
	make -C testing clean
	cd project;      rm -f Makefile;    
	cd requirements; rm -f Makefile;    
	cd design;       rm -f Makefile;    
	cd testing;      rm -f Makefile;    
	cd project;      rm -f include.tex; 
	cd requirements; rm -f include.tex; 
	cd design;       rm -f include.tex; 
	cd testing;      rm -f include.tex; 
	rm -f include.tex

