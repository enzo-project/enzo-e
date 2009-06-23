.PHONY: all links clean

all: links


.PHONY: test
test:
	echo $(PWD)/index-test.html
	(cd ../testing; cvs update; gen-www $(PWD)/index-test.html)

links:
	ln -sf ../tar .
	ln -sf ../doc/*.ps .
	ln -sf ../src/html html

clean:
	rm -f tar
