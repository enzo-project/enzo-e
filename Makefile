.PHONY: default
default: bin/enzo-p

bin/enzo-p:
	./build.sh bin/enzo-p

.PHONY: doc
doc:
	$(MAKE) -C src doc

.PHONY: diff
diff:
	./tools/org-diff.sh > diff.org
.PHONY: clean
clean:
	python scons.py -c


