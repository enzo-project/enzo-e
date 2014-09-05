.PHONY: default
default: enzop-8

.PHONY: enzop-%
enzop-%:
	python scons.py -j$* bin/enzo-p

.PHONY: doc
doc:
	$(MAKE) -C src doc

.PHONY: diff
diff:
	hg diff > out.diff
	./util/parse_diff.sh out.diff > diff.org
.PHONY: clean
clean:
	python scons.py -c


