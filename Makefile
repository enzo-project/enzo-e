#----------------------------------------------------------------------
.PHONY: bin/enzo-p
bin/enzo-p:
	./build.sh bin/enzo-p
#----------------------------------------------------------------------
.PHONY: help
help:
	@echo 
	@echo "make          Compile Enzo-P as ./bin/enzo-p"
	@echo "make test     Run regression tests"
	@echo "make doc      Generate doxygen documentation from source"
	@echo "make clean    Remove object and test files"
	@echo "make diff     Generate org-mode 'diff.org' file from 'hg diff' output"
#----------------------------------------------------------------------
.PHONY: doc
doc:
	$(MAKE) -C src doc
#----------------------------------------------------------------------
.PHONY: clean
clean:
	./build.sh clean
#----------------------------------------------------------------------
.PHONY: diff
diff:
	./tools/org-diff.sh > diff.org
#----------------------------------------------------------------------
.PHONY: log
log:
	hg log | ./tools/org-log.sh > log.org


