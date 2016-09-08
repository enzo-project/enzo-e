#----------------------------------------------------------------------
.PHONY: bin/enzo-p
bin/enzo-p:
	./build.sh bin/enzo-p
#----------------------------------------------------------------------
.PHONY: help
help:
	@echo 
	@echo "make [compile]  Compile Enzo-P as ./bin/enzo-p"
	@echo "make cccc       Compute code quality metrics in src/.cccc/cccc.html"
	@echo "make clean      Remove object and test files"
	@echo "make coverity   Compile enzo-p using the coverity static analysis tool"
	@echo "make diff       Generate org-mode 'diff.org' file from 'hg diff' output"
	@echo "make doc        Generate doxygen documentation from source in src-html"
	@echo "make log        Generate org-mode 'log.org' file from 'hg log' output"
	@echo "make reset      Clear any settings from an incomplete ./build.sh"
	@echo "make test       Run regression tests"
#----------------------------------------------------------------------
.PHONY: doc
doc:
	$(MAKE) -C src doc
#----------------------------------------------------------------------
.PHONY: clean
clean:
	./build.sh clean
#----------------------------------------------------------------------
.PHONY: reset
reset:
	./build.sh reset
#----------------------------------------------------------------------
.PHONY: test
test:
	./build.sh test
#----------------------------------------------------------------------
.PHONY: compile
compile:
	./build.sh compile
#----------------------------------------------------------------------
.PHONY: coverity
coverity:
	make clean
	rm -rf cov-int cov-int.tgz
	cov-build --dir cov-int $(MAKE)
	tar cfz cov-int.tgz cov-int
#----------------------------------------------------------------------
.PHONY: cccc
cccc:
	$(MAKE) -C src cccc
#----------------------------------------------------------------------
.PHONY: diff
diff:
	./tools/diff-org.sh > diff.org
#----------------------------------------------------------------------
.PHONY: log
log:
	hg log | ./tools/log-org.sh > log.org


