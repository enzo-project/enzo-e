AWK_ORG= ./tools/awk/diff-org.awk
#----------------------------------------------------------------------
.PHONY: bin/enzo-e
bin/enzo-e:
	(cd build; ../compile.sh)
default:
	(cd build-default; ./compile.sh)
random:
	(cd build-random; ./compile.sh)

#----------------------------------------------------------------------
.PHONY: help
help:
	@echo 
	@echo "make [compile]  Compile Enzo-E as ./bin/enzo-e"
	@echo "make cccc       Compute code quality metrics in src/.cccc/cccc.html"
	@echo "make clean      Remove object and test files"
	@echo "make coverity   Compile enzo-e using the coverity static analysis tool"
	@echo "make diff       Generate org-mode 'diff.org' file from 'git diff' output"
	@echo "make gdb        Generate org-mode 'gdb.org' from gdb 'where' output in gdb.out"
	@echo "make doc        Generate html documentation from Sphinx in doc/build/html"
	@echo "make dox        Generate doxygen html documentation from source in doc/dox-html"
#	@echo "make log        Generate org-mode 'log.org' file from 'git log' output"
	@echo "make reset      Clear any settings from an incomplete ./build.sh"
	@echo "make test       Run regression tests"
#----------------------------------------------------------------------
.PHONY: show-conflicts
show-conflicts:
	git diff --name-only --diff-filter=U | awk 'BEGIN{print "* [/] Conflicts";} {printf "** TODO [[file:%s][%s]]\n",$$1,$$1;}' > conflicts.org
#----------------------------------------------------------------------
.PHONY: doc
doc:
	$(MAKE) -C doc html
#----------------------------------------------------------------------
.PHONY: dox
dox:
	$(MAKE) -C src dox-html
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
	cd build; ../compile.sh; cd ..
#----------------------------------------------------------------------
.PHONY: errors
errors:
	./tools/awk/error-org.awk build/out.make >errors.org
#----------------------------------------------------------------------
.PHONY: errors-default
errors-default:
	./tools/awk/error-org.awk build-default/out.make >errors.org
#----------------------------------------------------------------------
.PHONY: errors-random
errors-random:
	./tools/awk/error-org.awk build-random/out.make >errors.org
#----------------------------------------------------------------------
.PHONY: warnings
warnings:
	./tools/awk/warning-org.awk build/out.make >warnings.org
#----------------------------------------------------------------------
.PHONY: coverity
coverity:
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
	echo '* TODO [/] git diff branch (' `git rev-parse --abbrev-ref HEAD` ')' > diff.org
	echo '** TODO [/] Unstaged diff'                     >> diff.org
	git diff -b |$(AWK_ORG)                >> diff.org
	echo '** TODO [/] Staged diff'                       >> diff.org
	git diff --cached HEAD -b | $(AWK_ORG) >> diff.org
#----------------------------------------------------------------------
.PHONY: gdb
gdb:
	awk -f ./tools/awk/gdb-org.awk <gdb.out > gdb.org
#----------------------------------------------------------------------
#.PHONY: log
#log:
#	git log | ./tools/log-org.sh > log.org


