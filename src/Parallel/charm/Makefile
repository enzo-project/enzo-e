.PHONY: run

run: test
	~/charm/charm/bin/charmrun +p4 ./test

test: test.cpp test.decl.h test.def.h
	~/charm/charm/bin/charmc -language charm++ test.cpp -o test
test.decl.h test.def.h: test.ci
	~/charm/charm/bin/charmc -language charm++ test.ci

.PHONY: clean
clean:
	rm -f test.decl.h test.def.h test *~


