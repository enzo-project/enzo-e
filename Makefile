.PHONY: charm-%
charm-%:
	python scons.py -j$* bin/charm/enzo-p type=charm

.PHONY: mpi-%
mpi-%:
	python scons.py -j$* bin/mpi/enzo-p type=mpi

.PHONY: clean
clean:
	scons -c

