.PHONY: charm
charm:
	python scons.py bin/charm/enzo-p type=charm

.PHONY: mpi
mpi:
	python scons.py bin/mpi/enzo-p type=mpi

.PHONY: clean
clean:
	scons -c

