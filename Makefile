.PHONY: default
default: charm-1

.PHONY: charm-%
charm-%:
	python scons.py -j$* bin/charm/enzo-p type=charm

.PHONY: mpi-%
mpi-%:
	python scons.py -j$* bin/mpi/enzo-p type=mpi

.PHONY: serial-%
serial-%:
	python scons.py -j$* bin/serial/enzo-p type=serial


.PHONY: clean
clean:
	python scons.py -c

