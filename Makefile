.PHONY: default
default: enzop-8

.PHONY: enzop-%
enzop-%:
	python scons.py -j$* bin/enzo-p


.PHONY: clean
clean:
	python scons.py -c


