.PHONY: dimension-%
dimension-%:
	ln -sf config-dimension-$*.incl config-dimension.incl
.PHONY: max-level-%
max-level-%:
	ln -sf config-max-level-$*.incl config-max-level.incl
.PHONY: stopping-cycle-%
stopping-cycle-%:
	ln -sf config-stopping-cycle-$*.incl config-stopping-cycle.incl
.PHONY: refresh-%
refresh-%:
	ln -sf config-refresh-$*.incl config-refresh.incl
.PHONY: output-%
output-%:
	ln -sf config-output-$*.incl config-output.incl

.PHONY: processors-%
processors-%:
	ln -sf config2-processors-$*.incl config2-processors.incl
	ln -sf config3-processors-$*.incl config3-processors.incl
.PHONY: root-size-%
root-size-%:
	ln -sf config2-root-size-$*.incl config2-root-size.incl
	ln -sf config3-root-size-$*.incl config3-root-size.incl

.PHONY: default
default:
	-@$(MAKE) dimension-3
	-@$(MAKE) processors-8
	-@$(MAKE) root-size-32
	-@$(MAKE) max-level-1
	-@$(MAKE) stopping-cycle-5
	-@$(MAKE) refresh-counter
	-@$(MAKE) output-1
.PHONY: show-config
show-config:
	-@grep FILE config-dimension.incl
	-@grep FILE config2-processors.incl
	-@grep FILE config3-processors.incl
	-@grep FILE config2-root-size.incl
	-@grep FILE config3-root-size.incl
	-@grep FILE config-max-level.incl
	-@grep FILE config-stopping-cycle.incl
	-@grep FILE config-refresh.incl
	-@grep FILE config-output.incl

.PHONY: clean
clean:
	-@rm -f config-dimension.incl
	-@rm -f config2-processors.incl
	-@rm -f config3-processors.incl
	-@rm -f config2-root-size.incl
	-@rm -f config3-root-size.incl
	-@rm -f config-max-level.incl
	-@rm -f config-stopping-cycle.incl
	-@rm -f config-output.incl


