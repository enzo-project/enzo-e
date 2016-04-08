.PHONY: dimension-%
dimension-%:
	ln -sf config/config-dimension-$*.incl config-dimension.incl
.PHONY: max-level-%
max-level-%:
	ln -sf config/config-max-level-$*.incl config-max-level.incl
.PHONY: stopping-cycle-%
stopping-cycle-%:
	ln -sf config/config-stopping-cycle-$*.incl config-stopping-cycle.incl
.PHONY: stopping-interval-%
stopping-interval-%:
	ln -sf config/config-stopping-interval-$*.incl config-stopping-interval.incl
.PHONY: adapt-interval-%
adapt-interval-%:
	ln -sf config/config-adapt-interval-$*.incl config-adapt-interval.incl

.PHONY: output-schedule-%
output-schedule-%:
	ln -sf config/config-output-schedule-$*.incl config-output-schedule.incl

.PHONY: processors-%
processors-%:
	ln -sf config/config2-processors-$*.incl config2-processors.incl
	ln -sf config/config3-processors-$*.incl config3-processors.incl
.PHONY: root-size-%
root-size-%:
	ln -sf config/config2-root-size-$*.incl config2-root-size.incl
	ln -sf config/config3-root-size-$*.incl config3-root-size.incl
.PHONY: trace-%
trace-%:
	ln -sf config/config-trace-$*.incl config-trace.incl

.PHONY: default
default:
	-@$(MAKE) dimension-3
	-@$(MAKE) processors-8
	-@$(MAKE) root-size-32
	-@$(MAKE) trace-off
	-@$(MAKE) max-level-3
	-@$(MAKE) stopping-cycle-5
	-@$(MAKE) stopping-interval-1
	-@$(MAKE) adapt-interval-1
	-@$(MAKE) output-schedule-0
.PHONY: show-config
show-config:
	-@grep FILE config-dimension.incl
	-@grep FILE config2-processors.incl
	-@grep FILE config3-processors.incl
	-@grep FILE config2-root-size.incl
	-@grep FILE config3-root-size.incl
	-@grep FILE config-max-level.incl
	-@grep FILE config-trace.incl
	-@grep FILE config-stopping-cycle.incl
	-@grep FILE config-stopping-interval.incl
	-@grep FILE config-adapt-interval.incl
	-@grep FILE config-output-schedule.incl

.PHONY: clean
clean:
	-@rm -f config-dimension.incl
	-@rm -f config2-processors.incl
	-@rm -f config3-processors.incl
	-@rm -f config2-root-size.incl
	-@rm -f config3-root-size.incl
	-@rm -f config-max-level.incl
	-@rm -f config-trace.incl
	-@rm -f config-stopping-cycle.incl
	-@rm -f config-stopping-interval.incl
	-@rm -f config-adapt-interval.incl
	-@rm -f config-output-schedule.incl


