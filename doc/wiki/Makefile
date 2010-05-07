# tcsh used for >&
SHELL = /bin/tcsh

TRAC_PATH = /usr/local/trac-cello
WIKI_FILES = \
	NotesComponentData \
	NotesComponentSchedule \
	NotesComponentTest \
	NotesComponentUser \

	NotesCodingStandards \
	NotesComponentAmr \
	NotesComponentAnalysis \
	NotesComponentArray \
	NotesComponentArrayBlock \
	NotesComponentArrayItBlock \
	NotesComponentArrayLayout \
	NotesComponentControl \
	NotesComponentDisk \
	NotesComponentDistribute \
	NotesComponentError \
	NotesComponentField \
	NotesComponentMemory \
	NotesComponentMethod \
	NotesComponentMonitor \
	NotesComponentParallel \
	NotesComponentParallelMpi \
	NotesComponentTask \
	NotesComponentParameters \
	NotesComponentParticles \
	NotesComponentPerformance \
	NotesComponentPortal \
	NotesComponentProblem \
	NotesComponentRegion \
	NotesComponents \
	NotesComponentSimulation \
	NotesDesignPatterns \
	NotesEnzoTestProblems \
	NotesParameters \
	NotesParametersCello \
	NotesParametersEnzo \
	NotesTestProblemImplosion \
	NotesTestProblems \
	NotesTestProblemShockPool \
	NotesTestProblemWavePool \
	NotesWishList \
	WikiStart

# 	NotesParametersEnzoPoll
help:
	-@echo
	-@echo "help:    Display this help summary"
	-@echo "write:   Write all wiki pages 'write/*' to the Trac site"
	-@echo "read:    Read all wiki pages from the Trac site to 'read/*'"
	-@echo "list:    List all wiki pages 'Notes*' on the Trac site"
	-@echo

.PHONY: read write-%
write: $(patsubst %,import-%,$(WIKI_FILES))
write-%:
	-@echo "importing import/$*"
	-@(trac-admin $(TRAC_PATH) wiki import $* write/$*)

.PHONY: read read-%
read: $(patsubst %,read-%,$(WIKI_FILES))
read-%: 
	-@echo "exporting read/$*"
	-@(if ( ! -e read ) mkdir read)
	-@rm -f read/$*
	-@(trac-admin $(TRAC_PATH) wiki export $* read/$*) >& /dev/null

.PHONY: list
list:
	-@(trac-admin $(TRAC_PATH) wiki list | awk '/Notes/ {print $$1}')

