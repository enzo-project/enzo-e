SHELL = /bin/tcsh
TRAC_PATH = /usr/local/trac-cello
WIKI_FILES = \
	NotesCodingStandards \
	NotesComponentAmr \
	NotesComponentAnalysis \
	NotesComponentArray \
	NotesComponentControl \
	NotesComponentDisk \
	NotesComponentError \
	NotesComponentField \
	NotesComponentMemory \
	NotesComponentMethod \
	NotesComponentMonitor \
	NotesComponentParallel \
	NotesComponentParameters \
	NotesComponentParticles \
	NotesComponentPerformance \
	NotesComponentPortal \
	NotesComponentProblem \
	NotesComponentRegion \
	NotesComponents \
	NotesComponentSimulation \
	NotesEnzoTestProblems \
	NotesParameters \
	NotesParametersEnzo \
	NotesParametersMethod \
	NotesTestProblems \
	NotesTestProblemShockPool \
	NotesTestProblemWavePool \
	NotesWishList \
	WikiStart

# 	NotesParametersEnzoPoll

write: $(patsubst %,import-%,$(WIKI_FILES))

write-%:
	-@echo "importing import/$*"
	-@(trac-admin $(TRAC_PATH) wiki import $* write/$*)

read: $(patsubst %,read-%,$(WIKI_FILES))
read-%: 
	-@echo "exporting read/$*"
	-@(trac-admin $(TRAC_PATH) wiki export $* read/$*) >& /dev/null

trac:
	-@mkdir trac
diff: $(patsubst %,diff-%,$(WIKI_FILES))
diff-%:
	-@(gmake export-$*)
	-diff $* trac/$*


