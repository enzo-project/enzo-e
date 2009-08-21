SHELL = /bin/tcsh
TRAC_PATH = /usr/local/trac-cello
WIKI_FILES = \
	NotesCodingStandards \
	NotesComponentAmr \
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
	NotesTestProblems \
	NotesTestProblemShockPool \
	NotesTestProblemWavePool \
	NotesWishList \
	WikiStart

# 	NotesParametersEnzoPoll

import: $(patsubst %,import-%,$(WIKI_FILES))

import-%: export-%
	-@(trac-admin $(TRAC_PATH) wiki import $* $*) >& /dev/null

export: trac $(patsubst %,export-%,$(WIKI_FILES))
export-%: 
	-@rm -f trac/$*
	-@(trac-admin $(TRAC_PATH) wiki export $* trac/$*) >& /dev/null

trac:
	-@mkdir trac
diff: $(patsubst %,diff-%,$(WIKI_FILES))
diff-%:
	-@(gmake export-$*)
	-diff $* trac/$*


