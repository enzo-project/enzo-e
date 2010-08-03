# tcsh used for >&
SHELL = /bin/tcsh

TRAC_PATH = /usr/local/trac-cello
WIKI_FILES = \
	NotesCodingStandards \
	NotesComponentControl \
	NotesComponentData \
	NotesComponentDisk \
	NotesComponentDistribute \
	NotesComponentError \
	NotesComponentField \
	NotesComponentMemory \
	NotesComponentMesh \
	NotesComponentMonitor \
	NotesComponentParallel \
	NotesComponentParameters \
	NotesComponentParticles \
	NotesComponentPerformance \
	NotesComponentPortal \
	NotesComponents \
	NotesComponentSchedule \
	NotesComponentSimulation \
	NotesComponentTask \
	NotesComponentTest \
	NotesComponentUser \
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
	NotesFieldsEnzo \
	NotesOpenIssues  \
	NotesTestProblemDoubleMach  \
	NotesTestProblemNoh  \
	NotesTestProblemSedovBlast  \
	WikiStart

# BOOK ORDER
#
#	WikiStart
#           NotesWishList
#           NotesParameters
#               NotesParametersEnzo
#               NotesParametersCello
#           NotesComponents
#               NotesComponentControl
#               NotesComponentData
#               NotesComponentDisk
#               NotesComponentDistribute
#               NotesComponentError
#               NotesComponentField
#               NotesComponentMemory
#               NotesComponentMesh
#               NotesComponentMethod
#               NotesComponentMonitor
#               NotesComponentParallel
#               NotesComponentParameters
#               NotesComponentParticles
#               NotesComponentPerformance
#               NotesComponentPortal
#               NotesComponentSchedule
#               NotesComponentSimulation
#               NotesComponentTask
#               NotesComponentTest
#               NotesComponentUser
#           NotesCodingStandards
#           NotesTestProblems
#               NotesTestProblemDoubleMach
#               NotesTestProblemImplosion
#               NotesTestProblemNoh
#               NotesTestProblemSedovBlast
#               NotesTestProblemShockPool
#               NotesTestProblemWavePool

# NOT IN BOOK

#       NotesOpenIssues
#       NotesDesignPatterns


# OBSOLETE
#
#       NotesFieldsEnzo
#	NotesComponentAmr
#	NotesComponentAnalysis
#	NotesComponentArrayBlock
#	NotesComponentArrayItBlock
#	NotesComponentArrayLayout
#	NotesComponentArray
#	NotesComponentMethod
#	NotesComponentParallelMpi
#	NotesComponentProblem
#	NotesComponentRegion
#	NotesEnzoTestProblems
# 	NotesParametersEnzoPoll

help:
	-@echo
	-@echo "help:    Display this help summary"
	-@echo "write:   Write all wiki pages 'write/*' to the Trac site"
	-@echo "read:    Read all wiki pages from the Trac site to 'read/*'"
	-@echo "list:    List all wiki pages 'Notes*' on the Trac site"
	-@echo

.PHONY: read write-%
write: $(patsubst %,write-%,$(WIKI_FILES))
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

