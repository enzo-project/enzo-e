.PHONY: run

CXXFLAGS = -g

CHARMRUN = ~/charm/charm/bin/charmrun +p4 
CHARMC   = ~/charm/charm/bin/charmc -language charm++  $(CXXFLAGS)


run: debug
	$(CHARMRUN) ./debug

debug: debug.o enzo_EnzoBlock.o mesh_Block.o
	$(CHARMC) debug.o enzo_EnzoBlock.o mesh_Block.o -o debug

debug.o: debug.decl.h debug.def.h
	$(CHARMC) -c debug.cpp 
enzo_EnzoBlock.o: debug.decl.h debug.def.h
	$(CHARMC) -c enzo_EnzoBlock.cpp
mesh_Block.o: debug.decl.h debug.def.h
	$(CHARMC) -c mesh_Block.cpp

debug.decl.h debug.def.h: debug.ci
	~/charm/charm/bin/charmc -language charm++ debug.ci

.PHONY: clean
clean:
	rm -f debug.decl.h debug.def.h debug *~ *.o charmrun


