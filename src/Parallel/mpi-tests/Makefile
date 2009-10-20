all: mpi-distance mpi-array
mpi-distance: mpi-distance.cpp
	mpiCC -o mpi-distance mpi-distance.cpp -I ../../../include -L../../../lib -lperformance
mpi-array: mpi-array.cpp
	mpiCC -o mpi-array mpi-array.cpp -I ../../../include -L../../../lib -lperformance -lparameters

clean:
	-@rm -f *.o *~ mpi-distance mpi-array
