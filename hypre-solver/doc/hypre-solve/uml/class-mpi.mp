input metauml;

beginfig(1);

   Class.mpi("Mpi")
	(
          "-comm_: MPI_Comm",
          "-np_: int",
          "-ip_: int "
        )
	(
          "+ip()",
          "+np()"
        );

	drawObjects(mpi);

endfig;
end

   