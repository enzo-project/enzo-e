input metauml;

beginfig(1);

   Class.hypre("Hypre")
    (
      "-stencil_: HYPRE_SStructStencil", 
      "+mpi: Mpi"
    )
    (
     "+init_hierarchy()",
     "+init_stencil()",
     "+init_graph()",
     "+init_matrix()",
     "+init_rhs()",
     "+init_solver()",
     "+solve()"
    );


	drawObjects(hypre);

endfig;
end

   