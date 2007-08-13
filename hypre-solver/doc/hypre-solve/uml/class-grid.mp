input metauml;

beginfig(1);

   Class.grid("Grid")
	(
         "-id_: int",
         "-id_parent_: int",
         "-ip_: int",
         "-xl_: Scalar[3]",
         "-xu_: Scalar[3]",
         "-n_: int[3]",
         "-discret_: Discret *",
         "-neighbors_: vector<Grid*>",
         "-children_: vector<Grid*>",
         "-hypre_grid_: HYPRE_SStructGrid",
         "-hypre_graph_: HYPRE_SStructGraph",
         "-level_: int",
         "-u_: Scalar*")
	(
         "+print()",
         "+write()",
         "+read()",
         "+check()",
         "+id()",
         "+id_parent()",
         "+set_child()",
         "+child()",
         "+num_children()",
         "+set_neighbor()",
         "+neighbor()",
         "+num_neighbors()",
         "+assert_neighbors()",
         "+hypre_grid()",
         "+hypre_graph()",
         "+set_level()",
         "+level()",
         "+num_unknowns()",
         "+upper_vertex()",
         "+lower_vertex()",
         "+unknown()",
         "+discret()",
         "+ip()",
         "+is_adjacent()",
         "+is_on_proc()",
         "+find_neighbor_indices()");

	drawObjects(grid);

endfig;
end

   