input metauml;

beginfig(1);

   Class.hierarchy("Hierarchy")
	("-dimension_: int",
	 "-levels_: vector<Level *>",
	 "-grids_: vector<Grid *>")
	("+insert_grid()",
	 "+init_levels()",
         "+set_dim()",
         "+print()",
         "+write()",
         "+check()",
         "+dimension()",
         "+level()",
         "+num_levels()",
         "+grid()",
         "+num_grids()",
         "+parent()",
         "+set_parent()");

	drawObjects(hierarchy);

endfig;
end

   