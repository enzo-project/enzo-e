input metauml;

beginfig(1);

   Class.level("Level")
	("-n_: int",
	 "-grids_: vector<Grid *>")
	("+insert_grid()",
	 "+print()",
         "+write()",
         "+check()",
         "+grid()",
         "+num_grids()");

	drawObjects(level);

endfig;
end

   