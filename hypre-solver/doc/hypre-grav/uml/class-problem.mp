input metauml;

beginfig(1);

   Class.problem("Problem")
	(
         "-spheres_: vector<Sphere *>",
         "-points_: vector<Point *>",
         "-hierarchy_: Hierarchy")
        (
         "+print()",
         "+write()",
         "+read()",
         "+dimension()",
         "+hierarchy()",
         "+sphere()",
         "+num_spheres()",
         "+point()",
         "+num_points()",
         "+grid()",
         "+num_grids()",
         "+level()",
         "+num_levels()");


	drawObjects(problem);

endfig;
end

   