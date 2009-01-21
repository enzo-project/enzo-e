input metauml;

beginfig(1);

   Class.grid("Grid")()();
   Class.gridlocal("GridLocal")()();
   Class.gridremote("GridRemote")()();

   centered_align_top(gridlocal, gridremote)(30, below(grid.s, 30));

   drawObjects(grid,gridlocal,gridremote);

   link(inheritance)(pathStepY(gridlocal.n, grid.s, 10));
   link(inheritance)(pathStepY(gridremote.n, grid.s, 10));

endfig;
end

   