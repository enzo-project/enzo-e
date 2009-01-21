input metauml;

beginfig(1);

   Class.hierarchy("Hierarchy")()();
   Class.level("Level")()();
   Class.grid("Grid")()();

   topToBottom(30)(hierarchy,level,grid);

   drawObjects(hierarchy,level,grid);

   link(composition)(level.n -- hierarchy.s);
   link(composition)(grid.n -- level.s);

   item(iAssoc)("1   ")(obj.ne = hierarchy.s);
   item(iAssoc)("*   ")(obj.se = level.n);
   item(iAssoc)("1   ")(obj.ne = level.s);
   item(iAssoc)("*   ")(obj.se = grid.n);

endfig;
end

   