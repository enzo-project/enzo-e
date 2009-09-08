prologues := 2
input metauml;
beginfig(1)
  Class.A("Amr")()();
  Class.B("Level")()();
  Class.Ca("Patch")()();
  Class.Cb("Tree")()();
  Class.Da("Layout")()();
  Class.Db("Box")()();
  Class.Dc("Node")()();
  B.n = A.s + (0,-20);
  Ca.n = B.s + (0,-20);
  Db.n = Ca.s + (0, -30);
  Cb.n = Ca.n + (50, 0);
  Da.n = Db.n + (-50,0);
  Dc.n = Db.n + (50,0);
  drawObjects(A,B,Ca,Cb,Da,Db,Dc);
  link(dashedLink)(Db.e -- Dc.w);
  link(aggregation)(B.n -- A.s);
  link(aggregation)(Ca.n -- B.s);
  link(aggregation)
	(pathStepY(Da.n, Ca.s, 10));
  link(aggregation)
	(pathStepY(Db.n, Ca.s, 10));
  link(aggregation)
	(pathStepY(Dc.n, Ca.s , 10));
  link(composition)(Dc.n -- Cb.s );
	
endfig;
end
