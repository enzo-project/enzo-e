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

  A.info.iName.iFont.scale := 1.5;
  B.info.iName.iFont.scale := 1.5;
  Ca.info.iName.iFont.scale := 1.5;
  Cb.info.iName.iFont.scale := 1.5;
  Da.info.iName.iFont.scale := 1.5;
  Db.info.iName.iFont.scale := 1.5;
  Dc.info.iName.iFont.scale := 1.5;

  B.n = A.s + (0,-30);
  Ca.n = B.s + (0,-30);
  Db.n = Ca.s + (0, -40);
  Cb.n = Ca.n + (60, 0);
  Da.n = Db.n + (-60,0);
  Dc.n = Db.n + (60,0);
  drawObjects(A,B,Ca,Cb,Da,Db,Dc);
  link(dashedLink)(Db.e -- Dc.w);
  link(aggregation)(B.n -- A.s);
  link(aggregation)(Ca.n -- B.s);
  link(aggregation)
	(pathStepY(Da.n, Ca.s, 15));
  link(aggregation)
	(pathStepY(Db.n, Ca.s, 15));
  link(aggregation)
	(pathStepY(Dc.n, Ca.s , 15));
  link(composition)(Dc.n -- Cb.s );
	
endfig;
end
