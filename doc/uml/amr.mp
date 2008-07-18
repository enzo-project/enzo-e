input metauml;
beginfig(1)
  Class.A("Hierarchy")()();
  Class.B("Level")()();
  Class.C("Patch")()();
  Class.Ca("Box")()();
  Class.Cb("Array")()();
  B.n = A.s + (0,-20);
  C.n = B.s + (0,-20);
  Ca.n = C.s + (-30,-30);
  Cb.n = C.s + (30,-30);
  drawObjects(A,B,C,Ca,Cb);
  link(aggregation)(B.n -- A.s);
  link(aggregation)(C.n -- B.s);
  link(composition)
	(pathStepY(Ca.n, C.s, 10));
  link(composition)
	(pathStepY(Cb.n, C.s, 10));
endfig;
end
