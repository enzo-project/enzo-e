prologues := 2
input metauml;
beginfig(1)
  Class.A("Method")()();
  Class.B("Field")()();
  Class.C("Task")()();
  Class.D("Control")()();
  Class.E("Parallel")()();

  A.info.iName.iFont.scale := 1.5;
  B.info.iName.iFont.scale := 1.5;
  C.info.iName.iFont.scale := 1.5;
  D.info.iName.iFont.scale := 1.5;
  E.info.iName.iFont.scale := 1.5;

  B.s = A.s + (80,0);
  C.s = A.s + (40,-60);
  D.s = C.s + (-40, -60);
  E.s = D.s + (80, 0);
  
  drawObjects(A,B,C,D,E);
  
  link(association)(pathStepY(C.n, A.s, 15));
  link(association)(pathStepY(C.n, B.s, 15));
  link(association)(pathStepY(D.n, C.s, 15));
  link(association)(pathStepY(E.n, C.s, 15));
	
endfig;
end
