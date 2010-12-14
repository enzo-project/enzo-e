prologues := 2
input metauml;
beginfig(1)
  Class.M("Mesh")()();
  Class.P("Patch")()();
  Class.T("Tree")()();
  Class.N("Node")()();
  Class.L("Layout")()();
  Class.B("Block")()();
  Class.F("Field")()();
  Class.R("Particles")()();

  T.n = M.n + (65, 0);

  P.n = M.s + (0, -30);
  N.n = T.s + (0, -30);

  L.n = P.n + (-60, 0);

  B.n = P.s + (0, -30);

  F.n = B.n + (-60, 15);
  R.n = B.n + (-60, -15);
  
  drawObjects(M,P,T,N,L,B,F,R);
  

  link(composition)(P.n -- M.s);
  link(composition)(B.n -- P.s);
  link(composition)(N.n -- T.s);
  link(aggregation)(T.w -- M.e);
  link(aggregation)(N.w -- P.e);
  link(association)(L.e -- P.w);

  link(association)(pathStepX(B.w, F.e, -10));
  link(association)(pathStepX(B.w, R.e, -10));

  item(iAssoc)("1")(obj.ne = M.s + (-5,0));
  item(iAssoc)("1")(obj.sw = M.e + (0,3));

  item(iAssoc)("1")(obj.ne = T.s + (-5,0));
  item(iAssoc)("1")(obj.se = T.w);

  item(iAssoc)("1..k")(obj.sw = L.e);

  item(iAssoc)("1")(obj.se = P.w);
  item(iAssoc)("1..*")(obj.se = P.n );
  item(iAssoc)("1")(obj.ne = P.s + (-5,0));
  item(iAssoc)("1")(obj.sw = P.e + (0,3));

  item(iAssoc)("1..*")(obj.se = N.n );
  item(iAssoc)("1")(obj.se = N.w);

  item(iAssoc)("1..*")(obj.se = B.n);

  
endfig;
end
