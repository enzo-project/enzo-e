prologues := 2
input metauml;
beginfig(1)
  Class.M("Mesh")()();
  Class.P("Patch")()();
  Class.T("Tree")()();
  Class.N("Node")()();
  Class.L("Layout")()();
  Class.B("Block")()();

  T.n = M.n + (65, 0);

  P.n = M.s + (0, -30);
  N.n = T.s + (0, -30);

  L.e = P.e + (-60, 0);

  B.n = P.s + (0, -30);
  
  drawObjects(M,P,T,N,L,B);
  
  link(composition)(P.n -- M.s);
  link(composition)(B.n -- P.s);
  link(composition)(N.n -- T.s);
  link(aggregation)(T.w -- M.e);
  link(aggregation)(N.w -- P.e);
  link(association)(L.e -- P.w);

  item(iAssoc)("1")(obj.ne = M.s + (-5,0));
  item(iAssoc)("1")(obj.sw = M.e + (0,3));

  item(iAssoc)("1")(obj.ne = T.s + (-5,0));
  item(iAssoc)("1")(obj.se = T.w);

  item(iAssoc)("1")(obj.sw = L.e);

  item(iAssoc)("1")(obj.se = P.w);
  item(iAssoc)("1..*")(obj.se = P.n );
  item(iAssoc)("1")(obj.ne = P.s + (-5,0));
  item(iAssoc)("1")(obj.sw = P.e + (0,3));

  item(iAssoc)("1..*")(obj.se = N.n );
  item(iAssoc)("1")(obj.se = N.w);

  item(iAssoc)("1..*")(obj.se = B.n);

	
  
endfig;
end
