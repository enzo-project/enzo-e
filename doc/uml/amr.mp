prologues := 2
input metauml;
beginfig(1)
  Class.M("Mesh")()();
  Class.P("Patch")()();
  Class.T("Tree")()();
  Class.N("Node")()();
  Class.L("Layout")()();
  Class.D("DataBlock")()();

  P.e = M.e + (0, -60);
  T.n = M.n + (60, 0);
  N.e = T.e + (0, -60);
  L.n = P.s + (-30, -40);
  D.n = P.s + (30, -40);
  
  drawObjects(M,P,T,N,L,D);
  
  link(composition)(P.n -- M.s);
  link(association)(M.e -- T.w);
  link(association)(P.e -- N.w);
  link(composition)(N.n -- T.s);
  link(aggregation)(pathStepY(L.n, P.s, 15));
  link(aggregation)(pathStepY(D.n, P.s, 15));

  item(iAssoc)("1")(obj.ne = M.s + (-5,0));
  item(iAssoc)("1")(obj.ne = T.s + (-5,0));
  item(iAssoc)("1..*")(obj.se = P.n );
  item(iAssoc)("1..*")(obj.se = N.n );
  item(iAssoc)("1")(obj.sw = M.e);
  item(iAssoc)("1")(obj.sw = P.e);
  item(iAssoc)("1")(obj.se = T.w);
  item(iAssoc)("1")(obj.se = N.w);
  item(iAssoc)("1")(obj.se = L.n);
  item(iAssoc)("1..*")(obj.se = D.n);

  

	
  
endfig;
end
