prologues := 2
input metauml;
beginfig(1)
  Class.M("Mesh")()();
  Class.P("Patch")()();
  Class.T("Tree")()();
  Class.N("Node")()();

  P.n = M.n + (80, 0);
  T.n = M.s + (0, -40);
  N.n = T.n + (80, 0);
  
  drawObjects(M,P,T,N);
  
  item(iAssoc)("1")(obj.sw = M.e + (10,0));
  item(iAssoc)("1")(obj.sw = T.e + (10,0));
  item(iAssoc)("1..*")(obj.se = P.w );
  item(iAssoc)("1..*")(obj.se = N.w );
  item(iAssoc)("1")(obj.nw = M.s);
  item(iAssoc)("1")(obj.nw = P.s);
  item(iAssoc)("1")(obj.sw = T.n);
  item(iAssoc)("1")(obj.sw = N.n);


 
  link(composition)(P.w -- M.e);
  link(association)(M.s -- T.n);
  link(association)(P.s -- N.n);
  link(composition)(N.w -- T.e);
	
  
endfig;
end
