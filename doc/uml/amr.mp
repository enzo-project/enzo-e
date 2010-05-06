prologues := 2
input metauml;
beginfig(1)
  Class.A("Amr")()();
  Class.P("Patch")()();
  Class.T("Tree")()();
  Class.N("Node")()();

  A.info.iName.iFont.scale := 1.5;
  P.info.iName.iFont.scale := 1.5;
  T.info.iName.iFont.scale := 1.5;
  N.info.iName.iFont.scale := 1.5;

  P.n = A.n + (80, 0);
  T.n = A.s + (0, -40);
  N.n = T.n + (80, 0);
  
  drawObjects(A,P,T,N);
  
  item(iAssoc)("1")(obj.sw = A.e + (10,0));
  item(iAssoc)("1")(obj.sw = T.e + (10,0));
  item(iAssoc)("1..*")(obj.se = P.w );
  item(iAssoc)("1..*")(obj.se = N.w );
  item(iAssoc)("1")(obj.nw = A.s);
  item(iAssoc)("1")(obj.nw = P.s);
  item(iAssoc)("1")(obj.sw = T.n);
  item(iAssoc)("1")(obj.sw = N.n);


 
  link(composition)(P.w -- A.e);
  link(association)(A.s -- T.n);
  link(association)(P.s -- N.n);
  link(composition)(N.w -- T.e);
	
  
endfig;
end
