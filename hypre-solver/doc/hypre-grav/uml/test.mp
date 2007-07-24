input metauml;
beginfig(1);
  Class.A("A")()();
  Class.B("B")()();
  B.e = A.w + (-20, 0);
  drawObjects(A, B);
  link(inheritance)(B.e -- A.w);
endfig;
end

