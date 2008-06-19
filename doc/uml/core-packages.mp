input metauml;
beginfig(1)
  Package.P("Problem")();
  Package.M("Method")();
  Package.C("Control")();
  Package.D("Data")();
  Package.I("IO")();
topToBottom(40)(P,D);
topToBottom(40)(M,I);
topToBottom(20)(M,C);


leftToRight(40)(P, M);


  drawObject(P);
  drawObject(M);
  drawObject(C);
  drawObject(D);
  drawObject(I);
endfig;
end
