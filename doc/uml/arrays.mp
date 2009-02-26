prologues := 2
input metauml;
beginfig(1)
  Interface.ABC("ArrayABC")();
  Class.A("Array")()();
  Interface.AD("ArrayDecorator")();
  Class.AB("ArrayBlocked")()();
  Class.AP("ArrayPadded")()();
  Class.AM("ArrayMpi")()();
  Class.AO("ArrayOmp")()();
   
  
  leftToRight.top(20)(A,AD);
  leftToRight(20)(AB,AP,AM,AO);

  ABC.s = .5[A.ne, AD.nw] + (0,30);
  AD.s  = .5[AP.ne, AM.nw] + (0,30);
  

  drawObjects(ABC,A,AD,AB,AP,AM,AO);

  stepY := 10;
  link(inheritance)(pathStepY(AD.n,ABC.s,stepY));
  link(inheritance)(pathStepY(A.n, ABC.s,stepY));

  link(inheritance)(pathStepY(AB.n,AD.s,stepY));
  link(inheritance)(pathStepY(AP.n,AD.s,stepY));
  link(inheritance)(pathStepY(AM.n,AD.s,stepY));
  link(inheritance)(pathStepY(AO.n,AD.s,stepY));

endfig;
end
