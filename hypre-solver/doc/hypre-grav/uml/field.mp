input metauml;

beginfig(1);

   ClassName.gridlocal("GridLocal");
   ClassName.field("Field");
   ClassName.array("Array");

   topToBottom(30)(gridlocal,field,array);

   drawObjects(gridlocal,field,array);

   link(composition)(field.n -- gridlocal.s);
   link(composition)(array.n -- field.s);

   item(iAssoc)("1   ")(obj.ne = gridlocal.s);
   item(iAssoc)("*   ")(obj.se = field.n);
   item(iAssoc)("1   ")(obj.ne = field.s);
   item(iAssoc)("*   ")(obj.se = array.n);

endfig;
end

   