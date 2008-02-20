input metauml;
beginfig(1)
  Class.A("Array")
	(
         "N_ int",
	 "n_ int [4]",
         "a_ Scalar *"
	)
	(
	"Array()",
	"Array(n0,n1,n2,n3)",
	"~Array()",
        "void copy(Array&)",
	"void resize(n0,n1,n2,n3)",
	"void size(n0,n1,n2,n3)",
	"int length()",
	"Scalar * values()",
	"Scalar & operator () (i0,i1,i2,i3)",
	);
  drawObject(A);
endfig;
end
