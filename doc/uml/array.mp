input metauml;
beginfig(1)
  Class.A("Array")
	(
         "N_ int",
	 "n_ int [3]",
         "a_ Scalar *"
	)
	(
	"Array()",
	"~Array()",
        "void copy()",
	"void resize()",
	"void size()",
	"int length()",
	"Scalar * values()",
	"Scalar & operator () ()",
	);
  drawObject(A);
endfig;
end
