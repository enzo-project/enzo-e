input metauml;

beginfig(1);

   Class.pnt("Point")
	(
          "-d_: static int",
          "-x_: Scalar *",
          "-m_: Scalar"
        )
	(
         "+print()",
         "+write()",
         "+read()",
         "+set_dim()"
        );

	drawObjects(pnt);

endfig;
end

   