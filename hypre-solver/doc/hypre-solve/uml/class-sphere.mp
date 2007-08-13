input metauml;

beginfig(1);

   Class.sphere("Sphere")
	("-c_: Scalar *",
         "-r_: Scalar",
         "-m_: Scalar",
         "-d_: static int")
        ("+set_dim()",
         "+check()",
         "+write()",
         "+print()",
         "+read()");

	drawObjects(sphere);

endfig;
end

   