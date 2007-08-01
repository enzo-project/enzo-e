input metauml;

beginfig(1);

   Class.discret("Discret")
   (
	"-neighbor_cell_: CELL_TYPE[3][2] *",
        "-n1_: int[3]",
        "-n2_: int[3]",
        "-n_: int[3]"
   )(
        "+neighbor_cell()",
        "+print()"
   );

	drawObjects(discret);

endfig;
end

   