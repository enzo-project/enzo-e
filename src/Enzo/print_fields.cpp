#include <stdio.h>
#include "cello_hydro.h"

void print_fields(int mx, int my)
{
  double b,bi,bw;

  for (int k=0; k<NumberOfBaryonFields; k++) {
    b = 0;
    bi = 0;
    bw = 0;
    for (int iy=0; iy<mx; iy++) {
      for (int ix=0; ix<my; ix++) {
	int i=ix + mx*iy;
	b += BaryonField[k][i];
	if (3<=ix && ix <mx+3 && 3<=iy && iy <my+3) bi += BaryonField[k][i];
	bw += double(i+1)*BaryonField[k][i];
      }
    }
    printf ("%d [%d %22.16g %22.16g %22.16g]\n",CycleNumber,k,b,bi,bw);
  }
}
