#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <map>

/**
 * @file      debug-matrix.cpp
 * @brief     Test program for debugging hypre matrices
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
 */

#define INDEX(i0,i1,i2,grid,N2)             ((i0)+(N2)*((i1)+(N2)*((i2)+(N2)*(grid))))
#define AINDEX(i0,i1,i2,k,grid,N2) ((k) + 7*((i0)+(N2)*((i1)+(N2)*((i2)+(N2)*(grid)))))

      

int main(int argc, char **argv)
{
  if (argc != 8) {
    fprintf (stderr,"Usage: %s N L offset procs dp0 dp1 dp2\n",argv[0]);
    exit(1);
  }

  // Usage: N levels type
  //   N = 4 (e.g.)
  //   levels = 1,2,3, etc.
  //   type = xyz   i0 i1 i2  j0 j1 j2  value
  //          ij    i j                 value
  //          none

  int iarg = 1;
  int N          = atoi(argv[iarg++]);
  int num_levels = atoi(argv[iarg++]);
  int is_offset  = atoi(argv[iarg++]);
  int np         = atoi(argv[iarg++]);
  int dp0        = atoi(argv[iarg++]);
  int dp1        = atoi(argv[iarg++]);
  int dp2        = atoi(argv[iarg++]);

  int N2 = N + 2;
  int N23 = N2*N2*N2;
  int Asize = N23*2*7;
  double * A = new double [Asize]; // 

  printf ("N = %d  L = %d  O = %d  P = %d\n",
	  N,num_levels,is_offset,np);

  FILE *fp;
  int ip;
  int i0,i1,i2;
  int j0,j1,j2;
  int i,j,k;
  double value;

  for (int level=0; level<num_levels; level++) {

    char filename[40];
    sprintf (filename,"A.%02d.00.00.00000",level);
    fp = fopen (filename,"r");

    int linenum = 0;
    while (fscanf (fp,"%d %d %d %d %d %lf",&ip,&i0,&i1,&i2,&k,&value) != EOF) {
      linenum++;

      i0 ++;
      i1 ++;
      i2 ++;

      i0 += dp0;
      i1 += dp1;
      i2 += dp2;

      j0 = i0;
      j1 = i1;
      j2 = i2;
      
      if (k==1) ++j0;
      if (k==2) --j0;
      if (k==3) ++j1;
      if (k==4) --j1;
      if (k==5) ++j2;
      if (k==6) --j2;

      i = INDEX(i0,i1,i2,level,N2);
      j = INDEX(j0,j1,j2,level,N2);

      if (value != 0) {
	if (! (0 <= i && i < Asize)) {
	  printf ("linenum=%d i=%d Asize=%d i0,i1,i2=%d,%d,%d k=%d level=%d N2=%d\n",
		  linenum,i,Asize,i0,i1,i2,k,level,N2);
	}
	if (! (0 <= j && j < Asize)) {
	  printf ("linenum=%d i=%d Asize=%d i0,i1,i2=%d,%d,%d k=%d level=%d N2=%d\n",
		  linenum,i,Asize,i0,i1,i2,k,level,N2);
	}
	if (! (AINDEX(i0,i1,i2,k,level,N2) >= 0)) {
	  printf ("linenum=%d i=%d Asize=%d i0,i1,i2=%d,%d,%d k=%d level=%d N2=%d\n",
		  linenum,i,Asize,i0,i1,i2,k,level,N2);
	}
	if (! (AINDEX(i0,i1,i2,k,level,N2) < Asize)) {
	  printf ("linenum=%d i=%d Asize=%d i0,i1,i2=%d,%d,%d k=%d level=%d N2=%d\n",
		  linenum, i,Asize,i0,i1,i2,k,level,N2);
	}
	assert (0 <= i && i < Asize);
	assert (0 <= j && j < Asize);
	assert (AINDEX(i0,i1,i2,k,level,N2) >= 0);
	assert (AINDEX(i0,i1,i2,k,level,N2) < Asize);
	A[AINDEX(i0,i1,i2,k,0,N2)] = value;
	printf ("%d  %d %d %g\n",level,i+1,j+1,value);
      }

    }
    fclose(fp);
  }

  typedef std::pair<const int,int> IJ_type;
  std::map<IJ_type,double> Agraph;

  fp = fopen ("A.UMatrix.00000","r");
  while (fscanf (fp,"%d %d %lf",&i,&j,&value) != EOF) {

    int ii=i;
    int jj=j;

    i0 = ii % N2;
    ii /= N2;
    i1 = ii % N2;
    ii /= N2;
    i2 = ii % N2;

    j0 = jj % N2;
    jj /= N2;
    j1 = jj % N2;
    jj /= N2;
    j2 = jj % N2;
      
    if (value != 0) {
      assert (0 <= i && i < Asize);
      assert (0 <= j && j < Asize);
      IJ_type p(i,j);
      Agraph[p] = value;
      printf ("12  %d %d %g\n",i+1,j+1,value);
    }
  }
  fclose(fp);

  printf ("Testing Stencil symmetry\n");

  // Check problem x-y-z symmetry for stencil values

  for (int level=0; level<num_levels; level++) {
    for (i0=0; i0<N2; i0++) {
      for (i1=0; i1<N2; i1++) {
	for (i2=0; i2<N2; i2++) {
	  for (int k=0; k<7; k++) {

	    int a1 = AINDEX(i0,i1,i2,k,level,N2);
	    int a2;

	    // Check diagonal
	    if (k==0) {
	      for (int k0=0; k0<2; k0++) {
		for (int k1=0; k1<2; k1++) {
		  for (int k2=0; k2<2; k2++) {
		    j0=(1-k0)*i0 + k0*(N2-i0-1);
		    j1=(1-k1)*i1 + k1*(N2-i1-1);
		    j2=(1-k2)*i2 + k2*(N2-i2-1);
		    a2 = AINDEX(j0,j1,j2,k,level,N2);
		    if (A[a1] != A[a2] && 
			(!is_offset || k0+k1+k2==0)) {
		      fprintf (stderr,"OOPS! level %d diagonal "
			       "(%d %d %d; %d) = %g  (%d %d %d; %d) = %g\n",
			       level,
			       i0,i1,i2,k,A[a1],
			       j0,j1,j2,k,A[a2]);
		    }
		  }
		}
	      }
	    }
	    // Check x+
	    if (k==1 || k==2) {
	      j0=N2-i0-1;
	      j1=i1;
	      j2=i2;
	      a2 = AINDEX(j0,j1,j2,3-k,level,N2);
	      if (A[a1] != A[a2] && !is_offset) {
		fprintf (stderr,"OOPS! level %d X "
			 "(%d %d %d; %d) = %g  (%d %d %d; %d) = %g\n",
			 level,
			 i0,i1,i2,k,A[a1],
			 j0,j1,j2,3-k,A[a2]);
	      }
	    }
	    if (k==3 || k==4) {
	      j0=i0;
	      j1=N2-i1-1;
	      j2=i2;
	      a2 = AINDEX(j0,j1,j2,7-k,level,N2);
	      if (A[a1] != A[a2] && !is_offset) {
		fprintf (stderr,"OOPS! level %d Y "
			 "(%d %d %d; %d) = %g  (%d %d %d; %d) = %g\n",
			 level,
			 i0,i1,i2,k,A[a1],
			 j0,j1,j2,7-k,A[a2]);
	      }
	    }
	    if (k==5 || k==6) {
	      j0=i0;
	      j1=i1;
	      j2=N2-i2-1;
	      a2 = AINDEX(j0,j1,j2,11-k,level,N2);
	      if (A[a1] != A[a2] && !is_offset) {
		fprintf (stderr,"OOPS! level %d Z "
			 "(%d %d %d; %d) = %g  (%d %d %d; %d) = %g\n",
			 level,
			 i0,i1,i2,k,A[a1],
			 j0,j1,j2,11-k,A[a2]);
	      }
	    }
	  }
	}
      }
    }
  }
  printf ("Testing Graph symmetry\n");

  
  std::map<IJ_type,double>::iterator nz;
  for (nz = Agraph.begin(); nz != Agraph.end(); nz++) {

    i = (nz->first).first;
    j = (nz->first).second;
    int ii=i;
    int jj=j;

    int ig = (ii > N23) ? 1 : 0;
    ii -= ig*N23;
    i0 = ii % N2;
    ii /= N2;
    i1 = ii % N2;
    ii /= N2;
    i2 = ii % N2;

    int jg = (jj > N23) ? 1 : 0;
    jj -= jg*N23;
    j0 = jj % N2;
    jj /= N2;
    j1 = jj % N2;
    jj /= N2;
    j2 = jj % N2;

    for (int k0=0; k0<2; k0++) {
      for (int k1=0; k1<2; k1++) {
	for (int k2=0; k2<2; k2++) {

	  int ib0=(1-k0)*i0 + k0*(N2-i0-1);
	  int ib1=(1-k1)*i1 + k1*(N2-i1-1);
	  int ib2=(1-k2)*i2 + k2*(N2-i2-1);

	  int jb0=(1-k0)*j0 + k0*(N2-j0-1);
	  int jb1=(1-k1)*j1 + k1*(N2-j1-1);
	  int jb2=(1-k2)*j2 + k2*(N2-j2-1);

	  int ib = ib0 + N2*(ib1 + N2*(ib2 + N2*ig));
	  int jb = jb0 + N2*(jb1 + N2*(jb2 + N2*jg));

	  if (nz->second != Agraph[IJ_type(ib,jb)]) {
	      printf ("(%d %d) (%d;%d %d %d) (%d;%d %d %d)  %g != "
		      "(%d %d) (%d;%d %d %d) (%d;%d %d %d)  %g\n",
		      i,j,
		      ig,i0,i1,i2,
		      jg,j0,j1,j2,
		      nz->second,
		      ib,jb,
		      ig,ib0,ib1,ib2,
		      jg,jb0,jb1,jb2,
		      Agraph[IJ_type(ib,jb)]);
	  }
	}
      }
    }
  }
  exit(0);
}
