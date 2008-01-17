
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <map>

#define INDEX(i0,i1,i2,grid,N2)             ((i0)+(N2)*((i1)+(N2)*((i2)+(N2)*(grid))))
#define AINDEX(i0,i1,i2,k,grid,N2) ((k) + 7*((i0)+(N2)*((i1)+(N2)*((i2)+(N2)*(grid)))))

      

main(int argc, char **argv)
{
  if (argc != 4) exit(1);

  // Usage: N mask type
  //   N = 4 (e.g.)
  //   mask = 1  coarse
  //          2  fine
  //          3  fine + coarse
  //   type = xyz   i0 i1 i2  j0 j1 j2  value
  //          ij    i j                 value
  //          none

  int N = atoi(argv[1]);
  int N2 = N + 2;
  int N23 = N2*N2*N2;
  int mask = atoi(argv[2]);
  int Asize = N23*2*7;
  double * A = new double [Asize]; // 
  char * type = argv[3];
  bool is_xyz   = strcmp(type,"xyz")==0;
  bool is_ij    = strcmp(type,"ij")==0;
  bool is_none  = strcmp(type,"none")==0;
  if (! is_xyz && ! is_ij & ! is_none) {
    printf ("type = %s  != [xyz|ij|none]\n");
    exit(1);
  }

  printf ("size = %d  mask = %d  type = %s\n",N,mask,type);

  FILE *fp;
  int i0,i1,i2;
  int j0,j1,j2;
  int i,j,k;
  double value;

  if (mask & 1 ) {
    fp = fopen ("A.00.00.00.00000","r");

    while (fscanf (fp,"%d %d %d %d %lf",&i0,&i1,&i2,&k,&value) != EOF) {

      i0 += N/2+1;
      i1 += N/2+1;
      i2 += N/2+1;

      j0 = i0;
      j1 = i1;
      j2 = i2;
      
      if (k==1) ++j0;
      if (k==2) --j0;
      if (k==3) ++j1;
      if (k==4) --j1;
      if (k==5) ++j2;
      if (k==6) --j2;

      i = INDEX(i0,i1,i2,0,N2);
      j = INDEX(j0,j1,j2,0,N2);

      if (value != 0) {
	assert (0 <= i && i < Asize);
	assert (0 <= j && j < Asize);
	A[AINDEX(i0,i1,i2,k,0,N2)] = value;
	if (is_ij)  printf ("0  %d %d %g\n",i+1,j+1,value);
	if (is_xyz) printf ("0  %d %d %d  %d %d %d  %g\n",
			    i0,i1,i2,j0,j1,j2,value);
      }

    }
    fclose(fp);
  }
  
  if (mask & 2 ) {
    fp = fopen ("A.01.00.00.00000","r");
    while (fscanf (fp,"%d %d %d %d %lf",&i0,&i1,&i2,&k,&value) != EOF) {


      i0 += N/2+1;
      i1 += N/2+1;
      i2 += N/2+1;

      j0 = i0;
      j1 = i1;
      j2 = i2;
      
      if (k==1) ++j0;
      if (k==2) --j0;
      if (k==3) ++j1;
      if (k==4) --j1;
      if (k==5) ++j2;
      if (k==6) --j2;

      i = INDEX(i0,i1,i2,1,N2);
      j = INDEX(j0,j1,j2,1,N2);

      if (value != 0) {
	assert (0 <= i && i < Asize);
	assert (0 <= j && j < Asize);
	A[AINDEX(i0,i1,i2,k,1,N2)] = value;
	if (is_ij)  printf ("1  %d %d %g\n",i+1,j+1,value);
	if (is_xyz) printf ("1  %d %d %d  %d %d %d  %g\n",
			    i0,i1,i2,j0,j1,j2,value);
      }

    }
    fclose(fp);
  }

  typedef std::pair<const int,int> IJ_type;
  std::map<IJ_type,double> Agraph;

  if (mask & 4 ) {
    fp = fopen ("A.UMatrix.00000","r");
    fscanf (fp,"%d %d %d %d",&i0,&i1,&i2,&k); // throw out matrix size
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
	if (is_ij)  printf ("12  %d %d %g\n",i+1,j+1,value);
	if (is_xyz) printf ("12  %d %d %d  %d %d %d  %g\n",
			    i0,i1,i2,j0,j1,j2,value);
      }
    }
    fclose(fp);
  }

  printf ("Testing Stencil symmetry\n");

  // Check problem x-y-z symmetry for stencil values

  for (int grid=0; grid<2; grid++) {
    for (i0=0; i0<N2; i0++) {
      for (i1=0; i1<N2; i1++) {
	for (i2=0; i2<N2; i2++) {
	  for (int k=0; k<7; k++) {

	    int a1 = AINDEX(i0,i1,i2,k,grid,N2);
	    int a2;

	    // Check diagonal
	    if (k==0) {
	      for (int k0=0; k0<2; k0++) {
		for (int k1=0; k1<2; k1++) {
		  for (int k2=0; k2<2; k2++) {
		    j0=(1-k0)*i0 + k0*(N2-i0-1);
		    j1=(1-k1)*i1 + k1*(N2-i1-1);
		    j2=(1-k2)*i2 + k2*(N2-i2-1);
		    a2 = AINDEX(j0,j1,j2,k,grid,N2);
		    if (A[a1] != A[a2]) {
		      fprintf (stderr,"OOPS! grid %d diagonal "
			       "(%d %d %d; %d) = %g  (%d %d %d; %d) = %g\n",
			       grid,
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
	      a2 = AINDEX(j0,j1,j2,3-k,grid,N2);
	      if (A[a1] != A[a2]) {
		fprintf (stderr,"OOPS! grid %d X "
			 "(%d %d %d; %d) = %g  (%d %d %d; %d) = %g\n",
			 grid,
			 i0,i1,i2,k,A[a1],
			 j0,j1,j2,3-k,A[a2]);
	      }
	    }
	    if (k==3 || k==4) {
	      j0=i0;
	      j1=N2-i1-1;
	      j2=i2;
	      a2 = AINDEX(j0,j1,j2,7-k,grid,N2);
	      if (A[a1] != A[a2]) {
		fprintf (stderr,"OOPS! grid %d Y "
			 "(%d %d %d; %d) = %g  (%d %d %d; %d) = %g\n",
			 grid,
			 i0,i1,i2,k,A[a1],
			 j0,j1,j2,7-k,A[a2]);
	      }
	    }
	    if (k==5 || k==6) {
	      j0=i0;
	      j1=i1;
	      j2=N2-i2-1;
	      a2 = AINDEX(j0,j1,j2,11-k,grid,N2);
	      if (A[a1] != A[a2]) {
		fprintf (stderr,"OOPS! grid %d Z "
			 "(%d %d %d; %d) = %g  (%d %d %d; %d) = %g\n",
			 grid,
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

}
