//345678901234567890123456789012345678901234567890123456789012345678901234567890
///-----------------------------------------------------------------------
/// Generate an input file for hypre-solve
///
/// Usage: 
///
///   Unigrid of global size N0**3 distributed along <np0,np1,np2>
///   processor grid
///
///       run unigrid N0 np0 np1 np2
///
///   Simple nested grids in L levels with global root-grid size N0**3
///   distributed along <np0,np1,np2> processor grid
///
///       run nested  N0 np0 np1 np2 L
///
///   (note that "run nested  N0 np0 np1 np2 0" is the same as
///              "run unigrid N0 np0 np1 np2")
///
///-----------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#define BOUNDARY "dirichlet"
#define DIR_MODE 0755
#define STRLEN   80
#define BOXSIZE 8e9

#define IS_OFFSET 1

//----------------------------------------------------------------------
void usage (char ** argv)
{
  // Print usage and exit abnormally

  printf ("\n"
	  "Usage: %s unigrid N0 np0 np1 np2\n"
	  "   or\n"
	  "       %s nested  N0 np0 np1 np2 L\n",argv[0],argv[0]);
  exit(1);
}

//----------------------------------------------------------------------

main(int argc, char **argv)
{

//-----------------------------------------------------------------------
// Parse command-line arguments
//-----------------------------------------------------------------------

  if (argc < 6 || argc > 7) {
    printf ("Number of arguments %d not between 5 and 6\n",argc-1);
    usage(argv);
  }

  char type[STRLEN];
  strcpy (type,argv[1]);
  int N0, np3[3],L;

  N0  = atoi(argv[2]);
  np3[0] = atoi(argv[3]);
  np3[1] = atoi(argv[4]);
  np3[2] = atoi(argv[5]);

  if (strcmp(type,"unigrid") == 0) 
    L   = 1;
  else if (strcmp(type,"nested") == 0)  
    L   = atoi(argv[6]);
  else {
    usage(argv);
  }

  int np = np3[0]*np3[1]*np3[2];

  //-----------------------------------------------------------------------
  // Generate input file
  //-----------------------------------------------------------------------

  // Get input, output, and directory names

  char dir[STRLEN],infile[STRLEN],outfile[STRLEN];
  sprintf (dir,    "%s.P%d%d%d.N%d.L%d",type,np3[0],np3[1],np3[2],N0,L);
  sprintf (infile,"in.%s",dir);
  FILE * fp = fopen (infile,"w");

  printf ("Running %s: ",dir);

  // Create new directory...

  if (mkdir (dir,DIR_MODE) != 0) {
    fprintf (stderr,"Cannot make directory %s.\n",dir);
    exit(1);
  }

  // And change path to the new directory

  char path[STRLEN];
  getcwd(path, STRLEN);
  strcat (path,"/");
  strcat (path,dir);
  chmod (path,DIR_MODE);

  // Copy required files to new directory

  link ("../hypre-solve",".");

  // Compute local grid sizes
  // WARNING: assumes # processors along each dimension evenly divides problem size

  int i,n3[3],m3[3];

  for (i=0; i<3; i++) {
    n3[i] = N0 / np3[i];  
    if ((n3[i] * np3[i]) != N0) {
      fprintf (stderr,"Problem size %d not evenly divisible by np3[%d] = %d",
	       N0,i,np3[i]);
      exit(1);
    }
  }

  // Loop though processors, generating a grid per processor in each level

  int ip3[3];
  double lp3[3],up3[3];
  int li3[3];
  int id_point[8];
  double point_pos[8][3];
  int level;
  int jp3[3]; // for determining parent id 

  double r = 1.0;
  for (level = 0; level < L; level++, r *= 0.5) {

    for (ip3[2] = 0; ip3[2] < np3[2]; ip3[2]++) {

      lp3[2] = (1-IS_OFFSET)*r*(-0.5*BOXSIZE + ip3[2] * BOXSIZE / np3[2]);
      up3[2] = (1-IS_OFFSET)*r*(-0.5*BOXSIZE + (ip3[2]+1) * BOXSIZE / np3[2]);
      up3[2] += r*IS_OFFSET*BOXSIZE/np3[2];
      li3[2] = int((1-IS_OFFSET)*0.5*N0*(1./r-1.) + ip3[2] * n3[2]);
      jp3[2] = (ip3[2])/2 + np3[2]/4;
      if (np3[2] == 2) jp3[2] = ip3[2];
   
      for (ip3[1] = 0; ip3[1] < np3[1]; ip3[1]++) {

	lp3[1] = (1-IS_OFFSET)*r*(-0.5*BOXSIZE + ip3[1] * BOXSIZE / np3[1]);
	up3[1] = (1-IS_OFFSET)*r*(-0.5*BOXSIZE + (ip3[1]+1) * BOXSIZE / np3[1]);
	up3[1] += r*IS_OFFSET*BOXSIZE/np3[1];
	li3[1] = int((1-IS_OFFSET)*0.5*N0*(1./r-1.) + ip3[1] * n3[1]);
	jp3[1] = (ip3[1])/2 + np3[1]/4;
	if (np3[1] == 2) jp3[1] = ip3[1];

	for (ip3[0] = 0; ip3[0] < np3[0]; ip3[0]++) {

	  lp3[0] = (1-IS_OFFSET)*r*(-0.5*BOXSIZE + ip3[0] * BOXSIZE / np3[0]);
	  up3[0] = (1-IS_OFFSET)*r*(-0.5*BOXSIZE + (ip3[0]+1) * BOXSIZE / np3[0]);
	  up3[0] += r*IS_OFFSET*BOXSIZE/np3[0];
	  li3[0] = int((1-IS_OFFSET)*0.5*N0*(1./r-1.) + ip3[0] * n3[0]);
	  jp3[0] = (ip3[0])/2 + np3[0]/4;
	  if (np3[0] == 2) jp3[0] = ip3[0];

	  int ip = ip3[0] + np3[0]*(ip3[1] + np3[1]*ip3[2]);

	  int id = ip + np*level;

	  int id_parent;
	  if (level == 0) {
	    id_parent = -1;
	  } else {
	    // 1  0                 0
	    // 2  01                01
	    // 4  0123              1122
	    // 8  01234567          22334455
	    // 16 0123456789012345  4455667788990011
	    
	    id_parent = jp3[0] + np3[0]*(jp3[1] + np3[1]*(jp3[2] + np3[2]*(level-1)));
	  }

	  fprintf (fp, "grid %d %d %d "
		   "%g %g %g  %g %g %g  "
		   "%d %d %d  %d %d %d\n",
		   id,id_parent,ip,
		   lp3[0],lp3[1],lp3[2], up3[0],up3[1],up3[2],
		   li3[0],li3[1],li3[2],  n3[0], n3[1], n3[2]);

	  // Store id for point.  Requires coarse-to-finer level loop

	  for (int k0=0; k0<2; k0++) {
	    double p0 = (k0-0.5)*BOXSIZE/(N0*pow(2.0,level-1))*0.5 
	      + 0.25*IS_OFFSET*BOXSIZE;
	    for (int k1=0; k1<2; k1++) {
	      double p1 = (k1-0.5)*BOXSIZE/(N0*pow(2.0,level-1))*0.5  
		+ 0.25*IS_OFFSET*BOXSIZE;
	      for (int k2=0; k2<2; k2++) {
		double p2 = (k2-0.5)*BOXSIZE/(N0*pow(2.0,level-1))*0.5  
		  + 0.25*IS_OFFSET*BOXSIZE;
		int k = k0 + 2*(k1 + 2*k2);
		if (lp3[0] < p0 && p0 < up3[0] &&
		    lp3[1] < p1 && p1 < up3[1] &&
		    lp3[2] < p2 && p2 < up3[2]) {
		  id_point[k] = id;
		  point_pos[k][0] = p0;
		  point_pos[k][1] = p1;
		  point_pos[k][2] = p2;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  // Create problem

  fprintf (fp, "dimension 3\n");
  if (IS_OFFSET) {
    fprintf (fp, "domain    3 0e9 0e9 0e9  8e9 8e9 8e9\n");
  } else {
    fprintf (fp, "domain    3 -4e9 -4e9 -4e9  4e9 4e9 4e9\n");
  }
  fprintf (fp, "boundary  " BOUNDARY "\n");
  fprintf (fp, "sphere    5.993985e27 6.378137e8 0.0 0.0 0.0\n");
  for (int k=0; k<8; k++) {
    fprintf (fp, "point     1e43   %g %g %g   %d\n",
	     point_pos[k][0],point_pos[k][1],point_pos[k][2],id_point[k]);
  }
  fprintf (fp, "discret constant\n");
  fprintf (fp, "solver fac\n");

  //========================================================================

  // Exit normally

  fclose (fp);
    
  exit(0);
}

