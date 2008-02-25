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

//----------------------------------------------------------------------
void usage (char ** argv)
{
  // Print usage and exit abnormally

  printf ("\n"
	  "Usage: %s N0 np0 np1 np2 num_levels offset\n",argv[0]);
  exit(1);
}

//----------------------------------------------------------------------

main(int argc, char **argv)
{

//-----------------------------------------------------------------------
// Parse command-line arguments
//-----------------------------------------------------------------------

  if (argc != 7) {
    printf ("Number of arguments %d is not 6\n",argc-1);
    usage(argv);
  }

  int N0, np3[3],num_levels,is_offset;

  N0         = atoi(argv[1]);
  np3[0]     = atoi(argv[2]);
  np3[1]     = atoi(argv[3]);
  np3[2]     = atoi(argv[4]);
  num_levels = atoi(argv[5]);
  is_offset  = atoi(argv[6]);

  int np = np3[0]*np3[1]*np3[2];

  //-----------------------------------------------------------------------
  // Generate input file
  //-----------------------------------------------------------------------

  // Get input, output, and directory names

  char dir[STRLEN],infile[STRLEN],outfile[STRLEN];
  sprintf (dir,    "N%d.P%d%d%d.L%d.O%d",N0,np3[0],np3[1],np3[2],num_levels,is_offset);
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

  int r = 1;
  for (level = 0; level < num_levels; level++, r *= 2) {

    for (ip3[2] = 0; ip3[2] < np3[2]; ip3[2]++) {

      // Centered box
      lp3[2] = (1-is_offset)*(-BOXSIZE/2 +  ip3[2]   *BOXSIZE / np3[2])/r;
      up3[2] = (1-is_offset)*(-BOXSIZE/2 + (ip3[2]+1)*BOXSIZE / np3[2])/r;

      // Offset box
      lp3[2] += is_offset*BOXSIZE* ip3[2]   /r/np3[2];
      up3[2] += is_offset*BOXSIZE*(ip3[2]+1)/r/np3[2];

      li3[2] = (1-is_offset)*N0*(r-1)/2 + ip3[2]*n3[2];

      jp3[2] = (ip3[2])/2 + np3[2]/4;

      if (np3[2] == 2) jp3[2] = ip3[2];
   
      for (ip3[1] = 0; ip3[1] < np3[1]; ip3[1]++) {

	// Centered box
	lp3[1] = (1-is_offset)*(-BOXSIZE/2 +  ip3[1]   *BOXSIZE / np3[1])/r;
	up3[1] = (1-is_offset)*(-BOXSIZE/2 + (ip3[1]+1)*BOXSIZE / np3[1])/r;

	// Offset box
	lp3[1] += is_offset*BOXSIZE* ip3[1]   /r/np3[1];
	up3[1] += is_offset*BOXSIZE*(ip3[1]+1)/r/np3[1];

	li3[1] = (1-is_offset)*N0*(r-1)/2 + ip3[1]*n3[1];
	jp3[1] = (ip3[1])/2 + np3[1]/4;

	if (np3[1] == 2) jp3[1] = ip3[1];

	for (ip3[0] = 0; ip3[0] < np3[0]; ip3[0]++) {

	  // Centered box
	  lp3[0] = (1-is_offset)*(-BOXSIZE/2 +  ip3[0]   *BOXSIZE / np3[0])/r;
	  up3[0] = (1-is_offset)*(-BOXSIZE/2 + (ip3[0]+1)*BOXSIZE / np3[0])/r;
	  // Offset box
	  lp3[0] += is_offset*BOXSIZE* ip3[0]   /r/np3[0];
	  up3[0] += is_offset*BOXSIZE*(ip3[0]+1)/r/np3[0];

	  li3[0] = (1-is_offset)*N0*(r-1)/2 + ip3[0]*n3[0];
	  jp3[0] = (ip3[0])/2 + np3[0]/4;

	  if (np3[0] == 2) jp3[0] = ip3[0];

	  int ip = ip3[0] + np3[0]*(ip3[1] + np3[1]*ip3[2]);

	  int id = ip + np*level;

	  int id_parent;
	  if (level == 0) {
	    id_parent = -1;
	  } else {
	    printf ("DEBUG %s:%d %d (%d %d %d) (%d %d %d)\n",__FILE__,__LINE__,
		    r,ip3[0],ip3[1],ip3[2],
		    np3[0],np3[1],np3[2]);
	    id_parent = ip3[0]/r + np3[0]*(ip3[1]/r + np3[1]*ip3[2]/r);
	    printf ("%d = %d/%d + %d*(%d/%d + %d*%d/%d\n",
		    id_parent,ip3[0],r,np3[0],ip3[1],r,np3[1],ip3[2],r);
	    //	    id_parent = jp3[0] + np3[0]*(jp3[1] + np3[1]*(jp3[2] + np3[2]*(level-1)));
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
	      + 0.25*is_offset*BOXSIZE;
	    for (int k1=0; k1<2; k1++) {
	      double p1 = (k1-0.5)*BOXSIZE/(N0*pow(2.0,level-1))*0.5  
		+ 0.25*is_offset*BOXSIZE;
	      for (int k2=0; k2<2; k2++) {
		double p2 = (k2-0.5)*BOXSIZE/(N0*pow(2.0,level-1))*0.5  
		  + 0.25*is_offset*BOXSIZE;
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
  if (is_offset) {
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

