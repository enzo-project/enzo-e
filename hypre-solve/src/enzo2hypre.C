#include <stdio.h>
#include <stdlib.h>

#include <string>

#define max(a,b) ((a)>(b) ? (a) : (b))

main(int argc, char ** argv)
{
  if (argc != 2) {
    fprintf (stderr,"Usage: %s <enzo-file>\n",argv[0]);
    exit(1);
  }

  //----------------------------------------------------------------------
  // Open output file
  //----------------------------------------------------------------------

  std::string outfile        = "in.hypre-solve";
  FILE *fpout = fopen (outfile.c_str(),"w");

  //----------------------------------------------------------------------
  // Open hierarchy_file to get number of grids "num_grids"
  //----------------------------------------------------------------------

  std::string hierarchy_file = std::string(argv[1]) + ".hierarchy";
  FILE *fpin  = fopen (hierarchy_file.c_str(), "r");

  char line[100];

  // Read input Enzo file to find number of grids "num_grids"

  int igrid,jgrid,num_grids=0;
  while (fgets(line,100,fpin) != NULL) {
    if (sscanf (line,"Grid = %d\n",&igrid)) num_grids = max(igrid,num_grids);
  }

  // Close hierarchy file
  fclose(fpin);

  //-----------------------------------------------------------------------
  // Allocate arrays
  //-----------------------------------------------------------------------

  int * processor = new int [num_grids+1];
  int * parent    = new int [num_grids+1];
  int * level     = new int [num_grids+1];
  double * xmin0  = new double [num_grids+1];
  double * xmin1  = new double [num_grids+1];
  double * xmin2  = new double [num_grids+1];
  double * xmax0  = new double [num_grids+1];
  double * xmax1  = new double [num_grids+1];
  double * xmax2  = new double [num_grids+1];
  int * imin0     = new int [num_grids+1];
  int * imin1     = new int [num_grids+1];
  int * imin2     = new int [num_grids+1];
  int * isize0     = new int [num_grids+1];
  int * isize1     = new int [num_grids+1];
  int * isize2     = new int [num_grids+1];
  int i;
  for (i=0; i<num_grids+1; i++) {
    processor[i] = 0;
    parent[i] = 0;
    level[i] = 0;
    xmin0[i] = 0;
    xmin1[i] = 0;
    xmin2[i] = 0;
    xmax0[i] = 0;
    xmax1[i] = 0;
    xmax2[i] = 0;
    imin0[i] = 0;
    imin1[i] = 0;
    imin2[i] = 0;
    isize0[i] = 0;
    isize1[i] = 0;
    isize2[i] = 0;
  }

  //-----------------------------------------------------------------------
  // Reread hierarchy file to compute arrays
  //-----------------------------------------------------------------------

  fpin  = fopen (hierarchy_file.c_str(), "r");

  printf ("Grids = %d\n",num_grids);

  char s [ 100 ];
  // Read input Enzo hierarchy file

  while (fgets(line,100,fpin) != NULL) {

    // Process grid number

    sscanf (line,"Grid = %d",&igrid); // get the current grid number

    // Process grid extent

    sscanf (line,"GridLeftEdge      = %lg %lg %lg",
	    &xmin0[igrid],&xmin1[igrid],&xmin2[igrid]);
    sscanf (line,"GridRightEdge      = %lg %lg %lg",
	    &xmax0[igrid],&xmax1[igrid],&xmax2[igrid]);

    // Process grid size

    sscanf (line,"GridDimension     = %d %d %d",
	    &isize0[igrid],&isize1[igrid],&isize2[igrid]);

    // Process pointer

    if (sscanf (line,"%s",s) && strcmp (s,"Pointer:")==0) {

      // Sibling grid
      if (sscanf (line+9,"Grid[%d]->NextGridThisLevel = %d",
		  &igrid,&jgrid)) {
	if (strstr (line,"NextGridThisLevel")) {
	  if (jgrid != 0) {
	    parent[jgrid] = parent[igrid];
	    level[jgrid]  = level[igrid];
	  }
	}
      }
      // Child grid
      if (sscanf (line+9,"Grid[%d]->NextGridNextLevel = %d",
		  &igrid,&jgrid)) {
	if (strstr (line,"NextGridNextLevel")) {
	  if (jgrid != 0) {
	    parent[jgrid] = igrid;
	    level[jgrid]  = level[igrid] + 1;
	  }
	}
      }
    }
  }

  // Close hierarchy file
  fclose(fpin);

  //-----------------------------------------------------------------------
  // Read processor map file to get processor-to-grid mapping
  //-----------------------------------------------------------------------

  std::string procmap_file = std::string(argv[1]) + ".procmap";
  fpin  = fopen (procmap_file.c_str(), "r");
  int iproc;
  char sgrid[20], sproc[100],dummy[20];
  while (fgets(line,100,fpin) != NULL) {
    sscanf (line,"%s %s %s",sgrid,sproc,dummy);
    sscanf (sgrid,"%d",&igrid);
    sscanf (strstr(sproc,"cpu")+3,"%d",&iproc);
    processor[igrid] = iproc;
  }

  //-----------------------------------------------------------------------
  // Compute imin given xmin, level, and grid size
  //-----------------------------------------------------------------------

  @@@

  //-----------------------------------------------------------------------
  // Print the arrays
  //-----------------------------------------------------------------------

  for (i=1; i<=num_grids; i++) {
    printf ("grid");
    printf (" %d",i);
    printf (" %d",parent[i]);
    printf (" %d",processor[i]); 
    printf (" %g %g %g",xmin0[i],xmin1[i],xmin2[i]);
    printf (" %g %g %g",xmax0[i],xmax1[i],xmax2[i]);
    printf (" %d %d %d",imin0[i],imin1[i],imin2[i]);
    printf (" %d %d %d",isize0[i],isize1[i],isize2[i]);
    printf ("\n");
  }

  printf ("%d\n",num_grids);

  // Deallocate arrays

  delete [] processor;
  delete [] parent;
  delete [] level;
  delete [] xmin0;
  delete [] xmin1;
  delete [] xmin2;
  delete [] xmax0;
  delete [] xmax1;
  delete [] xmax2;
  delete [] imin0;
  delete [] imin1;
  delete [] imin2;
  delete [] isize0;
  delete [] isize1;
  delete [] isize2;

  // Close output file
  fclose (fpout);
}
