#include <stdio.h>
#include <stdlib.h>

#include <string>

#include "hdf5.h"

/**
 * @file      enzo2hypre.cpp
 * @brief     Create a hypre-solve input file given an Enzo data dump
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
 */

//----------------------------------------------------------------------
// FUNCTIONS
//----------------------------------------------------------------------

void print_particles (std::string,FILE*fpout);

//----------------------------------------------------------------------
// MACROS
//----------------------------------------------------------------------

#define HDF_ERROR(MESSAGE) \
  { \
    fprintf (stderr,"%s:%d HDF5 ERROR: %s\n",__FILE__,__LINE__,MESSAGE); \
    exit(1); \
  }

#define max(a,b) ((a)>(b) ? (a) : (b))

//----------------------------------------------------------------------
// CONSTANTS
//----------------------------------------------------------------------

const int refinement_factor = 2;
const int max_levels        = 30;
const int max_line_length   = 100;

//----------------------------------------------------------------------
// MAIN
//----------------------------------------------------------------------

int main(int argc, char ** argv)
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
  // Write constant header
  //----------------------------------------------------------------------
  
  fprintf (fpout," dimension 3\n");
  fprintf (fpout," discret constant\n");
  fprintf (fpout," solver fac\n");
  fprintf (fpout," boundary periodic\n");
  fprintf (fpout," dump_hypre true\n");

  //----------------------------------------------------------------------
  // Open hierarchy file to get number of grids "num_grids"
  //----------------------------------------------------------------------

  std::string hierarchy_file = std::string(argv[1]) + ".hierarchy";

  FILE *fpin  = fopen (hierarchy_file.c_str(), "r");

  char line[max_line_length];

  int igrid,jgrid;
  int num_grids=0;

  while (fgets(line,max_line_length,fpin) != NULL) {
    if (sscanf (line,"Grid = %d\n",&igrid)) num_grids = max(igrid,num_grids);
  }

  fclose(fpin);   // Close hierarchy file

  //-----------------------------------------------------------------------
  // Allocate arrays
  //-----------------------------------------------------------------------

  // grid arrays

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
  int * ibegin0     = new int [num_grids+1];
  int * ibegin1     = new int [num_grids+1];
  int * ibegin2     = new int [num_grids+1];
  int * iend0     = new int [num_grids+1];
  int * iend1     = new int [num_grids+1];
  int * iend2     = new int [num_grids+1];
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
    ibegin0[i] = 0;
    ibegin1[i] = 0;
    ibegin2[i] = 0;
    iend0[i] = 0;
    iend1[i] = 0;
    iend2[i] = 0;
  }

  //-----------------------------------------------------------------------
  // Reread hierarchy file to compute arrays
  //-----------------------------------------------------------------------

  fpin  = fopen (hierarchy_file.c_str(), "r");

  char s [ max_line_length ];
  // Read input Enzo hierarchy file

  while (fgets(line,max_line_length,fpin) != NULL) {

    // Get grid number

    sscanf (line,"Grid = %d",&igrid); // get the current grid number

    // Get grid extent

    sscanf (line,"GridLeftEdge      = %lg %lg %lg",
	    &xmin0[igrid],&xmin1[igrid],&xmin2[igrid]);
    sscanf (line,"GridRightEdge      = %lg %lg %lg",
	    &xmax0[igrid],&xmax1[igrid],&xmax2[igrid]);

    // Get grid size

    // (note: GridDimension is not portable since different Enzo versions
    //  note: may be using different ghost zone counts)

    sscanf (line,"GridStartIndex    = %d %d %d",
	    &ibegin0[igrid],&ibegin1[igrid],&ibegin2[igrid]);
    sscanf (line,"GridEndIndex      = %d %d %d",
	    &iend0[igrid],&iend1[igrid],&iend2[igrid]);

    if (sscanf (line,"%s",s) && strcmp (s,"Pointer:")==0) {

      // Get sibling grid pointer

      if (sscanf (line+9,"Grid[%d]->NextGridThisLevel = %d",
		  &igrid,&jgrid)) {
	if (strstr (line,"NextGridThisLevel")) {
	  if (jgrid != 0) {
	    parent[jgrid] = parent[igrid];
	    level[jgrid]  = level[igrid];
	  }
	}
      }

      // Get child grid pointer

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

  fclose(fpin);   // Close hierarchy file

  //-----------------------------------------------------------------------
  // Read processor map file to get processor-to-grid mapping
  //-----------------------------------------------------------------------

  std::string procmap_file = std::string(argv[1]) + ".procmap";

  fpin  = fopen (procmap_file.c_str(), "r");

  int np = 0;

  bool have_procmap = (fpin != NULL);

  fprintf (fpout," enzo_density true\n");
  fprintf (fpout," enzo_procmap %s\n",have_procmap ? "true" : "false");
  fprintf (fpout," enzo_prefix  %s\n",argv[1]);

  if (have_procmap) {

    // If *.procmap file exists, set processor id's accordingly

    int iproc;
    char sgrid[20], sproc[100],dummy[20];
    while (fgets(line,max_line_length,fpin) != NULL) {
      sscanf (line,"%s %s %s",sgrid,sproc,dummy);
      sscanf (sgrid,"%d",&igrid);
      sscanf (strstr(sproc,"cpu")+3,"%d",&iproc);
      processor[igrid] = iproc;
      np = max(np,iproc);
    }
    np++;

    fclose(fpin);  // Close procmap file

  } else {

    // If *.procmap file does not exist, set all processor id's to 0

    np = 1;
    for (igrid = 1; igrid <= num_grids; igrid++) {
      processor[igrid] = 0;
    }

  }

  //-----------------------------------------------------------------------
  // Read restart file to get domain extents and top grid dimensions
  //-----------------------------------------------------------------------

  std::string restart_file = std::string(argv[1]);
  fpin  = fopen (restart_file.c_str(), "r");
  double dmin0,dmin1,dmin2;
  double dmax0,dmax1,dmax2;
  int n0,n1,n2;
  while (fgets(line,max_line_length,fpin) != NULL) {
    sscanf (line,"DomainLeftEdge         = %lg %lg %lg",&dmin0,&dmin1,&dmin2);
    sscanf (line,"DomainRightEdge        = %lg %lg %lg",&dmax0,&dmax1,&dmax2);
    sscanf (line,"TopGridDimensions   = %d %d %d",&n0,&n1,&n2);
  }

  fprintf (fpout, " domain 3 %g %g %g  %g %g %g\n",
	   dmin0,dmin1,dmin2,
	   dmax0,dmax1,dmax2);

  fclose(fpin);  // Close restart file

  //-----------------------------------------------------------------------
  // Compute imin given xmin, level, and grid size
  //-----------------------------------------------------------------------

  int levelpow[max_levels];
  levelpow[0] = 1;
  for (i=1; i<max_levels; i++) {
    levelpow[i] = levelpow[i-1]*refinement_factor;
  }
  
  for (i=1; i<=num_grids; i++) {
    int levelfactor = levelpow[level[i]];
    imin0[i] = int((xmin0[i] - dmin0) / (dmax0-dmin0) * n0 * levelfactor);
    imin1[i] = int((xmin1[i] - dmin1) / (dmax1-dmin1) * n1 * levelfactor);
    imin2[i] = int((xmin2[i] - dmin2) / (dmax2-dmin2) * n2 * levelfactor);
  }

  //-----------------------------------------------------------------------
  // Output grids
  //-----------------------------------------------------------------------

  for (i=1; i<=num_grids; i++) {
    fprintf (fpout,"grid");
    fprintf (fpout," %d",i-1);
    fprintf (fpout," %d",parent[i]-1);
    fprintf (fpout," %d",processor[i]); 
    fprintf (fpout," %g %g %g",xmin0[i],xmin1[i],xmin2[i]);
    fprintf (fpout," %g %g %g",xmax0[i],xmax1[i],xmax2[i]);
    fprintf (fpout," %d %d %d",imin0[i],imin1[i],imin2[i]);
    fprintf (fpout," %d %d %d",
	     iend0[i]-ibegin0[i]+1,
	     iend1[i]-ibegin1[i]+1,
	     iend2[i]-ibegin2[i]+1);

    fprintf (fpout,"\n");
  }

  //-----------------------------------------------------------------------
  // Output particles/points
  //-----------------------------------------------------------------------

  // particle arrays

  if (have_procmap) {

    // Loop over *.cpu* files

    for (int ip = 0; ip < np; ip++) {
      char cpufile[40];
      sprintf (cpufile,"%s.cpu%04d",std::string(argv[1]).c_str(),ip);
      print_particles (cpufile,fpout);
    }
  } else {

    // Loop over *.grid* files

    for (int igrid = 1; igrid <= num_grids; igrid++) {
      char gridfile[40];
      sprintf (gridfile,"%s.grid%04d",std::string(argv[1]).c_str(),igrid);
      print_particles (gridfile,fpout);
    }
  }
  //  for (

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
  delete [] ibegin0;
  delete [] ibegin1;
  delete [] ibegin2;
  delete [] iend0;
  delete [] iend1;
  delete [] iend2;

  fclose (fpout);  // Close output file

  exit (0);
}

void print_particles (std::string file_name, FILE * fpout)
{

  hid_t       file_id;

  // Open the file

  file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  if (file_id < 0) HDF_ERROR (file_name.c_str());

  hid_t group_id;

  // Open the group

  group_id = H5Gopen (file_id, "/");
  if (group_id < 0) HDF_ERROR ("H5Gopen");

  hsize_t num_groups;
  H5Gget_num_objs(group_id,&num_groups);

  herr_t  status;

  for (hsize_t igroup = 0; igroup < num_groups; igroup ++) {

    char group_name[20];
    H5Gget_objname_by_idx(group_id,igroup,group_name,20);

    std::string grid_group = group_name;

    // Open the particle datasets if any

    hid_t pmset_id;
    hid_t ppxset_id,ppyset_id,ppzset_id;

    pmset_id = H5Dopen(group_id, (grid_group + "/particle_mass").c_str());
    if (pmset_id < 0) break; // exit loop if no particles
    ppxset_id = H5Dopen(group_id, (grid_group + "/particle_position_x").c_str());
    ppyset_id = H5Dopen(group_id, (grid_group + "/particle_position_y").c_str());
    ppzset_id = H5Dopen(group_id, (grid_group + "/particle_position_z").c_str());

    // Get the dataset size

    hid_t pmspace_id;
    hid_t ppxspace_id,ppyspace_id,ppzspace_id;

    hsize_t num_particles = 0;
    pmspace_id  = H5Dget_space (pmset_id);
    if (pmspace_id < 0) HDF_ERROR ("H5Dopen");
    ppxspace_id = H5Dget_space (ppxset_id);
    if (ppxspace_id < 0) HDF_ERROR ("H5Dopen");
    ppyspace_id = H5Dget_space (ppyset_id);
    if (ppyspace_id < 0) HDF_ERROR ("H5Dopen");
    ppzspace_id = H5Dget_space (ppzset_id);
    if (ppzspace_id < 0) HDF_ERROR ("H5Dopen");

    // Get rank (should be same for all particle datasets)

    H5Sget_simple_extent_dims(ppxspace_id, &num_particles, NULL);

    double * pmass = new double [num_particles];
    double * ppos0 = new double [num_particles];
    double * ppos1 = new double [num_particles];
    double * ppos2 = new double [num_particles];
    
    int igrid = atoi(grid_group.c_str()+4);

    // Read the dataset

    status = H5Dread(pmset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
		     H5P_DEFAULT, pmass);
    if (status < 0) HDF_ERROR ("H5Dread");
    status = H5Dread(ppxset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
		     H5P_DEFAULT, ppos0);
    if (status < 0) HDF_ERROR ("H5Dread");
    status = H5Dread(ppyset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
		     H5P_DEFAULT, ppos1);
    if (status < 0) HDF_ERROR ("H5Dread");
    status = H5Dread(ppzset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
		     H5P_DEFAULT, ppos2);
    if (status < 0) HDF_ERROR ("H5Dread");

    for (hsize_t  i=0; i<num_particles; i++) {
      fprintf (fpout, "point");
      fprintf (fpout, " %g",pmass[i]);
      fprintf (fpout, " %g %g %g",ppos0[i],ppos1[i],ppos2[i]);
      fprintf (fpout, " %d",igrid-1);
      fprintf (fpout,"\n");
    }

    // Deallocate arrays 

    delete [] pmass;
    delete [] ppos0;
    delete [] ppos1;
    delete [] ppos2;

    // Close the dataset

    status = H5Dclose(pmset_id);
    if (status < 0) HDF_ERROR ("H5Dclose");
    status = H5Dclose(ppxset_id);
    if (status < 0) HDF_ERROR ("H5Dclose");
    status = H5Dclose(ppyset_id);
    if (status < 0) HDF_ERROR ("H5Dclose");
    status = H5Dclose(ppzset_id);
    if (status < 0) HDF_ERROR ("H5Dclose");

  }

  status = H5Gclose(group_id);  // Close the group

  status = H5Fclose(file_id);   // Close the file
}
