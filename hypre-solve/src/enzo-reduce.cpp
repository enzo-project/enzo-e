//----------------------------------------------------------------------
// 
// PROGRAM
//
//   enzo-reduce
//
// DESCRIPTION
//
//   Perform a reduction operation (projection, slice, etc.) along an
//   axis of an Enzo data dump
//
// USAGE
//
//   plot-amr [options]
//
// (*)  -a  --along     {x|y|z}             Axis along which to reduce
// (*)  -f  --field     <field-name>        Name of the field
// (*)  -g  --grid      <grid file>         Specify an Enzo grid file
// (*)  -h  --hierarchy <hierarchy dir>     Specify an Enzo hierarchy directory
// (*)  -p  --project                       Project values
// (*)  -s  --slice     [0:1]               Slice domain at the given point 
//
// EXAMPLES
//
// BUGS
//
//----------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>
#include <sys/stat.h>

#include <string>
#include <vector>

#include "hdf5.h"

#ifdef USE_MPI
#  include "mpi.h"
#endif

//----------------------------------------------------------------------
// USER DATATYPES
//----------------------------------------------------------------------

enum enum_amr_type { amr_type_unknown, amr_type_packed, amr_type_unpacked };

struct argument_struct {
  std::string reduce_axis;
  std::string field_name;
  std::string grid_file;
  std::string hierarchy_dir;
  std::string slice_axis;
  std::string operation;
};

struct grid_info_struct {
  int     index;
  int     level;
  double  left_edge[3];
  double  right_edge[3];
  int     start_index[3];
  int     stop_index[3];
};

struct hierarchy_info_struct {
  double  left_edge[3];
  double  right_edge[3];
  int     num_levels;
  int     start_index[3];
  int     stop_index[3];
};


//----------------------------------------------------------------------
// FUNCTION PROTOTYPES
//----------------------------------------------------------------------

void          read_enzo_grid   (std::string file_name);

enum_amr_type get_amr_type (std::string arg_hierarchy_dir);

bool          check_args (argument_struct *);

void          print_usage(char ** argv);

int           read_grid_info(std::string        hierarchy_name, 
			     grid_info_struct * grid_info);

void          get_hierarchy_info(hierarchy_info_struct * hierarchy_info,
				 grid_info_struct      * grid_info, 
				 int                     num_local_grids);


//----------------------------------------------------------------------
// CPP MACROS
//----------------------------------------------------------------------

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))

#define MAX_FILE_STRING_LENGTH 80

#define TRACE \
    printf ("%s:%d TRACE\n",__FILE__,__LINE__); fflush(stdout);

#define NOT_DONE \
    printf ("%s:%d Code not implemented: exiting\n",__FILE__,__LINE__); \
    fflush(stdout); \
    exit(1)();

#define ERROR(message) \
      fprintf (stderr, \
	       "%s:%d %s\n", \
	       __FILE__,__LINE__,message); \
      exit(1);

#define INDEX_TO_PROC(index) ((index) % mpi_size)

//====================================================================
// BEGIN main ()
//====================================================================

int mpi_size = 1, mpi_rank = 0;

int main(int argc, char **argv) 
{

  //====================================================================
  // PRINT HEADER
  //====================================================================

  //====================================================================
  // INITIALIZE MPI
  //====================================================================

#ifdef USE_MPI
  MPI_Init (&argc,&argv);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
#endif

  //====================================================================
  // PARSE COMMAND LINE
  //====================================================================

  const struct option long_options[] =
    {
      {"along",     required_argument, 0, 'a'},
      {"field",     required_argument, 0, 'f'},
      {"grid",      required_argument, 0, 'g'},
      {"hierarchy", required_argument, 0, 'h'},
      {"project",   no_argument,       0, 'p'},
      {"slice",     required_argument, 0, 's'},
      {0, 0, 0, 0}
    };

  //--------------------------------------------------
  // Parse arguments
  //--------------------------------------------------

  argument_struct arg;

  int c;
  int option;
  while ((c = getopt_long(argc, argv, "a:f:g:h:ps:", long_options, &option)) != -1) {
    switch (c) {
    case 0:
      // argument set in getopt_long()
      break;
    case 'a': arg.reduce_axis   = optarg; break;
    case 'f': arg.field_name    = optarg; break;
    case 'g': arg.grid_file     = optarg; break;
    case 'h': arg.hierarchy_dir = optarg; break;
    case 'p': arg.operation = "project";  break;
    case 's': arg.operation = "slice";	
      arg.slice_axis = optarg;
      break;
    case '?': print_usage(argv); break;
    default:
      fprintf (stderr,"%s:%d Unknown command-line parse error: exiting\n",
	       __FILE__,__LINE__);
      fprintf (stderr,"opterr = %d\n",opterr);
      fprintf (stderr,"optarg = %s\n",optarg);
      exit(1);
      break;
    }
  }

  //--------------------------------------------------
  // Verify arguments
  //--------------------------------------------------

  if (check_args(&arg)) { 

    printf ("\n");
    printf ("  Reduction axis      = %s\n", arg.reduce_axis.c_str());
    printf ("  Field name          = %s\n", arg.field_name.c_str());
    printf ("  Grid file           = %s\n", arg.grid_file.c_str());
    printf ("  Hierarchy directory = %s\n", arg.hierarchy_dir.c_str());
    printf ("  Slice axis          = %s\n", arg.slice_axis.c_str());
    printf ("  Operation           = %s\n", arg.operation.c_str() );
    printf ("\n");

  } else {

    print_usage(argv);

  }

  //====================================================================
  // DETERMINE LIST OF INPUT FILES
  //====================================================================

  bool      grid_specified = ! arg.grid_file.empty() ;
  bool hierarchy_specified = ! arg.hierarchy_dir.empty();
  

  std::vector <std::string> grid_file_list;

  if ((grid_specified)) {

    grid_file_list.push_back(arg.grid_file);

  } else if ((hierarchy_specified)) {

    chdir (arg.hierarchy_dir.c_str());
    // ERROR CHECK HERE


    // Read hierarchy file

    char file_name [MAX_FILE_STRING_LENGTH];

    enum_amr_type amr_type = get_amr_type (arg.hierarchy_dir);

    // Loop through files and add them to the list of files

    for (int i = mpi_rank; true; i += mpi_size) {

      if (amr_type == amr_type_unpacked) {
	sprintf (file_name, "%s.grid%04d",arg.hierarchy_dir.c_str(),i+1);
      } else if (amr_type == amr_type_packed) {
	sprintf (file_name, "%s.cpu%04d",arg.hierarchy_dir.c_str(),i);
      }

      // Check if file exists, and append to list of files if it does

      struct stat stFileInfo; 

      if (stat(file_name,&stFileInfo) == 0) {

	grid_file_list.push_back(file_name);

      } else {

	break;

      }

    }

    // ERROR CHECK ! read_successful
    
  }

  // Display list of files

  //  for (unsigned int i=0; i < grid_file_list.size(); i++) {
  //    printf ("%d %s\n",mpi_rank,grid_file_list[i].c_str());
  //  }

  // READ GRID INFO AND DETERMINE DOMAIN EXTENTS [left|right]_edge_global

  grid_info_struct *grid_info = new grid_info_struct[grid_file_list.size()];

  int num_local_grids = read_grid_info
    (arg.hierarchy_dir,grid_info);

  hierarchy_info_struct hierarchy_info;

  get_hierarchy_info(&hierarchy_info,grid_info,num_local_grids);

  //====================================================================
  // ACCUMULATE PROCESSOR-LOCAL REDUCTION
  //====================================================================

  
  //====================================================================
  // ACCUMULATE PROCESSOR-LOCAL REDUCTION
  //====================================================================

  // First loop through files to get size

  //====================================================================
  // REDUCE LOCAL REDUCTION TO GLOBAL ON ROOT
  //====================================================================

  //====================================================================
  // OUTPUT RESULT FROM ROOT
  //====================================================================

  //====================================================================
  // EXIT
  //====================================================================

#ifdef USE_MPI
  MPI_Finalize ();
#endif
  
  return 0;
}

//====================================================================
// END main ()
//====================================================================



//--------------------------------------------------------------------
// read_enzo_grid ()
//--------------------------------------------------------------------

// Read a dataset from an HDF5 file and return it in the given array

void read_enzo_grid 
(
 double   ** array, 
 std::string file_name, 
 std::string dataset_name
)

{

  hid_t       file_id;
  hid_t       dataset_id;
  hid_t       dataspace_id;
  herr_t      status;

  //--------------------------------------------------
  // Open the file
  //--------------------------------------------------
 
  file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  //--------------------------------------------------
  // Open the dataset
  //--------------------------------------------------

  dataset_id = H5Dopen(file_id, ((std::string) ("/" + dataset_name)).c_str());

  //--------------------------------------------------
  // Get the dataset size
  //--------------------------------------------------

  hsize_t n[3] = {0,0,0};

  dataspace_id = H5Dget_space (dataset_id);

  H5Sget_simple_extent_dims(dataspace_id, n, NULL);

  assert (n[0]!=0);
  assert (n[1]!=0);
  assert (n[2]!=0);

  //--------------------------------------------------
  // (Re)allocate storage for the dataset
  //--------------------------------------------------

  if (*array != NULL) delete [] *array;
  *array = new double [n[0]*n[1]*n[2]];

  //--------------------------------------------------
  // Read the dataset
  //--------------------------------------------------

  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
		   H5P_DEFAULT, *array);
    
  //--------------------------------------------------
  // Close the dataset
  //--------------------------------------------------

  status = H5Dclose(dataset_id);

  //--------------------------------------------------
  // Close the file
  //--------------------------------------------------

  status = H5Fclose(file_id);

}

int get_num_grids (std::string hierarchy_file)
{

  FILE *fp_hierarchy  = fopen (hierarchy_file.c_str(), "r");

  char line[MAX_FILE_STRING_LENGTH+1];

  int ind_grid;
  int num_grids=0;

  while (fgets(line, MAX_FILE_STRING_LENGTH,fp_hierarchy) != NULL) {
    if (sscanf (line,"Grid = %d\n",&ind_grid)) {
      num_grids = MAX(ind_grid,num_grids);
    }
  }

  fclose(fp_hierarchy);

  return num_grids;
}

//--------------------------------------------------------------------
// get_amr_type ()
//--------------------------------------------------------------------

// Determine whether the amr data is packed or unpacked

enum_amr_type get_amr_type (std::string hierarchy_dir)
{

  enum_amr_type amr_type;

  struct stat stFileInfo; 
  char file_name [MAX_FILE_STRING_LENGTH];

  amr_type = amr_type_unknown;

  sprintf (file_name, "%s.grid0001",hierarchy_dir.c_str());
  if (stat(file_name,&stFileInfo) == 0) {
    amr_type = amr_type_unpacked;
  }

  sprintf (file_name, "%s.cpu0000",hierarchy_dir.c_str());
  if (stat(file_name,&stFileInfo) == 0) {
    amr_type = amr_type_packed;
  }

  // exit if type cannot be determined

  if (amr_type == amr_type_unknown) {
    ERROR("Error in determining packed or unpacked)");
  }
  return amr_type;
}

//--------------------------------------------------------------------
// check_args()
//--------------------------------------------------------------------

bool check_args(argument_struct *arg_p)
{
  TRACE;
  bool args_ok = true;
  if (arg_p->reduce_axis == "") {
    printf ("argument '-a' required\n");
    args_ok = false;
  }

  if (arg_p->field_name  == "") {
    printf ("argument '-f' required\n");
    args_ok = false;
  }
	
  if ( (arg_p->grid_file == "" && arg_p->hierarchy_dir == "") ||
       (arg_p->grid_file != "" && arg_p->hierarchy_dir != "")) {
    printf ("argument '-g' or '-h' required\n");
    args_ok = false;
  }

  if ( arg_p->operation == "" ) {
    printf ("argument '-p' or '-s' required\n");
    args_ok = false;
  }

  return args_ok;

}

//--------------------------------------------------------------------
// print_usage()
//--------------------------------------------------------------------

void print_usage(char ** argv)
{
  printf ("\n\nUsage: %s [options]\n",argv[0]);
    
  printf ("\n");
  printf ("   -a (--along)      {x|y|z}             Axis along which to reduce\n");
  printf ("   -f (--field)      <field-name>        Name of the field\n");
  printf ("   -g (--grid)       <grid file>         Specify an Enzo grid file\n");
  printf ("   -h (--hierarchy)  <hierarchy dir>     Specify an Enzo hierarchy directory\n");
  printf ("   -p (--project)                        Project values\n");
  printf ("   -s (--slice)      [0:1]               Slice domain at the given point\n\n");
  printf ("\n");
  printf ("Required options: af[g|h][p|s]\n");
  printf ("\n");

  exit (1);

}


//-----------------------------------------------------------------------
// Read hierarchy file 
//-----------------------------------------------------------------------

int read_grid_info(std::string        hierarchy_name, 
		   grid_info_struct * grid_info)
{
  TRACE;
  std::string hierarchy_file = hierarchy_name + ".hierarchy";

  FILE *fp_hierarchy  = fopen (hierarchy_file.c_str(), "r");

  char line [ MAX_FILE_STRING_LENGTH ];

  // Read input Enzo hierarchy file

  int grid_index = 1;

  int index;

  while (fgets(line,MAX_FILE_STRING_LENGTH, fp_hierarchy) != NULL) {

    // Get grid number

    sscanf (line,"Grid = %d", &grid_index);

    if (INDEX_TO_PROC(grid_index-1) == mpi_rank) {

      int index = (grid_index - 1) / mpi_size;  // Update index

      // Get grid extent

      sscanf (line,"GridLeftEdge      = %lg %lg %lg",
	      &grid_info[index].left_edge[0],
	      &grid_info[index].left_edge[1],
	      &grid_info[index].left_edge[2]);
	    
      sscanf (line,"GridRightEdge      = %lg %lg %lg",
	      &grid_info[index].right_edge[0],
	      &grid_info[index].right_edge[1],
	      &grid_info[index].right_edge[2]);

      // Get grid size

      // (NOTE: Using GridDimension is not portable since different Enzo
      // versions may be using different ghost zone counts)

      sscanf (line,"GridStartIndex    = %d %d %d",
	      &grid_info[index].start_index[0],
	      &grid_info[index].start_index[1],
	      &grid_info[index].start_index[2]);

      sscanf (line,"GridEndIndex      = %d %d %d",
	      &grid_info[index].stop_index[0],
	      &grid_info[index].stop_index[1],
	      &grid_info[index].stop_index[2]);

    }
  }

  fclose(fp_hierarchy);   // Close hierarchy file

  int num_local_grids = index;

  // Determine grid levels

  double root_cell_width = 
    (grid_info[0].right_edge[0] - grid_info[0].left_edge[0]) / 
    (grid_info[0].stop_index[0] - grid_info[0].start_index[0] + 1) ;

  for (int i = 0; i < num_local_grids; i++) {

    double this_cell_width = 
      (grid_info[i].right_edge[0] - grid_info[i].left_edge[0]) /
      (grid_info[i].stop_index[0] - grid_info[i].start_index[0] + 1);

    grid_info[i].level = log (root_cell_width / this_cell_width) / log(2.0);

  }
  return index;
}

void get_hierarchy_info (hierarchy_info_struct * hierarchy_info,
			 grid_info_struct      * grid_info, 
			 int                     num_local_grids)
{
  double left_edge_local[3]  = {1e30,1e30,1e30};
  double right_edge_local[3] = {-1e30,-1e30,-1e30};

  for (int index=0; index<num_local_grids; index++) {
    
    // Update left_edge_local[] and right_edge_local[]
    
    left_edge_local[0]  = MIN( left_edge_local[0],grid_info[index].left_edge[0]);
    left_edge_local[1]  = MIN( left_edge_local[1],grid_info[index].left_edge[1]);
    left_edge_local[2]  = MIN( left_edge_local[2],grid_info[index].left_edge[2]);
    right_edge_local[0] = MAX(right_edge_local[0],grid_info[index].right_edge[0]);
    right_edge_local[1] = MAX(right_edge_local[1],grid_info[index].right_edge[1]);
    right_edge_local[2] = MAX(right_edge_local[2],grid_info[index].right_edge[2]);
  }

  // Ensure that it works without USE_MPI, or with USE_MPI but not called with mpirun

  hierarchy_info->left_edge[0] = left_edge_local[0];
  hierarchy_info->left_edge[1] = left_edge_local[1];
  hierarchy_info->left_edge[2] = left_edge_local[2];
  hierarchy_info->right_edge[0] = right_edge_local[0];
  hierarchy_info->right_edge[1] = right_edge_local[1];
  hierarchy_info->right_edge[2] = right_edge_local[2];

#ifdef USE_MPI  
  MPI_Allreduce (left_edge_local,hierarchy_info->left_edge,3,
		 MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce (right_edge_local,hierarchy_info->right_edge,3,
		 MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif

}
