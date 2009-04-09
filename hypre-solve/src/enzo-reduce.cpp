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

//----------------------------------------------------------------------
// FUNCTION PROTOTYPES
//----------------------------------------------------------------------

void read_enzo_grid   (std::string file_name);
int get_num_grids     (std::string arg_hierarchy_dir);
enum_amr_type get_amr_type (std::string arg_hierarchy_dir);

//----------------------------------------------------------------------
// CPP MACROS
//----------------------------------------------------------------------

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

//====================================================================
// BEGIN main ()
//====================================================================

int main(int argc, char **argv) 
{

  //====================================================================
  // PRINT HEADER
  //====================================================================

  //====================================================================
  // INITIALIZE MPI
  //====================================================================

  int mpi_size = 1, mpi_rank = 0;

#ifdef USE_MPI
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
#endif

  //====================================================================
  // PARSE COMMAND LINE
  //====================================================================

  std::string arg_reduce_axis   = "";
  std::string arg_field_name    = "";
  std::string arg_grid_file     = "";
  std::string arg_hierarchy_dir = "";
  std::string arg_slice_axis    = "";
  std::string arg_operation     = "";

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

  int c;
  int option;
  while ((c = getopt_long(argc, argv, "a:f:g:h:ps:", long_options, &option)) != -1) {
    switch (c) {
    case 0:
      // argument set in getopt_long()
      break;
    case 'a': arg_reduce_axis   = optarg; break;
    case 'f': arg_field_name    = optarg; break;
    case 'g': arg_grid_file     = optarg; break;
    case 'h': arg_hierarchy_dir = optarg; break;
    case 'p': arg_operation = "project";  break;
    case 's': arg_operation = "slice";	
      arg_slice_axis = optarg;
      break;
    case '?':	printf ("Hmmm...\n");     break;
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

  if (arg_reduce_axis == "") {

    printf ("\n\nUsage: %s [options]\n",argv[0]);
    
    printf ("\n");
    printf ("    --along (-a)     {x|y|z}             Axis along which to reduce\n");
    printf ("    --field (-f)     <field-name>        Name of the field\n");
    printf ("    --grid (-g)      <grid file>         Specify an Enzo grid file\n");
    printf ("    --hierarchy (-h) <hierarchy dir>     Specify an Enzo hierarchy directory\n");
    printf ("    --project (-p)                       Project values\n");
    printf ("    --slice (-s)     [0:1]               Slice domain at the given point\n\n");
    exit(1);

  } else {

    printf ("\n");
    printf ("  Reduction axis      = %s\n", arg_reduce_axis.c_str());
    printf ("  Field name          = %s\n", arg_field_name.c_str());
    printf ("  Grid file           = %s\n", arg_grid_file.c_str());
    printf ("  Hierarchy directory = %s\n", arg_hierarchy_dir.c_str());
    printf ("  Slice axis          = %s\n", arg_slice_axis.c_str());
    printf ("  Operation           = %s\n", arg_operation.c_str() );
    printf ("\n");

  }

  //====================================================================
  // DETERMINE LIST OF INPUT FILES
  //====================================================================

  bool      grid_specified = ! arg_grid_file.empty() ;
  bool hierarchy_specified = ! arg_hierarchy_dir.empty();
  
  
  std::vector <std::string> grid_file_list;

  if ( ( ! grid_specified && ! hierarchy_specified) ||
       (   grid_specified &&   hierarchy_specified) ) {

    ERROR("Exactly one of --grid (-g) and --hierarchy (-h) "
	  "must be specified");
    
  } else if ((grid_specified)) {

    grid_file_list.push_back(arg_grid_file);

  } else if ((hierarchy_specified)) {


    chdir (arg_hierarchy_dir.c_str());
    // ERROR CHECK HERE

    int num_grids = get_num_grids(arg_hierarchy_dir);


    // Read hierarchy file

    char file_name [MAX_FILE_STRING_LENGTH];

    enum_amr_type amr_type = get_amr_type (arg_hierarchy_dir);

    // Loop through files and add them to the list of files

    TRACE;

    for (int i = mpi_rank; true; i += mpi_size) {

      if (amr_type == amr_type_unpacked) {
	sprintf (file_name, "%s.grid%04d",arg_hierarchy_dir.c_str(),i+1);
      } else if (amr_type == amr_type_packed) {
	sprintf (file_name, "%s.cpu%04d",arg_hierarchy_dir.c_str(),i);
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

  for (unsigned int i=0; i < grid_file_list.size(); i++) {
    printf ("%s\n",grid_file_list[i].c_str());
  }

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

//--------------------------------------------------------------------
// get_num_grids ()
//--------------------------------------------------------------------

// Open the hierarchy file and determine the maximum grid number

int get_num_grids (std::string arg_hierarchy_dir)
{

  // Determine number of grids 'num_grids'

  std::string hierarchy_file = arg_hierarchy_dir + ".hierarchy";

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

enum_amr_type get_amr_type (std::string arg_hierarchy_dir)
{

  enum_amr_type amr_type;

  struct stat stFileInfo; 
  char file_name [MAX_FILE_STRING_LENGTH];

  amr_type = amr_type_unknown;

  sprintf (file_name, "%s.grid0001",arg_hierarchy_dir.c_str());
  if (stat(file_name,&stFileInfo) == 0) {
    amr_type = amr_type_unpacked;
  }

  sprintf (file_name, "%s.cpu0000",arg_hierarchy_dir.c_str());
  if (stat(file_name,&stFileInfo) == 0) {
    amr_type = amr_type_packed;
  }

  // exit if type cannot be determined

  if (amr_type == amr_type_unknown) {
    ERROR("Error in determining packed or unpacked)");
  }
  return amr_type;
}
