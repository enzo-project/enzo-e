//----------------------------------------------------------------------
// 
// PROGRAM
//
//   enzo-reduce
//
// DESCRIPTION
//
//   Perform a reduction operation (projection, slice, etc.) of an Enzo
//   data dump
//
// USAGE
//
//   plot-amr [options]
//
// (*)  -a  --along     {x|y|z}             Axis along which to reduce
// (*)  -f  --field     <field-name>        Name of the field
// (*)  -g  --grid      <grid file>         Specify an Enzo grid file
// (*)  -h  --hierarchy <hierarchy dir>     Specify an Enzo hierarchy directory
// ( )  -l  --lower     {x|y|z} [0:1]       Start at given axis coordinate
// (*)  -p  --project                       Project values
// (*)  -s  --slice     [0:1]               Slice domain at the given point 
// ( )  -u  --upper     {x|y|z} [0:1]       End at given axis coordinate 
//
// EXAMPLES
//
// BUGS
//
//----------------------------------------------------------------------

#include <stdio.h>
#include <assert.h>
#include <getopt.h>

#include <string>

#include "hdf5.h"

//----------------------------------------------------------------------

void read_enzo_grid (std::string file_name);

//----------------------------------------------------------------------
// COMMAND LINE PARAMETERS
//----------------------------------------------------------------------

std::string arg_reduce_axis   = "";
std::string arg_field_name    = "";
std::string arg_grid_file     = "";
std::string arg_hierarchy_dir = "";
std::string arg_slice_axis    = "";
std::string arg_operation     = "";

//----------------------------------------------------------------------

#define TRACE printf ("%s:%d TRACE\n",__FILE__,__LINE__); fflush(stdout);

//====================================================================
// BEGIN main ()
//====================================================================

int main(int argc, char **argv) 
{

  //====================================================================
  // PRINT HEADER
  //====================================================================
  //====================================================================
  // PARSE COMMAND LINE
  //====================================================================

  const struct option long_options[] =
    {
      {"along",     required_argument, 0, 'a'},
      {"field",     required_argument, 0, 'f'},
      {"grid",      required_argument, 0, 'g'},
      {"hierarchy", required_argument, 0, 'h'},
//       {"lower",     required_argument, 0, 'l'},
      {"project",   no_argument,       0, 'p'},
      {"slice",     required_argument, 0, 's'},
//       {"upper",     required_argument, 0, 'u'},
      {0, 0, 0, 0}
    };

  //--------------------------------------------------
  // Parse arguments
  //--------------------------------------------------

  int c;
  int option;
  while ((c = getopt_long(argc, argv, "a:f:g:h:ps:", long_options, &option)) != -1) {
    switch (c) 
      {
      case 0:
	// argument set in getopt_long()
	break;
      case 'a':
	arg_reduce_axis = optarg;
	break;
      case 'f':
	arg_field_name = optarg;
	break;
      case 'g':
	arg_grid_file = optarg;
	break;
      case 'h':
	arg_hierarchy_dir = optarg;
	break;
//       case 'l': // NOT IMPLEMENTED: requires multiple arguments
// 	break;
      case 'p':
	arg_operation = "project";
	break;
      case 's':
	arg_operation = "slice";
	arg_slice_axis = optarg;
	break;
//       case 'u': // NOT IMPLEMENTED: requires multiple arguments
// 	break;
      case '?':
	printf ("Hmmm...\n");
	break;
      default:
	fprintf (stderr,"Unknown command-line parse error: aborting\n");
	printf ("opterr = %d\n",opterr);
	printf ("optarg = %s\n",optarg);
	abort();
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
    printf ("    --lower (-l)     {x|y|z} [0:1]       Start at given axis coordinate\n");
    printf ("    --project (-p)                       Project values\n");
    printf ("    --slice (-s)     [0:1]               Slice domain at the given point \n");
    printf ("    --upper (-u)     {x|y|z} [0:1]       End at given axis coordinate \n\n");
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
  // INPUT ENZO HIERARCHY OR GRID
  //====================================================================

  //====================================================================
  // ACCUMULATE PROCESSOR-LOCAL REDUCTION
  //====================================================================

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



//====================================================================
// BEGIN read_enzo_grid ()
//====================================================================

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

//====================================================================
// BEGIN read_enzo_grid ()
//====================================================================
