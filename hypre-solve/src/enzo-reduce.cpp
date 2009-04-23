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
// ( )  -l  --level     <finest level>      Specify finest grid level
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

enum enum_operation_type { operation_type_unknown, 
			   operation_type_project, 
			   operation_type_slice };

struct argument_struct {
  std::string reduce_axis;
  std::string field_name;
  std::string grid_file;
  std::string hierarchy_dir;
  std::string slice_point;
  std::string operation;
};

struct grid_info_struct {
  int         index;
  std::string file_name;
  std::string group_name;
  int         level;
  double      cell_width;
  double      left_edge[3];
  double      right_edge[3];
  int         start_index[3];
  int         stop_index[3];
  int         start_index_fine[3];
  int         stop_index_fine[3];
  int         ratio_to_fine;
  int         parent;
};

struct hierarchy_info_struct {
  int     num_grids_local;
  int     max_level;
  double  fine_cell_width;
  double  left_edge[3];
  double  right_edge[3];
  int     start_index[3];
  int     stop_index[3];
};


//----------------------------------------------------------------------
// FUNCTION PROTOTYPES
//----------------------------------------------------------------------

void read_enzo_grid   
(    double          ** grid_array,
     grid_info_struct * grid_info,
     std::string        field_name  );

enum_amr_type get_amr_type 
(    std::string        arg_hierarchy_dir );

bool check_args (argument_struct *);

void print_usage(char ** argv);

void read_grids 
(    std::string        hierarchy_name, 
     grid_info_struct * grid_info,
     int                num_grids_local );

void init_hierarchy 
(    hierarchy_info_struct * hierarchy_info,
     grid_info_struct      * grid_info, 
     int                     num_grids_local );

int get_num_grids (std::string hierarchy_name);

void print_hierarchy_info(hierarchy_info_struct * hierarchy_info);

void print_grid_info(grid_info_struct * grid_info);

void write_array 
(    std::string        file_name,
     std::string        dataset_name,
     double           * array,
     int                nx,
     int                ny );

bool is_contained_in (grid_info_struct * grid_info,
		      int grid_1,
		      int grid_2);

//----------------------------------------------------------------------
// CPP MACROS
//----------------------------------------------------------------------

#define BIG_INT    2147483647
#define BIG_DOUBLE 1e30

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

#define INDEX(ix,iy,iz,nx,ny) ((ix) + (nx)*((iy) + (ny)*(iz)))

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
  // Initialize default arguments
  //--------------------------------------------------

  argument_struct arg;

  arg.reduce_axis = "x";
  arg.operation   = "project";
  arg.field_name  = "Density";

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
    case 'a': arg.reduce_axis   = optarg; break;
    case 'f': arg.field_name    = optarg; break;
    case 'g': arg.grid_file     = optarg; break;
    case 'h': arg.hierarchy_dir = optarg; break;
    case 'p': arg.operation = "project";  break;
    case 's': arg.operation = "slice";	
      arg.slice_point = optarg;
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

  if ( ! check_args(&arg) ) { 

    print_usage(argv);

  }

  //====================================================================
  // DETERMINE LIST OF INPUT FILES
  //====================================================================

  bool      grid_specified = ! arg.grid_file.empty() ;
  bool hierarchy_specified = ! arg.hierarchy_dir.empty();
  

  std::vector <std::string> grid_file_list;
  enum_amr_type amr_type = amr_type_unknown;

  if ((grid_specified)) {

    grid_file_list.push_back(arg.grid_file);

  } else if ((hierarchy_specified)) {

    chdir (arg.hierarchy_dir.c_str());
    // ERROR CHECK HERE


    // Read hierarchy file

    char file_name [MAX_FILE_STRING_LENGTH];

    amr_type = get_amr_type (arg.hierarchy_dir);

    // Loop through files and add them to the list of files

    if (amr_type == amr_type_unpacked) {

      for (int i = mpi_rank; true; i += mpi_size) {

	sprintf (file_name, "%s.grid%04d",arg.hierarchy_dir.c_str(),i+1);
	
	// Check if file exists, and append to list of files if it does

	struct stat stFileInfo; 

	if (stat(file_name,&stFileInfo) == 0) {

	  grid_file_list.push_back(file_name);

	} else {

	  break;

	}

      }

    } else if (amr_type == amr_type_packed) {

	sprintf (file_name, "%s.cpu%04d",arg.hierarchy_dir.c_str(),mpi_rank);
    }
    // ERROR CHECK ! read_successful
    
  }

  // Display list of files

  // READ GRID INFO AND DETERMINE DOMAIN EXTENTS [left|right]_edge_global

  int num_grids_local = get_num_grids (arg.hierarchy_dir);
  grid_info_struct *grid_info = new grid_info_struct[num_grids_local];
  hierarchy_info_struct hierarchy_info;

  read_grids (arg.hierarchy_dir,grid_info,num_grids_local);

  init_hierarchy (&hierarchy_info,grid_info,num_grids_local);


  // DEBUG
  print_hierarchy_info(&hierarchy_info);

  // DEBUG
  for (int i=0; i<num_grids_local; i++) {
    print_grid_info(&grid_info[i]);
  }


  //====================================================================
  // ACCUMULATE PROCESSOR-LOCAL REDUCTION
  //====================================================================

  // determine value array size

  int scale = 1;
  for (int i=0; i<hierarchy_info.max_level; i++) scale *= 2;

  int n3[3]; // Dimensions of the grid hierarchy at the finest level
  n3[0] =(hierarchy_info.stop_index[0] - hierarchy_info.start_index[0] + 1);
  n3[1] =(hierarchy_info.stop_index[1] - hierarchy_info.start_index[1] + 1);
  n3[2] =(hierarchy_info.stop_index[2] - hierarchy_info.start_index[2] + 1);

  int axis;
  int axis_char = (arg.reduce_axis.c_str())[0];
  switch (axis_char) {
  case 'x':  case 'X':  case '0':    axis = 0;    break;
  case 'y':  case 'Y':  case '1':    axis = 1;    break;
  case 'z':  case 'Z':  case '2':    axis = 2;    break;
  default:
    fprintf (stderr,"%s:%d Illegal reduce axis '%s'\n",
	     __FILE__,__LINE__,arg.reduce_axis.c_str());
    exit(1);
    break;
  }

  int ir = axis;           // reduction axis
  int ix = (axis + 1) % 3; // x-axis of array
  int iy = (axis + 2) % 3; // y-axis of array

  int nx = n3[ix];
  int ny = n3[iy];

  int n = nx*ny;

  printf ("Full size = %d %d  %d\n",nx,ny,n);
  
  // Allocate array_local[index]

  double * array_local = new double [n];

  double * grid_array = NULL;
  double * parent_array = NULL;

  if (amr_type == amr_type_packed) {
    fprintf (stderr,"%s:%d amr_type == amr_type_packed is not supported yet\n",
	     __FILE__,__LINE__);
    exit(1);
  }

  assert (amr_type == amr_type_unpacked);

  //====================================================================
  // ACCUMULATE PROCESSOR-LOCAL REDUCTION
  //====================================================================

  enum_operation_type operation_type;
  if (arg.operation == "project") operation_type = operation_type_project;
  if (arg.operation == "slice")   operation_type = operation_type_slice;

  for (int i=0; i<num_grids_local; i++) {

    read_enzo_grid (& grid_array, 
		    &grid_info[i],
		    arg.field_name);

    // And the grid's parent if it has one

    bool has_parent = grid_info[i].level > 0;

    if (has_parent) {
      read_enzo_grid (& parent_array, 
		      &grid_info[grid_info[i].parent],
		      arg.field_name);
    }

    // Loop over grid cells
    
    int    i3[3];   // current index of grid cell on grid level
    int    if3[3];  // current index of grid cell on finest level
    int    ip3[3];  // current index of parent cell on finest level

    int index_array, index_grid, index_parent;

    // first and last indices of the grid on the fine level

    int grid_index_start[3], grid_index_stop[3];

      grid_index_start[0] = grid_info[i].start_index_fine[0];
      grid_index_start[1] = grid_info[i].start_index_fine[1];
      grid_index_start[2] = grid_info[i].start_index_fine[2];
      grid_index_stop[0]  = grid_info[i].stop_index_fine[0];
      grid_index_stop[1]  = grid_info[i].stop_index_fine[1];
      grid_index_stop[2]  = grid_info[i].stop_index_fine[2];

    int ip = grid_info[i].parent;

    // first and last indices of the parent on the fine level



    int parent_index_start[3],parent_index_stop[3];

    if (has_parent) {
      parent_index_start[0] = grid_info[ip].start_index_fine[0];
      parent_index_start[1] = grid_info[ip].start_index_fine[1];
      parent_index_start[2] = grid_info[ip].start_index_fine[2];
      parent_index_stop[0]  = grid_info[ip].stop_index_fine[0];
      parent_index_stop[1]  = grid_info[ip].stop_index_fine[1];
      parent_index_stop[2]  = grid_info[ip].stop_index_fine[2];
    }

    // Determine ratio between grid level and finest level cell sizes

    int grid_subcells = 1; 
    for (int j=hierarchy_info.max_level; j>grid_info[i].level; j--) {
      grid_subcells*=2;
    }

    // Determine ratio between parent level and finest level cell sizes

    int parent_subcells = 1;
    if (has_parent) {
      for (int j=hierarchy_info.max_level; j>grid_info[ip].level; j--) {
	parent_subcells*=2;
      }
    }

    int grid_size[3] = {
      grid_info[i].stop_index[0] - grid_info[i].start_index[0] + 1,
      grid_info[i].stop_index[1] - grid_info[i].start_index[1] + 1,
      grid_info[i].stop_index[2] - grid_info[i].start_index[2] + 1
    };
    int parent_size[3] = {
      grid_info[ip].stop_index[0] - grid_info[ip].start_index[0] + 1,
      grid_info[ip].stop_index[1] - grid_info[ip].start_index[1] + 1,
      grid_info[ip].stop_index[2] - grid_info[ip].start_index[2] + 1
    };


    printf ("%d grid %d of %d   %4d:%4d %4d:%4d %4d:%4d \n",
	    mpi_rank,i,hierarchy_info.num_grids_local,
	    grid_index_start[0],grid_index_stop[0],
	    grid_index_start[1],grid_index_stop[1],
	    grid_index_start[2],grid_index_stop[2] );

    // Loop over fine grid cells on array covered by the grid

    for (if3[iy] =  grid_index_start[iy]; 
	 if3[iy] <= grid_index_stop[iy];
	 if3[iy]++) {

      i3[iy]  = (if3[iy] - grid_index_start[iy])   / grid_subcells;
      ip3[iy] = (if3[iy] - parent_index_start[iy]) / parent_subcells;

      for (if3[ix] = grid_index_start[ix]; 
	   if3[ix] <= grid_index_stop[ix]; 
	   if3[ix]++) {

	i3[ix]  = (if3[ix] - grid_index_start[ix])  / grid_subcells;
	ip3[ix] = (if3[ix] - parent_index_start[ix]) / parent_subcells;

	double value = 0;

	if ( operation_type == operation_type_project ) {

	  for (if3[ir] = grid_index_start[ir]; 
	       if3[ir] <= grid_index_stop[ir];
	       if3[ir]++) {

	    i3[ir]  = (if3[ir] - grid_index_start[ir] ) / grid_subcells;
	    ip3[ir] = (if3[ir] - parent_index_start[ir]) / parent_subcells;

	    // Add value for grid along reduction axis

	    index_grid = INDEX(i3[0],i3[1],i3[2],grid_size[0],grid_size[1]);
	    value +=  grid_array [index_grid];

	    // Subtract value for parent along reduction axis
 
	    if (has_parent) {
	      index_parent = INDEX(ip3[0],ip3[1],ip3[2],parent_size[0],parent_size[1]);
	      value -=  parent_array [index_parent];
	    }

	  }

	  // Update the array with the grid reduction value

	  index_array = if3[ix] + nx*if3[iy];
	  assert (0 <= index_array && 
		  (index_array  < nx*ny));
	  array_local[index_array] += value;

	} else if ( operation_type == operation_type_slice ) {
	  assert(0);
	}
      }
    }
  }

  //====================================================================
  // REDUCE LOCAL REDUCTION TO GLOBAL ON ROOT
  //====================================================================



  //====================================================================
  // POST-PROCESS 
  //====================================================================

  for (int i=0; i<n; i++) {
    if (array_local[i] > 0) {
      array_local[i] = log(array_local[i]);
    }	else {
      array_local[i] = 0;
    }
  };

  //====================================================================
  // OUTPUT RESULT FROM ROOT
  //====================================================================

  char file_name [MAX_FILE_STRING_LENGTH];

  sprintf (file_name,"enzo-reduce.out.%d.h5",mpi_rank);
  write_array (file_name,arg.field_name,array_local,nx,ny);

  //====================================================================
  // EXIT
  //====================================================================

  delete [] grid_info;
  delete [] array_local;
  

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
 double          ** array, 
 grid_info_struct * grid_info,
 std::string        dataset_name
 )

{

  hid_t       file_id;
  hid_t       dataset_id;
  hid_t       dataspace_id;
  herr_t      status;

  //--------------------------------------------------
  // Open the file
  //--------------------------------------------------

  file_id = H5Fopen(grid_info->file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

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
// write_array ()
//--------------------------------------------------------------------

void write_array
(
 std::string        file_name,
 std::string        dataset_name,
 double           * array,
 int                nx,
 int                ny
 )

{

  herr_t      status;

  //--------------------------------------------------
  // Create the file
  //--------------------------------------------------

  hid_t file_id = H5Fcreate ( file_name.c_str(), 
			      H5F_ACC_TRUNC, 
			      H5P_DEFAULT, 
			      H5P_DEFAULT );

  if (file_id < 0) {
    fprintf (stderr,"%s:%d %d Error creating file %s\n",
	     __FILE__,__LINE__,mpi_rank,file_name.c_str());
    exit(1);
  }

  //--------------------------------------------------
  // Create the dataspace
  //--------------------------------------------------

  int d = 2;
  hsize_t n2[2] = {nx,ny};
  hid_t space_id = H5Screate_simple(d, n2, 0);

  //--------------------------------------------------
  // Create the dataset
  //--------------------------------------------------

  std::string name = "/" + dataset_name;

  hid_t type_id = H5T_NATIVE_DOUBLE;
  
  hid_t dataset_id = H5Dcreate(file_id, name.c_str(), type_id, space_id, H5P_DEFAULT);

  if (dataset_id < 0) {
    fprintf (stderr,"%s:%d %d Error creating dataset %s in file %s\n",
	     __FILE__,__LINE__,mpi_rank,name.c_str(),file_name.c_str());
    exit(1);
  }

  //--------------------------------------------------
  // Write array to the dataset
  //--------------------------------------------------

  status = H5Dwrite (dataset_id,type_id,space_id,H5S_ALL,H5P_DEFAULT,array);

  //--------------------------------------------------
  // Close the file
  //--------------------------------------------------

  status = H5Fclose(file_id);

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
  if (mpi_rank == 0) {
    printf ("\n\nUsage: %s [options]\n",argv[0]);
    
    printf ("\n");
    printf ("   -a (--along)      {x|y|z}             Axis along which to reduce [default -a x]\n");
    printf ("   -f (--field)      <field-name>        Name of the field [default Density]\n");
    printf ("   -g (--grid)       <grid file>         Specify an Enzo grid file\n");
    printf ("   -h (--hierarchy)  <hierarchy dir>     Specify an Enzo hierarchy directory\n");
    printf ("   -p (--project)                        Project values [default]\n");
    printf ("   -s (--slice)      [0:1]               Slice domain at the given point\n\n");
    printf ("\n");
    printf ("Required options: af[g|h][p|s]\n");
    printf ("\n");
  }
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  exit (1);

}


//-----------------------------------------------------------------------
// Read hierarchy file 
//-----------------------------------------------------------------------

void read_grids (std::string        hierarchy_name, 
		 grid_info_struct * grid_info,
		 int                num_grids_local)
{
  std::string hierarchy_file = hierarchy_name + ".hierarchy";

  FILE *fp_hierarchy  = fopen (hierarchy_file.c_str(), "r");

  char line [ MAX_FILE_STRING_LENGTH ];

  // Read input Enzo hierarchy file

  int grid_index = 1;

  int index = 0;

  while (fgets(line,MAX_FILE_STRING_LENGTH, fp_hierarchy) != NULL) {

    // Get grid number

    sscanf (line,"Grid = %d", &grid_index);

    if (INDEX_TO_PROC(grid_index-1) == mpi_rank) {

      // Update grid array index

      index = (grid_index - 1) / mpi_size;

      // Get grid index

      grid_info[index].index = grid_index;

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

  // Determine grid levels using mesh widths
  // ASSUMES AT LEAST ONE GRID AND THAT IT IS A ROOT-LEVEL GRID

  double root_cell_width = 
    (grid_info[0].right_edge[0] - grid_info[0].left_edge[0]) / 
    (grid_info[0].stop_index[0] - grid_info[0].start_index[0] + 1) ;

  char file_name[MAX_FILE_STRING_LENGTH];

  for (int i = 0; i < num_grids_local; i++) {

    double local_cell_width = 
      (grid_info[i].right_edge[0] - grid_info[i].left_edge[0]) /
      (grid_info[i].stop_index[0] - grid_info[i].start_index[0] + 1);

    // Determine cell width

    grid_info[i].cell_width = local_cell_width;

    // Determine level (0.5 required since conversion from float to int truncates)

    grid_info[i].level = 0.5+log(root_cell_width/grid_info[i].cell_width)/log(2.0);

    // Determine filenme

    sprintf (file_name,"%s.grid%04d",hierarchy_name.c_str(), i*mpi_size+mpi_rank+1);
    grid_info[i].file_name = file_name;

  }

}

//-----------------------------------------------------------------------

void init_hierarchy (hierarchy_info_struct * hierarchy_info,
		     grid_info_struct      * grid_info, 
		     int                     num_grids_local)
{
  double left_edge_local[3]  = {BIG_DOUBLE,BIG_DOUBLE,BIG_DOUBLE};
  double right_edge_local[3] = {-BIG_DOUBLE,-BIG_DOUBLE,-BIG_DOUBLE};

  hierarchy_info->num_grids_local = num_grids_local;

  for (int index=0; index<num_grids_local; index++) {
    
    // Update left_edge_local[] and right_edge_local[]

    for (int i=0; i<3; i++) {
      left_edge_local[i]  = MIN(left_edge_local[i], grid_info[index].left_edge[i]);
      right_edge_local[i] = MAX(right_edge_local[i],grid_info[index].right_edge[i]);
    }

  }

  // Ensure that it works without USE_MPI, or with USE_MPI but not called with mpirun

  for (int i=0; i<3; i++) {
    hierarchy_info->left_edge[i] = left_edge_local[i];
    hierarchy_info->right_edge[i] = right_edge_local[i];
  }

#ifdef USE_MPI  
  if (mpi_size > 1) {
    MPI_Allreduce (left_edge_local,hierarchy_info->left_edge,3,
		   MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce (right_edge_local,hierarchy_info->right_edge,3,
		   MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  }
#endif

  // Determine maximum level in hierarchy

  hierarchy_info->max_level = 0;

  for (int i = 0; i < num_grids_local; i++) {

    hierarchy_info->max_level = MAX(hierarchy_info->max_level,grid_info[i].level);

  }

  // Determine finest mesh width (assumes first grid is coarse)

  double fine_to_coarse_ratio = 1;
  for (int i=0; i<hierarchy_info->max_level; i++) fine_to_coarse_ratio *= 0.5;
  hierarchy_info->fine_cell_width = fine_to_coarse_ratio * grid_info[0].cell_width;

  // grid.ratio_to_fine
  // grid.start_index_fine
  // grid.stop_index_fine

  for (int i=0; i<num_grids_local; i++) {

    // Determine ratio_to_fine grid cell ratio between level and finest

    grid_info[i].ratio_to_fine = 1;
    for (int j=grid_info[i].level; j<hierarchy_info->max_level; j++) {
      grid_info[i].ratio_to_fine *= 2;
    }

    // Compute start_index_fine and stop_index_fine

    for (int j=0; j<3; j++) {
    grid_info[i].start_index_fine[j] 
      = ((grid_info[i].left_edge[j] - hierarchy_info->left_edge[j]) 
	    / hierarchy_info->fine_cell_width);
    grid_info[i].stop_index_fine[j] 
      = ((grid_info[i].right_edge[j] - hierarchy_info->left_edge[j]) 
	    / hierarchy_info->fine_cell_width - 1);
    }
  }

  // hierarchy.start_index
  // hierarchy.stop_index

  int    start_index_local[3] = {BIG_INT,BIG_INT,BIG_INT};
  int    stop_index_local[3] = {-BIG_INT,-BIG_INT,-BIG_INT};

  for (int index=0; index<num_grids_local; index++) {
    for (int i=0; i<3; i++) {
      start_index_local[i] = 
	MIN(start_index_local[i],grid_info[index].start_index_fine[i]);
      stop_index_local[i] = 
	MAX(stop_index_local[i],grid_info[index].stop_index_fine[i]);
    }
  }

  for (int i=0; i<3; i++) {
    hierarchy_info->start_index[i] = start_index_local[i];
    hierarchy_info->stop_index[i] = stop_index_local[i];
  }

#ifdef USE_MPI  
  if (mpi_size > 1) {
    MPI_Allreduce (start_index_local,hierarchy_info->start_index,3,
		   MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce (stop_index_local,hierarchy_info->stop_index,3,
		   MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  }
#endif

//   // determine parents 

  if (mpi_size != 1) {
    fprintf (stderr,"%s:%d Parent determination is not parallelized\n",
	     __FILE__,__LINE__);
    exit(1);
  }

  for (int grid_1=0; grid_1 < num_grids_local; grid_1++) {
    // find parent of grid_1 
    grid_info[grid_1].parent = -1;
    int level_parent = -1;
    for (int grid_2=0; grid_2 < num_grids_local; grid_2++) {
      if (grid_info[grid_1].level >= grid_info[grid_2].level + 1 &&
	  is_contained_in (grid_info,grid_1,grid_2) ) {
	if (grid_info[grid_2].level > level_parent) {
	  level_parent = grid_info[grid_2].level;
	  grid_info[grid_1].parent = grid_2;
	  printf ("grid %d parent %d\n",grid_1+1,grid_2+1);
	}
      }
    }
  }

  bool err = false;
  for (int grid_1=0; grid_1 < num_grids_local; grid_1++) {
    if ((grid_info[grid_1].level == 0 && grid_info[grid_1].parent != -1) ||
	(grid_info[grid_1].level != 0 && grid_info[grid_1].parent == -1) ) {
      err = true;
      printf ("Grid %d parent is incorrect\n",grid_1+1);
    }
  }

  if (err) exit(1);
  

}

//-----------------------------------------------------------------------

int get_num_grids (std::string hierarchy_name)
{
  std::string hierarchy_file = hierarchy_name + ".hierarchy";
  FILE *fp_hierarchy  = fopen (hierarchy_file.c_str(), "r");
  
  char line [ MAX_FILE_STRING_LENGTH ];

  // Read input Enzo hierarchy file

  int grid_index = 1;

  int num_local_grids = 0;

  while (fgets(line,MAX_FILE_STRING_LENGTH, fp_hierarchy) != NULL) {
    sscanf (line,"Grid = %d", &grid_index);
    if (INDEX_TO_PROC(grid_index-1) == mpi_rank) {
      num_local_grids = MAX(num_local_grids,(grid_index-1)/mpi_size+1);
    }
  }

  fclose(fp_hierarchy);   // Close hierarchy file

  return num_local_grids;
}

//-----------------------------------------------------------------------

void print_hierarchy_info (hierarchy_info_struct * hierarchy_info)
{
  printf ("%d Hierarchy info\n",mpi_rank);
  printf ("%d  num_grids_local  = %d\n",
	  mpi_rank,hierarchy_info->num_grids_local);
  printf ("%d  max_level        = %d\n",
	  mpi_rank,hierarchy_info->max_level);
  printf ("%d  fine_cell_width  = %g\n",
	  mpi_rank,hierarchy_info->fine_cell_width);
  printf ("%d  left_edge        = %g %g %g\n",
	  mpi_rank,
	  hierarchy_info->left_edge[0],
	  hierarchy_info->left_edge[1],
	  hierarchy_info->left_edge[2]);
  printf ("%d  right_edge       = %g %g %g\n",
	  mpi_rank,
	  hierarchy_info->right_edge[0],
	  hierarchy_info->right_edge[1],
	  hierarchy_info->right_edge[2]);
  printf ("%d  start_index      = %d %d %d\n",
	  mpi_rank,
	  hierarchy_info->start_index[0],
	  hierarchy_info->start_index[1],
	  hierarchy_info->start_index[2]);
  printf ("%d  stop_index       = %d %d %d\n",
	  mpi_rank,
	  hierarchy_info->stop_index[0],
	  hierarchy_info->stop_index[1],
	  hierarchy_info->stop_index[2]);

  printf ("\n");
}

//-----------------------------------------------------------------------

void print_grid_info (grid_info_struct * grid_info)
{
  printf ("%d  index         = %d\n",mpi_rank, grid_info->index);
  printf ("%d  file_name     = %s\n",mpi_rank, grid_info->file_name.c_str());
  printf ("%d  group_name    = %s\n",mpi_rank, grid_info->group_name.c_str());
  printf ("%d  level         = %d\n",mpi_rank, grid_info->level);
  printf ("%d  cell_width    = %g\n",mpi_rank, grid_info->cell_width);
  printf ("%d  left_edge     = %g %g %g\n",mpi_rank, grid_info->left_edge[0],grid_info->left_edge[1],grid_info->left_edge[2]);
  printf ("%d  right_edge    = %g %g %g\n",mpi_rank, grid_info->right_edge[0],grid_info->right_edge[1],grid_info->right_edge[2]);
  printf ("%d  start_index   = %d %d %d\n",mpi_rank, grid_info->start_index[0],grid_info->start_index[1],grid_info->start_index[2]);
  printf ("%d  stop_index    = %d %d %d\n",mpi_rank, grid_info->stop_index[0],grid_info->stop_index[1],grid_info->stop_index[2]);
  printf ("%d  start_index_fine = %d %d %d\n",
	  mpi_rank,
	  grid_info->start_index_fine[0],
	  grid_info->start_index_fine[1],
	  grid_info->start_index_fine[2]);
  printf ("%d  stop_index_fine = %d %d %d\n",
	  mpi_rank,
	  grid_info->stop_index_fine[0],
	  grid_info->stop_index_fine[1],
	  grid_info->stop_index_fine[2]);
  printf ("%d  ratio_to_fine = %d\n",	  mpi_rank,  grid_info->ratio_to_fine);
  printf ("%d  parent        = %d\n",	  mpi_rank,  grid_info->parent);
	  
  printf ("\n");
}

bool is_contained_in (grid_info_struct * grid_info,
		      int grid_1,
		      int grid_2)
// Return true if grid 1 is contained in grid 2
{
  return (grid_info[grid_2].left_edge[0]  <= grid_info[grid_1].left_edge[0] &&
	  grid_info[grid_2].left_edge[1]  <= grid_info[grid_1].left_edge[1] &&
	  grid_info[grid_2].left_edge[2]  <= grid_info[grid_1].left_edge[2] &&
	  grid_info[grid_1].right_edge[0] <= grid_info[grid_2].right_edge[0] &&
	  grid_info[grid_1].right_edge[1] <= grid_info[grid_2].right_edge[1] &&
	  grid_info[grid_1].right_edge[2] <= grid_info[grid_2].right_edge[2]);
      
}
