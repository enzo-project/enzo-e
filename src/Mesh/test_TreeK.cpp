// See LICENSE_CELLO file for license and copyright information

/// @file     test_TreeK.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-28
/// @brief    Test program for Tree2K and Tree3K classes

#include "mesh_tree.hpp"
#include "test.hpp"

#define index(ix,iy,iz,n) ((ix) + (n)*((iy) + (n)*(iz)))

const bool debug    = false;
const bool geomview = false;

const int  cell_size = 1;
const int  line_width = 1;
const int  gray_threshold = 127;


//----------------------------------------------------------------------

int * create_image_array (const char * filename, int * nx, int * ny, int max_levels);
//int * create_level_array3 (int * n3, int max_levels);
int * create_sphere_array (int * n, int max_levels);
void write_image(std::string filename, float * image, 
		 int nx, int ny, int nz=0 );

void create_tree (int * level_array, int nx, int ny, int nz, int k,  int d, 
		   std::string name, int max_level);
void print_usage(int, char**);
//----------------------------------------------------------------------


#ifdef CONFIG_USE_CHARM
#   include "main.decl.h"
#endif

PARALLEL_MAIN_BEGIN

{

  PARALLEL_INIT;

  unit_init();

  // Parse command line

  if (PARALLEL_ARGC != 5 && PARALLEL_ARGC != 4) {
    print_usage(PARALLEL_ARGC,PARALLEL_ARGV);
  }

  // Check arguments

  int dimension  = atoi(PARALLEL_ARGV[1]);
  int refinement = atoi(PARALLEL_ARGV[2]);
  int max_level  = atoi(PARALLEL_ARGV[3]);
  const char * png_file = 0;
  if (PARALLEL_ARGC == 5) png_file = PARALLEL_ARGV[4];

  if (dimension != 2 && 
      dimension != 3) print_usage(PARALLEL_ARGC,PARALLEL_ARGV);
  
  if (refinement != 2 && 
      refinement != 4 &&
      refinement != 8 &&
      refinement != 16) print_usage(PARALLEL_ARGC,PARALLEL_ARGV);
  
  if (! (0 < max_level && max_level <= 14)) print_usage(PARALLEL_ARGC,PARALLEL_ARGV);

  char filename[80];
  sprintf (filename,"TreeK-D=%d-R=%d-L=%d",dimension,refinement,max_level);

  int nx,ny,nz;
  int * level_array;

  if (dimension == 2) {

    nx = 1;
    ny = 1;
    nz = 1;
    // Read png image file
    level_array = create_image_array(png_file,&nx,&ny,max_level);
    // Change back to parent directory

  } else {

    level_array = create_sphere_array(&nx,max_level);
    ny = nx;
    nz = nx;
  }

  create_tree (level_array, nx, ny, nz, refinement, dimension, filename,max_level);

  delete [] level_array;

  Memory * memory = Memory::instance();
  memory->print();

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

//----------------------------------------------------------------------
void print_usage(int argc, char **argv)
{
  fprintf (stderr,"\n");
  fprintf (stderr,"Usage: %s <dimension> <refinement> <levels> [file]\n",argv[0]);
  fprintf (stderr,"\n");
  fprintf (stderr,"   where \n");
  fprintf (stderr,"\n");
  fprintf (stderr,"         <dimension>  = [2|3]\n");
  fprintf (stderr,"         <refinement> = [2|4|8|16]\n");
  fprintf (stderr,"         [file]       = filename if dimension is 2\n");
  fprintf (stderr,"\n");
  exit(1);
}
//----------------------------------------------------------------------

// read in the gimp-generated image data into a level array
// values are set to [0:max_levels)

int * create_image_array (const char * pngfile, 
			  int * nx, int * ny, int max_levels)
{

  pngwriter png;

  png.readfromfile(pngfile);

  int width  = png.getwidth();
  int height = png.getheight();

  int size = (width > height) ? width: height;

  int * level_array = new int [size*size];

  for (int i=0; i<size*size; i++) level_array[i] = 0;

  *nx = size;
  *ny = size;

  int ix0 = (size - width) / 2;
  int iy0 = (size - height) / 2;

  int max = 0;
  for (int iy=0; iy<height; iy++) {
    for (int ix=0; ix<width; ix++) {
      int pixel = png.read(ix+1,iy+1);
      int i = (ix+ix0) + size*(iy+iy0);
      level_array[i] = max_levels * 1.0*pixel / 256;
      if (level_array[i] > max) max = level_array[i];
    }
  }
  png.close();
  return level_array;
}

//----------------------------------------------------------------------

int * create_sphere_array (int * n3, int max_levels)
{

  *n3 = int(exp(max_levels * log (2.0)) + 0.5) * 2;
  int n = *n3;

  PARALLEL_PRINTF ("%d %d\n",max_levels,n);
  int * level_array = new int [n*n*n];

  for (int i=0; i<n*n*n; i++) level_array[i] = 0;

  const double R = 0.3;  // radius
  double R2 = R*R;
  
  double x,y,z;

  for (int iz=0; iz<n/2; iz++) {
    z = double(iz) / n - 0.5;
    for (int iy=0; iy<n/2; iy++) {
      y = double(iy) / n - 0.5;
      for (int ix=0; ix<n/2; ix++) {
	x = double(ix) / n - 0.5;
	double r2 = x*x + y*y + z*z;
	double v = r2 < R2 ? max_levels : 0;

	level_array[index(     ix,     iy,     iz,n)] = v;
	level_array[index(n-ix-1,     iy,     iz,n)] = v;
	level_array[index(     ix,n-iy-1,     iz,n)] = v;
	level_array[index(n-ix-1,n-iy-1,     iz,n)] = v;
	level_array[index(     ix,     iy,n-iz-1,n)] = v;
	level_array[index(n-ix-1,     iy,n-iz-1,n)] = v;
	level_array[index(     ix,n-iy-1,n-iz-1,n)] = v;
	level_array[index(n-ix-1,n-iy-1,n-iz-1,n)] = v;
      }
    }
  }
  return level_array;
}

//----------------------------------------------------------------------

void write_image(std::string filename, float * image, int nx, int ny, int nz)
{
  if (nx > 8194 || ny > 8194 || nz > 8194) {
    PARALLEL_PRINTF ("%s:%d (nx,ny,nz) = (%d,%d,%d)\n",__FILE__,__LINE__,nx,ny,nz);
    exit(1);
  }
  // Write HDF5 file

  // Hdf5 hdf5;
  // hdf5.file_open((filename+".hdf5").c_str(),"w");
  // PARALLEL_PRINTF ("write_image %d %d %d\n",nx,ny,nz);
  // hdf5.dataset_open_write ("tree_image",nx,ny,nz);
  // hdf5.write(image);
  // hdf5.dataset_close ();
  // hdf5.file_close();

  // Write PNG image

  float min=image[0];
  float max=image[0];
  for (int i=0; i<nx*ny*nz; i++) {
    min = MIN(min,image[i]);
    max = MAX(max,image[i]);
  }

  Monitor * monitor = Monitor::instance();

  monitor->image ((filename+".png").c_str(),
		  nx,ny,
		  image,
		  nx,ny,1,
		  nx,ny,1,
		  0,0,0,
		  axis_z,reduce_sum,
		  min,max);

}

//----------------------------------------------------------------------

void create_tree 
(
 int * level_array, 
 int nx, int ny, int nz,
 int k,  int d, 
 std::string name,
 int max_level
 )
{

  TreeK * tree = 0;
  if (d==2) tree = new Tree2K(k);
  if (d==3) tree = new Tree3K(k);

  float mem_per_node;

  Memory * memory = Memory::instance();
  memory->reset();

  PARALLEL_PRINTF ("--------------------------------------------------\n");
  PARALLEL_PRINTF ("k=%d d=%d\n",k,d);
  PARALLEL_PRINTF ("--------------------------------------------------\n");

  int full_nodes = true;

  //--------------------------------------------------
  // Refine the tree
  //--------------------------------------------------

  PARALLEL_PRINTF ("\nINITIAL TREE\n");

  memory->set_active(true);
  tree->refine(level_array,nx,ny,nz,max_level,full_nodes);
  memory->print();
  memory->set_active(false);

  mem_per_node = (float) memory->bytes(0) / tree->num_nodes();
  PARALLEL_PRINTF ("nodes      = %d\n",tree->num_nodes());
  PARALLEL_PRINTF ("levels     = %d\n",tree->levels());
  PARALLEL_PRINTF ("bytes/node = %g\n",mem_per_node);

  //--------------------------------------------------
  // Balance the tree
  //--------------------------------------------------

  // Determine image size
  float * image;
  int image_size = cell_size + 2*line_width;
  for (int i=0; i<tree->levels(); i++) {
    image_size = 2*image_size - line_width;
  }

  PARALLEL_PRINTF ("\nBALANCED TREE\n");

  memory->set_active(true);
  tree->balance(full_nodes);
  memory->print();
  memory->set_active(false);

  if (d==2) {
    image = tree->create_image(image_size,line_width);
    write_image(name + "-0",image,image_size,image_size,1);
    delete [] image;
  } else {
    image = tree->create_image(image_size,line_width,0);
    write_image(name + "-x" + "-0",image,image_size,image_size,1);
    delete [] image;
    image = tree->create_image(image_size,line_width,1);
    write_image(name + "-y" + "-0",image,image_size,image_size,1);
    delete [] image;
    image = tree->create_image(image_size,line_width,2);
    write_image(name + "-z" + "-0",image,image_size,image_size,1);
    delete [] image;
  }

  mem_per_node = (float) memory->bytes(0) / tree->num_nodes();
  PARALLEL_PRINTF ("nodes      = %d\n",tree->num_nodes());
  PARALLEL_PRINTF ("levels     = %d\n",tree->levels());
  PARALLEL_PRINTF ("bytes/node = %g\n",mem_per_node);

  //--------------------------------------------------
  // Coalesce Patches In the tree
  //--------------------------------------------------

  PARALLEL_PRINTF ("\nCOALESCED TREE\n");

  memory->set_active(true);
  tree->optimize();
  memory->print();
  memory->set_active(false);

  if (d==2) {
    image = tree->create_image(image_size,line_width);
    write_image(name + "-1",image,image_size,image_size,1);
    delete [] image;
  } else {
    image = tree->create_image(image_size,line_width,0);
    write_image(name + "-x" + "-1",image,image_size,image_size,1);
    delete [] image;
    image = tree->create_image(image_size,line_width,1);
    write_image(name + "-y" + "-1",image,image_size,image_size,1);
    delete [] image;
    image = tree->create_image(image_size,line_width,2);
    write_image(name + "-z" + "-1",image,image_size,image_size,1);
    delete [] image;
  }

  if (geomview) tree->geomview(name + ".gv");

  PARALLEL_PRINTF ("nodes      = %d\n",tree->num_nodes());
  PARALLEL_PRINTF ("levels     = %d\n",tree->levels());
  PARALLEL_PRINTF ("bytes/node = %g\n",mem_per_node);

  delete tree;

}

#ifdef CONFIG_USE_CHARM
#   include "main.def.h"
#endif
