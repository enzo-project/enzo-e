/* #define DEBUG */
//----------------------------------------------------------------------
// 
//   Reads an Enzo data dump, and generates a *.png plot
//
//   Usage: plot-amr <color map #> <dim-1> <dim-2> <dim-3> <dataset> <file[s]> 
//   
//          where
//
//        1.  <color map #> :  Number identifying the color table
//    2,3.4.  [dim-i]       :  How to handle ith dimension
//              "x", "y"          plot on the given axis
//              "+"               sum (projection)
//              "<"               minimum
//              ">"               maximum
//        5.  <dataset>     :  Dataset name
//        6.  <file>        :  Enzo hierarchy file name
//
//   Examples:
//
//          grid-plot 0 x y + Density data0000
//
//----------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "assert.h"
#include "hdf5.h"
#include "pngwriter.h"
#include "colormap.h"

//======================================================================

const bool use_log = true;   // True if plotting log(value) not value

double min (double x, double y) { (x < y) ? x : y; }
double max (double x, double y) { (x > y) ? x : y; }

//======================================================================

#define trace printf ("TRACE %s:%d\n",__FILE__,__LINE__); fflush(stdout);

//======================================================================

void image_clear (double *image, int nx, int ny)
{  
  int ix,iy;
  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      image[ix + nx*iy] = 0;
    }
  }
}

//----------------------------------------------------------------------

void image_assemble (double      *image,     
		     int         *m, 
		     double      *dset_data, 
		     hsize_t     *n,
		     std::string operation)
{
  int i0,i1,i2;
  
  for (i2 = 0; i2 < n[2]; i2++) {
    for (i1 = 0; i1 < n[1]; i1++) {
      for (i0 = 0; i0 < n[0]; i0++) {
	int iv = i0 + n[0]*(i1 + n[1]*i2);
	int ip = i0*m[0] + i1*m[1] + i2*m[2];

	double v  = dset_data[iv];
	double &p = image[ip];

	if (use_log) {
	  if (v != 0.0) {
	    v = log(1.0+fabs(v));
	  }
	}
	   
	if (operation == "min") p = min(p,v);
	if (operation == "max") p = max(p,v);
	if (operation == "sum") p += v;
      }
    }
  }
}

//----------------------------------------------------------------------

double image_rmin (double *image, int nx, int ny)
{
  int i,ix,iy;
  double rmin = image[0];
  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      i = ix + nx*iy;
      if (rmin > image[i]) rmin = image[i];
    }
  }
  return rmin;
}

//----------------------------------------------------------------------

double image_rmax (double *image, int nx, int ny)
{
  int i,ix,iy;
  double rmax = image[0];
  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      i = ix + nx*iy;
      if (rmax < image[i]) rmax = image[i];
    }
  }
  return rmax;
}

//----------------------------------------------------------------------

void image_create (double      *image,
		   int         nx,
		   int         ny, 
		   double      *red,
		   double      *green,
		   double      *blue,
		   std::string filename)

{

  int ix,iy,i;

  pngwriter png (nx,ny,0,filename.c_str());

  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      i = ix + nx*iy;
      png.plot(ix+1,iy+1,red[i],green[i],blue[i]);
    }
  }
  png.close();
}

//----------------------------------------------------------------------

void image_init (int         &nx,
		 int         &ny,
		 int         *m,
		 hsize_t     *n,
		 std::string *dim)
{
  // x-axis values m[?] and nx
  int i;
  for (i=0; i<3; i++) {
    m[i]=0;
    if (dim[i] == "x") {
      m[i] = 1;
      nx = n[i];
    }
  }
  // Determine y-axis values m[?] and ny
  for (i=0; i<3; i++) {
    if (dim[i] == "y") {
      m[i] = nx;
      ny = n[i];
    }
  }
}

//----------------------------------------------------------------------

void generate_image_colormap (double *red,
			      double *green,
			      double *blue,
			      double *image,
			      int nx,
			      int ny,
			      double rmin,
			      double rmax,
			      int i_colormap)
{

  int ix,iy,i,index;

  for (ix = 0; ix < nx; ix++) {
    for (iy = 0; iy < ny; iy++) {
      i = ix + nx*iy;
      if (rmax > rmin) {
	index = int(((image[i] - rmin)/(rmax-rmin))*255);
      } else {
	// Uniform field--set to 127 to avoid dividing by 0
	index = 127;
      }
      if (! (0 <= index && index <= 255)) {
	printf ("OOPS!  index = %d\n",index);
	printf ("image[i] = %g\n",image[i]);
	printf ("rmin,rmax = %g,%g",rmin,rmax);
      }
      assert (0 <= index && index <= 255);
      red[i]   = colormap[i_colormap][index][0] / 255.0;
      green[i] = colormap[i_colormap][index][1] / 255.0;
      blue[i]  = colormap[i_colormap][index][2] / 255.0;
      assert (0 <= red[i]   && red[i]   <= 1);
      assert (0 <= green[i] && green[i] <= 1);
      assert (0 <= blue[i]  && blue[i]  <= 1);
    }
  }
}

//======================================================================

int main(int argc, char **argv) 
{

  int red_map[256], green_map[256], blue_map[256];

  // Check argument count

  if (argc < 7) {
    fprintf (stderr, "Usage: %s <color map #> <dim-1> <dim-2> <dim-3> <dataset> <file[s]>\n",argv[0]);
    exit(1);
  }

  // Input arguments

  int iarg=1;

  // Argument 1. Colormap number

  int i_colormap = atoi(argv[iarg++]);

  for (int line=0; line<256; line++) {
    red_map  [line] = colormap[i_colormap][line][0];
    green_map[line] = colormap[i_colormap][line][1];
    blue_map [line] = colormap[i_colormap][line][2];
  }

  // Arguments 2.3.4. Dimension operations

  std::string dim[3];
  dim[0] = argv[iarg++];
  dim[1] = argv[iarg++];
  dim[2] = argv[iarg++];

  // Argument 5. Dataset name

  std::string dataset_name (argv[iarg++]);

  // Argument(s) 6... File names

  int iarg0 = iarg;

  double rmin_all;
  double rmax_all;

  bool first_file = true;

  //----------------------------------------------------------------------

  while (iarg < argc) {

    std::string file_name    (argv[iarg++]);

    int i0,i1,i2;
   
    hid_t       file_id;

    // Open the file
 
    hid_t       dataset_id;
    hid_t       dataspace_id;
    hid_t       datatype_id;

    file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    // Open the dataset

    dataset_id = H5Dopen(file_id, ((std::string) ("/" + dataset_name)).c_str());

    // Get the dataset size

    hsize_t n[3] = {0,0,0};

    dataspace_id = H5Dget_space (dataset_id);

    int rank = H5Sget_simple_extent_ndims(dataspace_id);

    H5Sget_simple_extent_dims(dataspace_id, n, NULL);

    // Read the dataset

    herr_t      status;

    if (n[0]==0) n[0]=1;
    if (n[1]==0) n[1]=1;
    if (n[2]==0) n[2]=1;

    double       *dset_data = new double [n[0]*n[1]*n[2]];

    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
		     H5P_DEFAULT, dset_data);
    
    // Close the dataset

    status = H5Dclose(dataset_id);

    // Close the file

    status = H5Fclose(file_id);

    // Generate the plot

    int nx=0, ny=0;
    double x0,x1,y0,y1;
    int m[3];

    // Determine nx, ny, and m[] given n and dim[]

    image_init (nx,ny,m,n,dim);

    // Determine operation
    std::string operation = "";
    int i;
    for (i=0; i<3; i++) {
      if (dim[i] == "+") operation = "sum";
      if (dim[i] == "<") operation = "min";
      if (dim[i] == ">") operation = "max";
    }

    // Clear image

    double *image = new double[nx*ny];
    image_clear (image,nx,ny);

    // Assemble values

    image_assemble (image,m,dset_data,n,operation);

    // Determine rmin and rmax values for colormap scaling

    double rmin = image_rmin (image,nx,ny);
    double rmax = image_rmax (image,nx,ny);

    if (first_file) {
      first_file = false;
      rmin_all = rmin;
      rmax_all = rmax;
    } else {
      rmin_all = (rmin < rmin_all) ? rmin : rmin_all;
      rmax_all = (rmax > rmax_all) ? rmax : rmax_all;
    }

    delete [] image;
  }

  iarg = iarg0;

  while (iarg < argc) {
    
    std::string file_name    (argv[iarg++]);

    int i0,i1,i2;
   
    hid_t       file_id;

    // Open the file
 
    hid_t       dataset_id;
    hid_t       dataspace_id;
    hid_t       datatype_id;

    file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    // Open the dataset

    dataset_id = H5Dopen(file_id, ((std::string) ("/" + dataset_name)).c_str());

    // Get the dataset size

    hsize_t n[3] = {0,0,0};

    dataspace_id = H5Dget_space (dataset_id);

    int rank = H5Sget_simple_extent_ndims(dataspace_id);

    H5Sget_simple_extent_dims(dataspace_id, n, NULL);

    // Read the dataset

    herr_t      status;

    if (n[0]==0) n[0]=1;
    if (n[1]==0) n[1]=1;
    if (n[2]==0) n[2]=1;

    double       *dset_data = new double [n[0]*n[1]*n[2]];

    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
		     H5P_DEFAULT, dset_data);
    
    // Close the dataset

    status = H5Dclose(dataset_id);

    // Close the file

    status = H5Fclose(file_id);

    // Generate the plot

    int nx=0, ny=0;
    double x0,x1,y0,y1;
    int m[3];

    // Determine nx, ny, and m[] given n and dim[]

    image_init (nx,ny,m,n,dim);

    // Determine operation
    std::string operation = "";
    int i;
    for (i=0; i<3; i++) {
      if (dim[i] == "+") operation = "sum";
      if (dim[i] == "<") operation = "min";
      if (dim[i] == ">") operation = "max";
    }

    // Clear image

    double *image = new double[nx*ny];
    image_clear (image,nx,ny);

    // Assemble values

    image_assemble (image,m,dset_data,n,operation);

    // Convert values to colors

    double *red = new double[nx*ny];
    double *green = new double[nx*ny];
    double *blue = new double[nx*ny];
    int index;

    generate_image_colormap (red,green,blue,image,nx,ny,rmin_all,rmax_all,i_colormap);

    // Create the image

    std::string pngfile_name = file_name + ".png";
    image_create (image,nx,ny,red,green,blue,pngfile_name);

    delete [] dset_data;
    delete [] image;
    delete [] red;
    delete [] green;
    delete [] blue;
  }
}
