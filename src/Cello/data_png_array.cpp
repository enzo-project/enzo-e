#include <stdio.h>
#include <limits>
#include "pngwriter.h"

// #define DEBUG_PNG

// #define PLOT
void png_array (const char * filename,
		float * array,
		int gx,int gy,int gz,
		int mx,int my,int mz,
		const char * file, int line,
		int axis = 2,
		int px=0, int py=0,
		double scale = 1.0
		)
{
  //  double colormap[4][3] = {{0.0, 0.0, 0.0},
  //			   {1.0, 0.0, 0.0},
  //			   {1.0, 1.0, 0.0},
  //			   {1.0, 1.0, 1.0}};
  // const int nc = 4;

  double colormap[2][3] = {{ 0.0, 0.0, 0.0},
			   { 1.0, 1.0, 1.0} };
  const int nc = 2;

#ifdef PLOT  
  pngwriter * png;
  if (px==0) px = 1;
  if (py==0) py = 1;
  px = (mx-2*gx)*px;
  py = (my-2*gy)*py;
  png = new pngwriter(px, py,0,filename);
#endif
  
  double min=std::numeric_limits<double>::max();
  double max=-std::numeric_limits<double>::max();

  double sum_real=0.0;
  double sum_abs=0.0;
  double sum_ghost=0.0;
  double sum2_real=0.0;
  double sum2_ghost=0.0;
  for (int iy=0; iy<my; iy++) {
    bool ly = (gy <= iy && iy < my-gy);
    for (int ix=0; ix<mx; ix++) {
      bool lx = (gx <= ix && ix < mx-gx);
      double sum = 0.0;
      for (int iz=0; iz<mz; iz++) {
	bool lz = (gz <= iz && iz < mz-gz);
	int i=ix + mx*(iy + my*iz);
	if (lx&&ly&&lz) sum += scale*array[i];
	if (lx&&ly&&lz) sum_real+=array[i];
	if (lx&&ly&&lz) sum_abs+=fabs(array[i]);
	if (lx&&ly&&lz) sum2_real+=array[i]*array[i];
	sum_ghost+=array[i];
	sum2_ghost+=array[i]*array[i];
      }
      if (lx&&ly) {
	min = std::min(sum,min);
	max = std::max(sum,max);
      }
    }
  }
#ifdef DEBUG_PNG  
  printf ("DEBUG_PNG FILE %s: %s:%d\n",filename,file,line);
  printf ("DEBUG_PNG %s: min %20.15g max %20.15g\n",filename,min,max);
  printf ("DEBUG_PNG %s: 0xyz %20.15g : %20.15g %20.15g %20.15g\n",filename,
	  array[(gx)+mx*((gy) + my*(gz))],
	  array[(gx+1)+mx*((gy) + my*(gz))],
	  array[(gx)+mx*((gy+1) + my*(gz))],
	  array[(gx)+mx*((gy) + my*(gz+1))]);
  printf ("DEBUG_PNG %s: sum_real %20.15g sum_ghost %20.15g\n",filename,sum_real,sum_ghost);
  printf ("DEBUG_PNG %s: sum_abs %20.15g\n",filename,sum_abs);
  printf ("DEBUG_PNG %s: sum2_real %20.15g sum2_ghost %20.15g\n",filename,sum2_real,sum2_ghost);
#endif
  int nx=mx-2*gx;
  int ny=my-2*gy;
#ifdef PLOT  
  for (int ky=0; ky<py; ky++) {
    int iy = 1.0*ky*ny/py+gy;
    for (int kx=0; kx<px; kx++) {
      int ix = 1.0*kx*nx/px+gx;
      double sum = 0.0;
      for (int iz=gz; iz<mz-gz; iz++) {
	int i=ix + mx*(iy + my*iz);
	sum += scale*array[i];
      }
      size_t k = (nc-1)*(sum-min) / (max-min);
      if (k > nc-2) k = nc-2;
      double lo = min +  k   *(max-min)/(nc-1);
      double hi = min + (k+1)*(max-min)/(nc-1);
      double ratio = (sum - lo) / (hi-lo);
      double r = (1-ratio)*colormap[k][0] + ratio*colormap[k+1][0];
      double g = (1-ratio)*colormap[k][1] + ratio*colormap[k+1][1];
      double b = (1-ratio)*colormap[k][2] + ratio*colormap[k+1][2];
      png->plot(kx+1,ky+1,r,g,b);
    }
  }
  png->close();
  
  delete png;
#endif
}
