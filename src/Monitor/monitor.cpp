//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      monitor.cpp
 * @brief     Routines for simple output of text, plots, and graphs
 * @author    James Bordner
 * @date      Thu Feb 21 12:45:56 PST 2009
 * @bug       
 * @note      
 *
 * DESCRIPTION 
 * 
 *    Routines for simple output of text, plots, and graphs
 *
 * PACKAGES
 *
 *    none
 * 
 * INCLUDES
 *  
 *    pngwriter.h
 *
 * PUBLIC FUNCTIONS
 *  
 *    
 *
 * PRIVATE FUCTIONS
 *  
 *    
 *
 * $Id$
 *
 *********************************************************************
 */

#include "cello.h"

#include "pngwriter.h"

#include "monitor.hpp" 

void Monitor::plot_png 
(std::string name, 
 Scalar * array, 
 int nx, int ny,
 Scalar min, Scalar max, 
 const int * map, 
 int map_length)
/**
 *********************************************************************
 *
 * @param  name       File name
 * @param  array      Array of values to plot
 * @param  nx,ny      Size of the array
 * @param  min,max    Bounds for color map
 * @param  map        Color map [r0, g0, b0, r1, g1, b1, ...]
 * @param  map_length Length of color map / 3
 *
 * Plot an array as a png file
 *
 *********************************************************************
 */
{
  pngwriter png (nx,ny,0,name.c_str());
  double h = (max-min) / (map_length-1);
  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      int i = ix + nx*iy;
      // find lower index (0:map_length-2) into color map
      int index = (array[i] - min) / h;
      if (index > map_length - 2) index = map_length-2;
      double ratio = (array[i] - min - index*h) / h;
      if (index < 0 || index > map_length -2) {
	printf ("%s:%d out of range index = %d\n",__FILE__,__LINE__,index);
	printf ("min=%g max=%g value=%g map_length=%d\n",
		min,max,array[i],map_length);
	exit(1);
      }
      if (ratio < 0 || ratio > 1) {
	printf ("%s:%d out of range ratio = %g\n",__FILE__,__LINE__,ratio);
	printf ("min=%g max=%g value=%g map_length=%d\n",
		min,max,array[i],map_length);
	exit(1);
      }
      double r = (1-r)*map[3*index+0] + r*map[3*(index+1)+0];
      double g = (1-r)*map[3*index+1] + r*map[3*(index+1)+1];
      double b = (1-r)*map[3*index+2] + r*map[3*(index+1)+2];
      png.plot(ix+1,iy+1,r,g,b);
      
    }
  }
  png.close();
}

