// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_Soup.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-09-10
/// @brief    Definition of the Soup class for "alphabet soup" test problem

#include "enzo.hpp"

const int EnzoInitialSoup::position_[] = 
  {
      0,1,2,4,5,6,7,9,10,12,13,14,15,16,17,19,20,21,22,24,25,27,28,29,30,31
  };

//======================================================================

EnzoInitialSoup::EnzoInitialSoup
(int cycle,
 double time,
 std::string filename,
 int rank,
 bool rotate,
 int    nax,    int nay,    int naz,
 double dpx, double dpy, double dpz,
 double dsx, double dsy, double dsz,
 double density,
 double pressure_in, double pressure_out) throw ()
  : Initial(cycle,time),
    file_name_(filename),
    rank_(rank),
    rotate_(rotate),
    png_(NULL),
    density_(density),
    pressure_in_(pressure_in),
    pressure_out_(pressure_out)
{
  array_[0] = nax;
  array_[1] = nay;
  array_[2] = naz;
  d_pos_[0] = dpx;
  d_pos_[1] = dpy;
  d_pos_[2] = dpz;
  d_size_[0] = dsx;
  d_size_[1] = dsy;
  d_size_[2] = dsz;

  png_ = new pngwriter;
    
  png_->readfromfile(file_name_.c_str());

  int nx = png_->getwidth();
  int ny = png_->getheight();

  ASSERT3 ("EnzoInitialSoup::EnzoInitialSoup()",
	   "Error reading input PNG file %s: size must be %d x %d",
	   file_name_.c_str(),SOUP_IMAGE_NX,SOUP_IMAGE_NY,
	   (nx == SOUP_IMAGE_NX) &&
	   (ny == SOUP_IMAGE_NY));

  // seed random number generator to ensure all processors generate the
  // same letter_ array
  srand(31415);
  
  letter_ = new char [nax*nay*naz];
  for (int kz=0; kz<naz; kz++) {
    for (int ky=0; ky<nay; ky++) {
      for (int kx=0; kx<nax; kx++) {
	int k = kx + nax*(ky + nay*kz);
	letter_[k] = 'A' + 26*(rand()/(RAND_MAX+1.0));
      }
    }
  }
   
}


//----------------------------------------------------------------------

EnzoInitialSoup::~EnzoInitialSoup() throw ()
{
  delete png_;
  png_ = NULL;
  delete [] letter_;
  letter_ = NULL;
}

//----------------------------------------------------------------------

void EnzoInitialSoup::enforce_block
( Block            * block, 
  const FieldDescr * field_descr,
  const ParticleDescr * particle_descr,
  const Hierarchy  * hierarchy
  ) throw()
{
  Field field = block->data()->field();

  enzo_float * d  = (enzo_float *) field.values("density");
  enzo_float * te = (enzo_float *) field.values ("total_energy");
  enzo_float * vx = (enzo_float *) field.values("velocity_x");
  enzo_float * vy = (enzo_float *) field.values("velocity_y");
  enzo_float * vz = (enzo_float *) field.values("velocity_z");

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  double xmb,ymb,zmb;
  double xpb,ypb,zpb;
  block->data()->lower(&xmb,&ymb,&zmb);
  block->data()->upper(&xpb,&ypb,&zpb);

  double hx,hy,hz;
  field.cell_width(xmb,xpb,&hx,
		   ymb,ypb,&hy,
		   zmb,zpb,&hz);

  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  int mx = nx + 2*gx;
  int my = ny + 2*gy;
  int mz = nz + 2*gz;

  // array of explosions

  double xmd,ymd,zmd;
  hierarchy->lower(&xmd,&ymd,&zmd);
  double xpd,ypd,zpd;
  hierarchy->upper(&xpd,&ypd,&zpd);

  double hxa = (xpd-xmd) / array_[0];
  double hya = (ypd-ymd) / array_[1];
  double hza = (zpd-zmd) / array_[2];

  // bounds of possible explosions intersecting this Block
  
  int kxm = MAX((int)floor(xmb/(xpd-xmd)*array_[0])-1,0);
  int kym = MAX((int)floor(ymb/(ypd-ymd)*array_[1])-1,0);
  int kzm = MAX((int)floor(zmb/(zpd-zmd)*array_[2])-1,0);
  int kxp = MIN((int) ceil(xpb/(xpd-xmd)*array_[0])+1,array_[0]);
  int kyp = MIN((int) ceil(ypb/(ypd-ymd)*array_[1])+1,array_[1]);
  int kzp = MIN((int) ceil(zpb/(zpd-zmd)*array_[2])+1,array_[2]);

  for (int kz=kzm; kz<kzp; kz++) {
    double cz = hza*(0.5+kz);
    for (int ky=kym; ky<kyp; ky++) {
      double cy = hya*(0.5+ky);
      for (int kx=kxm; kx<kxp; kx++) {
	const double cx = hxa*(0.5+kx);
	const int k = kx + array_[0]*(ky + array_[1]*kz);
	char letter = letter_[k];
	int lx=position_[letter-'A'] % 8;
	int ly=position_[letter-'A'] / 8;
	CkPrintf ("%d %d %d %c  %d %d\n",kx,ky,kz,letter,lx,ly);
	fflush(stdout);
	// (explosion center cx,cy,cz)

	for (int iz=0; iz<mz; iz++) {
	  double z = zmb + (iz - gz + 0.5)*hz - cz;
	  for (int iy=0; iy<my; iy++) {
	    double y = ymb + (iy - gy + 0.5)*hy - cy;
	    for (int ix=0; ix<mx; ix++) {
	      double x = xmb + (ix - gx + 0.5)*hx - cx;
	      // map (x,y,z) to pixel in letter_[lmx:lpx,lmy:lpy]
	      // letter center = (cx,cy,cz)
	      // scaling factor for letter_:
	      //    (radius/block width) * (block cells/letter-cells)
	      //    [ fraction of block
	      //      covered by letter ]
	      // (x - cx)
	      // for given block index (ix,iy,iz)
	      //    map (x,y,z) to place in letter_[lmx:lpx,lmy:lpy]
	      //    letter_[lmx:lpx,lmy:lpy] extents == sphere with radius radius_
	      
	    }
	  }
	}
      }
    }
  }

}
