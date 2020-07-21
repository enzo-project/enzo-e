// See LICENSE_CELLO file for license and copyright information

/// @file     data_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-10-15
/// @brief    [\ref Data] Implementation of the FaceFluxes class

#include "data.hpp"

// #define DEBUG_REFRESH

// #define TRACE CkPrintf ("%d TRACE %s:%d\n",CkMyPe(),__FILE__,__LINE__); fflush(stdout);
#define TRACE /*...*/

FaceFluxes::FaceFluxes (Face face, int index_field,
                        int nx, int ny, int nz,
                        int level, double dt,
                        int cx, int cy, int cz)
  : face_(face),
    level_block_(level),
    level_neighbor_(level),
    dt_block_(dt),
    dt_neighbor_(dt),
    fluxes_(),
    index_field_(index_field),
    nx_(nx),ny_(ny),nz_(nz),
    cx_(cx),cy_(cy),cz_(cz)
{
  int fx,fy,fz; // face
  face.get_face(&fx,&fy,&fz);
  // set face size to 1 when orthogonal to face
  // (should be 0 for axes greater than rank)
  if (fx && nx) nx_ = 1;
  if (fy && ny) ny_ = 1;
  if (fz && nz) nz_ = 1;
  ASSERT3 ("FaceFluxes::FaceFluxes",
          "Block dimensions %d %d %d must be even or 1",
          nx_,ny_,nz_,
           (( nx_==1 || ((nx_&1)==0) ) &&
            ( ny_==1 || ((ny_&1)==0) ) &&
            ( nz_==1 || ((nz_&1)==0) )));
}

//----------------------------------------------------------------------

void FaceFluxes::pup (PUP::er &p)
{
  TRACEPUP;

  p | face_;

  p | level_block_;
  p | level_neighbor_;
  p | dt_block_;
  p | dt_neighbor_;

  p | fluxes_;

  p | index_field_;

  p | nx_;
  p | ny_;
  p | nz_;

  p | cx_;
  p | cy_;
  p | cz_;
}

//----------------------------------------------------------------------

void FaceFluxes::get_size
(int * nx, int * ny, int * nz, int rank, int level_neighbor) const
{
  int ix,iy,iz;
  face_.get_face(&ix,&iy,&iz);
  int mx,my,mz;
  get_dimensions (&mx,&my,&mz);
  (*nx) = mx;
  (*ny) = (rank >= 2) ? my : 1;
  (*nz) = (rank >= 3) ? mz : 1;
}

//----------------------------------------------------------------------

void FaceFluxes::set_flux_array ( std::vector<cello_float> array,
                                  int dx, int dy, int dz)
{
  int mx,my,mz;
  get_dimensions(&mx,&my,&mz);
  ASSERT5("FaceFluxes::set_flux_array",
          "Input array size %lu is smaller than required %d = %d * %d *%d",
          array.size(),mx*my*mz,mx,my,mz,
          (array.size() >= mx*my*mz) );
  ASSERT4("FaceFluxes::set_flux_array",
          "Input array size %lu is too small for strides (dx,dy,dz)=(%d,%d,%d)",
          array.size(),dx,dy,dz,
          (array.size() >= (mx-1)*dx+(my-1)*dy+(mz-1)*dz) );

  ASSERT5("FaceFluxes::set_flux_array",
          "Flux array size %lu is smaller than required %d = %d * %d *%d",
          fluxes_.size(),mx*my*mz,mx,my,mz,
          (fluxes_.size() >= mx*my*mz) );

  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
        int i_f = ix + mx*(iy+my*iz);
        int i_a = ix*dx + iy*dy + iz*dz;
        fluxes_[i_f] = array[i_a];
      }
    }
  }
}
  
//----------------------------------------------------------------------

std::vector<cello_float> & FaceFluxes::flux_array (int * dx, int * dy, int *dz)
{
  int mx,my,mz;
  get_dimensions(&mx,&my,&mz);
  if (dx) (*dx) = 1;
  if (dy) (*dy) = mx;
  if (dz) (*dz) = mx*my;

  return fluxes_;
}
  
//----------------------------------------------------------------------

float ratio_cell_volume
(const FaceFluxes & ff_1, const FaceFluxes & ff_2, const int rank)
{
  const int L_1 = ff_1.level_neighbor();
  const int L_2 = ff_2.level_neighbor();
  if (L_1 <= L_2) {
    const int shift = (1 << (L_2 - L_1)*rank);
    return 1.0 * shift;
  } else if (L_1 > L_2) {
    const int shift = (1 << (L_1 - L_2)*rank);
    return 1.0 / shift;
  }
}

//----------------------------------------------------------------------

float ratio_time_step (const FaceFluxes & ff_1, const FaceFluxes & ff_2)
{
  return ff_1.time_step_block() / ff_2.time_step_block();
}
  
//----------------------------------------------------------------------

void FaceFluxes::coarsen (int cx, int cy, int cz, int rank)
{
  std::vector<cello_float> fluxes_fine = fluxes_;
  std::fill(fluxes_.begin(),fluxes_.end(),0.0);

  int mxf,myf,mzf;
  get_dimensions(&mxf,&myf,&mzf);

  const int mxc = (nx_ > 1) ? (nx_/2 + cx_) : 1;
  const int myc = (ny_ > 1) ? (ny_/2 + cy_) : 1;
  const int mzc = (nz_ > 1) ? (nz_/2 + cz_) : 1;

  if (rank == 2) {

    const int dxc = 1;
    const int dyc = mxc;
    
    const int ixc0 = cx*mxc*dxc;
    const int iyc0 = cy*myc*dyc;
    
    // Ratio of flux areas
    const cello_float ratio = 0.5;
    if (mxf == 1) {

      const int ixf = 0;
      const int ic0 = iyc0;
      for (int iyf=0; iyf<myf; iyf++) {
        int i_f = ixf + mxf*iyf;
        int iycm = dyc*(iyf>>1);
        int iycp = dyc*((iyf+cy_)>>1);
        fluxes_[ic0+iycm] += ratio*fluxes_fine[i_f];
        fluxes_[ic0+iycp] += ratio*fluxes_fine[i_f];
      }

    } else if (myf == 1) {

      const int iyf = 0;
      const int ic0 = ixc0;
      for (int ixf=0; ixf<mxf; ixf++) {
        int i_f = ixf + mxf*iyf;
        int ixcm = dxc*(ixf>>1);
        int ixcp = dxc*((ixf+cx_)>>1);
        fluxes_[ic0+ixcm] += ratio*fluxes_fine[i_f];
        fluxes_[ic0+ixcp] += ratio*fluxes_fine[i_f];
      }
    }

  } else if (rank == 3) {
    
    const int dxc = 1;
    const int dyc = mxc;
    const int dzc = mxc*myc;
    
    const int ixc0 = cx*mxc;
    const int iyc0 = cy*myc;
    const int izc0 = cz*mzc;

    // Ratio of flux areas
    const cello_float ratio = 0.25;
    
    if (mxf == 1) {
      
      const int ixf = 0;
      const int ic0 =iyc0 + izc0;
      for (int izf=0; izf<mzf; izf++) {
        int izcm = dzc*(izf>>1);
        int izcp = dzc*((izf+cz_)>>1);
        for (int iyf=0; iyf<myf; iyf++) {
          int i_f = ixf + mxf*(iyf + myf*izf);
          int iycm = dyc*(iyf>>1);
          int iycp = dyc*((iyf+cy_)>>1);

          fluxes_[ic0+iycm+izcm] += ratio*fluxes_fine[i_f];
          fluxes_[ic0+iycm+izcp] += ratio*fluxes_fine[i_f];
          fluxes_[ic0+iycp+izcm] += ratio*fluxes_fine[i_f];
          fluxes_[ic0+iycp+izcp] += ratio*fluxes_fine[i_f];
        }
      }
    } else if (myf == 1) {

      const int iyf=0;
      const int ic0 =ixc0 + izc0;
      for (int izf=0; izf<mzf; izf++) {
        int izcm = dzc*(izf>>1);
        int izcp = dzc*((izf+cz_)>>1);
        for (int ixf=0; ixf<mxf; ixf++) {
          int ixcm = dxc*(ixf>>1);
          int ixcp = dzc*((ixf+cx_)>>1);
          int i_f = ixf + mxf*(iyf + myf*izf);
          fluxes_[ic0+izcm+ixcm] += ratio*fluxes_fine[i_f];
          fluxes_[ic0+izcp+ixcm] += ratio*fluxes_fine[i_f];
          fluxes_[ic0+izcm+ixcp] += ratio*fluxes_fine[i_f];
          fluxes_[ic0+izcp+ixcp] += ratio*fluxes_fine[i_f];
        }
      }
    } else if (mzf == 1) {

      const int izf=0;
      const int ic0 =ixc0 + iyc0;
      for (int iyf=0; iyf<myf; iyf++) {
        int iycm = dyc*(iyf>>1);
        int iycp = dyc*((iyf+cy_)>>1);
        for (int ixf=0; ixf<mxf; ixf++) {
          int ixcm = dxc*(ixf>>1);
          int ixcp = dxc*((ixf+cx_)>>1);
          int i_f = ixf + mxf*(iyf + myf*izf);
          fluxes_[ic0+ixcm+iycm] += ratio*fluxes_fine[i_f];
          fluxes_[ic0+ixcp+iycm] += ratio*fluxes_fine[i_f];
          fluxes_[ic0+ixcm+iycp] += ratio*fluxes_fine[i_f];
          fluxes_[ic0+ixcp+iycp] += ratio*fluxes_fine[i_f];
        }
      }
    }

  } else {

    ERROR1 ("FaceFluxes::coarsen()",
            "Unsupported rank = %d",
            rank);

  }
  
  --level_neighbor_;
}
  
//----------------------------------------------------------------------

void FaceFluxes::accumulate
(const FaceFluxes & ff_2, int cx,int cy, int cz, int rank)
{
  FaceFluxes & ff_1 = *this;

  int nx1,ny1,nz1;
  int mx1,my1,mz1;
  
  ff_1.get_size(&nx1,&ny1,&nz1,rank,ff_2.level_block());
  ff_1.get_dimensions(&mx1,&my1,&mz1);
  
  int nx2,ny2,nz2;
  int mx2,my2,mz2;
  ff_2.get_size(&nx2,&ny2,&nz2,rank,ff_1.level_block());
  ff_2.get_dimensions(&mx2,&my2,&mz2);

  ASSERT8 ("FaceFluxes::operator +=()",
          "Flux array sizes level %d (%d %d %d) and level %d (%d %d %d) must be conforming",
           ff_1.level_block(),nx1,ny1,nz1,
           ff_2.level_block(),nx2,ny2,nz2,
           ((nx1==nx2) && (ny1==ny2) && (nz1==nz2)));

  for (int iz=0; iz<nz1; iz++) {
    for (int iy=0; iy<ny1; iy++) {
      for (int ix=0; ix<nx1; ix++) {
        int i1=ix + mx1*(iy + my1*iz);
        int i2=ix + mx2*(iy + my2*iz);
        ff_1.fluxes_[i1] += ff_2.fluxes_[i2];
      }
    }
  }
 
}
  
//----------------------------------------------------------------------

FaceFluxes & FaceFluxes::operator *= (double weight)
{
  int mx,my,mz;
  get_dimensions(&mx,&my,&mz);
  const int m=mx*my*mz;

  for (int i=0; i<m; i++) fluxes_[i] *= weight;
  
  return *this;
}
  
//======================================================================

int FaceFluxes::data_size () const
{
  int size = 0;

#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH FaceFluxes::data_size()\n",CkMyPe());
  fflush(stdout);
#endif  
  size += face_.data_size();
  size += sizeof(int);    // std::vector<cello_float> fluxes_;
  int n = fluxes_.size();
  size += n*sizeof(cello_float);
  size += sizeof(int);    // int level_block_;
  size += sizeof(int);    // int level_neighbor_;
  size += sizeof(double);    // double dt_block_;
  size += sizeof(double);    // double dt_neighbor_;
  size += sizeof(int);    // int index_field_;
  size += 3*sizeof(int);  // int nx_,ny_,nz_;
  size += 3*sizeof(int);  // int cx_,cy_,cz_;
#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH FaceFluxes::data_size() %d\n",CkMyPe(),size);
#endif  

  return size;
}

//----------------------------------------------------------------------

char * FaceFluxes::save_data (char * buffer) const
{
#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH FaceFluxes::save_data()\n",CkMyPe());
  fflush(stdout);
#endif  
  union {
    int  * pi;
    char * pc;
    double * pd;
    cello_float * pcf;
  };

  pc = (char *) buffer;

  pc = face_.save_data(pc);
  int n = fluxes_.size();
  (*pi++) = n;
  for (int i=0; i<n; i++) {
    (*pcf++) = fluxes_[i];
  }
  (*pi++) = level_block_;
  (*pi++) = level_neighbor_;
  (*pd++) = dt_block_;
  (*pd++) = dt_neighbor_;
  (*pi++) = index_field_;
  (*pi++) = nx_;
  (*pi++) = ny_;
  (*pi++) = nz_;
  (*pi++) = cx_;
  (*pi++) = cy_;
  (*pi++) = cz_;

  
  ASSERT2("FaceFluxes::save_data()",
	  "Buffer has size %ld but expecting size %d",
	  (pc-buffer),data_size(),
	  ((pc-buffer) == data_size()));

#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH FaceFluxes::save_data() %d\n",CkMyPe(),pc-buffer);
#endif  
  return pc;
}

//----------------------------------------------------------------------

char * FaceFluxes::load_data (char * buffer)
{
  TRACE;
#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH FaceFluxes::load_data()\n",CkMyPe());
  fflush(stdout);
#endif  
  union {
    int  * pi;
    char * pc;
    double * pd;
    cello_float * pcf;
  };

  pc = (char *) buffer;

  TRACE;
  pc = face_.load_data(pc);
  int n = (*pi++);
  TRACE;
  fluxes_.resize(n);
  for (int i=0; i<n; i++) {
    fluxes_[i] = (*pcf++);
  }
  TRACE;
  level_block_ = (*pi++);
  level_neighbor_ = (*pi++);
  dt_block_ = (*pd++);
  dt_neighbor_ = (*pd++);
  index_field_ = (*pi++);
  nx_ = (*pi++);
  ny_ = (*pi++);
  nz_ = (*pi++);
  cx_ = (*pi++);
  cy_ = (*pi++);
  cz_ = (*pi++);
  TRACE;

#ifdef DEBUG_REFRESH
  CkPrintf ("%d DEBUG_REFRESH FaceFluxes::load_data() %d\n",CkMyPe(),pc-buffer);
#endif  
  ASSERT2("FaceFluxes::load_data()",
	  "Buffer has size %ld but expecting size %d",
	  (pc-buffer),data_size(),
	  ((pc-buffer) == data_size()));
#ifdef DEBUG_REFRESH
  print("BLOCK","FaceFluxes::load_data()");
#endif  

  TRACE;
  return pc;
}

void FaceFluxes::print(const std::string block_name,const char * msg) const
{
  const int ip=CkMyPe();
  int ix,iy,iz;
  face().get_face(&ix,&iy,&iz);
  int rx,ry,rz;
  face().get_normal(&rx,&ry,&rz);
  int mx,my,mz;
  get_dimensions(&mx,&my,&mz);
  const int m=mx*my*mz;
  double sum=0.0;
  for (int i=0; i<m; i++) sum += std::abs(fluxes_[i]);
  CkPrintf ("%d FFP: %s FaceFluxes %s %p size %d level %d field %d face %d %d %d avg %g\n",
            ip,block_name.c_str(),msg,(void*)this,m,level_block_,index_field_,ix,iy,iz,sum/m);
  //  CkPrintf ("%d FFP:       index_field  %d\n",ip,index_field_);
  //  CkPrintf ("%d FFP:          ix,iy,iz  %d %d %d\n",ip,ix,iy,iz);
  //  CkPrintf ("%d FFP:          rx,ry,rz  %d %d %d\n",ip,rx,ry,rz);
  //  CkPrintf ("%d FFP:          nx,ny,nz  %d %d %d\n",ip,nx_,ny_,nz_);
  //  CkPrintf ("%d FFP:          cx,cy,cz  %d %d %d\n",ip,cx_,cy_,cz_);
  //  CkPrintf ("%d FFP:    block level dt  %d %g\n",ip,level_block_,dt_block_);
  //  CkPrintf ("%d FFP: neighbor level dt  %d %g\n",
  //            ip,level_neighbor_,dt_neighbor_);
  //  CkPrintf ("%d FFP:        avg fluxes  %g\n",ip,sum/m);
}
