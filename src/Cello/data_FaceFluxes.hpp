// See LICENSE_CELLO file for license and copyright information

/// @file     data_FaceFluxes.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-10-15
/// @brief    [\ref Data] Declaration of the FaceFluxes class

#ifndef DATA_FACE_FLUXES_HPP
#define DATA_FACE_FLUXES_HPP

class Face;

class FaceFluxes {

  /// @class    FaceFluxes
  /// @ingroup  Data
  /// @brief    [\ref Data] 

public: // interface

  /// Create an ininitialized FaceFluxes object.  Called when unpacking
  FaceFluxes ()
  { }
  
  /// Create a FaceFluxes object for the given face, field, block
  /// size, mesh width and time interval.  Optionally include
  /// centering adjustment (0 <= cx,cy,cz <= 1) for facet-, edge-, or
  /// corner-located field values.
  FaceFluxes (Face face, int index_field,
              int nx, int ny, int nz,
              int level, double dt,
              int cx=0, int cy=0, int cz=0);

  /// CHARM++ Pack / Unpack method
  void pup (PUP::er &p);

  /// Allocate the flux array and clear values
  void allocate ()
  {
    int mx,my,mz;
    get_dimensions (&mx,&my,&mz);
    fluxes_.resize (mx*my*mz);
    clear();
  }

  /// Deallocate the flux array.
  void deallocate()
  { fluxes_.resize(0); }

  /// Clear flux array values
  void clear()
  { std::fill (fluxes_.begin(),fluxes_.end(),0.0); }

  /// Return the face associated with the FaceFluxes.
  Face face () const
  { return face_; }

  /// Return the mesh level associated with the neighbor block (may
  /// change if fluxes are coarsened)
  int level_neighbor () const
  { return level_neighbor_; }

  /// Return the mesh level associated with the block on which the
  /// fluxes were originally defined
  int level_block () const
  { return level_block_; }
  
  /// Return the time step associated with the neighbor block (may
  /// change if multiple time steps accumulated)
  double time_step_neighbor () const
  { return dt_neighbor_; }
    
  /// Return the time step associated with the block on which the
  /// fluxes were originally defined
  double time_step_block () const
  { return dt_block_; }
  
  /// Return loop limits on this face relative to given neighbor
  /// face fluxes.  Must be conforming.
  void get_start(int * ix0, int * iy0, int * iz0,
                 int cx, int cy, int cz,
                 const FaceFluxes & ff) const
  {
    if (level_block() >= ff.level_block()) {
      // neighbor block is not in finer level
      (*ix0) = 0;
      (*iy0) = 0;
      (*iz0) = 0;
    } else {
      // neighbor block is finer, use neighbor child array for start
      int ix,iy,iz;
      face_.get_face(&ix,&iy,&iz);
      (*ix0) = (ix != 0) ? 0 : cx*nx_/2;
      (*iy0) = (iy != 0) ? 0 : cy*ny_/2;
      (*iz0) = (iz != 0) ? 0 : cz*nz_/2;
    }
  }

  void get_size(int * nx, int * ny, int * nz,
                const FaceFluxes & ff) const
  {
    int ix,iy,iz;
    face_.get_face(&ix,&iy,&iz);
    int mx,my,mz;
    get_dimensions (&mx,&my,&mz);
    if (level_block() >= ff.level_block()) {
      // neighbor block is not in finer level
      (*nx) = mx;
      (*ny) = my;
      (*nz) = mz;
    } else {
      // neighbor block is finer, size is halved
      (*nx) = (ix != 0) ? 1 : (mx+cx_) / 2;
      (*ny) = (iy != 0) ? 1 : (my+cy_) / 2;
      (*nz) = (iz != 0) ? 1 : (mz+cz_) / 2;
    }
  }
  
  /// Return the array dimensions of the flux array, including
  /// adjustments for centering.  Indexing is ix + mx*(iy + my*iz).
  void get_dimensions (int *mx, int *my, int *mz) const
  {
    if (mx) (*mx) = (nx_==0 || nx_==1) ?  1 : (nx_ + cx_);
    if (my) (*my) = (ny_==0 || ny_==1) ?  1 : (ny_ + cy_);
    if (mz) (*mz) = (nz_==0 || nz_==1) ?  1 : (nz_ + cz_);
  }
   
  /// Copy flux values from an array to the FluxFaces. Array element
  /// array[ix*dx + iy*dy + iz*dz] should correspond to flux value
  /// (ix,iy,iz), where (0,0,0) <= (ix,iy,iz) < (mx,my,mz).
  void set_flux_array ( std::vector<double> array, int dx, int dy, int dz);
  
  /// Return the array of fluxes and associated strides (dx,dy,dz)
  /// such that the (ix,iy,iz) flux value is fluxes[ix*dx + iy*dy +
  /// iz*dz], where (0,0,0) <= (ix,iy,iz) < (mx,my,mz).
  std::vector<double> & flux_array (int * dx=0, int * dy=0, int * dz=0);
  
  /// Return the ratio of volume element resolutions h(ff_1) / h(ff_2)
  /// = {0.5, 1.0, 2.0} along each dimension between stored fluxes in
  /// two FaceFluxes objects. FaceFluxes are assumed to be associated
  /// with the same face. Must be 1.0 to compute sum or difference.
  friend float ratio_cell_volume
  (const FaceFluxes & ff_1, const FaceFluxes & ff_2, int rank);

  /// Return the ratio of time steps dt(ff_1) / dt(ff_2) = {0.5, 1.0,
  /// 2.0} of fluxes between two FaceFluxes objects. FaceFluxes are
  /// assumed to be associated with the same face. Ratio must be 1.0
  /// to compute difference.
  friend float ratio_time_step
  (const FaceFluxes & ff_1, const FaceFluxes & ff_2);
  
  /// Coarsen a FaceFluxes object by reducing dimensions by two along
  /// each face dimension, and summing fine elements contained in each
  /// coarse flux element. Updates level_neighbor_ accordingly. Used for
  /// coarsening fine-level fluxes to match coarse level fluxes.
  void coarsen ();
  
  /// Add FaceFluxes object to this one. FaceFluxes are assumed to be
  /// associated with the same face. Used for accumulating fluxes with
  /// finer time steps until they match the coarser time step. Updates
  /// dt_neighbor accordingly. Assumes spacially-conforming FaceFlux
  /// objects: FaceFluxes must be associated with the same face, and
  void accumulate (const FaceFluxes & face_fluxes, int cx, int cy, int cz);
  
  /// Scale the fluxes array by a scalar constant.
  FaceFluxes & operator *= (double weight);
  
  //--------------------------------------------------

  /// Return the number of bytes required to serialize the data object
  int data_size () const;

  /// Serialize the object into the provided empty memory buffer.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * save_data (char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * load_data (char * buffer);

private: // attributes

  // NOTE: change pup() method whenever attributes change

  // Face for which the fluxes are defined
  Face face_;

  // Array of fluxes ix + mx*(iy+my*iz)
  std::vector<double> fluxes_;

  // Original mesh refinement level
  int level_block_;

  // Mesh refinement level associated with neighbor block
  int level_neighbor_;
  
  // Original time step 
  double dt_block_;

  // Time step associated with neighbor block
  double dt_neighbor_;

  // Index of the conserved field the fluxes are associated with
  int index_field_;

  // Size of the flux array, excluding centering adjustments; 1 for
  // axes orthogonal to face; 0 for axes greater than rank.
  int nx_,ny_,nz_;

  // Centering adjustment to flux array, 0 for centered
  int cx_,cy_,cz_;

};

#endif /* DATA_FACE_FLUXES_HPP */

