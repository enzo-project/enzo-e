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

  /// Create an ininitialized FaceFluxes object. Called when unpacking
  FaceFluxes ()
  { }
  
  /// Create a FaceFluxes object for the given face, field, and block
  /// size. Optionally include centering adjustment (0 <= cx,cy,cz <=
  /// 1) for facet-, edge-, or corner-located field values
  FaceFluxes (Face face, int index_field,
              int nx, int ny, int nz,
              int cx=0, int cy=0, int cz=0);

  /// CHARM++ Pack / Unpack method
  void pup (PUP::er &p);

  /// Allocate the flux array and initialize values to 0.0
  void allocate ()
  {
    int mx,my,mz;
    get_size (&mx,&my,&mz);
    fluxes_.resize (mx*my*mz);
    clear();
  }

  /// Deallocate the flux array
  void deallocate()
  { fluxes_.clear(); }

  /// Set flux array values to 0.0
  void clear()
  { std::fill (fluxes_.begin(),fluxes_.end(),0.0); }

  /// Return the face associated with the FaceFluxes.
  Face face () const
  { return face_; }

  /// Return the associated Field index
  int index_field () const
  { return index_field_; }
  
    
  /// Return the array dimensions of the flux array, including any
  /// adjustments for centering. Indexing is ix + mx*(iy + my*iz).
  void get_size (int *mx, int *my, int *mz) const
  {
    if (mx) (*mx) = (nx_==0 || nx_==1) ?  1 : (nx_ + cx_);
    if (my) (*my) = (ny_==0 || ny_==1) ?  1 : (ny_ + cy_);
    if (mz) (*mz) = (nz_==0 || nz_==1) ?  1 : (nz_ + cz_);
  }
   
  /// Copy flux values from an array to the FluxFaces flux
  /// array. Array element array[ix*dx + iy*dy + iz*dz] should
  /// correspond to flux value (ix,iy,iz), where (0,0,0) <= (ix,iy,iz)
  /// < (mx,my,mz).
  void set_flux_array
  ( std::vector<cello_float> array, int dx, int dy, int dz);
  
  /// Return the array of fluxes and associated strides (dx,dy,dz)
  /// such that the (ix,iy,iz) flux value is fluxes[ix*dx + iy*dy +
  /// iz*dz], where (0,0,0) <= (ix,iy,iz) < (mx,my,mz).
  std::vector<cello_float> & flux_array (int * dx=0, int * dy=0, int * dz=0);
  
  /// Used for coarsening fine-level fluxes to match coarse level
  /// fluxes. Arguments (cx,cy,cz) specify the child indices of the
  /// block within its parent (not to be confused with centering
  /// (cx_,cy_,cz_); flux array size is kept the same, with offset
  /// determined by child indices
  void coarsen (int cx, int cy, int cz, int rank);
  
  /// Add FaceFluxes object to this one. Used for accumulating fluxes
  /// with finer time steps until they match the coarser time
  /// step. Assumes spacially-conforming FaceFlux objects
  void accumulate
  (const FaceFluxes & face_fluxes, int cx, int cy, int cz, int rank);
  
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
  std::vector<cello_float> fluxes_;

  // Index of the conserved field the fluxes are associated with
  int index_field_;

  // Size of the flux array, excluding centering adjustments; 1 for
  // axes orthogonal to face; 0 for axes greater than rank.
  int nx_,ny_,nz_;

  // Centering adjustment to flux array, 0 for centered
  int cx_,cy_,cz_;

};

#endif /* DATA_FACE_FLUXES_HPP */

