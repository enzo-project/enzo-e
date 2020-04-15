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

  /// Create a FaceFluxes object for the given face, field,
  /// block size, mesh width and time interval.
  FaceFluxes (Face face, int index_field,
              int nx, int ny, int nz,
              double hx, double hy, double hz, double dt);

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Set ghost zones if any
  void set_ghost(int gx, int gy, int gz);

  /// Set centering if non-centered
  void set_centering(int cx, int cy, int cz);

  /// Allocate the flux array.
  void allocate ();

  /// Deallocate the flux array.
  void deallocate();

  /// Clear flux array values
  void clear();

  /// Return the face associated with the FaceFluxes.
  Face face () const;

  /// Return the volume size of the fluxes. May differ from Block cell
  /// width if coarsening operations have been performed.
  void get_element_size (double *hx, double *hy, double * hz) const;

  /// Return loop limits on this face relative to neighbor face
  void get_limits(int * ixl, int * ixu,
                  int * iyl=0, int * iyu=0,
                  int * izl=0, int * izu=0) const;

  /// Return the time step dt of the accumulated fluxes. May differ
  /// from Block time step if multiple time steps have been
  /// accumulated.
  double time_step () const;
  
  /// Return the array dimensions, including adjustments for ghost or
  /// centering. One or more of mx,my,mz will be 1.  Stored block size
  /// (nx_,ny_,nz_) used unless provided explicitly by (nx,ny,nz)
  void get_dimensions (int *mx, int *my, int *mz) const;
   
  /// Copy flux values from an array to the FluxFaces. Array element
  /// array[ix*dx + iy*dy + iz*dz] should correspond to flux value
  /// (ix,iy,iz), where (0,0,0) <= (ix,iy,iz) < (mx,my,mz).
  void set_fluxes ( std::vector<double> array, int dx, int dy, int dz);
  
  /// Return the array of fluxes and associated strides (dx,dy,dz)
  /// such that the (ix,iy,iz) flux value is fluxes[ix*dx + iy*dy +
  /// iz*dz], where (0,0,0) <= (ix,iy,iz) < (mx,my,mz).
  std::vector<double> & get_fluxes (int * dx, int * dy, int * dz);
  
  /// Return the ratio of volume element resolutions h(ff_1) / h(ff_2)
  /// = {0.5, 1.0, 2.0} along each dimension between stored fluxes in
  /// two FaceFluxes objects. FaceFluxes are assumed to be associated
  /// with the same face. Must be 1.0 to compute sum or difference.
  friend float ratio_cell_width
  (const FaceFluxes & ff_1, const FaceFluxes & ff_2);

  /// Return the ratio of time steps dt(ff_1) / dt(ff_2) = {0.5, 1.0,
  /// 2.0} of fluxes between two FaceFluxes objects. FaceFluxes are
  /// assumed to be associated with the same face. Ratio must be 1.0
  /// to compute difference.
  friend float ratio_time_step
  (const FaceFluxes & ff_1, const FaceFluxes & ff_2);
  
  /// Coarsen a FaceFluxes object by reducing dimensions by two along
  /// each face dimension, and summing fine elements contained in each
  /// coarse flux element. Updates volume element resolution
  /// (hx,hy,hz) accordingly. Used for coarsening fine-level fluxes to
  /// match coarse level fluxes.
  void coarsen ();
  
  /// Add FaceFluxes object to this one. FaceFluxes are assumed to be
  /// associated with the same face. Used for accumulating fluxes with
  /// finer time steps until they match the coarser time step. Updates
  /// dt accordingly. Assumes spacially-conforming FaceFlux objects:
  /// FaceFluxes must be associated with the same face, and ratio of
  /// cell_widths must be 1.0
  FaceFluxes & operator += (const FaceFluxes & face_fluxes);
  
  /// Scale the fluxes array by a scalar constant.
  FaceFluxes & operator *= (double weight);
  
  /// Return a FaceFluxes object containing the difference (ff_1 -
  /// ff_2) between the two FaceFluxes. Used to compute flux
  /// correction factors. Assumes fully-conforming FaceFlux objects:
  /// FaceFluxes must be associated with the same face, and ratios of
  /// both cell_widths and time_steps must be 1.0.
  friend FaceFluxes operator -
  ( const FaceFluxes & ff_1, const FaceFluxes & ff_2 );
  
private: // functions

  /// Update dimensions (mx_,my_,mz_) after any of n?_, g?_, c?_ change
  void update_dimensions_()
  {
    const int rank = face_.rank();
    mx_ = (rank >= 1 && nx_ > 1) ? (nx_ + 2*gx_ + cx_) : 1;
    my_ = (rank >= 1 && ny_ > 1) ? (ny_ + 2*gy_ + cy_) : 1;
    mz_ = (rank >= 1 && nz_ > 1) ? (nz_ + 2*gz_ + cz_) : 1;
  }
  
private: // attributes

  // NOTE: change pup() function whenever attributes change

  // Face for which the fluxes are defined
  Face face_;

  // Index of the conserved field the fluxes are associated with
  int index_field_;

  // Dimensions of the fluxes array
  int mx_,my_,mz_;

  // Size of the block face
  int nx_,ny_,nz_;

  // Centering adjustment to flux array; same as FieldDescr
  int cx_,cy_,cz_;

  // Ghost zone adjustment to flux array
  int gx_,gy_,gz_;

  // Array of fluxes
  std::vector<double> fluxes_;

  // Physical cell size associated with the fluxes
  double hx_,hy_,hz_;

  // Time step associated with the fluxes
  double dt_;
};

#endif /* DATA_FACE_FLUXES_HPP */

