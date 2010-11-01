#ifndef PATCH_HPP
#define PATCH_HPP

#include "parallel.def"

#include "test_parallel_jacobi.decl.h"

class Patch : public CBase_Patch
{
  
public:

  /// Create a Patch object
  Patch(int patch_count, 
	int patch_size, 
	int cycles_max, 
	CProxy_Main main_proxy) ;

  /// Delete a Patch object
  ~Patch() ;

  /// Migrate a Patch object
  Patch(CkMigrateMessage *) ;

public: // parallel methods

  /// Initial patch advance, ending with receive_()
  void p_evolve();

  /// Perform one time step of the computation, ending with main::p_next()
  void compute_();

  /// [FORK] Calls p_receive() for each face
  void receive_ ();

  /// [JOIN] Receive a single face, the last continues with compute_()
  void p_receive(int axis, int face, int n, double * buffer);

private: // functions

  /// Allocate values and buffers, called from constructor
  void allocate_();

  /// Set initial conditions, called from constructor
  void initialize_();

  /// Clear boundary values
  bool clear_boundary_ (int axis, int face, double * buffer);

  /// Copy values from a face to a buffer
  void face_to_buffer_ (int axis, int face, double * buffer);

  /// Copy values from a buffer to ghosts
  void buffer_to_ghost_(int axis, int face, double * buffer);

  /// Display values to stdout
  void print_ ();

  /// Write data to an HDF5 file
  void store_ ();

  /// Return the local sum-of-squares
  double norm_ ();

//   /// Return the Patch ID
//   int id_();

private: // attributes

  CProxy_Main main_proxy_;

  /// Number of blocks along each axis
  int nbx_,nby_,nbz_;
  /// Depth of ghost zone layer along each block face
  int ngx_,ngy_,ngz_;  
  /// Number of values along each block axis
  int nvx_,nvy_,nvz_;
  /// Allocated array size along each block axis
  int nax_,nay_,naz_;
  /// Lower indices for block values
  int ixl_,iyl_,izl_;
  /// Upper indices for block values
  int ixu_,iyu_,izu_;

  double * values_;
  double * buffer_[3][2];


  int cycle_store_;
  int cycle_values_;
  int cycle_max_;

  Counter receives_;

};

#endif /* PATCH_HPP */
