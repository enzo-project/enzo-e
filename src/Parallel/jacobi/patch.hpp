#ifndef PATCH_HPP
#define PATCH_HPP

#include "parallel.def"
#include "counter.hpp"

#include "test_jacobi.decl.h"

class Patch : public CBase_Patch
{
  
public:

  Patch(int patch_count, int patch_size, CProxy_Main main_proxy) ;
  ~Patch() ;

  Patch(CkMigrateMessage *) ;

public: // entry methods

  void p_evolve();
  void p_receive(int axis, int face, int n, double * buffer);

private: // functions

  void initialize_();

  void allocate_();
  void compute_();

  bool clear_boundary_ (int axis, int face, double * buffer);
  void face_to_buffer_ (int axis, int face, double * buffer);
  void buffer_to_ghost_(int axis, int face, double * buffer);
  void print_ ();
  void store_ ();
  double norm_ ();
  int id_();

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
  int ilvx_,ilvy_,ilvz_;
  /// Upper indices for block values
  int iuvx_,iuvy_,iuvz_;
  /// Lower indices for block ghosts
  int ilgx_,ilgy_,ilgz_;
  /// Upper indices for block ghosts
  int iugx_,iugy_,iugz_;

  double * values_;
  double * buffer_[3][2];


  int cycle_store_;
  int cycle_values_;
  int cycle_ghosts_[6];

  Counter receives_;

};

#endif /* PATCH_HPP */
