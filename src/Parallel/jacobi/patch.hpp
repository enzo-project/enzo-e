#ifndef PATCH_HPP
#define PATCH_HPP

#include "parallel.def"

enum face_type {
  face_xm,
  face_xp,
  face_ym,
  face_yp,
  face_zm,
  face_zp
};

#include "test_jacobi.decl.h"

class Patch : public CBase_Patch
{
  
public:

  Patch() ;
  ~Patch() ;

  Patch(CkMigrateMessage *) ;

public: // entry methods

  void p_evolve(int patch_count, int patch_size);
  void p_receive(int face, int n, double * buffer);

private: // functions

  void initialize_();
  double initial_(double x, double y, double z);

  void allocate_(int n);

  void face_to_buffer_ (int face, double * buffer);
  void buffer_to_ghost_(int face, double * buffer);
  void print_ ();


private: // attributes

  int block_count_;
  int block_size_;

  double * values_;
  double * ghosts_[6];

  int cycle_values_;
  int cycle_ghosts_[6];

  int count_receive_;

};

#endif /* PATCH_HPP */
