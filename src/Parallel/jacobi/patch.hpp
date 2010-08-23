#ifndef PATCH_HPP
#define PATCH_HPP

#include "parallel.def"

enum neighbor_type {
  neighbor_xm,
  neighbor_xp,
  neighbor_ym,
  neighbor_yp,
  neighbor_zm,
  neighbor_zp
};

#include "test_jacobi.decl.h"

class Patch : public CBase_Patch
{
  
public:

  Patch() ;
  ~Patch() ;

  Patch(CkMigrateMessage *) ;

  void advance(int n);
  void allocate(int n);
  void neighbors (CProxy_Patch xm, CProxy_Patch xp,
		  CProxy_Patch ym, CProxy_Patch yp,
		  CProxy_Patch zm, CProxy_Patch zp);

  void receive(int n, double * values);
  void initialize();

  void receive(int type, int n, double * buffer);
  void prepare_(int type, double * buffer);
  void unpack_(int type, double * buffer);

private:

  double initial_(double x, double y, double z);

  int n_;

  double * values_;
  double * ghosts_[6];

  int cycle_values_;
  int cycle_ghosts_;

  CProxy_Patch neighbor_[6];
};

#endif /* PATCH_HPP */
