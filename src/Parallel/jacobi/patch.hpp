#ifndef PATCH_HPP
#define PATCH_HPP

#include "parallel.def"
#include PARALLEL_CHARM_INCLUDE(test_jacobi.decl.h)
PARALLEL_CLASS_DECL(Patch) 
{
 public:
  Patch() ;
  Patch(CkMigrateMessage *) ;
  void advance() ;
  void set_bounds(double xm, double xp,
		  double ym, double yp,
		  double zm, double zp );
private:
  double lower_[3];      // lower corner
  double upper_[3];      // upper corner
  
};

#endif /* PATCH_HPP */
