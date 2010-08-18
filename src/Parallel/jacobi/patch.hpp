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
};

#endif /* PATCH_HPP */
