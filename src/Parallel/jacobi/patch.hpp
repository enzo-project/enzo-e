#ifndef PATCH_HPP
#define PATCH_HPP

#include "test_jacobi.decl.h"

class Patch : public CBase_Patch {
 public:
  Patch() ;
  Patch(CkMigrateMessage *) ;
  void advance() ;
};

#endif /* PATCH_HPP */
