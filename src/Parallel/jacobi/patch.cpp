#include <stdio.h>

#include "patch.hpp"
#include "parallel.def"

Patch::Patch() 
{
  PARALLEL_PRINTF ("Patch()\n"); 
  lower_[0] = 0.0;
  lower_[0] = 0.0;
  lower_[0] = 0.0;
  upper_[0] = 1.0;
  upper_[0] = 1.0;
  upper_[0] = 1.0;
}
Patch::Patch(CkMigrateMessage *) 
{
  PARALLEL_PRINTF ("Patch(migrate)\n"); 
}
void Patch::advance() 
{
  PARALLEL_PRINTF ("advance()\n"); 
}
void Patch::set_bounds
(
 double xm, double xp,
 double ym, double yp,
 double zm, double zp 
 )
{
}
