#include "patch.hpp"
#include "parallel.def"

Patch::Patch() 
{
  PARALLEL_PRINTF ("Patch()\n"); 
}
Patch::Patch(CkMigrateMessage *) 
{
  PARALLEL_PRINTF ("Patch(migrate)\n"); 
}
void Patch::advance() 
{
  PARALLEL_PRINTF ("advance()\n"); 
}
