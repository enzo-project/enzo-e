#include "test_jacobi.decl.h"
#include "patch.hpp"

Patch::Patch() 
{
  CkPrintf ("Patch()\n"); 
}
Patch::Patch(CkMigrateMessage *) 
{
  CkPrintf ("Patch(migrate)\n"); 
}
void Patch::advance() 
{
  CkPrintf ("advance()\n"); 
}
