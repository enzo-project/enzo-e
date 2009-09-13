#include <stdio.h>
#include "cello.h"
#include "tree4.h"

main()
{
  const int width = 4;
  const int height = 4;
  const bool array[] = {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};

  Tree4 tree(array,width, height);
  printf ("sizeof (Tree4) = %d\n",sizeof(Tree4));
	    
  
}
