#include <stdio.h>
#include <stdlib.h>

main(int argc, char ** argv)
{
  if (argc != 2) {
    fprintf (stderr,"Usage: %s <enzo-file>\n",argv[0]);
    exit(1);
  }
  
}
