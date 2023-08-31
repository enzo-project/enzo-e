
#include <sys/types.h>
#include <time.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include "cnames.h"
#include "enzo_defines.hpp"

int curr_rand_seed = 0;

void initrandomx(int* x)
{
  srand(*x);
}

void initrandom()
{
  srand((int)time(0));
}

double randomnr()
{
  return ((double)rand_r(&curr_rand_seed))/RAND_MAX;
}

int crandseed()
{
  return curr_rand_seed;
}
void FORTRAN_NAME(setcrandseed)(int * rand_seed)
{
  curr_rand_seed = *rand_seed;
}

void seedrang( int* ss )
{
   int fd = open("/dev/urandom",O_RDONLY);
   int mask = 0x0001FFFF;
   int nr;
   nr = read(fd,ss,sizeof(int));
   *ss &= mask;
   close(fd);
}
