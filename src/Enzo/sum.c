#include <stdio.h>

int main(int argc, char ** argv)
{
  double x,s=0;
  int i = 0;
  int n = atoi(argv[1]);
  while (scanf ("%lf",&x)!= EOF) {
    s += x;
    if (++i == n) {
      printf ("%18.14f\n",s/n);
      s = 0.0;
      i = 0;
    }
  }
}
