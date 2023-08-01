#ifdef __cplusplus
#include <cstdio>
extern "C" {
#include <grackle.h>
}
#else
#include <stdio.h>
#include <grackle.h>
#endif
int main(int argc, char* argv[]){
  grackle_version gversion = get_grackle_version();
#ifdef __cplusplus
  std::puts(gversion.version);
#else
  puts(gversion.version);
#endif
  return 0;
}