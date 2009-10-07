#include "parameters.hpp"

main(int argc, char **argv)
{
  /*  yydebug=1; */
  cello_parameters_read(stdin);
  cello_parameters_print();
}
