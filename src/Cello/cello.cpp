#include "cello.hpp"
#include "error.hpp"
//----------------------------------------------------------------------

namespace cello {

  // @@@ KEEP IN SYNCH WITH precision_enum in cello.hpp
  const char * precision_name[7] = {
    "unknown",
    "default",
    "single",
    "double",
    "extended80",
    "extended96",
    "quadruple"
  };

  // @@@ KEEP IN SYNCH WITH type_enum in cello.hpp
  const char * type_name[NUM_TYPES] = {
    "unknown",     // unknown type
    "default",    // "default" floating-point precision, e.g. enzo_float
    "single",
    "double",
    "extended80",
    "extended96",
    "quadruple",
    "int8",
    "int16",
    "int32",
    "int64"
  };

  // @@@ KEEP IN SYNCH WITH type_enum in cello.hpp
  const int type_bytes[NUM_TYPES] = {
    0,
    0, // default
    4, // single
    8, // double
    10, // extended80
    12, // extended96
    16, // quadruple
    1, // int8
    2, // int16
    4, // int32
    8 // int64
  };

  bool type_is_int(int type) {
    return (type == type_int8 || 
	    type == type_int16 || 
	    type == type_int32 || 
	    type == type_int64);
  }

  bool type_is_float(int type) {
    return (type == type_single || 
	    type == type_double || 
	    type == type_quadruple ||
	    type == type_extended80 ||
	    type == type_extended96);
  }

  //----------------------------------------------------------------------

  int sizeof_precision(precision_type precision)
  {
    int size = 0;
    switch (precision) {
    case precision_unknown:
      break;
    case precision_default:
      size = sizeof_precision(default_precision);
      break;
    case precision_single:
      size = 4;
      break;
    case precision_double:
      size = 8;
      break;
    case precision_extended80:
      size = 10;
      break;
    case precision_extended96:
      size = 12;
      break;
    case precision_quadruple:
      size = 16;
      break;
    default:
      break;
    }
    return size;
  }

  //----------------------------------------------------------------------

  int is_precision_supported(precision_type precision)
  {
    int is_supported = false;
    switch (precision) {
    case precision_unknown:
      is_supported = false;
      break;
    case precision_default:
      is_supported = is_precision_supported(default_precision);
      break;
    case precision_single:
      is_supported = (sizeof(float)==4);
      break;
    case precision_double:
      is_supported = (sizeof(double)==8);
      break;
    case precision_extended80:
      is_supported = (sizeof(long double)==10);
      break;
    case precision_extended96:
      is_supported = (sizeof(long double)==12);
      break;
    case precision_quadruple:
      is_supported = (sizeof(long double)==16);
      break;
    default:
      is_supported = 0;
      break;
    }
    return is_supported;
  }

  //----------------------------------------------------------------------

  double machine_epsilon (precision_type precision)
  {
    double epsilon = 0.0;
    switch (precision) {
    case precision_unknown:
      epsilon = 0.0;
      break;
    case precision_default:
      epsilon = machine_epsilon(default_precision);
      break;
    case precision_single: // 32-bit
      epsilon =  5.96e-08; // http://en.wikipedia.org/wiki/Machine_epsilon
      break;
    case precision_double: // 64-bit
      epsilon =  1.11e-16; // http://en.wikipedia.org/wiki/Machine_epsilon
      break;
    default: // 9.63e-35 for 128-bit, but compliers vary on extended
	     // precision support
      WARNING ("cello::machine_epsilon", 
	       "Machine epsilon for extended precision unknown;"
	       " assuming double");
      epsilon =  1.11e-16; // http://en.wikipedia.org/wiki/Machine_epsilon;
      break;
    }
    return epsilon;
  }

  //----------------------------------------------------------------------

  void backtrace(const char * msg)
  {
    int j, nptrs;
#define SIZE 100
    void *buffer[100];
    char **strings;

    nptrs = ::backtrace(buffer, SIZE);

    /* The call backtrace_symbols_fd(buffer, nptrs, STDOUT_FILENO)
       would produce similar output to the following: */

    strings = backtrace_symbols(buffer, nptrs);
    if (strings == NULL) {
      perror("backtrace_symbols");
      exit(EXIT_FAILURE);
    }

    for (j = 0; j < nptrs; j++)
      printf("%s\n", strings[j]);

    free(strings);
  }

  int index_static()
  {
    return CkMyPe() % CONFIG_NODE_SIZE;
  }

}
