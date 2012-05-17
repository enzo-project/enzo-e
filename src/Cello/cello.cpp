#include "cello.hpp"
#include "error.hpp"
//----------------------------------------------------------------------

namespace cello {

  // @@@ KEEP IN SYNCH WITH precision_enum in cello_precision.hpp
  const char * precision_name[7] = {
    "unknown",
    "default",
    "single",
    "double",
    "extended80",
    "extended96",
    "quadruple"
  };

  //----------------------------------------------------------------------

  int sizeof_precision(enum precision_type precision)
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

}
