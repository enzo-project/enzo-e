
#include "cello.hpp"

//----------------------------------------------------------------------

namespace cello {

  const char * precision_name[8] = {
    "unknown",
    "default",
    "single",
    "double",
    "extended80",
    "extended96",
    "quadruple"
  };

  //----------------------------------------------------------------------

  int sizeof_precision(enum precision_enum precision)
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

  int is_precision_supported(enum precision_enum precision)
  {
    switch (precision) {
    case precision_unknown:
      return 0;
      break;
    case precision_default:
      return is_precision_supported(default_precision);
      break;
    case precision_single:
      return (sizeof(float)==4);
      break;
    case precision_double:
      return (sizeof(double)==8);
      break;
    case precision_extended80:
      return (sizeof(long double)==10);
      break;
    case precision_extended96:
      return (sizeof(long double)==12);
      break;
    case precision_quadruple:
      return (sizeof(long double)==16);
      break;
    default:
      return 0;
      break;
    }
  }
}
