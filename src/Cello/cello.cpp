#include "cello.hpp"
#include "error.hpp"
#include "charm_simulation.hpp"
#include "simulation.hpp"
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
  bool type_is_valid(int type) {
    return (0 <= type && type < NUM_TYPES) &&
      (type_bytes[type] > 0);
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

  //---------------------------------------------------------------------- 

  Simulation * simulation()
  {
    return proxy_simulation.ckLocalBranch();
  }

  //---------------------------------------------------------------------- 
  ScalarDescr * scalar_descr_double()
  {
    static ScalarDescr * scalar_descr_double_[CONFIG_NODE_SIZE] = {nullptr};
    const int in = cello::index_static();
    if (scalar_descr_double_[in] == nullptr) {
      scalar_descr_double_[in] = simulation() ?
        simulation()->scalar_descr_double() : new ScalarDescr;
    }
    return scalar_descr_double_[in];
  }
  ScalarDescr * scalar_descr_long_double()
  {
    static ScalarDescr * scalar_descr_long_double_[CONFIG_NODE_SIZE] = {nullptr};
    const int in = cello::index_static();
    if (scalar_descr_long_double_[in] == nullptr) {
      scalar_descr_long_double_[in] = simulation() ?
        simulation()->scalar_descr_long_double() : new ScalarDescr;
    }
    return scalar_descr_long_double_[in];
  }
  ScalarDescr * scalar_descr_int()
  {
    static ScalarDescr * scalar_descr_int_[CONFIG_NODE_SIZE] = {nullptr};
    const int in = cello::index_static();
    if (scalar_descr_int_[in] == nullptr) {
      scalar_descr_int_[in] = simulation() ?
        simulation()->scalar_descr_int() : new ScalarDescr;
    }
    return scalar_descr_int_[in];
  }
  ScalarDescr * scalar_descr_sync()
  {
    static ScalarDescr * scalar_descr_sync_[CONFIG_NODE_SIZE] = {nullptr};
    const int in = cello::index_static();
    if (scalar_descr_sync_[in] == nullptr) {
      scalar_descr_sync_[in] = simulation() ?
        simulation()->scalar_descr_sync() : new ScalarDescr;
    }
    return scalar_descr_sync_[in];
  }
  ScalarDescr * scalar_descr_void()
  {
    static ScalarDescr * scalar_descr_void_[CONFIG_NODE_SIZE] = {nullptr};
    const int in = cello::index_static();
    if (scalar_descr_void_[in] == nullptr) {
      scalar_descr_void_[in] = simulation() ?
        simulation()->scalar_descr_void() : new ScalarDescr;
    }
    return scalar_descr_void_[in];
  }
  //---------------------------------------------------------------------- 

  CProxy_Block block_array()
  {
    return cello::simulation()->hierarchy()->block_array();
  }

  //---------------------------------------------------------------------- 

  const Config * config()
  {
    return simulation() ? simulation()->config() : NULL;
  }

  //---------------------------------------------------------------------- 

  const Parameters * parameters()
  {
    return simulation() ? simulation()->parameters() : NULL;
  }

  //---------------------------------------------------------------------- 

  Problem * problem()
  {
    return simulation() ? simulation()->problem() : NULL;
  }

  //---------------------------------------------------------------------- 

  Hierarchy * hierarchy()
  {
    return simulation() ? simulation()->hierarchy() : NULL;
  }

  //---------------------------------------------------------------------- 

  FieldDescr * field_descr()
  {
    static FieldDescr * field_descr_[CONFIG_NODE_SIZE] = {nullptr};
    const int in = cello::index_static();
    if (field_descr_[in] == nullptr) {
      field_descr_[in] = simulation() ?
        simulation()->field_descr() : new FieldDescr;
    }
    return field_descr_[in];
  }

  //---------------------------------------------------------------------- 

  Monitor * monitor()
  {
    return simulation() ? simulation()->monitor() : Monitor::instance();
  }

  //---------------------------------------------------------------------- 

  ParticleDescr * particle_descr()
  {
    static ParticleDescr * particle_descr_[CONFIG_NODE_SIZE] = {nullptr};
    const int in = cello::index_static();
    if (particle_descr_[in] == nullptr) {
      particle_descr_[in] = simulation() ?
        simulation()->particle_descr() : new ParticleDescr;
    }
    return particle_descr_[in];
  }

  //---------------------------------------------------------------------- 

  Output * output(int index)
  {
    return problem() ? problem()->output(index) : NULL;
  }

  //---------------------------------------------------------------------- 

  Solver * solver(int index)
  {
    return problem() ? problem()->solver(index) : NULL;
  }

  //---------------------------------------------------------------------- 

  Units * units()
  {
    return problem() ? problem()->units() : NULL;
  }

  //----------------------------------------------------------------------

  int rank()
  {
    return simulation() ? simulation()->rank() : 1;
  }

  //----------------------------------------------------------------------

  int num_children()
  {
    int r = rank();
    return (r==1) ? 2 : ( (r==2) ? 4 : 8 );
  }
  

}
