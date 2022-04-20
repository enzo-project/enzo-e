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

  int digits_max(int precision)
  {
    switch (precision) {
    case precision_single:
      return std::numeric_limits<float>::digits10;
    case precision_double:
      return std::numeric_limits<double>::digits10;
    case precision_quadruple:
      return std::numeric_limits<long double>::digits10;
    default:
      return 0;
    }
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
    return cello::hierarchy()->block_array();
  }

  //---------------------------------------------------------------------- 

  const Config * config()
  {
    return simulation() ? simulation()->config() : nullptr;
  }

  //---------------------------------------------------------------------- 

  const Parameters * parameters()
  {
    return simulation() ? simulation()->parameters() : nullptr;
  }

  //---------------------------------------------------------------------- 

  Problem * problem()
  {
    return simulation() ? simulation()->problem() : nullptr;
  }

  //---------------------------------------------------------------------- 

  Boundary * boundary(int i)
  {
    return problem() ? problem()->boundary(i) : nullptr;
  }

  //---------------------------------------------------------------------- 

  Hierarchy * hierarchy()
  {
    return simulation() ? simulation()->hierarchy() : nullptr;
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

  Grouping * field_groups()
  { return field_descr()->groups(); }
    
  //---------------------------------------------------------------------- 

  int define_field (std::string field_name, int cx, int cy, int cz)
  {
    FieldDescr * field_descr = cello::field_descr();
    Config   * config  = (Config *) cello::config();
    if( ! field_descr->is_field( field_name )){
      const int id_field = field_descr->insert_permanent( field_name );
 
      field_descr->set_precision(id_field, config->field_precision);
 
      if ( cx != 0 || cy != 0 || cz != 0 ) {
        field_descr->set_centering(id_field, cx, cy, cz);
      }
    }
    return field_descr->field_id(field_name);
  }
 
  //---------------------------------------------------------------------- 
 
  int define_field_in_group
  (std::string field_name, std::string group_name,
   int cx, int cy, int cz)
  {
    int out = define_field (field_name,cx,cy,cz);
    field_groups()->add(field_name,group_name);
    return out;
  }
 
  //---------------------------------------------------------------------- 
  void finalize_fields ()
  {
    FieldDescr * field_descr = cello::field_descr();
    Config   * config  = (Config *) cello::config();
    field_descr->reset_history(config->field_history);
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

  Grouping * particle_groups()
  { return particle_descr()->groups(); }
    
  //---------------------------------------------------------------------- 

  Output * output(int index)
  {
    return problem() ? problem()->output(index) : nullptr;
  }

  //---------------------------------------------------------------------- 

  Solver * solver(int index)
  {
    return problem() ? problem()->solver(index) : nullptr;
  }

  //---------------------------------------------------------------------- 

  Units * units()
  {
    return problem() ? problem()->units() : nullptr;
  }

  //----------------------------------------------------------------------

  Refresh * refresh(int ir)
  {
    return simulation() ? &simulation()->refresh_list(ir) : nullptr;
  }

  //----------------------------------------------------------------------

  int rank()
  {
    return simulation() ? simulation()->rank() : 1;
  }

  //----------------------------------------------------------------------

  int num_children()
  { return 1 << rank(); }

  //----------------------------------------------------------------------

  int num_children(int r)
  { return (1 << r); }

  //----------------------------------------------------------------------

  size_t num_blocks_process()
  {
    return hierarchy() ? hierarchy()->num_blocks() : 0;
  }
  
  //----------------------------------------------------------------------

  double relative_cell_volume (int level)
  {
    return (1.0/pow(1.0*num_children(),1.0*level));
  }
  
  //----------------------------------------------------------------------

  std::string expand_name
  (
   const std::vector<std::string> * file_format,
   int counter,
   Block * block
   )
  {
    if (file_format->size()==0) return "";
  
    const std::string & name = (*file_format)[0];
  
    const int MAX_BUFFER = 255;

    char buffer[MAX_BUFFER+1];
    char buffer_new[MAX_BUFFER+1];

    // Error check no \% in file name

    ASSERT1 ("cello::expand_name",
             "File name %s cannot contain '\\%%'",
             name.c_str(),
             name.find("\\%") == std::string::npos);

    // Error check variable count equals format conversion specifier count

    std::string rest = name;
    size_t count = 0;
    size_t pos = 0;
    size_t len;
    while ((pos = rest.find("%")) != std::string::npos) {
      count ++;
      len = rest.size();
      rest = rest.substr(pos+1,len-pos-1);
    }

    ASSERT3 ("cello::expand_name",
             "The number of format conversion specifiers %lu "
             "associated with file name %s "
             "must equal the number of variables %lu",
             count, name.c_str(),file_format->size()-1,
             file_format->size()-1 == count);

    // loop through file_format[] from the right and replace 
    // format strings with variable values

    std::string left  = name;
    std::string middle = "";
    std::string right = "";

    for (size_t i=0; i<file_format->size()-1; i++) {

      // visit variables from right to left
      const std::string & arg = (*file_format)[file_format->size() - 1 - i];

      size_t pos = left.rfind("%");
      size_t len = left.size();

      middle = left.substr(pos,len-pos);
      left  = left.substr(0,pos);

      strncpy (buffer, middle.c_str(),MAX_BUFFER);
      if      (arg == "cycle") { sprintf (buffer_new,buffer, block->cycle()); }
      else if (arg == "time")  { sprintf (buffer_new,buffer, block->time()); }
      else if (arg == "count") { sprintf (buffer_new,buffer, counter); }
      else if (arg == "proc")  { sprintf (buffer_new,buffer, CkMyPe()); }
      else if (arg == "flipflop")  { sprintf (buffer_new,buffer, counter % 2); }
      else 
        {
          ERROR3("cello::expand_name",
                 "Unknown file variable #%d '%s' for file '%s'",
                 int(i),arg.c_str(),name.c_str());
        }

      right = std::string(buffer_new) + right;

    }

    return left + right;
  }

  //----------------------------------------------------------------------

  std::string directory
  (
   const std::vector<std::string> * path_format,
   int counter,
   Block * block
   )
  {
    std::string dir = ".";
    std::string name_dir = expand_name(path_format,counter, block);

    // Create subdirectory if any
    if (name_dir != "") {
      dir = name_dir;
      boost::filesystem::path directory(name_dir);
      if (! boost::filesystem::is_directory(directory)) {

        boost::filesystem::create_directory(directory);

        ASSERT1 ("cello::directory()",
                 "Error creating directory %s",
                 name_dir.c_str(),
                 boost::filesystem::is_directory(directory));
      }
    }

    return dir;
  }

  //----------------------------------------------------------------------
}
