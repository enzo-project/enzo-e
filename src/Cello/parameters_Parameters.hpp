// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Parameters.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul  9 15:44:21 PDT 2009
/// @brief    [\ref Parameters] Declaration for the Parameters class

#ifndef PARAMETERS_PARAMETERS_HPP
#define PARAMETERS_PARAMETERS_HPP

/// @def      MAX_PARAMETER_FILE_WIDTH
/// @brief    Maximum allowed width of a line in a parameter file
#define MAX_PARAMETER_FILE_WIDTH 255

class Parameters {

  /// @class    Parameters
  /// @ingroup  Parameters
  /// @brief    [\ref Parameters] Read in a parameter file and access
  /// parameter values

public: // interface


  /// Create an empty Parameters object
  Parameters(Monitor * monitor = 0) throw();

  /// Create a new Parameters object and read parameters from the given file
  Parameters(const char * file_name, Monitor * monitor = 0) throw();

  /// Copy constructor
  Parameters(const Parameters & parameters) throw();

  /// Assignment operator

  Parameters & operator= (const Parameters & parameters) throw();

  /// Delete a Parameters object (singleton design pattern)
  ~Parameters();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    WARNING("Parameters::pup","current_group_ not pup'ed");
    //    PUParray(p,current_group_,MAX_GROUP_DEPTH);
    p | current_group_depth_;
    WARNING("Parameters::pup","parameter_map_ not pup'ed");
    //    p | parameter_map_;
    p | *parameter_tree_;
    p | *monitor_; 

  }
#endif

  /// Read in parameters from a file
  void read (const char * file_name);
  /// Write parameters to a file
  void write (const char * file_name);


  /// Return the integer-valued parameter
  int value_integer (std::string , int deflt = 0) throw();
  /// Return the floating-point valued parameter
  double value_float (std::string, double deflt = 0.0) throw();
  /// Return the logical-valued parameter
  bool value_logical (std::string , bool deflt = false) throw();
  /// Return the string-valued parameter
  const char * value_string ( std::string , const char * deflt = "") throw();

  /// Return the length of the list parameter
  int list_length (std::string parameter);
  /// Access an integer list element
  int list_value_integer (int , std::string , int deflt = 0) throw();
  /// Access a floating point list element
  double list_value_float (int , std::string , double deflt = 0.0) throw();
  /// Access a logical list element
  bool list_value_logical (int ,std::string , bool deflt = false) throw();
  /// Access a string list element
  const char * list_value_string (int,std::string,const char * d= "") throw();

  /// Assign a value to the integer-valued parameter
  void set_integer ( std::string parameter, int value ) throw();
  /// Assign a value to the floating-point valued parameter
  void set_float ( std::string parameter, double value ) 
    throw();
  /// Assign a value to the logical-valued parameter
  void set_logical ( std::string parameter, bool value ) 
    throw();
  /// Assign a value to the string-valued parameter
  void set_string ( std::string parameter, const char * value ) 
    throw();

  /// Return the length of the list parameter
  void set_list_length (std::string parameter, int length);
  /// Assign a value to an integer-valued list element parameter
  void set_list_integer (int , std::string , int value) throw();
  /// Assign a value to the floating-point list-element parameter
  void set_list_float (int , std::string , double value) throw();
  ///  Assign a value to the logical-valued list-element parameter
  void set_list_logical (int ,std::string , bool value) throw();
  /// Assign a value to the string-valued list-element parameter
  void set_list_string (int,std::string,const char * value) throw();
  void done_set_list ( std::string parameter ) throw()
  { monitor_write_(parameter); }

  /// Evaluate the floating-point valued parameter expression
  void evaluate_float 
  (
   std::string parameter,
   int         n, 
   double    * result, 
   double    * deflt,
   double    * x, 
   double    * y, 
   double    * z, 
   double    * t)
    throw();

  /// Evaluate the logical-valued parameter expression
  void evaluate_logical 
  (
   std::string parameter,
   int         n, 
   bool      * result, 
   bool      * deflt,
   double    * x, 
   double    * y, 
   double    * z, 
   double    * t)
    throw();

  /// Evaluate the floating-point valued list element expression
  void list_evaluate_float 
  (
   int ,
   std::string parameter,
   int         n, 
   double    * result, 
   double    * deflt,
   double    * x, 
   double    * y, 
   double    * z, 
   double    * t
   )    
    throw();

  /// Evaluate the logical-valued list element expression
  void list_evaluate_logical
  (
   int ,
   std::string parameter,
   int         n, 
   bool      * result, 
   bool      * deflt,
   double    * x, 
   double    * y, 
   double    * z, 
   double    * t
   )    
    throw();

  //--------------------------------------------------
  // PARAMETER GROUPS

  /// Return the ith group in the grouping
  std::string group(int i) const throw();

  /// Return the depth of groups
  int group_depth() const throw();

  /// Return the number of subgroups in the current grouping
  int group_count() const throw();

  /// Push the group onto the current grouping
  void group_push(std::string) throw();

  /// Replace the top-level group in the grouping with the specified group;
  /// equivalent to a pop() then push(group)
  void group_replace(std::string group) throw() {
    group_pop();
    group_push(group);
  }

  /// pop the head from the current grouping, checking that it matches
  /// the optionally provided group name
  void group_pop(std::string group = "") throw();

  /// Set the ith group to specified group, clearing all higher groups
  void group_set(int, std::string) throw();

  /// Clear all groups
  void group_clear() throw();
  
  /// Return the type of the given parameter
  parameter_type type(std::string) throw();

  /// Return the type of the given parameter
  parameter_type list_type(int, std::string) throw();

  //--------------------------------------------------

private: // functions

  /// Read in the next line of the input file
  int readline_ (FILE* fp, char * buffer, int n) throw();

  /// Return the full parameter name Group:group:group:...:parameter
  std::string parameter_name_(std::string parameter)
  {
    bool is_full_parameter = (parameter.find(":") != std::string::npos);

    std::string parameter_name;

    if (! is_full_parameter) {
      parameter_name = "";
      for (int i=0; current_group_[i] != 0 && i<MAX_GROUP_DEPTH; i++) {
	parameter_name = parameter_name + current_group_[i] + ":";
      }
      parameter_name = parameter_name + parameter;
    } else {
      parameter_name = parameter;
    }
    return parameter_name;
  }

  /// Return the Param pointer for the specified parameter
  Param * parameter_ (std::string parameter)
  {
    return parameter_map_[parameter_name_(parameter)];
  };

  /// Return the Param pointer for the specified list parameter element
  Param * list_element_ (std::string parameter, int index) throw();

  /// Display through monitor that a parameter was accessed
  void monitor_access_ (std::string parameter,
			std::string deflt_string,
			int index=-1) throw();

  /// Display through monitor that a parameter was assigned a value
  void monitor_write_ (std::string parameter) throw();

  /// Create a new parameter with the given grouping and name
  void new_param_ ( std::string full_parameter,
		    Param * param ) throw();

  size_t extract_groups_( const std::string parameter, std::string * group);

  //--------------------------------------------------

private: // attributes

  /// Stack of current grouping
  char * current_group_[MAX_GROUP_DEPTH];

  /// Top of the current_group_ stack
  int current_group_depth_;

  /// Map parameter name to Param object
  std::map<std::string, Param *>  parameter_map_;

  /// Parameters represented as a tree with groups as internal nodes
  ParamNode                     * parameter_tree_;

  /// Monitor object for parameters
  Monitor * monitor_; 

};

//----------------------------------------------------------------------

extern "C" { 
  /// C function for reading parameters from a file
  struct param_struct * cello_parameters_read(const char *, FILE *);
  /// C function for printing parameters to stdout
  void cello_parameters_print();
}

#endif /* PARAMETERS_PARAMETERS_HPP */

