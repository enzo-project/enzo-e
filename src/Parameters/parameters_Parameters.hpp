// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PARAMETERS_PARAMETERS_HPP
#define PARAMETERS_PARAMETERS_HPP

/// @file     parameters_Parameters.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul  9 15:44:21 PDT 2009
/// @todo     Add output support for individual parameters, e.g. for Monitor output
/// @todo     set_list(), etc functions for lists and expressions
/// @todo     set_value() using string values for all types
/// @todo     assert_required() to make given parameter required (don't exit since more than one required parameter may be missing)
/// @todo     Add "check()" function to check individual parameters, or all
/// @todo     Move Group and subgroup to parameter lists at end with "" default
/// @todo     Convert current_group_ to std::stack<std::string>
/// @bug      "x - 0.5" broken since intepreted as "x (-0.5)"; workaround "x - (0.5)"
/// @bug      Evaluating logical expressions with parentheses seg-faults (fixed?)
/// @brief    [\ref Parameters] Declaration for the Parameters class

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

  /// Get single instance of the Parameters object
//   static Parameters * instance() throw ()
//   { return & instance_; }

  /// Read in parameters from a file
  void read (const char * file_name);

  /// Write parameters to a file
  void write (const char * file_name);

  // /// Return the parameter value of specified type
  // void value (std::string, parameter_enum type, 
  // 	      void * value, 
  // 	      void * deflt = 0);

  /// Return the integer-valued parameter
  int value_integer (std::string , int deflt = 0) 
    throw(ExceptionParametersBadType);

  void set_integer ( std::string parameter, int value ) 
    throw(ExceptionParametersBadType);

  /// Return the scalar-valued parameter
  double value_scalar (std::string, double deflt = 0.0) 
    throw(ExceptionParametersBadType);

  void set_scalar ( std::string parameter, double value ) 
    throw(ExceptionParametersBadType);

  /// Return the logical-valued parameter
  bool value_logical (std::string , bool deflt = false) 
    throw(ExceptionParametersBadType);

  void set_logical ( std::string parameter, bool value ) 
    throw(ExceptionParametersBadType);

  /// Return the string-valued parameter
  const char * value_string ( std::string , const char * deflt = "") 
    throw(ExceptionParametersBadType);

  void set_string ( std::string parameter, const char * value ) 
    throw(ExceptionParametersBadType);

  /// Evaluate the scalar-valued parameter expression
  void evaluate_scalar 
  (
   std::string parameter,
   int         n, 
   double    * result, 
   double    * deflt,
   double    * x, 
   double    * y, 
   double    * z, 
   double    * t)
    throw(ExceptionParametersBadType);

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
    throw(ExceptionParametersBadType);

  /// Return the length of the list parameter
  int list_length (std::string parameter);

  /// Access an integer list element
  int list_value_integer (int , std::string , int deflt = 0)    
    throw(ExceptionParametersBadType);

  /// Access a scalar list element
  double list_value_scalar (int , std::string , double deflt = 0.0)    
    throw(ExceptionParametersBadType);

  /// Access a logical list element
  bool list_value_logical (int ,std::string , bool deflt = false)    
    throw(ExceptionParametersBadType);

  /// Access a string list element
  const char * list_value_string (int ,std::string , const char * deflt = "")    
    throw(ExceptionParametersBadType);

  /// Evaluate the scalar-valued list element expression
  void list_evaluate_scalar 
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
    throw(ExceptionParametersBadType);

  /// Evaluate the scalar-valued list element expression
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
    throw(ExceptionParametersBadType);

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
  void set_group(int, std::string) throw();

  //--------------------------------------------------

  /// Return the type of the given parameter
  parameter_enum type(std::string) throw();

  /// Return the type of the given parameter
  parameter_enum list_type(int, std::string) throw();

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

  void monitor_access_ (std::string parameter,
		      std::string deflt_string,
		      int index=-1) throw();
  void monitor_write_ (std::string parameter) throw();

  /// Create a new parameter with the given grouping
  void new_param_ ( char * group[],
		    std::string parameter,
		    Param * param ) throw();

private: // attributes

  /// Single instance of the Parameters object (singleton design pattern)
//   static Parameters instance_;

  /// Stack of current grouping
  char * current_group_[MAX_GROUP_DEPTH];

  /// Top of the current_group_ stack
  int current_group_depth_;

  std::map<std::string, Param *>  parameter_map_;
  ParamNode                     * parameter_tree_;

  /// Monitor object for parameters
  Monitor * monitor_; 

};


extern "C" { 
  /// C function for reading parameters from a file
  struct param_struct * cello_parameters_read(FILE *);
  /// C function for printing parameters to stdout
  void cello_parameters_print();
}

#endif /* PARAMETERS_PARAMETERS_HPP */

