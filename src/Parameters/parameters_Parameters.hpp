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
  Parameters() throw();

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
  void read (FILE * file_pointer);

  /// Write parameters to a file
  void write (FILE * file_pointer);

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

  /// Return the number of groups
  int group_count() throw ();

  /// Return the name of the ith group
  std::string group(int i) throw ();

  /// Return the number of subgroups in the current group
  int subgroup_count() throw ();

  /// Return the ith subgroup in the current group
  std::string subgroup(int i) throw ();

  /// Set the current group.  Clears current subgroup
  void set_current_group  (std::string group, std::string subgroup = "") throw ()
  { 
    current_group_    = group;
    current_subgroup_ = subgroup;
  };

  /// Set the current subgroup
  void set_current_subgroup  (std::string subgroup) throw ()
  {
    current_subgroup_ = subgroup;
  }

  /// Get the current group
  std::string current_group () throw ()
  { return current_group_; };

  /// Get the current subgroup
  std::string current_subgroup () throw ()
  { return current_subgroup_; };

  /// Return the type of the given parameter
  parameter_enum type(std::string) throw();

  /// Return the type of the given parameter
  parameter_enum list_type(int, std::string) throw();

private: // functions

  /// Read in the next line of the input file
  int readline_ (FILE* fp, char * buffer, int n) throw();

  /// Return the Param pointer for the specified parameter
  Param * parameter_ (std::string parameter)
  {
    std::string p = current_group_ + ":" + current_subgroup_ + ":" + parameter;
    return parameter_map_[p];
  };

  /// Return the Param pointer for the specified list parameter element
  Param * list_element_ (std::string parameter, int index) throw();

  void monitor_read_ (std::string parameter,
		      std::string deflt_string,
		      int index=-1) throw();
  void monitor_write_ (std::string parameter) throw();

  void new_param_ ( std::string group,
		    std::string subgroup,
		    std::string parameter,
		    Param * param ) throw();

private: // attributes

  /// Single instance of the Parameters object (singleton design pattern)
//   static Parameters instance_;

  /// Current Group
  std::string current_group_;

  /// Current subgroup
  std::string current_subgroup_;

  std::map<std::string, Param *>  parameter_map_;
  ParamNode                     * parameter_tree_;

};


extern "C" { 
  /// C function for reading parameters from a file
  struct param_struct * cello_parameters_read(FILE *);
  /// C function for printing parameters to stdout
  void cello_parameters_print();
}

#endif /* PARAMETERS_PARAMETERS_HPP */

