// $Id: parameters.hpp 1272 2010-03-08 20:07:32Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef PARAMETERS_PARAMETERS_HPP
#define PARAMETERS_PARAMETERS_HPP

/// @file     parameters_parameters.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Jul  9 15:44:21 PDT 2009
/// @todo     Add output support for individual parameters, e.g. for Monitor output
/// @brief    Declaration for the Parameters class

/// @def      MAX_PARAMETER_FILE_WIDTH
/// @brief    Maximum allowed width of a line in a parameter file
#define MAX_PARAMETER_FILE_WIDTH 255

#include "parameters_Param.hpp"

class Parameters {

  /// @class    Parameters
  /// @ingroup  Parameters
  /// @brief    Read in a parameter file and access parameter values

public: // interface

  /// Get single instance of the Parameters object
  static Parameters * instance() throw ()
  { return & instance_; }

  /// Read in parameters from a file
  void read (FILE * file_pointer);

  /// Write parameters to a file
  void write (FILE * file_pointer);

  /// Return the integer-valued parameter
  int value_integer (std::string , int deflt = 0) 
    throw(ExceptionParametersBadType);

  /// Return the scalar-valued parameter
  double value_scalar (std::string, double deflt = 0.0) 
    throw(ExceptionParametersBadType);

  /// Return the logical-valued parameter
  bool value_logical (std::string , bool deflt = false) 
    throw(ExceptionParametersBadType);

  /// Return the string-valued parameter
  std::string value_string ( std::string , std::string deflt = "") 
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
  std::string list_value_string (int ,std::string , std::string deflt = "")    
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

  /// Get number of groups
  int group_count() throw ();

  /// Get number of groups
  std::string group(int i) throw ();

  /// Get number of subgroups
  int subgroup_count() throw ();

  /// Get number of subgroups
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

private: // functions

  /// Create an empty Parameters object (singleton design pattern)
  Parameters() throw();

  /// Copy constructor
  Parameters(const Parameters & parameters) throw();

  /// Assignment operator

  Parameters & operator= (const Parameters & parameters) throw();

  /// Delete a Parameters object (singleton design pattern)
  ~Parameters();

  /// Read in the next line of the input file
  int readline_ (FILE* fp, char * buffer, int n) throw();

  /// Return the Param pointer for the specified parameter
  Param * parameter_ (std::string parameter)
  {
    std::string p = current_group_ + ":" + current_subgroup_ + ":" + parameter;
    return parameter_map_[p];
  };

  /// Return the Param pointer for the specified list parameter element
  Param * list_element_ (std::string parameter, int index)
  {
    Param * list = parameter_(parameter);
    Param * param = NULL;
    int list_length = list->value_list_->size();
    if (list != NULL && 0 <= index && index < list_length ) {
      param =  (*(list->value_list_))[index];
    }
    return param;
  }

  void monitor_log (std::string parameter) throw();

private: // attributes

  /// Single instance of the Parameters object (singleton design pattern)
  static Parameters instance_;

  std::string current_group_;
  std::string current_subgroup_;

  std::map<std::string, Param *>  parameter_map_;
  ParamNode                     * parameter_tree_;

};


extern "C" { 
  /// C function for reading parameters from a file
  struct param_type * cello_parameters_read(FILE *);
  /// C function for printing parameters to stdout
  void cello_parameters_print();
}

#endif /* PARAMETERS_PARAMETERS_HPP */

