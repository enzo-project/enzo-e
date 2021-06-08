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
class Monitor;
class Parameters {

  /// @class    Parameters
  /// @ingroup  Parameters
  /// @brief    [\ref Parameters] Read in a parameter file and access
  /// parameter values
  friend Expression; // This is only a friend class so that it can call the
                     // construct_or_rebuild_Expression_ method

public: // interface

  /// Create an empty Parameters object
  Parameters(Monitor * monitor = 0) throw();

  /// Create a new Parameters object and read parameters from the given file
  Parameters(const char * file_name, 
	     Monitor * monitor = 0) throw();

  /// Copy constructor
  Parameters(const Parameters & parameters) throw();

  /// Assignment operator

  Parameters & operator= (const Parameters & parameters) throw();

  /// Delete a Parameters object (singleton design pattern)
  ~Parameters();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Read in parameters from a file
  void read (const char * file_name);
  /// Write parameters to a file
  void write (const char * file_name, int write_type = param_write_cello);
  void write (FILE * fp, int write_type = param_write_cello);

  int value (std::string s, int deflt) throw()
  { return value_integer(s,deflt); }

  double value (std::string s, double deflt) throw()
  { return value_float(s,deflt); }

  bool value (std::string s, bool deflt) throw()
  { return value_logical(s,deflt); }

  std::string value (std::string s, std::string deflt) throw()
  { return value_string(s,deflt); }

  int value (int i,std::string s, int deflt) throw()
  { return list_value_integer(i,s,deflt); }

  double value (int i,std::string s, double deflt) throw()
  { return list_value_float(i,s,deflt); }

  bool value (int i,std::string s, bool deflt) throw()
  { return list_value_logical(i,s,deflt); }

  std::string value (int i,std::string s, const char * deflt) throw()
  { return list_value_string(i,s,deflt); }
  std::string value (int i,std::string s, std::string deflt) throw()
  { return list_value_string(i,s,deflt); }

  /// Return the integer-valued parameter
  int value_integer (std::string , int deflt = 0) throw();
  /// Return the floating-point valued parameter
  double value_float (std::string, double deflt = 0.0) throw();
  /// Return the logical-valued parameter
  bool value_logical (std::string , bool deflt = false) throw();
  /// Return the string-valued parameter
  std::string value_string ( std::string , std::string deflt = "") throw();
  /// Return the Expression-valued parameter
  Expression value_Expression ( std::string ) throw();

  /// Return the length of the list parameter
  int list_length (std::string parameter);
  /// Access an integer list element
  int list_value_integer (int , std::string , int deflt = 0) throw();
  /// Access a floating point list element
  double list_value_float (int , std::string , double deflt = 0.0) throw();
  /// Access a logical list element
  bool list_value_logical (int ,std::string , bool deflt = false) throw();
  /// Access a string list element
  std::string list_value_string (int,std::string, std::string d= "") throw();
  /// Return an Expression list element
  Expression list_value_Expression ( int, std::string ) throw();

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

  /// Returns a vector holding the names of all leaf parameters in the current
  /// group
  std::vector<std::string> leaf_parameter_names() const throw();

  /// Return the full name of the parameter including group
  std::string full_name (std::string parameter) const throw();

  /// Return the type of the given parameter
  parameter_type type(std::string) throw();

  /// Return the type of the given parameter
  parameter_type list_type(int, std::string) throw();

  /// Set whether to output 
  void set_monitor (bool lmonitor) { lmonitor_ = lmonitor; };

  void check();

  /// Return the Param pointer for the specified parameter
  Param * param (std::string parameter, int index = -1);

  //--------------------------------------------------

private: // functions

  /// Read in the next line of the input file
  int readline_ (FILE* fp, char * buffer, int n) throw();

  /// Return the full parameter name Group:group:group:...:parameter
  std::string parameter_name_(std::string parameter) const throw()
  {
    bool is_full_parameter = (parameter.find(":") != std::string::npos);
    return (is_full_parameter ? parameter : this->full_name (parameter));
  }

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

  /// Return the pointer to the const Param for the specified parameter
  ///
  /// @param[in]  parameter  Name of the parameter
  /// @param[in]  index      List index of the parameter. If the parameter is
  ///     not in a list, then this should be -1.
  /// @returns This returns a pair of pointers to const Param objects. The
  ///     first entry is a pointer to the const Param for the specified
  ///     parameter (or a nullptr if it can't be found). When ``index > -1``,
  ///     the second entry is a pointer to the constant ``Param`` object
  ///     representing a list that holds the specified parameter (this is a
  ///     nullptr if it can't be found or ``index == -1``)
  std::pair<const Param*, const Param*> const_param_(std::string parameter,
						     int index = -1) const;

  /// This is used to construct or rebuild Expression objects
  Expression construct_or_rebuild_Expression_
  (const std::string & parameter_name, int param_index = -1) const throw();

  //--------------------------------------------------

private: // attributes

  /// Stack of current grouping
  std::vector <std::string> current_group_;

  /// Map parameter name to Param object
  std::map<std::string, Param *>  parameter_map_;

  /// Parameters represented as a tree with groups as internal nodes
  ParamNode                     * parameter_tree_;

  /// Monitor object for parameters
  Monitor * monitor_; 

  /// Whether monitor should output accessed parameters when requested
  bool lmonitor_;
};

//----------------------------------------------------------------------

extern "C" { 
  /// C function for reading parameters from a file
  struct param_struct * cello_parameters_read(const char *, FILE *);
  /// C function for printing parameters to stdout
  void cello_parameters_print();
}

extern Parameters g_parameters;

//----------------------------------------------------------------------

#endif /* PARAMETERS_PARAMETERS_HPP */

