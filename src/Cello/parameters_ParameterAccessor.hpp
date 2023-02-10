// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_ParameterAccessor.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Jan 17 2023
/// @brief    [\ref Parameters] Declaration for the ParameterAccessor class

#ifndef PARAMETERS_PARAMETER_ACCESSOR_HPP
#define PARAMETERS_PARAMETER_ACCESSOR_HPP

class ParameterAccessor {

  /// @class    ParameterAccessor
  /// @ingroup  Parameters
  /// @brief    [\ref Parameters] Acts as an interface for accessing parameters
  ///
  /// A lightweight wrapper around a Parameters object that provides access
  /// access to parameters stored by that object, which share a common root
  /// "parameter_path" (defined below). Instances of this class are intended
  /// to be used when initializing objects in other software components. The
  /// two guiding principles include:
  ///     1. restricting access to parameters within the associated root path
  ///        This is to discourage design of objects that are configured by
  ///        parameters scattered throughout the parameter file. In rare cases
  ///        (e.g. deprecating a parameter), exceptions need to be made. Thus,
  ///        an "escape-hatch" is provided to directly access the wrapped
  ///        Parameters object.
  ///     2. providing an explicit/transparent way to alias the root-path's
  ///        name. This can be useful when multiple instances of the same class
  ///        (e.g. a Method class), but each instance has different
  ///        configurations.
  ///
  /// @par "parameter_path"
  /// Cello uses a "hierarchical" parameter file: the parameters themselves are
  /// essentially leaf nodes in a tree-hierarchy of "groups" (like how a file
  /// is organized in a directory hierarchy). A parameter_path unambiguously
  /// identifes a parameter in this hierarchy. A given parameter_path lists the
  /// names of ancestor "groups", separated by colons, and lists the name
  /// of the parameter at the end. An example parameter_path is
  /// ``"Method:null:dt"``
  ///
  /// @note
  /// alternative names for this class include ParameterClient, ParameterView,
  /// ParameterTerminal... They reflect the fact that the class provides access
  /// to the parameters held in a centralized Parameters object
public:

  /// construct a new ParameterAccessor object
  ParameterAccessor(Parameters &p, const std::string& root_parameter_path);

  // default implementations of copy and move constructors
  ParameterAccessor(const ParameterAccessor&) = default;
  ParameterAccessor(ParameterAccessor&&) = default;

  // the following operations are deleted (because the class holds a reference
  // and a const member as attributes)
  ParameterAccessor() = delete;
  ParameterAccessor& operator=(const ParameterAccessor&) = delete;
  ParameterAccessor& operator=(ParameterAccessor&&) = delete;

  const std::string& get_root_parpath() const noexcept
  {return root_parameter_path;}

  void set_temporary_alias_root_path(const std::string& path);


  int value (std::string s, int deflt) noexcept;
  double value (std::string s, double deflt) noexcept;
  bool value (std::string s, bool deflt) throw();
  std::string value (std::string s, std::string deflt) throw();

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
  

  /// Only use in case of emergencies
  Parameters& wrapped_Parameters_ref() noexcept { return wrapped_p; }

private:

  /// the wrapped Parameters object
  ///
  /// the Parameters object is implicitly assumed to outlive the instance
  /// holding this reference
  Parameters &wrapped_p;

  /// The associated root parameter path.
  ///
  /// This will never have a trailing colon or be empty.
  ///
  /// An invariant of this class is that this will NOT change. If we're ever
  /// tempted to allow this attribute to change, we should prefer creation of a
  /// new ParameterAccessor instance (since they are light)
  const std::string root_parameter_path;

  /// This is a user-configurable attribute that serves as an alias in their
  /// specified paths for the root_parameter_path.
  ///
  /// An empty string means that it's unset. This never has a trailing colon
  std::string alias_root_path;
};

#endif /* PARAMETERS_PARAMETER_ACCESSOR_HPP */
