// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_ParameterGroup.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Jan 17 2023
/// @brief    [\ref Parameters] Declaration for the ParameterGroup class

#ifndef PARAMETERS_PARAMETER_ACCESSOR_HPP
#define PARAMETERS_PARAMETER_ACCESSOR_HPP

class ParameterGroup {

  /// @class    ParameterGroup
  /// @ingroup  Parameters
  /// @brief    [\ref Parameters] Acts as an interface for accessing parameters
  ///           within a parameter group.
  ///
  /// A lightweight wrapper around a Parameters object that provides methods
  /// for accessing parameters in a parameter group, that are stored within the
  /// wrapped object.
  ///
  /// All parameters in a group, all share a common root "parameter_path"
  /// (defined below). The common root "parameter_path" is specified when an
  /// instance of this object is initialized. When a string is passed to one of
  /// the accessor methods, the string is internally appended to the end of the
  /// root "parameter_path" and the result represents the full name of the
  /// queried parameter.
  ///
  /// Instances of this class are intended to be used when initializing objects
  /// in other software components (an instance of this class should be passed
  /// to a factory-method or constructor).
  ///
  /// This class was designed with a few guiding principles:
  ///     1. restricting access to parameters within the associated root path
  ///        This is to discourage design of objects that are configured by
  ///        parameters scattered throughout the parameter file. In rare cases
  ///        (e.g. deprecating a parameter), exceptions need to be made. Thus,
  ///        an "escape-hatch" is provided to directly access the wrapped
  ///        Parameters object.
  ///     2. providing an explicit way to signal that all parameters accessed
  ///        through a given instance of this class share a single (queryable &
  ///        immutable) root-path without requiring the root path to be
  ///        explicitly specified as part of every parameter-path.
  ///        - This capability is very useful when multiple instances of the
  ///          same class (e.g. a Method class) needs to be initialized but
  ///          each instance has different configurations.
  ///        - In that scenario, each instance of the class is typically
  ///          initialized from parameters where the last section(s) of the
  ///          parameter-paths are unchanged, but the root path differs.
  ///     3. We do NOT allow the common root-path to be mutated. Thus, if
  ///        an instance is passed to a helper function, you can always be
  ///        confident that the root-path is unchanged by the helper function
  ///        (without needing to check the implementation of that function).
  ///        If we are ever tempted to mutate the root-path, we should just
  ///        initialize a new ParameterGroup (since they are lightweight)
  ///
  /// @par "parameter_path"
  /// Cello uses a "hierarchical" parameter file: the parameters themselves are
  /// essentially leaf nodes in a tree-hierarchy of "groups" (like how a file
  /// is organized in a directory hierarchy). A parameter_path unambiguously
  /// identifes a parameter in this hierarchy. A given parameter_path lists the
  /// names of ancestor "groups", separated by colons, and lists the name
  /// of the parameter at the end. An example parameter_path is
  /// ``"Method:null:dt"``

public:

  /// construct a new ParameterGroup object
  ParameterGroup(Parameters &p, const std::string& root_parameter_path)
    : wrapped_p_(p),
      root_parameter_path_(root_parameter_path)
  {
    ASSERT("ParameterGroup::ParameterGroup",
           "root_parameter_path must not have a trailing colon",
           root_parameter_path.back() != ':');
  }

  // default implementations of copy and move constructors
  ParameterGroup(const ParameterGroup&) = default;
  ParameterGroup(ParameterGroup&&) = default;

  // the following operations are deleted (because the class holds a reference
  // and a const member as attributes)
  ParameterGroup() = delete;
  ParameterGroup& operator=(const ParameterGroup&) = delete;
  ParameterGroup& operator=(ParameterGroup&&) = delete;

  /// query the path for the represented group of parameters
  ///
  /// For the sake of clarity, we highlight a few examples:
  /// - if ``*this`` represents the group containing the parameters with paths
  ///   "foo:bar:par_1" and "foo:bar:par_2", then this method would return
  ///   ``"foo:bar``".
  /// - if ``*this`` represents the group containing the parameters with paths
  ///   "foo:par_A" and "foo:par_B", then this method would return ``"foo"``
  /// - if ``*this`` represents the group containing the parameters with paths
  ///   "par_alpha" and "par_beta", then this method would return ``""``
  ///   (Note: by convention, parameters are never actually put in this group)
  const std::string& get_group_path() const noexcept
  {return root_parameter_path_;}

  int value (std::string s, int deflt) noexcept
  { return value_integer(s,deflt); }
  double value (std::string s, double deflt) noexcept
  { return value_float(s,deflt); }
  bool value (std::string s, bool deflt) noexcept
  { return value_logical(s, deflt); }
  std::string value (std::string s, std::string deflt) noexcept
  { return value_string(s, deflt); }

  int value (int i,std::string s, int deflt) noexcept
  { return list_value_integer(i,s,deflt); }
  double value (int i,std::string s, double deflt) noexcept
  { return list_value_float(i,s,deflt); }
  bool value (int i,std::string s, bool deflt) noexcept
  { return list_value_logical(i,s,deflt); }
  std::string value (int i,std::string s, const char * deflt) noexcept
  { return list_value_string(i,s,deflt); }
  std::string value (int i,std::string s, std::string deflt) noexcept
  { return list_value_string(i,s,deflt); }

  /// Return the type of the given parameter
  parameter_type type(std::string param) noexcept;

  /// Return the Param pointer for the specified parameter
  Param * param (std::string parameter);

  /// Return the integer-valued parameter
  int value_integer (std::string s, int deflt = 0) noexcept;
  /// Return the floating-point valued parameter
  double value_float (std::string s, double deflt = 0.0) noexcept;
  /// Return the logical-valued parameter
  bool value_logical (std::string s, bool deflt = false) noexcept;
  /// Return the string-valued parameter
  std::string value_string ( std::string s, std::string deflt = "") noexcept;

  /// Return the length of the list parameter
  int list_length (std::string parameter);
  /// Access an integer list element
  int list_value_integer (int i, std::string s, int deflt = 0) noexcept;
  /// Access a floating point list element
  double list_value_float (int i, std::string s, double deflt = 0.0) noexcept;
  /// Access a logical list element
  bool list_value_logical (int i, std::string s, bool deflt = false) noexcept;
  /// Access a string list element
  std::string list_value_string (int, std::string, std::string d="") noexcept;

  /// Return the full name of the parameter (including the root parameter path)
  std::string full_name(const std::string& parameter) const noexcept
  { return root_parameter_path_ + ":" + parameter; }

  /// Returns a vector holding the names of all leaf parameters that share the
  /// root parameter path encapsulated by this object
  std::vector<std::string> leaf_parameter_names() const noexcept
  { return wrapped_p_.leaf_parameter_names(root_parameter_path_); }

private:

  std::vector<std::string> pop_wrapped_p_groups_()
  {
    const int n = wrapped_p_.group_depth();
    std::vector<std::string> grps(n);
    for (int i = 0; i < n; i++) { grps[i] = wrapped_p_.group(i); }
    wrapped_p_.group_clear();
    return grps;
  }

  void restore_wrapped_p_groups_(const std::vector<std::string>& groups)
  {
    wrapped_p_.group_clear();
    for (const std::string& grp : groups) { wrapped_p_.group_push(grp); }
  }

private: // attributes
  /// the wrapped Parameters object
  ///
  /// the Parameters object is implicitly assumed to outlive the instance
  /// holding this reference
  Parameters &wrapped_p_;

  /// The associated root parameter path. All parameters in a given group share
  /// this path
  ///
  /// This will never have a trailing colon.
  ///
  /// An invariant of this class is that this will NOT change. If we're ever
  /// tempted to allow this attribute to change, we should prefer creation of a
  /// new ParameterGroup instance (since they are light)
  const std::string root_parameter_path_;
};

#endif /* PARAMETERS_PARAMETER_ACCESSOR_HPP */
