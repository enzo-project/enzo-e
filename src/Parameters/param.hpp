// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PARAM_HPP
#define PARAM_HPP

/// @file     param.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Sun Oct 11 14:55:25 PDT 2009
/// @todo     Rename print_expression to match cello_parameters_print_list
/// @brief    Interface for the Param class

#include <vector>
#include <string>

//----------------------------------------------------------------------

/// @fun print_expression
/// @brief Print a parameter expression
extern "C" { 
  void print_expression(struct node_expr * node,
		      FILE * fp = stdout);
}

/// @fun cello_parameters_print_list
/// @brief Print a parameter list
extern "C" { 
  void cello_parameters_print_list(struct param_type * head, int level);
}

//----------------------------------------------------------------------

class Param {

  /// @class    Param
  /// @ingroup  Parameter
  /// @brief    Represents and evaluates various parameter types and expressions

  friend class Parameters;

public: // interface

  /// Initialize a Param object
  Param () {};

  /// Delete a Param object
  ~Param () { dealloc_(); };

  /// Evaluate a scalar expression given vectos x,y,z,t
  void evaluate_scalar  
  ( struct node_expr * node, 
    int                n, 
    double *           result, 
    double *           x, 
    double *           y, 
    double *           z, 
    double *           t);

  /// Evaluate a logical expression given vectos x,y,z,t
  void evaluate_logical  
  ( struct node_expr * node, 
    int                n, 
    bool *             result, 
    double *           x, 
    double *           y, 
    double *           z, 
    double *           t);

  /// Set the parameter type and value
  void set(struct param_type * param);

  /// Write the parameter to the file
  void write(FILE * file_pointer,
	     std::string parameter);

  /// Return whether the parameter is an integer
  bool is_integer()      { return type_ == type_integer_; };

  /// Return whether the parameter is a scalar
  bool is_scalar()       { return type_ == type_scalar_; };

  /// Return whether the parameter is a logical
  bool is_logical()      { return type_ == type_logical_; };

  /// Return whether the parameter is a string
  bool is_string()       { return type_ == type_string_; };

  /// Return whether the parameter is a list
  bool is_list()         { return type_ == type_list_; };

  /// Return whether the parameter is a scalar expression
  bool is_scalar_expr()  { return type_ == type_scalar_expr_; };

  /// Return whether the parameter is a logical expression
  bool is_logical_expr() { return type_ == type_logical_expr_; };

private: // attributes

  /// @enum type_param
  /// @brief Parameter type 
  enum type_param {
    type_unknown_,
    type_integer_,
    type_scalar_,
    type_logical_,
    type_string_,
    type_list_,
    type_scalar_expr_,
    type_logical_expr_
  };

  /// Parameter type
  enum type_param type_;

  /// Type definition for a list of parameters
  typedef std::vector<class Param *> list_type;

  /// Value of the parameter, with type depending on type_
  union {
    int                value_integer_; 
    double             value_scalar_; 
    bool               value_logical_; 
    char *             value_string_;
    list_type        * value_list_;
    struct node_expr * value_expr_;
  };

private: // functions

  /// Set an integer parameter
  void set_integer_ (int value)    
  { 
    type_ = type_integer_; 
    value_integer_ = value; 
  };

  /// Set a scalar parameter (note that scalar parameters are doubles)
  void set_scalar_  (double value) 
  { 
    type_ = type_scalar_; 
    value_scalar_ = value; 
  };

  /// Set a logical parameter
  void set_logical_ (int value)    
  { 
    type_ = type_logical_; 
    value_logical_ = (value != 0); 
  };

  /// Set a string parameter (note that string is aliased)
  void set_string_ (char * value) 
  { 
    type_ = type_string_; 
    value_string_ = value;
  };

  /// Set a list parameter given a list of parameter values
  ///
  /// Lists are bounded by "sentinel" types: the list values are
  /// taken to be between the first sentinel and the next sentinel
  void set_list_ (struct param_type * value) 
  { 
    type_ = type_list_; 
    value_list_ = new list_type;
    value = value->next; // Skip sentinel
    while (value->type != enum_parameter_sentinel) {
      Param * param = new Param;
      param->set(value);
      (*value_list_).push_back(param);
      value = value->next;
    }
  };

  /// Set a scalar expression parameter
  void set_scalar_expr_ (struct node_expr * value)
  { 
    type_ = type_scalar_expr_;
    value_expr_ = value; 
  };

  /// Set a logical expression parameter
  void set_logical_expr_ (struct node_expr * value)
  { 
    type_ = type_logical_expr_;
    value_expr_ = value; 
  };

  /// Deallocate the parameter
  void dealloc_();

  /// Deallocate a string
  void dealloc_string_() { free (value_string_); } 

  /// Deallocate a list of parameters
  void dealloc_list_     (list_type *value_list_);

  /// Deallocate an expression parameter
  void dealloc_node_expr_ (struct node_expr * p);

  ///
  void write_scalar_expr_(FILE * file_pointer,
			  struct node_expr * value_expr_);
};

//----------------------------------------------------------------------

#endif
