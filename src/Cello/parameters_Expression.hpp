// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_Expression.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri Jun 4 2021
/// @brief    [\ref Parameters] Declaration for the Expression class


#ifndef PARAMETERS_EXPRESSIONS_HPP
#define PARAMETERS_EXPRESSIONS_HPP
class Expression {

  /// @class    Expression
  /// @ingroup  Parameters
  /// @brief    [\ref Parameters] Represents scalar or logical expressions
  ///
  /// To conserve memory, this acts as a glorified pointer to struct node_expr,
  /// whose lifetime and serialization is managed by a Param object.
  ///
  /// Properly initialized versions of this class can only be configured by a
  /// Parameters object (and through the pup routine during deserialization).

  friend Parameters; // This is a friend since it calls the main constructor

public:

  /// Default Constructor
  Expression() throw()
    : value_expr_(nullptr),
      is_scalar_expression_(false),
      param_name_(""),
      param_index_(-1)
  { }

  /// Default Copy and Move constructors
  Expression(const Expression &) = default;
  Expression(Expression&&) = default;
  Expression & operator= (const Expression & expr) = default;
  Expression & operator= (Expression &&expr) = default;

  // These expressions could probably just inspect the underlying pointer
  bool is_scalar_expression() const noexcept {return is_scalar_expression_;}
  bool is_logical_expression() const noexcept {return !is_scalar_expression_;}

  /// Indicates whether the object has been initialized
  bool initialized() const noexcept{ return value_expr_ != nullptr; } 

  /// Evaluate a floating-point expression given vectors x,y,z,t
  void evaluate_float
  ( int                      n,
    double *                 result,
    double *                 x,
    double *                 y,
    double *                 z,
    double                   t,
    const struct node_expr * node = nullptr ) const noexcept;

  /// Evaluate a logical expression given vectors x,y,z,t
  void evaluate_logical
  ( int                      n,
    bool *                   result,
    double *                 x,
    double *                 y,
    double *                 z,
    double                   t,
    const struct node_expr * node = nullptr) const noexcept;

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

private: // functions

  /// Main Constructor.
  ///
  /// This only get's called by Parameters::construct_or_rebuild_Expression_
  Expression(const struct node_expr * value_expr,
             bool scalar_expression,
             const std::string & param_name,
             int param_index = -1) noexcept
    : value_expr_(value_expr),
      is_scalar_expression_(scalar_expression),
      param_name_(param_name),
      param_index_(param_index)
  { }

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Pointer to the top level node_expr object where the parameter was defined.
  /// (The lifetime of this pointer is managed by a Param object)
  ///
  /// NOTE: This could theoretically become a dangling pointer (this is highly
  /// unlikely in practice). It would theoretically be safer to use a
  /// shared_ptr or manage a locally cloned hierarchy of node_expr pointers.
  const struct node_expr * value_expr_;

  /// Indicates whether the encapsulated expression is a scalar expression or
  /// a logical expression
  bool is_scalar_expression_;

  /// Name of the parameter where the encapsulated expression was defined.
  /// (Necessary for restoring this object after a restart).
  std::string param_name_;

  /// List index corresponding to the parameter location where the expression
  /// was defined. This is -1 if the parameter was not a list. (Necessary for
  /// restoring this object after a restart).
  int param_index_;

};

#endif /* PARAMETERS_EXPRESSIONS_HPP */
