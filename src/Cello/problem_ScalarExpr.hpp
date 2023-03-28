// See LICENSE_CELLO file for license and copyright information

/// @file     problem_ScalarExpr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-03-28
/// @brief    [\ref Problem] Declaration of the ScalarExpr class

#ifndef PROBLEM_SCALAR_EXPR_HPP
#define PROBLEM_SCALAR_EXPR_HPP

class ScalarExpr {

  /// @class    ScalarExpr
  /// @ingroup  Problem
  /// @brief    [\ref Problem]
  ///
  /// This class effectively acts as a union of a ``double`` and a pointer to a
  /// ``Param`` instance (that holds a holds a floating-point expression). In
  /// the latter case, the ``Param`` pointer's lifetime is not managed by this
  /// object.

public: // interface

  /// Default Constructor
  ///
  /// The resulting instance is identical to an instance that was constructed
  /// from a ``Param`` instance storing a floating point value of ``0.0``.
  ScalarExpr() throw()
  : param_(nullptr),
    value_(0)
  { }

  /// construct an instance from a non-null pointer to a ``Param`` instance.
  ScalarExpr(Param * param) throw();

  /// destructor
  ///
  /// Earlier versions of this method tried to ``delete`` param_ when it's not
  /// a ``nullptr``. This was technically unsound because the global
  /// ``Parameter`` object would subsequently tried to deallocate the same
  /// ``Param`` object. In practice, this wasn't ever a problem, because this
  /// destructor wasn't ever called
  ~ScalarExpr() throw() { }

  /// We delete the copy constructor/assignment-operator because the historic
  /// implementation has been broken.
  ///
  /// Historically, it sought to clone the ``param_`` instance (when it wasn't a
  /// ``nullptr``) using ``Param``'s copy constructor. However, ``Param``'s
  /// copy constructor just provides a warning saying it isn't implemented.
  ///
  /// @note
  /// It may be worth revisiting these operations in the future. However, if we
  /// reintroduce the older behavior will require us to revisit what the
  /// primary constructor and destructor do
  ScalarExpr(const ScalarExpr & scalar_expr) throw() = delete;
  ScalarExpr & operator= (const ScalarExpr & scalar_expr) throw() = delete;

  /// implement the default move constructor and move assignment operations
  ScalarExpr(ScalarExpr &&) = default;
  ScalarExpr & operator= (ScalarExpr&&) = default;

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;

    // NOTE: change this function whenever attributes change

    bool is_null = (param_ == nullptr);
    p | is_null;
    if (!is_null) {
      ERROR("ScalarExpr::pup()",
            "not implemented when param_ isn't a nullptr");
    }

    p | value_;
  }

  /// Evaluate mask at a point
  double evaluate (double t, double x, double y, double z,
		   std::shared_ptr<Mask> mask = nullptr, double deflt = 0) const;

  /// Return mask values in an array
  template <class T>
  void evaluate (T *value, double t,
		 int ndx, int nx, double * x,
		 int ndy, int ny, double * y,
		 int ndz, int nz, double * z, 
		 std::shared_ptr<Mask> mask = nullptr, T * deflt = 0) const;

  /// Return mask values in an array
  template <class T>
  void evaluate (T *value, double t,
		 int ndx, int nx, double * x,
		 int ndy, int ny, double * y,
		 int ndz, int nz, double * z) const
  {
    evaluate(value,t,ndx,nx,x,ndy,ny,y,ndz,nz,z,0,0);
  }


  /// method used for debugging
  std::string expr_to_string() const throw();

private: // attributes

  /// Value if param_ type is precision_float_expr
  Param * param_;

  /// Value if param_ type is precision_float
  double value_;

};
#endif /* PROBLEM_SCALAR_EXPR_HPP */

