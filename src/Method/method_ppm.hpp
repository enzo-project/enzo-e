#ifndef METHOD_PPM_HPP
#define METHOD_PPM_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      method_ppm.hpp
 * @brief     Defines the PPM Enzo hydrodynamics Method
 * @author    James Bordner
 * @date      Sat Jul 18 13:45:13 PDT 2009
 *
 * Defines the PPM Enzo hydrodynamics Method
 *
 * $Id$
 *
 *********************************************************************
 */

//  *     d) iflatten = 0 -> no flattening
//  *     e) iflatten = 1 -> The simple flattening scheme described in eqs
//  *        A1 & A2 of CW84.
//  *     f) iflatten = 2 -> Flattening needed for Lagrangean hydrodynamics
//  *        (note: not L+remap).  Stuck in for completeness, but not tested.
//  *        This is described in equations. A.4 - A.6.
//  *     g) iflatten = 3 -> Flattening recomended for multidimensional
//  *        calculations.  Described in CW84 A7-A10.

enum flatten_type {
  flatten_none,
  flatten_simple,
  flatten_lagrangean,
  flatten_multidimension
};

class MethodPpm : public Method {

/** 
 *********************************************************************
 *
 * @class     MethodPpm
 * @brief     PPM Enzo hydrodynamics Method
 * @ingroup   Method
 *
 * PPM Enzo hydrodynamics Method
 *
 *********************************************************************
 */

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

public:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------


  /// Create a new PPM method
  MethodPpm() throw();

  /// Copy a PPM method
  MethodPpm(const MethodPpm &) throw();

  /// Initialize the PPM method
  void initialize() throw(); 

  /// Apply the method
  void apply() throw(); 

private:

  void solve_hydro_equations_() throw ();

  flatten_type flatten_type_;

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

};

#endif /* METHOD_PPM_HPP */

