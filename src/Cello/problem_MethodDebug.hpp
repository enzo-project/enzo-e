// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodDebug.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2020-03-19
/// @brief    [\ref Problem] Declaration for the MethodDebug class

#ifndef PROBLEM_METHOD_DEBUG_HPP
#define PROBLEM_METHOD_DEBUG_HPP

class MethodDebug : public Method
{
  /// @class    MethodDebug
  /// @ingroup  MethodDebug
  /// @brief    [\ref MethodDebug] Declaration of MethodDebug
  ///
  /// Correct fluxes of conserved fields at AMR refinement level jumps

public: // interface

  /// Create a new MethodDebug
  MethodDebug
  (int num_fields,
   bool l_print,
   bool l_coarse,
   bool l_ghost) throw();

  /// Destructor
  virtual ~MethodDebug() throw()
  {};

  /// Charm++ PUP::able declarations
  PUPable_decl(MethodDebug);

  /// Charm++ PUP::able migration constructor
  MethodDebug (CkMigrateMessage *m)
    : Method(m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;
    Method::pup(p);
    p | field_sum_;
    p | field_min_;
    p | field_max_;
    p | field_count_;
    p | l_print_;
    p | l_coarse_;
    p | l_ghost_;
  };

  void compute_continue_sum_fields ( Block * block, CkReductionMsg * msg) throw();

public: // virtual functions

  /// Apply the method to advance a block one timestep 

  virtual void compute ( Block * block) throw();

  /// Return the name of this MethodDebug
  virtual std::string name () throw () { return "debug"; }

protected: // attributes

  std::vector<long double> field_sum_;
  std::vector<long double> field_min_;
  std::vector<long double> field_max_;
  std::vector<long double> field_count_;
  bool l_print_;
  bool l_coarse_;
  bool l_ghost_;
};


#endif /* PROBLEM_METHOD_DEBUG_HPP */
