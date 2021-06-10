// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodFluxCorrect.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2020-03-19
/// @brief    [\ref Problem] Declaration for the MethodFluxCorrect class

#ifndef PROBLEM_METHOD_FLUX_CORRECT_HPP
#define PROBLEM_METHOD_FLUX_CORRECT_HPP

class MethodFluxCorrect : public Method
{
  /// @class    MethodFluxCorrect
  /// @ingroup  MethodFluxCorrect
  /// @brief    [\ref MethodFluxCorrect] Declaration of MethodFluxCorrect
  ///
  /// Correct fluxes of conserved fields at AMR refinement level jumps

public: // interface

  /// Create a new MethodFluxCorrect
  MethodFluxCorrect
  (const std::string group, bool enable,
   const std::vector<std::string>& min_digits_fields,
   const std::vector<double>& min_digits_values) throw();

  /// Destructor
  virtual ~MethodFluxCorrect() throw()
  {};

  /// Charm++ PUP::able declarations
  PUPable_decl(MethodFluxCorrect);

  /// Charm++ PUP::able migration constructor
  MethodFluxCorrect (CkMigrateMessage *m)
    : Method(m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;
    Method::pup(p);
    p | ir_pre_;
    p | group_;
    p | enable_;
    p | min_digits_map_;
    p | field_sum_;
    p | field_sum_0_;
    // don't pup scratch_
  };

  void compute_continue_refresh ( Block * block) throw();
  void compute_continue_sum_fields ( Block * block, CkReductionMsg * msg) throw();

public: // virtual functions

  /// Apply the method to advance a block one timestep 

  virtual void compute ( Block * block) throw();

  /// Return the name of this MethodFluxCorrect
  virtual std::string name () throw ()
  { return "flux_correct"; }

protected: // functions

  void flux_correct_ (Block * block);
  
protected: // attributes

  /// Refresh id
  int ir_pre_;
  /// Field group to apply flux-correction to
  std::string group_;

  /// Whether to actually perform the flux-correction.  Setting
  /// to false still computes conserved values and fails if below
  /// min_digits
  bool enable_;

  /// Used for testing, write FAIL if number of conserved digits falls below
  /// this number for the specified field(s). This is empty by default (which
  /// effectively deactivates this checking).
  std::map<std::string,double> min_digits_map_;

  std::vector<long double> field_sum_;
  std::vector<long double> field_sum_0_;

  /// scratch space for performing the flux correction
  std::vector<cello_float> scratch_;
};


#endif /* PROBLEM_METHOD_FLUXCORRECT_HPP */
