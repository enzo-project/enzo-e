// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Refine.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-08-14
/// @brief    [\ref Mesh] Declaration of the Refine Base class
///

#ifndef MESH_REFINE_HPP
#define MESH_REFINE_HPP

class Refine : public PUP::able 
{

  /// @class    Refine
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  Refine(double min_refine,
	 double max_coarsen,
	 int    max_level,
	 bool include_ghosts,
	 std::string output) throw()
    : min_refine_(min_refine),
      max_coarsen_(max_coarsen),
      max_level_(max_level),
      include_ghosts_(include_ghosts),
      output_(output),
      schedule_(NULL)
  {};

  /// CHARM++ PUP::able declaration
  PUPable_decl(Refine);

  /// CHARM++ migration constructor for PUP::able

  Refine (CkMigrateMessage *m)
    : PUP::able(m),
      min_refine_(0.0),
      max_coarsen_(0.0),
      max_level_(0),
      include_ghosts_(false),
      output_(""),
      schedule_(NULL)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Evaluate the refinement criteria, updating the refinement field
  virtual int apply (Block * block) throw ()
  {
    ERROR("Refine::apply",
	 "'Abstract' virtual function Refine::apply() should not be called");
    return 0;
  };

  /// Return the name of the refinement criteria
  virtual std::string name () const { return "unknown"; }

  /// Clear the output field to the default coarsen (-1)
  void * initialize_output_(FieldData * field_data);

  /// Return the Schedule object pointer
  Schedule * schedule() throw() 
  { return schedule_; }

  /// Set schedule
  void set_schedule (Schedule * schedule) throw();

protected: // functions

  /// Don't refine if already at max_level_
  void adjust_for_level_ (int * adapt_result, int level) const throw ()
  {
    if (level >= max_level_ && *adapt_result == adapt_refine) {
      *adapt_result = adapt_same;
    }
  }

protected:

  /// Minimum allowed value before refinement kicks in
  double min_refine_;

  /// Maximum allowed value before coarsening is allowed
  double max_coarsen_;

  /// Maximum level allowed to refine.  May be different than global
  /// maximum for the simulation
  int max_level_;

  /// Whether to include ghost zones when evaluating refinement criteria
  bool include_ghosts_;

  /// Field name to write refinement field to (-1 coarsen 0 same +1 refine)
  std::string output_;
  
  /// Schedule for refinement; NULL if none
  Schedule * schedule_;

};

#endif /* MESH_REFINE_HPP */

