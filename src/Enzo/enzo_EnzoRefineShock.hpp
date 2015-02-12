// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRefineShock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Jul 21 16:23:57 PDT 2014
/// @brief    [\ref Enzo] Declaration of the EnzoRefineShock class
///

#ifndef ENZO_ENZO_REFINE_SHOCK_HPP
#define ENZO_ENZO_REFINE_SHOCK_HPP

class EnzoRefineShock : public Refine {

  /// @class    EnzoRefineShock
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] 

public: // interface

  /// Constructor
  EnzoRefineShock
  (const FieldDescr * field_descr,
   double pressure_min_refine,
   double pressure_max_coarsen,
   double energy_ratio_min_refine,
   double energy_ratio_max_coarsen,
   double gamma,
   int comoving_coordinates,
   std::string output) throw();

  /// default constructor
  // EnzoRefineShock () throw() : Refine() {};

  PUPable_decl(EnzoRefineShock);

  EnzoRefineShock(CkMigrateMessage *m) : Refine (m) {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change

    Refine::pup(p);

    TRACEPUP;

    p | pressure_min_refine_;
    p | pressure_max_coarsen_;
    p | energy_ratio_min_refine_;
    p | energy_ratio_max_coarsen_;
    p | gamma_;
    p | comoving_coordinates_;
  }

  /// Evaluate the refinement criteria, updating the refinement field
  virtual int apply (Block * block,
		     const FieldDescr * field_descr) throw();

  virtual std::string name () const { return "shock"; };

private: // functions

  template <class T>
  void evaluate_block_( const T ** v3,
			const T * te,
			const T * de,
			const T * p,
			T * output,
			int ndx, int ndy, int ndz,
			int nx, int ny, int nz,
			int gx, int gy, int gz,
			bool *any_refine,
			bool *all_coarsen, 
			int rank);

private: // attributes

  /// Refine when pressure becomes greater than this somewhere
  double pressure_min_refine_;

  /// Coarsen when all pressure is less than this everywhere
  double pressure_max_coarsen_;

  /// Refine when the energy ratio becomes greater than this somewhere
  double energy_ratio_min_refine_;

  /// Coarsen when all the energy ratio is less than this everywhere
  double energy_ratio_max_coarsen_;

  /// Gamma
  double gamma_;

  /// Comoving coordinates
  int comoving_coordinates_;
};

#endif /* ENZO_ENZO_REFINE_SHOCK_HPP */

