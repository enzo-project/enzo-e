// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemann.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of the Riemann Solver abstract base
/// class. This class should be subclassed to implement various riemann solvers.

#ifndef ENZO_ENZO_RIEMANN_HPP
#define ENZO_ENZO_RIEMANN_HPP



class EnzoRiemann : public PUP::able
{
  /// @class    EnzoRiemann
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate approximate Riemann Solvers

public: // interface

  /// Factory method for constructing the EnzoRiemann object. (The signature
  /// this method must be modified as additional physics gets implemented
  static EnzoRiemann* construct_riemann(std::string solver,
					const EnzoFieldConditions cond);

  /// Constructor
  EnzoRiemann(EnzoFieldConditions cond,
	      std::vector<std::string> &extra_passive_groups,
	      FluxFunctor** flux_funcs, int n_funcs);

  /// Virtual destructor
  virtual ~EnzoRiemann();

  /// CHARM++ PUP::able declaration
  PUPable_abstract(EnzoRiemann);

  /// CHARM++ migration constructor for PUP::able
  EnzoRiemann (CkMigrateMessage *m)
    : PUP::able(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  {
    PUP::able::pup(p);
  }

  /// Computes the Riemann Fluxes for each conserved field along a given
  /// dimension, dim
  /// @param block holds data to be processed
  /// @param priml_group,primr_group holds field names where reconstructed
  ///  left/right face-centered primitives are stored. The relevant fields
  ///  should be formally defined as cell-centered (to allow for reuse). During
  ///  the calculation, they are treated as face-centered (without having
  ///  values on the exterior faces of the block). Consequentially there will
  ///  be some unused space at the end of the arrays.
  /// @param flux_group holds field names where the calculated fluxes will be
  ///  stored. The relevant fields should be face-centered along the dimension
  ///  of the calculation (without having values on the exterior faces of the
  ///  block)
  /// @param cons_group,cons_group holds field names where reconstructed
  ///  left/right face-centered conserved quantities are stored. The relevant
  ///  fields are expected to have the same properties as priml_group and
  ///  primr_group, respectively (with respect to the field cell-centering).
  ///  The values in these fields should be computed from the fields in
  ///  priml_group and primr_group ahead of time.
  /// @param dim Dimension along which to reconstruct interface values. Values
  ///  of 0, 1, and 2 correspond to the x, y, and z directions, respectively.
  /// @param eos Instance of the fluid's EnzoEquationOfState object
  virtual void solve (Block *block, Grouping &priml_group,
		      Grouping &primr_group, Grouping &flux_group,
		      Grouping &consl_group, Grouping &consr_group, int dim,
		      EnzoEquationOfState *eos) = 0;

};

#endif /* ENZO_ENZO_RIEMANN_HPP */
