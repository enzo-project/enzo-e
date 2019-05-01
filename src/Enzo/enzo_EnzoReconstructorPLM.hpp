// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoReconstructor.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed May 1 2019
/// @brief    [\ref Enzo] Implementation of Enzo's Piecewise Linear
///           Reconstruction

#ifndef ENZO_ENZO_RECONSTRUCTOR_PLM_HPP
#define ENZO_ENZO_RECONSTRUCTOR_PLM_HPP

class EnzoReconstructorPLM : public EnzoReconstructor
{
  /// @class    EnzoReconstructorPLM
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates piecwise linear reconstruction of
  ///           primitives at interface

public: // interface

  /// Create a new EnzoReconstructorPLM
  EnzoReconstructorPLM(std::vector<std::string> group_names)
    : EnzoReconstructor(group_names)
  { }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoReconstructorPLM);

  /// CHARM++ migration constructor for PUP::able
  EnzoReconstructorPLM (CkMigrateMessage *m)
    : EnzoReconstructor(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    EnzoReconstructor::pup(p);
  };

  void reconstruct_interface (Block *block, Grouping &prim_group,
			      Grouping &priml_group, Grouping &primr_group,
			      int dim, EnzoEquationOfState *eos);

  int get_staling_rate()
  {
    return 2;
  }

};

#endif /* ENZO_ENZO_RECONSTRUCTOR_PLM_HPP */
