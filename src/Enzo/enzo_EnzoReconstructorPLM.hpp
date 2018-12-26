// Right now, the primitive left and right fields are large enough that they
// can be treated as ordinary face-centered centered fields (with N+1 values)
// (Strictly speaking, only N-1 fields are needed along the dimension of
//  reconstruction)
// For the face-centered quantities, the ith element corresponds to the value
// of the face at i-1/2
// Additionally, reconstruction is performed so that fluxes can be computed for
// every cell (except the outermost faces)

// If this reconstructor is used N times for a single time step, then there
// must be at least N+1 ghost zones

// Still Need to implement limiters and primitive value floors (to ensure
// values don't
// get too small).



#ifndef ENZO_ENZO_RECONSTRUCTOR_PLM_HPP
#define ENZO_ENZO_RECONSTRUCTOR_PLM_HPP
class EnzoReconstructorPLM : public EnzoReconstructor
{
  /// @class    EnzoReconstructor
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates piecwise linear reconstruction of
  //            primitives at interface

public: // interface

  /// Create a new EnzoReconstructorPLM
  EnzoReconstructorPLM()
    : EnzoReconstructor()
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
    PUP::able::pup(p);
  };

  // Reconstructs the interface values
  // priml and primr are formally defined as corner-centered. However all
  // operations assume that they are face-centered along only 1 dimension.
  // This amounts to having some extra space at the end of the array
  void reconstruct_interface (Block *block, Grouping &prim_group,
			      Grouping &priml_group, Grouping &primr_group,
			      int dim, EnzoEquationOfState *eos);
};

#endif /* ENZO_ENZO_RECONSTRUCTOR_PLM_HPP */
