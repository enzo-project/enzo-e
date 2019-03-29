#ifndef ENZO_ENZO_RECONSTRUCTOR_NN_HPP
#define ENZO_ENZO_RECONSTRUCTOR_NN_HPP
class EnzoReconstructorNN : public EnzoReconstructor
{
  /// @class    EnzoReconstructorNN
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates nearest neighbor reconstruction of
  //            primitives at interface

public: // interface

  /// Create a new EnzoReconstructorNN
  EnzoReconstructorNN(std::vector<std::string> group_names)
    : EnzoReconstructor(group_names)
  { }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoReconstructorNN);

  /// CHARM++ migration constructor for PUP::able
  EnzoReconstructorNN (CkMigrateMessage *m)
    : EnzoReconstructor(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    EnzoReconstructor::pup(p);
  };

  // Reconstructs the interface values
  // priml and primr are formally defined as corner-centered. However all
  // operations assume that they are face-centered along only 1 dimension.
  // This amounts to having some extra space at the end of the array
  void reconstruct_interface (Block *block, Grouping &prim_group,
			      Grouping &priml_group, Grouping &primr_group,
			      int dim, EnzoEquationOfState *eos);
};

#endif /* ENZO_ENZO_RECONSTRUCTOR_NN_HPP */
