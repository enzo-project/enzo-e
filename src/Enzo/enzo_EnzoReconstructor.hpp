
#ifndef ENZO_ENZO_RECONSTRUCTOR_HPP
#define ENZO_ENZO_RECONSTRUCTOR_HPP

#include <pup_stl.h>

class EnzoReconstructor : public PUP::able
{
  /// @class    EnzoReconstructor
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates reconstruction of primitives at
  //            interface

public: // interface

  // Factory method for constructing EnzoReconstructor
  // The signature of this method may need to be modified
  static EnzoReconstructor* construct_reconstructor(std::string name,
						    const EnzoFieldConditions cond);

  /// Create a new EnzoReconstructor
  EnzoReconstructor(std::vector<std::string> group_names) throw()
    : group_names_(group_names)
  { }

  /// Virtual destructor
  virtual ~EnzoReconstructor()
  {  }

  /// CHARM++ PUP::able declaration
  PUPable_abstract(EnzoReconstructor);

  /// CHARM++ migration constructor for PUP::able
  EnzoReconstructor (CkMigrateMessage *m)
    : PUP::able(m)  
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    PUP::able::pup(p);
    p|group_names_;
  }

  // Reconstructs the interface values
  // priml and primr are formally defined as corner-centered. However all
  // operations assume that they are face-centered along only 1 dimension.
  // This amounts to having some extra space at the end of the array
  virtual void reconstruct_interface (Block *block, Grouping &prim_group,
				      Grouping &priml_group,
				      Grouping &primr_group, int dim,
				      EnzoEquationOfState *eos)=0;
protected:
  std::vector<std::string> group_names_;
};

#endif /* ENZO_ENZO_RECONSTRUCTOR_HPP */
