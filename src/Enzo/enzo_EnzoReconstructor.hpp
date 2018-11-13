
#ifndef ENZO_ENZO_RECONSTRUCTOR_HPP
#define ENZO_ENZO_RECONSTRUCTOR_HPP
class EnzoReconstructor : public PUP::able
{
  /// @class    EnzoReconstructor
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates reconstruction of primitives at
  //            interface

public: // interface

  /// Create a new EnzoReconstructor
  EnzoReconstructor() throw()
  {}

  /// CHARM++ PUP::able declaration
  //PUPable_abstract(EnzoReconstructor);

  /// CHARM++ migration constructor for PUP::able
  //EnzoReconstructor (CkMigrateMessage *m)
  //  : PUP::able(m)
  //{  }

  /// Virtual destructor
  virtual ~EnzoReconstructor()
  {  }

  /// CHARM++ Pack / Unpack function
  //void pup (PUP::er &p)
  //{
  //  PUP::able::pup(p);
  //}

  // Reconstructs the interface values
  // (May want to avoid reconstructing longitudinal values)
  virtual void reconstruct_interface (Block *block,
				      const std::vector<int> &prim_ids,
				      const std::vector<int> &priml_ids,
				      const std::vector<int> &primr_ids)=0;
};

#endif /* ENZO_ENZO_RECONSTRUCTOR_HPP */
