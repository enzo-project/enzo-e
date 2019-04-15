// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineRegion.hpp
/// @author   Andrew Emerick (aemerick11@gmail.com)
/// @date     2019-04-13
/// @brief    [\ref Mesh] Declaration of the RefineRegion class
///

#ifndef MESH_REFINE_REGION_HPP
#define MESH_REFINE_REGION_HPP

class RefineRegion : public Refine {

  /// @class     RefineRegion
  /// @ingroup   Mesh
  /// @brief     [\ref Mesh]
  ///
  /// RefineRegion allows for the creation of rectangular
  /// regions with separate minimum and maximum refinement levels.
  /// Additional refinemnet criteria operate within these bounds.
  /// Static refinement regions can be defined by setting min and max
  /// level to be the same. Regions do not have to be contiguous or
  /// nested. The magic happens in control_adapt

public: // interface

   /// Constructor
   RefineRegion() throw();

   PUPable_decl(RefineRegion);

   RefineRegion(CkMigrateMessage *m) : Refine (m) {}

   /// CHARM++ Pack / Unpack function
   inline void pup (PUP::er &p)
   {
     TRACEPUP;
     // NOTE: change this function whenever attributes change
     Refine::pup(p);
   }

   /// Evaluate the refinement criteria, update the refinement field
   virtual int apply (Block * block) throw();

   virtual std::string name () const { return "region"; };

private: // functions

};

#endif /* MESH_REFINE_REGION_HPP */
