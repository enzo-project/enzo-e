// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialShockTube.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed April 17 2019
/// @brief    [\ref Enzo] Initialization routine for RJ2a shocktube problem
///           mentioned by Stone et al 2008 and that is detailed by in Ryu &
///           Jones (1995). This is primarily to test the VL+CT MHD integrator.

#ifndef ENZO_ENZO_INITIAL_SHOCK_TUBE_HPP
#define ENZO_ENZO_INITIAL_SHOCK_TUBE_HPP

class EnzoInitialShockTube : public Initial {
  /// @class    EnzoInitialShockTube
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initializer for the axis-aligned Shock Tube test
  /// problem for the VLCT method from Ryu and Jones (95) and Stone et al. (08).
  /// Assumes an adiabatic ideal gas, with gamma = 5/3

public: // interface

  /// Constructor
  EnzoInitialShockTube(int cycle, double time, std::string aligned_ax)
    : Initial(cycle, time),  aligned_ax_(aligned_ax)
  {
    ASSERT("EnzoInitialShockTube",
	   "Invalid aligned_ax value specified (must 0, 1, or 2).",
	   aligned_ax_ == "x" || aligned_ax_ == "y" || aligned_ax_ == "z");
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialShockTube);

  /// CHARM++ migration constructor
  EnzoInitialShockTube(CkMigrateMessage *m)
    : Initial (m), aligned_ax_("x")
  {  }

  /// Destructor
  virtual ~EnzoInitialShockTube() throw()
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block
  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

protected: // functions

  /// Allocate and determine the slices along the ax_aligned_ dimension in for
  /// block that corresponds to the left and right states. If one of the slices
  /// is not present in the current block, then the corresponding slice is not
  /// allocated
  void prep_aligned_slices_(Block *block, ESlice *l_slice, ESlice *r_slice);

  /// sets the subarray of arr given by the slice along ax_aligned_ and the
  /// full extents of the other dimensions to val.
  void initializer_helper_(ESlice &slice, enzo_float val, EFlt3DArray &arr);

private: // attributes

  /// indicates the axis along which the problem is initialized: "x", "y", "z"
  std::string aligned_ax_;

};

#endif /* ENZO_ENZO_INITIAL_SHOCK_TUBE_HPP */
