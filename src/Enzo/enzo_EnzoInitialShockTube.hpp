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
  EnzoInitialShockTube(double gamma, int cycle, double time,
		       std::string setup_name, std::string aligned_ax_name);

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialShockTube);

  /// CHARM++ migration constructor
  EnzoInitialShockTube(CkMigrateMessage *m)
    : Initial (m), gamma_(0.0), setup_name_(""), aligned_ax_(0)
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
  void prep_aligned_slices_(Block *block, CSlice **l_slice, CSlice **r_slice);

  /// sets the subarray of arr given by the slice along ax_aligned_ and the
  /// full extents of the other dimensions to val.
  void initializer_helper_(CSlice &slice, enzo_float val, EFlt3DArray &arr);

private: // attributes

  /// adiabatic index
  double gamma_;

  /// indicates the type of shock tube to set up
  std::string setup_name_;

  /// indicates the axis along which the problem is initialized: "x", "y", "z"
  int aligned_ax_;
};

#endif /* ENZO_ENZO_INITIAL_SHOCK_TUBE_HPP */
