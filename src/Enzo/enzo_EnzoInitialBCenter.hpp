/// See LICENSE_CELLO file for license and copyright information
/// @file     enzo_EnzoInitialBCenter.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed Jul 3 2019
/// @brief    [\ref Problem] VLCT Bfield initialization method

#ifndef ENZO_INITIAL_B_CENTER_HPP
#define ENZO_INITIAL_B_CENTER_HPP

// The primary purpose of this class is to initialize the
// cell-centered B-fields:
//   - given the a function for the vector potential in the parameter file
//     (this also initializes face-centered B-fields)
//   - given pre-initialized face-centered B-fields
//
// Additionally, it provides the static method initialize_bfield_center for
// other Initial classes to call directly.

class EnzoInitialBCenter : public Initial {

  /// @class    EnzoInitialPpmlTest
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initializes cell-centered B-field from a Vector
  /// potential (face-centered B-fields are also initialized) or from a
  /// pre-initialized face-centered B-fields

public: // interface

  /// Constructor
  EnzoInitialBCenter (Parameters * parameters,
		      int cycle, double time,
		      bool update_etot) throw ();

  /// Destructor
  ~EnzoInitialBCenter()
  {
    for (int i = 0; i < 3; i++){
      if (values_[i] != nullptr) {
	delete values_[i];
	values_[i] = nullptr;
      }
    }
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialBCenter);

  /// CHARM++ migration constructor
  EnzoInitialBCenter(CkMigrateMessage *m)
    : Initial (m),
      values_(),
      update_etot_(false)
  {  }


  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    // NOTE: update this member function whenever class attributes change
    Initial::pup(p);
    TRACEPUP;

    // This would be simplified if Value were declared as PUPable
    // First we pup which of the components exist
    bool has_values[3];
    for (int i=0; i<3; i++){
      has_values[i] = (values_[i] != nullptr);
      p|(has_values[i]);
    }
    // Now we pup the components that exist
    for (int i=0; i<3; i++){
      if (has_values[i]){
	if (p.isUnpacking()) {values_[i] = new Value;}
	p|*(values_[i]);
      } else{
	values_[i] = nullptr;
      }
    }
    p|update_etot_;
  }

  /// static method that initializes the cell-centered bfield from previously
  /// initialized face-centered bfields
  static void initialize_bfield_center( Block * block );

  /// static method that initializes the cell-centered (and face-centered)
  /// magnetic field given a Vector potential.
  ///
  /// Currently this assumes that the problem is 3D
  ///
  /// If the mesh has dimensions mx,my,mz (including ghost zones), then
  /// the arrays of vector potential values should have the following shapes:
  ///  Ax: (mz+1,my+1,  mx)
  ///  Ay: (mz+1,  my,mx+1)
  ///  Az: (  mz,my+1,mx+1)
  static void initialize_bfield_interface( Block * block,
					   CelloArray<double,3> &Ax,
					   CelloArray<double,3> &Ay,
					   CelloArray<double,3> &Az);

public: // virtual functions

  /// Initialize a Block
  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

protected: // attributes

  /// Each value will hold a pointer to a component of the vector potential
  Value* values_[3];

  /// Whether or not to update the total specific energy using the new
  /// calculated Bfield
  bool update_etot_;
};

#endif /* ENZO_INITIAL_B_CENTER_HPP */
