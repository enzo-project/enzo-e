// See LICENSE_CELLO file for license and copyright information

/// @file     data_FieldFluxes.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-02-20
/// @brief    [\ref Data] Declaration of the FieldFluxes class

#ifndef DATA_FIELD_FLUXES_HPP
#define DATA_FIELD_FLUXES_HPP

template <class T>
class FieldFluxes {

  /// @class    FieldFluxes
  /// @ingroup  Data
  /// @brief    [\ref Data] 

public: // interface

  /// Create a FieldFluxes object given a Field
  FieldFluxes (Field field)
  {
  }
  
  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Convert from Cello FieldFluxes object to ENZO fluxes arrays
  void save_to_arrays
  (
   Field field,
   T ** LeftFluxes[3],
   T ** RightFluxes[3]
   )
  {
  }
  
  /// Convert from ENZO fluxes arrays to Cello FieldFluxes object
  void load_from_arrays
  (
   Field field,
   T ** LeftFluxes[3],
   T ** RightFluxes[3]
   )
  {
  }

  /// Update fluxes from finer Block
  void update (Field field, FieldFace * field_face);

  /// Finalize Fluxes
  void finalize (Field field);

  /// Add to boundary fluxes
  void add_to_boundary (Field field, FieldFace * field_face);

  
private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  std::vector < std::vector<T *> > fluxes_[3][2]; // [index_field][axis][face]

};

#endif /* DATA_FIELD_FLUXES_HPP */

