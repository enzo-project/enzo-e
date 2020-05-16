// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoCenteredFieldRegistry.hpp
/// @author   Matthew W. Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs March 28 2019
/// @brief    [\ref Enzo] Declaration of the EnzoCenteredFieldRegistry class.
///
/// This class serves to track all cell-centered MHD quantities that are not
/// passively advected scalars. The class also tracks the group names of
/// passively advected scalars. The API supports constructing Groupings of
/// fields (with names matching the the names of registered names) and
/// determining which quantities are conserved, specific, or other.

#ifndef ENZO_ENZO_CENTERED_FIELD_REGISTRY_HPP
#define ENZO_ENZO_CENTERED_FIELD_REGISTRY_HPP


#include <type_traits> // std::is_arithmetic
#include <string>      // std::string, std::to_string

#include <cstdio> // For help with nicely printing values


//----------------------------------------------------------------------
/// @var      passive_group_names
/// @brief    vector of passively advected group names
///
/// This should be updated as new groups are added
const std::vector<std::string> passive_group_names = {"colour"};


//----------------------------------------------------------------------
/// @def      FieldCat
/// @brief    Categorizes a field as conserved, specific, or other
///
/// The categories are:
///  - conserved - the field is a conserved fluid property (e.g. density)
///  - specific - the field is a specific primitive fluid property (multiply by
///    density to get conserved, e.g. velocity) primitive
///  - other - anything else (e.g. temperature, pressure)
enum class FieldCat{conserved, specific, other};


//----------------------------------------------------------------------
/// @def      FIELD_TABLE
/// @brief    xmacro table of cell-centered registered (non-passive) fields
///
/// Table of fluid properties represented by cell-centered fields. This table
/// must be updated as new physics gets added. The columns of the table are:
///   1. Quantity name.
///   2. Mathematical classificiation (SCALAR or VECTOR)
///   3. classification of the quantity as conserved, specific or other
///   4. Advection classification - if T, then in some context, the quantity is
///      an actively advected quantity (i.e. a Riemann Flux is computed for the
///      term) - otherwise it should be F (quantities that always
///      act as source terms should not be included)
///
/// The entries in this table determine the name of the field to represent the
/// quantity:
///    - A SCALAR quantity is represented by a field with the name in column 1.
///    - A VECTOR quantity is represented by 3 fields. The fields are named
///      {column_1}_x, {column_1}_y, {column_2}_z
///
/// As additional physics are added this table should be updated (e.g.
/// cr_energy, cr_flux).
#define FIELD_TABLE                                                          \
  /* Hydro & MHD Fields directly related to Advection */                     \
  ENTRY(                density, SCALAR, FieldCat::conserved, T)             \
  ENTRY(               velocity, VECTOR,  FieldCat::specific, T)             \
  ENTRY(           total_energy, SCALAR,  FieldCat::specific, T)             \
  ENTRY(                 bfield, VECTOR, FieldCat::conserved, T)             \
  /* only an actively advected quantity under the dual energy formalism */   \
  ENTRY(        internal_energy, SCALAR,  FieldCat::specific, T)             \
  /* Derived Fields */                                                       \
  ENTRY(               pressure, SCALAR,     FieldCat::other, F)             \
  ENTRY(            temperature, SCALAR,     FieldCat::other, F)             \
  ENTRY(           cooling_time, SCALAR,     FieldCat::other, F)             \
  ENTRY(           acceleration, VECTOR,  FieldCat::specific, F)             \
  /* Grackle Related Fields */                                               \
  ENTRY(  specific_heating_rate, SCALAR,     FieldCat::other, F)             \
  ENTRY(volumetric_heating_rate, SCALAR,     FieldCat::other, F)             \
  ENTRY(             HI_density, SCALAR, FieldCat::conserved, F)             \
  ENTRY(            HII_density, SCALAR, FieldCat::conserved, F)             \
  ENTRY(            HeI_density, SCALAR, FieldCat::conserved, F)             \
  ENTRY(           HeII_density, SCALAR, FieldCat::conserved, F)             \
  ENTRY(          HeIII_density, SCALAR, FieldCat::conserved, F)             \
  ENTRY(              e_density, SCALAR, FieldCat::conserved, F)             \
  ENTRY(             HM_density, SCALAR, FieldCat::conserved, F)             \
  ENTRY(            H2I_density, SCALAR, FieldCat::conserved, F)             \
  ENTRY(           H2II_density, SCALAR, FieldCat::conserved, F)             \
  ENTRY(             DI_density, SCALAR, FieldCat::conserved, F)             \
  ENTRY(            DII_density, SCALAR, FieldCat::conserved, F)             \
  ENTRY(            HDI_density, SCALAR, FieldCat::conserved, F)             \
  ENTRY(          metal_density, SCALAR, FieldCat::conserved, F)             \
  /* Field for Gravity */                                                    \
  ENTRY(              potential, SCALAR,  FieldCat::specific, F)             \
  /* Fields for Turbulence */                                                \
  ENTRY(                driving, VECTOR,  FieldCat::specific, F)

//----------------------------------------------------------------------

/// @typedef FT_row
/// @brief The type of the of the tuple holding the FIELD_TABLE in memory
typedef std::tuple<std::string, FieldCat, bool> FT_row;

//----------------------------------------------------------------------

class EnzoCenteredFieldRegistry
{
  /// @class    EnzoCenteredFieldRegistry
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Serves as a registry for non-passively advected
  ///           cell-centered fields
  ///
  /// Specifically, the registered fields are given by the entries of
  /// FIELD_TABLE. The name of the field(s) that correspond to an entry in the
  /// table are dependent on the mathematical type of the quantity and on the
  /// name given in column 1.
  ///  - A SCALAR simply corresponds to a field with a name given in column 1.
  ///  - A VECTOR corresponds to 3 fields. The fields are named {column_1}_x,
  ///    {column_1}_y, {column_2}_z
  /// To register a new field, add a new entry to FIELD_TABLE

public:

  EnzoCenteredFieldRegistry();

  // some utility methods - some of these are called in more important methods
  // These are not particularly well optimized

  /// Returns a vector of registered quantities.
  std::vector<std::string> get_registered_quantities() const
  { return table_keys_; }

  /// Returns the vector of registered field names
  std::vector<std::string> get_registered_fields() const;

  /// Checks that that the quantity names are in FIELD_TABLE. If not then,
  /// raises an error.
  void check_known_quantity_names(const std::vector<std::string> names) const;

  /// provides the quantity properties listed in FIELD_TABLE (if present)
  ///
  /// returns true when successful (i.e. quantity is actually included in the
  /// table) and false when unsuccesful
  bool quantity_properties(std::string name, bool* vector_quantity = 0,
			   FieldCat* category = 0,
			   bool* actively_advected = 0) const;

  /// Checks is the name of an actively advected vector component
  bool is_actively_advected_vector_component(std::string name,
					     bool ijk_suffix) const noexcept;

  /// Determine the actively advected quantity associated with the given name
  /// If there is not an associated quantity, `""` is returned.
  ///
  /// There are 2 ways that a name can be associated with a quantity
  ///   1. They can exactly match the name of the quantity.
  ///   2. For VECTOR quantities, they can be composed of the quantity name
  ///      followed by a 2 character suffix. If ijk_suffix is true, then the
  ///      suffixes are {'_i', '_j', '_k'}. Otherwise they are
  ///      {'_x', '_y', '_z'}.
  std::string get_actively_advected_quantity_name
  (std::string name, bool ijk_suffix) const noexcept;

  // more important methods:
  
  /// Constructs a Grouping of fields and yields a pointer to it from a vector
  /// of quantity names. The quantity names must match entries of FIELD_TABLE.
  ///
  /// The names of the names of the fields in each group depend on the the
  /// leading_prefix argument and the quantity type indicated in FIELD_TABLE
  /// The rules for determining the field names are:
  /// - SCALAR quantity: leading_prefix + group_name
  /// - Vector conserved quantity:
  ///   - leading_prefix + group_name + "_x"
  ///   - leading_prefix + group_name + "_y"
  ///   - leading_prefix + group_name + "_z"
  Grouping* build_grouping(const std::vector<std::string> quantity_names,
			   const std::string leading_prefix = "") const;
  
  /// Returns a vector of passive scalar group names
  ///
  /// To register new names add entry to the static constant vector variable
  /// called passive_group_names
  static std::vector<std::string> passive_scalar_group_names()
  { return passive_group_names; }

private:

  /// Helper method of build_grouping. Adds a fields to the grouping
  void add_group_fields_(Grouping *grouping, const std::string group_name,
			 const std::string quantity_type, FieldCat category,
			 const std::string leading_prefix) const;


  /// Helper method of prepare_lut. Determines conserved_start, conserved_stop,
  /// specific_start, specific_stop, other_start, other_stop, nfields
  void prepare_lut_(const std::vector<std::string> quantity_names,
		    int &conserved_start, int &conserved_stop,
		    int &specific_start, int &specific_stop,
		    int &other_start, int &other_stop, int &nfields,
		    const std::vector<std::string> flagged_quantities) const;

private: // attributes

  /// representation of FIELD_TABLE in memory
  std::map<std::string, FT_row> field_table_;

  /// vector of keys for field_table_;
  std::vector<std::string> table_keys_;
};

#endif /* ENZO_ENZO_CENTERED_FIELD_REGISTRY_HPP */
