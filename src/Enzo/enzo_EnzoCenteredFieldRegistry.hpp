// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoCenteredFieldRegistry.hpp
/// @author   Matthew W. Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs March 28 2019
/// @brief    [\ref Enzo] Declaration of the EnzoCenteredFieldRegistry class.
///
/// FIELD_TABLE is a table (defined as a macro) that lists the names and
/// properties of all physical quantities that the code represented as fields.
/// It also tracks which fields play an active role in MHD integration (It is
/// used in the definition of Riemann Solvers).
///
/// EnzoCenteredFieldRegistry is effectively a glorified namespace that
/// group functions that provide an API for accessing information listed
/// in FIELD_TABLE at runtime. It also
///  - supports the use of this information to construct a Groupings of fields
///    (with names matching the registered names quantites)
///  - provides a list of groupings of passively advected scalars.

#ifndef ENZO_ENZO_CENTERED_FIELD_REGISTRY_HPP
#define ENZO_ENZO_CENTERED_FIELD_REGISTRY_HPP


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

/// @typedef ft_row
/// @brief Holds the data associated with a single row of FIELD_TABLE
struct ft_row{
  bool vector_quantity;
  FieldCat category;
  bool actively_advected;
};

//----------------------------------------------------------------------

/// @typedef ft_map
/// @brief The map type used to hold FIELD_TABLE in memory
typedef std::map<std::string, ft_row> ft_map;

//----------------------------------------------------------------------

class EnzoCenteredFieldRegistry
{
  /// @class    EnzoCenteredFieldRegistry
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Serves as a registry for non-passively advected
  ///           cell-centered fields
  ///
  /// The registered fields are given by the entries of FIELD_TABLE. The name
  /// of the field(s) that correspond to an entry in the table are dependent on
  /// the mathematical type of the quantity and on the name given in column 1.
  ///  - A SCALAR simply corresponds to a field with a name given in column 1.
  ///  - A VECTOR corresponds to 3 fields. The fields are named {column_1}_x,
  ///    {column_1}_y, {column_2}_z
  /// To register a new field, add a new entry to FIELD_TABLE

public:

  EnzoCenteredFieldRegistry() = delete;

  /// Returns a vector of passive scalar group names
  ///
  /// To register new names add entry to the static constant vector variable
  /// called passive_group_names
  static std::vector<std::string> passive_scalar_group_names()
  { return passive_group_names; }

  /// Returns a vector of registered quantities.
  ///
  /// @param enumerate_components Determines whether the individual vector
  ///     components are listed or the just the vector component is listed.
  ///     When true, the function effectively returns a vector containing all
  ///     registered fields. Otherwise, the returned vector containing the
  ///     names of all quantities listed in FIELD_TABLE.
  static std::vector<std::string> get_registered_quantities
  (bool enumerate_components);

  /// provides the quantity properties listed in FIELD_TABLE (if present)
  ///
  /// returns true when successful (i.e. quantity is actually included in the
  /// table) and false when unsuccesful
  static bool quantity_properties(const std::string &name,
                                  bool* vector_quantity = 0,
                                  FieldCat* category = 0,
                                  bool* actively_advected = 0) noexcept;

  /// Determine the actively advected quantity associated with the given name
  /// If there is not an associated quantity, `""` is returned.
  ///
  /// There are 2 ways that a name can be associated with a quantity
  ///   1. They can exactly match the name of the quantity.
  ///   2. For VECTOR quantities, they can be composed of the quantity name
  ///      followed by a 2 character suffix. If ijk_suffix is true, then the
  ///      suffixes are {'_i', '_j', '_k'}. Otherwise they are
  ///      {'_x', '_y', '_z'}.
  static std::string get_actively_advected_quantity_name
  (std::string name, bool ijk_suffix) noexcept;
  
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
  static Grouping* build_grouping(const std::vector<std::string> quan_names,
                                  const std::string leading_prefix = "");

private: // attributes

  /// representation of FIELD_TABLE in memory
  static const ft_map field_table_;
};

#endif /* ENZO_ENZO_CENTERED_FIELD_REGISTRY_HPP */
