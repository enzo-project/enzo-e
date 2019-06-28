// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoCenteredFieldRegistry.hpp
/// @author   Matthew W. Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs March 28 2019
/// @brief    [\ref Enzo] Implementation of the EnzoAdvectionFieldLUT and
///           declaration of the EnzoCenteredFieldRegistry class.
///
///
/// This class serves to track all cell-centered MHD quantities that are not
/// passively advected scalars. The class also tracks the group names of
/// passively advected scalars. The API supports constructing Groupings of
/// fields (with names matching the the names of registered names), determining
/// which quantities are conserved, specific, or other, preparing instances of
/// EnzoAdvectionFieldLUT to serve as lookup tables in EnzoRiemannImpl for and
/// to construct an array of Field Arrays appropriate for a given instance of
/// EnzoAdvectionFieldLUT

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
  ENTRY(        internal_energy, SCALAR,  FieldCat::specific, T)             \
  ENTRY(           total_energy, SCALAR,  FieldCat::specific, T)             \
  ENTRY(                 bfield, VECTOR, FieldCat::conserved, T)             \
  /* Derived Fields */                                                       \
  ENTRY(               pressure, SCALAR,     FieldCat::other, F)             \
  ENTRY(            temperature, SCALAR,     FieldCat::other, F)             \
  ENTRY(           cooling_time, SCALAR,     FieldCat::other, F)             \
  ENTRY(           acceleration, VECTOR,  FieldCat::specific, F)             \
  /* Grackle Related Fields */                                               \
  ENTRY(  specific_heating_rate, SCALAR,     FieldCat::other, F)             \
  ENTRY(volumetric_heating_rate, SCALAR,     FieldCat::other, F)             \
  /* Field for Gravity */                                                    \
  ENTRY(              potential, SCALAR,  FieldCat::specific, F)


// Yields a combined token
#define COMBINE(prefix, suffix) prefix##suffix
// Yields a stringized version of the combined token
#define STRINGIZE_COMBINE(prefix,suffix) STR_C_1(COMBINE(prefix,suffix))
#define STR_C_1(s) STR_C_2(s)
#define STR_C_2(s) #s

//----------------------------------------------------------------------
/// @def      ADVEC_STRUCT_MEMBER
/// @brief    Macro that helps define advection field struct members based off
///           the entries in FIELD_TABLE
///
/// This includes 4 macros:
/// STRUCT_MEMBER_T_SCALAR, STRUCT_MEMBER_T_VECTOR,
/// STRUCT_MEMBER_F_SCALAR, STRUCT_MEMBER_F_VECTOR
/// which have specialized behavior based on whether the quantity is advective
/// and if its a SCALAR or VECTOR quantities 
#define ADVEC_STRUCT_MEMBER_T_SCALAR(name,type) type name
#define ADVEC_STRUCT_MEMBER_T_VECTOR(name,type)                               \
  type COMBINE(name,_i);   type COMBINE(name,_j);  type COMBINE(name,_k)
#define ADVEC_STRUCT_MEMBER_F_SCALAR(name,type) /* ... */
#define ADVEC_STRUCT_MEMBER_F_VECTOR(name,type) /* ... */


//----------------------------------------------------------------------
/// @def      ADVEC_STRUCT_UNARY_FUNC
/// @brief    Macro that applies a unary function to members of a struct named
///           for advection related quantities in FIELD_TABLE
#define ADVEC_STRUCT_UNARY_FUNC_T_SCALAR(func, struc, name)                   \
  func(#name, struc.name)
#define ADVEC_STRUCT_UNARY_FUNC_T_VECTOR(func, struc, name)                   \
  func(STRINGIZE_COMBINE(name,_i), struc.COMBINE(name,_i));                 \
  func(STRINGIZE_COMBINE(name,_j), struc.COMBINE(name,_j));                 \
  func(STRINGIZE_COMBINE(name,_k), struc.COMBINE(name,_k))
#define ADVEC_STRUCT_UNARY_FUNC_F_SCALAR(func, struc, name) /* ... */
#define ADVEC_STRUCT_UNARY_FUNC_F_VECTOR(func, struc, name) /* ... */


//----------------------------------------------------------------------
/// template helper function that applies unary function on the members of
/// struct that have been named for advection related quantities in FIELD_TABLE
/// @param obj Reference to a class/struct with members named after advection
/// quantities in FIELD_TABLE
/// @param fn Unary function or that accepts the name of the struct member and
/// the value of the member as arguments
template<class AdvectStruct, class Function>
void unary_advec_struct_for_each_(AdvectStruct &obj, Function fn){
  #define ENTRY(name, math_type, category, if_advection)                      \
      ADVEC_STRUCT_UNARY_FUNC_##if_advection ## _ ## math_type (fn, obj, name);
  FIELD_TABLE
  #undef ENTRY
}

class EnzoAdvectionFieldLUT{

  /// @class    EnzoAdvectionFieldLUT
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Serves as a lookup table for fluid advection
  ///           calculations
  ///
  /// This is supposed to be a bare bones class - exactly like a C-struct with
  /// members named for different advection quantites listed in FIELD_TABLE
  /// that each hold integer values. The mathematical type of the quantity
  /// determines the name of the member. A SCALAR simply corresponds to a
  /// member with a name matching the value in column 1. A VECTOR corresponds
  /// members named {column_1}_i, {column_1}_j, {column_2}_k.
  ///
  /// Basically instances of this object are supposed to be used in advection
  /// based calculations (like computing Riemann fluxes) to serve as a look up
  /// table. Basically one would use EnzoCenteredFieldRegistry to help
  /// construct an array of fields and an instance of this class. Then the
  /// value indicated by a given member of the instance would indicate the
  /// index of the array where the corresponding field value would be stored.
  /// (E.g. a density value would be stored at index obj.density).
  ///
  /// Members that don't have correspond to a value are set to -1.

public: // interface
  
  /// Default constructor
  EnzoAdvectionFieldLUT() = default;

  /// Destructor
  ~EnzoAdvectionFieldLUT()
  {  }
  
public: // attributes

  // Add attributes named for entries in LUT
  #define ENTRY(name, math_type, category, if_advection)                      \
    ADVEC_STRUCT_MEMBER_ ## if_advection ## _ ## math_type(name, int);
  FIELD_TABLE
  #undef ENTRY
};


//----------------------------------------------------------------------
/// utility function that prints the names of the members and the values of an
/// instance of EnzoAdvectionFieldLUT.
///
/// This is provided to assist with debugging, especially gdb encounters
/// problems with instances of EnzoAdvectionFieldLUT (presumably, because the
/// members are specified with macros
void print_lut(const EnzoAdvectionFieldLUT lut);

//----------------------------------------------------------------------
/// utility function for debugging that prints the entries in array that are
/// labelled by an instance of EnzoAdvectionFieldLUT.
void print_looked_up_vals(const EnzoAdvectionFieldLUT lut,
			  const enzo_float* array);

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

  // some utility methods - some of these are called in more important methods
  // These are not particularly well optimized

  /// Returns a vector of registered quantities
  std::vector<std::string> get_registered_quantities() const;

  /// Returns the vector of registered field names
  std::vector<std::string> get_registered_fields() const;
  
  /// Checks that that the quantity names are in FIELD_TABLE. If not then,
  /// raises an error.
  void check_known_quantity_names(const std::vector<std::string> names) const;

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

  /// Prepares an instance of EnzoAdvectionFieldLUT for a specified set of
  /// advection related quantities.
  ///
  /// It also yields the start and stop indices to iterate over all non-flagged
  /// FieldCat::conserved quantities, FieldCat::specific quantities and
  /// FieldCat::other quanties. Additionally, it determines the size of the
  /// array necessary all of the relevant field values (including the flagged
  /// quantites). Note that quantity names are used as group names in other
  /// contexts.
  ///
  /// @param quantity_names Vector of the quantity names that should be
  ///   included in the table. These names should match the names of entries
  ///   in FIELD_TABLE. For example it might include density, velocity
  ///   (not velocity_x), etc.
  /// @param conserved_start,conserved_stop References values that will be
  ///   modified to indicate the start and stop indices for iterating over all
  ///   unflagged conserved quantities.
  /// @param specific_start,specific_stop References values that will be
  ///   modified to indicate the start and stop indices for iterating over all
  ///   unflagged specific quantities.
  /// @param other_start,other_stop References values that will be modified to
  ///   indicate the start and stop indices for iterating over all unflagged
  ///   other quantities.
  /// @param nfields Number of fields that the array will have to hold all
  ///   specified quantities
  /// @param flagged_quantities Optional argument. Vector of flagged quantity
  ///   names. A flagged group name will not be included in the yielded
  ///   iteration ranges (An example of where this might be desirable that when
  ///   applying fluxes in an integrator using Constrained Transport, we don't
  ///   want to add the Bfield fluxes). The name of the flagged quantity must
  ///   also be included in quantity names. If this is not passed, then no
  ///   quantites are flagged
  EnzoAdvectionFieldLUT prepare_advection_lut
  (const std::vector<std::string> quantity_names,
   int &conserved_start, int &conserved_stop, int &specific_start,
   int &specific_stop, int &other_start, int &other_stop, int &nfields,
   const std::vector<std::string> flagged_quantities) const;


  EnzoAdvectionFieldLUT prepare_advection_lut
  (const std::vector<std::string> quantity_names,
   int &conserved_start, int &conserved_stop, int &specific_start,
   int &specific_stop, int &other_start, int &other_stop, int &nfields) const
  {
    std::vector<std::string> flagged_quantities;
    return prepare_advection_lut(quantity_names,
				 conserved_start, conserved_stop,
				 specific_start, specific_stop,
				 other_start, other_stop, nfields,
				 flagged_quantities);
  }
  
  /// Returns a vector of passive scalar group names
  ///
  /// To register new names add entry to the static constant vector variable
  /// called passive_group_names
  std::vector<std::string> passive_scalar_group_names() const
  {
    return passive_group_names;
  }

  /// Constructs an array of instances of EFlt3DArray where each instance holds
  /// the field data of the corresponding field in lut.
  ///
  /// @param block The mesh block from which data will be loaded
  /// @param lut,nfields A look up table and the size of the output array
  ///   determined by a previous call to of the prepare_advection_lut method.
  ///   These which fields are loaded and the indices at which the field data
  ///   is stored.
  /// @param grouping Reference to a grouping that with groups named after the
  ///   relevant quantities specified in the generation of lut. The data is
  ///   actually loaded from the fields in block with names matching the
  ///   relevant field names specified in this object.
  /// @param dim Optional integer specifying which dimension is the ith
  ///   direction. This is used for mapping the i,j,k vector components listed
  ///   in lut to the x,y,z field components. Values of 0, 1, and 2 correspond
  ///   the ith direction pointing parallel to the x, y, and z directions,
  ///   respectively. Note that each of the fields in grouping are assumed to
  ///   be face-centered along this dimension (excluding the exterior faces of
  ///   the mesh). This allows for appropriate loading of reconstructed fields.
  /// @param stale_depth indicates the current stale_depth for the loaded
  ///   quanties.
  ///
  /// 
  EFlt3DArray* load_array_of_fields(Block *block,
				    const EnzoAdvectionFieldLUT lut,
				    const int nfields, Grouping &grouping,
				    const int dim,
				    const int stale_depth) const;

  /// @overload
  ///
  /// The fields are assumed not reconstructed in this case (the stored shape
  /// of the fields is accurate). The components i,j,k of a vector listed by
  /// the lut always map to x,y,z
  EFlt3DArray* load_array_of_fields(Block *block,
				    const EnzoAdvectionFieldLUT lut,
				    const int nfields, Grouping &grouping,
				    const int stale_depth) const
  {
    return load_array_of_fields(block, lut, nfields, grouping,
				-1, stale_depth);
  }

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
};

#endif /* ENZO_ENZO_CENTERED_FIELD_REGISTRY_HPP */
