// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoCenteredFieldRegistry.hpp
/// @author   Matthew W. Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs March 28 2019
/// @brief    [\ref Enzo] Implementation of the EnzoFieldConditions class,
///           declaration of the field_lut struct and declaration of the
///           EnzoCenteredFieldRegistry class.
///
///
/// This class serves to track all cell-centered MHD fields that are not
/// passively advected scalars and are not specific to a single integration
/// method. The class also tracks the group names of passively advected
/// scalars. The API supports constructing Groupings of conserved and primitive
/// quantities, returning vectors of conserved, primitive, or passively
/// advected group names, preparing instances of field_lut to serve as lookup
/// tables in EnzoRiemann for primitive or conserved quantities, and to
/// construct an array of Field Arrays appropriate for an instance of field_lut.
/// of primitive

#ifndef ENZO_ENZO_CENTERED_FIELD_REGISTRY_HPP
#define ENZO_ENZO_CENTERED_FIELD_REGISTRY_HPP

struct EnzoFieldConditions : public PUP::able
{
  /// @class    EnzoFieldConditions
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] PUPable struct with members that are used to
  ///           identify which fields are required from the registry. Additional
  ///           should be added in the future, like dual_energy_formalism or
  ///           cosmic rays as new physics is added.
  
  bool hydro;
  bool MHD;

  EnzoFieldConditions() = default;

  /// CHARM++ migration constructor for PUP::able
  EnzoFieldConditions (CkMigrateMessage *m)
    : PUP::able(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    PUP::able::pup(p);
    p|hydro;
    p|MHD;
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoFieldConditions);
};


/// Enum that classifies a quantity as conserved, primitive, or both
enum class FieldCat{
  conserved = 1,
  primitive = 2,
  overlap = 3
};

/// Checks to see if one instance of the FieldCat overlaps with the category of
/// another instance
inline bool category_overlap_(FieldCat a, FieldCat b){  
  return (static_cast<int>(a) & static_cast<int>(b)) > 0;
}


// list of passively advected group names
const std::vector<std::string> passive_group_names = {"colors", "species"};

// Table of fluid properties represented by cell-centered fields. This table
// must be updated as new physics gets added. The columns of the table are:
//   1. Quantity name.
//   2. member of EnzoFieldConditions that indicates if the quantity is needed
//   3. classification of the quantity as primitive, conserved, or both 
//   4. quantity type (SCALAR or VECTOR)
//
// The entries in this table determine the name of the field to represent the
// quantity:
//    - A SCALAR quantity is represented by a field with the name in column 1.
//    - A VECTOR quantity is represented by 3 fields. The fields are named
//      {column_1}_x, {column_1}_y, {column_2}_z
//
// As additional physics are added this table should be updated (e.g.
// internal_energy, cr_energy, cr_flux).

#define FIELD_TABLE                                                   \
  ENTRY(     density,        hydro,     FieldCat::overlap, SCALAR)    \
  ENTRY(    momentum,        hydro,   FieldCat::conserved, VECTOR)    \
  ENTRY(    velocity,        hydro,   FieldCat::primitive, VECTOR)    \
  ENTRY(total_energy,        hydro,   FieldCat::conserved, SCALAR)    \
  ENTRY(    pressure,        hydro,   FieldCat::primitive, SCALAR)    \
  ENTRY(      bfield,          MHD,     FieldCat::overlap, VECTOR)

  
#define VECTOR_COMPONENT(prefix, suffix) prefix##suffix

#define LUT_MEMBER_SCALAR(name) int name
#define LUT_MEMBER_VECTOR(name) int VECTOR_COMPONENT(name,_i);                \
                                int VECTOR_COMPONENT(name,_j);                \
				int VECTOR_COMPONENT(name,_k)

// A struct that serves as a fast lookup table for EnzoRiemann. For every
// scalar in FIELD_TABLE, field_lut has a member of the same name. For every
// vector in FIELD_TABLE, field_lut has 3 members, named {col_1}_i, {col_1}_j,
// and {col_1}_k.
//
// When in use, each relevant member get's dynamically assigned an consecutive
// int starting from 0 which corresponds to a the quantity of the same name's
// location in an array. All members corresponding to non-relevant quantities
// are set to -1.
typedef struct {
  #define ENTRY(name, req_name, category, quantity_type)                      \
    LUT_MEMBER_##quantity_type(name);
  FIELD_TABLE
  #undef ENTRY
  } field_lut;

#define PRINT_SCALAR(lut,name) CkPrintf("%s: %d\n", #name, lut.name)
#define PRINT_VECTOR(lut,name)                                                \
  CkPrintf("%s: %d\n", (#name "_i"), lut.VECTOR_COMPONENT(name,_i));          \
  CkPrintf("%s: %d\n", (#name "_j"), lut.VECTOR_COMPONENT(name,_j));          \
  CkPrintf("%s: %d\n", (#name "_k"), lut.VECTOR_COMPONENT(name,_k))

inline void print_lut(const field_lut lut)
{
  #define ENTRY(name, req_name, category, quantity_type)                      \
    PRINT_##quantity_type( lut, name);
  FIELD_TABLE
  #undef ENTRY
 }


class EnzoCenteredFieldRegistry
{
  /// @class    EnzoCenteredFieldRegistry
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Serves as a registry for centered fields
  ///           representing non-passively advected primitive or conserved
  ///           quantities

public:

  /// Constructs a Grouping of conserved quantites. The names of each group
  /// matches the first column from FIELD_TABLE. The rules for field names are:
  ///    - Scalar conserved quantity:
  ///        - leading_prefix + (overlap_cat_prefix +) {col_1}
  ///    - Vector conserved quantity:
  ///        - leading_prefix + (overlap_cat_prefix +) {col_1} + "_x"
  ///        - leading_prefix + (overlap_cat_prefix +) {col_1} + "_y"
  ///        - leading_prefix + (overlap_cat_prefix +) {col_1} + "_z"
  /// overlap_cat_prefix is only used if a quantity is both conserved and
  /// primitive
  Grouping* build_cons_grouping(const EnzoFieldConditions cond,
				const std::string leading_prefix = "",
				const std::string overlap_cat_prefix =""){
    return build_grouping_(cond, leading_prefix, overlap_cat_prefix,
			   FieldCat::conserved);
  }

  /// Constructs a Grouping of primitive quantites. The names of each group
  /// matches the first column from FIELD_TABLE. The rules for field names are
  /// the same as for build_cons_grouping
  Grouping* build_prim_grouping(const EnzoFieldConditions cond,
				const std::string leading_prefix ="",
				const std::string overlap_cat_prefix =""){
    return build_grouping_(cond, leading_prefix, overlap_cat_prefix,
			   FieldCat::primitive);
  }

  /// Returns the vector of conserved group names
  std::vector<std::string> cons_group_names(const EnzoFieldConditions cond,
					    bool include_passive = false){
    return group_names_(cond, include_passive, FieldCat::conserved);
  }

  /// Returns the vector of primitive group names
  std::vector<std::string> prim_group_names(const EnzoFieldConditions cond,
					    bool include_passive = false){
    return group_names_(cond, include_passive, FieldCat::primitive);
  }

  /// Returns the vector of passive scalar group names
  std::vector<std::string> passive_scalar_group_names()
  {
    return passive_group_names;
  }

  /// Prepares an instance of field_lut for conserved quantities and determines
  /// number of entries needed in an array to hold this information
  field_lut prepare_cons_lut(const EnzoFieldConditions cond, int* nfields){
    return prepare_lut_(cond, nfields, FieldCat::conserved);
  }

  /// Prepares an instance of field_lut for primitive quantities and determines
  /// number of entries needed in an array to hold this information
  field_lut prepare_prim_lut(const EnzoFieldConditions cond, int* nfields){
    return prepare_lut_(cond, nfields, FieldCat::primitive);
  }

  /// Constructs an array of instances of EFlt3DArray where each instance holds
  /// the field data of the corresponding field in lut.
  EFlt3DArray* load_array_of_fields(Block *block, const field_lut lut,
				    const int nfields, Grouping &grouping,
				    const int dim);

private:
  Grouping* build_grouping_(const EnzoFieldConditions cond,
			    const std::string leading_prefix,
			    const std::string overlap_cat_prefix,
			    FieldCat target_cat);

  // Helper method of build_grouping_
  void add_group_fields_(Grouping *grouping, const std::string group_name,
			 const std::string quantity_type, FieldCat category,
			 const std::string leading_prefix,
			 const std::string overlap_cat_prefix);

  std::vector<std::string> group_names_(const EnzoFieldConditions cond,
					bool include_passive,
					FieldCat target_cat);
  
  field_lut prepare_lut_(const EnzoFieldConditions cond, int *nfields,
			 FieldCat target_cat);

  // Helper method of load_array_of_fields
  int load_array_of_fields_(EFlt3DArray* arr, int cur_index, Grouping &grouping,
			    std::string group_name, std::string quantity_type,
			    int dim, EnzoPermutedCoordinates coord,
			    EnzoFieldArrayFactory array_factory);
};

#endif /* ENZO_ENZO_CENTERED_FIELD_REGISTRY_HPP */
