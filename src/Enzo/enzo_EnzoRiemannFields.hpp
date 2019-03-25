#ifndef ENZO_ENZO_RIEMANN_FIELDS_HPP
#define ENZO_ENZO_RIEMANN_FIELDS_HPP

// Presently, this where we register all non-passive scalar fields used by the
// Riemann Solver. (Any cell-centered field with a non-trivial flux OR is
// involved in the calculation of the wave speeds
//
// The eventual goal is to make this into a registry for all such fields and
// use it in other areas of the code.


// Scoped Enum to classify quantities as conserved, primitive, or both
enum class FieldCat{
  conserved = 1,
  primitive = 2,
  overlap = 3
};

inline bool category_overlap_(FieldCat a, FieldCat b){  
  return (static_cast<int>(a) & static_cast<int>(b)) > 0;
}


// class used to indicate conditions employed by the integrator. In the
// future, additional members may be added like dual_energy_formalism or
// cosmic_rays
// (it's literally a struct that can be pupped
struct EnzoFieldConditions : public PUP::able
{
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



// Table of fluid properties represented by cell-centered fields.
// The columns of the table are:
//   1. Quantity name.
//   2. member of EnzoFieldConditions that indicates when the quantity is used
//   3. classification of the quantity as primitive, conserved, or both 
//   4. quantity type (SCALAR or VECTOR)
//
// As additional non-passive sacalars are added (e.g. internal energy,
// cosmic ray energy, cosmic ray flux), the table should be extended

#define FIELD_TABLE                                                   \
  ENTRY(     density,        hydro,     FieldCat::overlap, SCALAR)    \
  ENTRY(    momentum,        hydro,   FieldCat::conserved, VECTOR)    \
  ENTRY(    velocity,        hydro,   FieldCat::primitive, VECTOR)    \
  ENTRY(total_energy,        hydro,   FieldCat::conserved, SCALAR)    \
  ENTRY(    pressure,        hydro,   FieldCat::primitive, SCALAR)    \
  ENTRY(      bfield,          MHD,     FieldCat::overlap, VECTOR)

// The current convention (within EnzoMethodVlct) is that the main fields
// associated with each quantity are: 
//    - the same as the group name if the quantity is a scalar
//    - the same as the group name + "_x", "_y", or "_z" if the quantity is a
//      vector.


// The Riemann Solver makes extensive use of the following derivatives of
// FIELD_TABLE:
//
//   - The Riemann Solver makes use of two sets of arrays.
//
//       1. The first set is made up of arrays of instances of EFlt3DArray.
//          Within a given array, each instance of EFlt3DArray encapsulates the
//          data for a different field. The solver uses:
//            A. arrays of left/right reconstructed primitive fields
//            B. arrays of left/right reconstructed conserved fields
//            C. array of fields where the calculated Riemann Flux is stored
//
//       2. The second set is made up of arrays of instances of enzo_float.
//          These arrays serve as temporary storage buffers which store 
//          quantities for a single cell interface. The solver uses
//            A. arrays of left/right reconstructed primitives quantities
//            B. arrays of left/right reconstructed conserved quantities
//            C. arrars of left/right fluxes
//
//     Each of the above arrays only include quantities that are required for
//     the calculations of wave speeds or have non-trivial flux calculations.
//     Passively advected scalars are not included in these arrays (their flux
//     is computed separately).
//
//   - The general code-flow of the RiemannSolver is:
//
//       a For a given cell-interface on a grid, the reconstructed primitive
//         conserved quantites from arrays 1A and 1B into the temporary arrays
//         2A and 2B.
//
//       b The standard MHD fluxes are computed at that location and saved into
//         the arrays left/right fluxes.
//
//       c Optional functions are also applied to compute additional left/right
//         fluxes. Pointers to these functions are specified upon construction
//         of the solver and are used for non-standard fluxes (e.g. cosmic ray
//         energy/fluxes)
//
//       d This step is implemented in a virtual function implemented by a
//         subclass. The wavespeeds at the current interface is computed.
//         Then for each conserved (non-passive scalar) field, the array of
//         Riemann Flux (1C) at the current interface, is set equal to the
//         flux computed from the left/right reconstructed conserved values
//         (2B) and left/right fluxes (2C). This is achieved by iterating over
//         the indices of each array simultaneously.
//
//   - Use of the following field_lut struct:
//
//       - The calculation of fluxes and wave speeds requires random access of
//         specific fields. We also need to be able to iterate over the entries
//         of multiple arrays simultaneously (e.g. to accomplish part d, above)
//         Unlike Enzo, we wanted to avoid statically declaring which indices
//         correspond to which fields (adding additional sets of fields, like
//         internal energy and cosmic ray energy/fluxes becomes harder)
//
//       - We settled on using the field_lut struct as a lookup table. The
//         struct has members named for every quantity listed in FIELD_TABLE
//           - For a SCALAR, the member name directly matches the name in
//             column 1
//           - For a VECTOR, there are 3 members: {name}_i, {name}_j, {name}_k
//             ({name} cooresponds to the name appearing in column 1)
//         Each struct contains members named for all quantities in the table
//         (it includes conserved AND primitive quantites).
//
//       - Example: If we have an array of reconstructed primitives, wl, and
//         an instance of field_lut, prim_lut, that stores the indices of
//         primitives, then wl[prim_lut.density] and wl[prim_lut.pressure]
//         indicates the entries reserved for density and pressure (an
//         instance of field_lut storing indices for conserved quantities
//         would also indicate the index where density - since density is BOTH
//         conserved AND primitive)
//
//       - Given an instance of EnzoFieldConditions, the prepare_conserved_lut
//         function yields an initialized instance of field_lut and the
//         the length necessary for an array to hold values related to
//         conserved quantities. All members of field_lut corresponding to
//         conserved quantities are set equal to consectuive integers starting
//         from 0 (All members that don't correspond to conserved quantites
//         are set to -1).
//
//       - prepare_primitive_lut is analogous to prepare_conserved_lut except
//         it applies to primitive quantities
//     
//       - load_array_of_fields constructs an array of instances of EFlt3DArray
//         that correspond to reconstructed primitive fields, reconstructed
//         conserved fields OR flux fields. The function requires a pointer to
//         an instance of Block, the relevant initialized field_lut, the
//         length of the output field, an instance of grouping (group names
//         must match the relevant quantity names) and the direction along
//         which we are computing fluxes

#define VECTOR_COMPONENT(prefix, suffix) prefix##suffix

#define LUT_MEMBER_SCALAR(name) int name
#define LUT_MEMBER_VECTOR(name) int VECTOR_COMPONENT(name,_i);                \
                                int VECTOR_COMPONENT(name,_j);                \
				int VECTOR_COMPONENT(name,_k)

typedef struct {
  #define ENTRY(name, req_name, category, quantity_type)                      \
    LUT_MEMBER_##quantity_type(name);
  FIELD_TABLE
  #undef ENTRY
  } field_lut;

// For every scalar in FIELD_TABLE, field_lut has a member of the same name.
// For every vector in FIELD_TABLE, field_lut has 3 members. The members use
// the name from a table as a prefix, and have _i, _j, and _k as suffixes.


#define PREPARE_LUT_SCALAR(name, req_name, in_category)                       \
  out.name = (cond.req_name && in_category) ? i++ : -1
#define PREPARE_LUT_VECTOR(name, req_name, in_category)                       \
  out.VECTOR_COMPONENT(name,_i) = (cond.req_name && in_category) ? i++ : -1;  \
  out.VECTOR_COMPONENT(name,_j) = (cond.req_name && in_category) ? i++ : -1;  \
  out.VECTOR_COMPONENT(name,_k) = (cond.req_name && in_category) ? i++ : -1

inline field_lut prepare_lut_(EnzoFieldConditions cond, int *nfields,
			      FieldCat target_cat){
  field_lut out;
  int i = 0;
  #define ENTRY(name, req_name, category, quantity_type)                      \
    PREPARE_LUT_##quantity_type(name, req_name,                               \
                                category_overlap_(category, target_cat));
  FIELD_TABLE
  #undef ENTRY
  *nfields = i;
  return out;
}

// prepares lookup table for conserved quantites. (The members of field_lut
// that correspond to conserved quantites are set to an index of an array - the
// any other ).
inline field_lut prepare_conserved_lut(EnzoFieldConditions cond, int *nfields)
{  return prepare_lut_(cond, nfields, FieldCat::conserved); }


inline field_lut prepare_primitive_lut(EnzoFieldConditions cond, int *nfields)
{  return prepare_lut_(cond, nfields, FieldCat::primitive); }


#define USE_QUANTITY_SCALAR(lut, name) (lut.name > -1)
#define USE_QUANTITY_VECTOR(lut, name) (lut.VECTOR_COMPONENT(name,_i) > -1)

inline int load_array_of_fields_(EFlt3DArray* arr, int cur_index,
				 Grouping &grouping, std::string group_name,
				 std::string quantity_type, int dim,
				 EnzoPermutedCoordinates coord,
				 EnzoFieldArrayFactory array_factory)
{
  // Sanity Check:
  int group_size = grouping.size(group_name);

  ASSERT("load_fluid_fields_",
	 "Groups of fields don't have the correct number of fields.",
	 (quantity_type == "VECTOR" && group_size == 3) ||
	 (quantity_type == "SCALAR" && group_size == 1));

  if (quantity_type == "VECTOR"){
    arr[cur_index] = array_factory.reconstructed_field(grouping, group_name,
						       coord.i_axis(), dim);
    arr[cur_index+1] = array_factory.reconstructed_field(grouping,
							 group_name,
							 coord.j_axis(),
							 dim);
    arr[cur_index+2] = array_factory.reconstructed_field(grouping,
							 group_name,
							 coord.k_axis(),
							 dim);
    return cur_index+3;
  } else {

    arr[cur_index] = array_factory.reconstructed_field(grouping, group_name,
						       0, dim);
    return cur_index+1;
  }
}

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


inline EFlt3DArray* load_array_of_fields(Block *block, const field_lut lut,
					 const int nfields, Grouping &grouping,
					 const int dim)
{
  EFlt3DArray* arr = new EFlt3DArray[nfields];
  int cur_index = 0;
  EnzoPermutedCoordinates coord(dim);
  EnzoFieldArrayFactory array_factory(block);

  #define ENTRY(name, req_name, category, quantity_type)                      \
    if (USE_QUANTITY_##quantity_type(lut, name)) {		              \
      std::string group_name = std::string( #name );                          \
      std::string quantity_type = std::string( #quantity_type );              \
      cur_index = load_array_of_fields_(arr, cur_index, grouping, group_name, \
                                        quantity_type, dim, coord,            \
				        array_factory);		              \
  }
  FIELD_TABLE
  #undef ENTRY
  return arr;
}

#endif /* ENZO_ENZO_RIEMANN_FIELDS_HPP */
