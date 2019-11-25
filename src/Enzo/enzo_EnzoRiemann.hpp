// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemann.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of the Riemann Solver abstract base
/// class. This class should be subclassed to implement various riemann solvers.

#ifndef ENZO_ENZO_RIEMANN_HPP
#define ENZO_ENZO_RIEMANN_HPP



class EnzoRiemann : public PUP::able
{
  /// @class    EnzoRiemann
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate approximate Riemann Solvers

public: // interface

  /// Factory method for constructing the EnzoRiemann object. (The signature
  /// may need to be modified as additional physics get added)
  ///
  /// @param integrable_groups A vector of integrable quantities (listed as
  ///     advected quantities in FIELD_TABLE). These are used as group names in
  ///     the Grouping objects that store field names. In effect this is used
  ///     to register the quantities operated on by the Riemann Solver
  /// @param passive_groups A vector with the names of the groups of passively
  ///     advected scalars that may be included. (If a group is listed here but
  ///     the Grouping object doesn't actually provide any fields in the group,
  ///     no problems are caused)
  /// @param solver The name of the Riemann solver to use. Valid names include
  ///     "hll", "hlle", and "hlld"
  static EnzoRiemann* construct_riemann
    (std::vector<std::string> integrable_groups,
     std::vector<std::string> passive_groups, std::string solver);

  EnzoRiemann() throw()
  {}

  /// Virtual destructor
  virtual ~EnzoRiemann()
  {}

  /// CHARM++ PUP::able declaration
  PUPable_abstract(EnzoRiemann);

  /// CHARM++ migration constructor for PUP::able
  EnzoRiemann (CkMigrateMessage *m)
    : PUP::able(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    PUP::able::pup(p);
  }

  /// Computes the Riemann Fluxes for each conserved field along a given
  /// dimension, dim
  /// @param block holds data to be processed
  /// @param priml_group,primr_group holds field names where the left/right
  ///     reconstructed face-centered integrable primitives are stored. The
  ///     relevant fields should be formally defined as cell-centered (to allow
  ///     for reuse). During the calculation, they are treated as face-centered
  ///     (without having values on the exterior faces of the block). As a
  ///     result, there is some unused space at the end of the arrays.
  /// @param pressure_name_l,pressure_name_r are the names of the fields
  ///     storing the left/right pressure values that have already been
  ///     computed from the reconstructed values. The face-centering is expected
  ///     to match fields contained by priml_group, primr_group
  /// @param flux_group holds field names where the calculated fluxes will be
  ///     stored. The relevant fields should be face-centered along the
  ///     dimension of the calculation (without having values on the exterior
  ///     faces of the block)
  /// @param dim Dimension along which to compute Riemann fluxes. Values of 0,
  ///     1, and 2 correspond to the x, y, and z directions, respectively.
  /// @param eos Instance of the fluid's EnzoEquationOfState object
  /// @param stale_depth indicates the number of field entries from the
  ///     outermost field value that the region including "stale" values (need
  ///     to be refreshed) extends over (0 means there are no "stale" values).
  /// @param interface_velocity_name indicates the name of field where the
  ///     value of the interface velocity along dimension `dim` should be
  ///     stored. This quantity is used to compute the internal energy source
  ///     term (needed under the dual energy formalism). If the value is `""`
  ///     (the default) then the interface velocity is not stored.
  virtual void solve (Block *block,
		      Grouping &priml_group, Grouping &primr_group, 
		      std::string pressure_name_l, std::string pressure_name_r,
		      Grouping &flux_group, int dim, EnzoEquationOfState *eos,
		      int stale_depth,
		      std::string interface_velocity_name = "") const = 0;

  /// computes the (adiabatic) fast magnetosonic speed along dimension i at a
  /// single location (vector names in lut are denoted as being along the i, j
  /// and k directions)
  ///
  /// @param prim_vals array containing integrable primitives
  /// @param lut An instance EnzoAdvectionFieldLUT that serves as a look up
  ///     table. It has members named for each relevent field (although
  ///     vectors have the suffixes _i, _j, _k instead of _x, _y, _z). The
  ///     value of a given member corresponds to the index in prim_vals where
  ///     the value of the field is stored.
  /// @param pressure The value of the thermal pressure. (This is never
  ///     considered an integrable primitive)
  /// @param gamma The adiabatic index
  static enzo_float fast_magnetosonic_speed_(const enzo_float prim_vals[],
					     const EnzoAdvectionFieldLUT lut,
					     enzo_float pressure,
					     const enzo_float gamma) throw()
  {
    enzo_float bi = prim_vals[lut.bfield_i];
    enzo_float bj = prim_vals[lut.bfield_j];
    enzo_float bk = prim_vals[lut.bfield_k];

    // TODO: optimize calc of cs2 to omit sqrt and pow
    enzo_float cs2 = std::pow(sound_speed_(prim_vals, lut, pressure, gamma),2);
    enzo_float B2 = (bi*bi + bj*bj + bk *bk);
    if (B2 == 0){
      return std::sqrt(cs2);
    }
    enzo_float va2 = B2/prim_vals[lut.density];
    // TODO: replace va2 * cos2 with va2_cos2 = bi*bi/prim_vals[lut.density]
    enzo_float cos2 = bi*bi / B2;
    return std::sqrt(0.5*(va2+cs2+std::sqrt(std::pow(cs2+va2,2) -
					    4.*cs2*va2*cos2)));
  }

  /// @overload
  ///
  /// This function should be called when we to fix cos2 to some predetermined
  /// value
  static enzo_float fast_magnetosonic_speed_(const enzo_float prim_vals[],
					     const EnzoAdvectionFieldLUT lut,
					     enzo_float pressure,
					     const enzo_float gamma,
					     const enzo_float cos2) throw()
  {
    enzo_float bi = prim_vals[lut.bfield_i];
    enzo_float bj = prim_vals[lut.bfield_j];
    enzo_float bk = prim_vals[lut.bfield_k];

    enzo_float cs2 = std::pow(sound_speed_(prim_vals, lut, pressure, gamma),2);
    enzo_float B2 = (bi*bi + bj*bj + bk *bk);
    enzo_float va2 = B2/prim_vals[lut.density];
    return std::sqrt(0.5*(va2+cs2+std::sqrt(std::pow(cs2+va2,2) -
					    4.*cs2*va2*cos2)));
  }

  /// computes the (adiabatic) magnetic pressure at a single location
  ///
  /// @param prim_vals array containing integrable primitives
  /// @param lut An instance EnzoAdvectionFieldLUT that serves as a look up
  ///     table. It has members named for each relevent field (although
  ///     vectors have the suffixes _i, _j, _k instead of _x, _y, _z). The
  ///     value of a given member corresponds to the index in prim_vals where
  ///     the value of the field is stored.
  static enzo_float mag_pressure_(const enzo_float prim_vals[],
				  const EnzoAdvectionFieldLUT lut) throw()
  {
    enzo_float bi = prim_vals[lut.bfield_i];
    enzo_float bj = prim_vals[lut.bfield_j];
    enzo_float bk = prim_vals[lut.bfield_k];
    return 0.5 * (bi*bi + bj*bj + bk *bk);
  }

  /// computes the (adiabatic) sound speed at a single location
  ///
  /// @param prim_vals array containing integrable primitives
  /// @param lut An instance EnzoAdvectionFieldLUT that serves as a look up
  ///     table. It has members named for each relevent field (although
  ///     vectors have the suffixes _i, _j, _k instead of _x, _y, _z). The
  ///     value of a given member corresponds to the index in prim_vals where
  ///     the value of the field is stored.
  /// @param pressure The value of the thermal pressure. (This is never
  ///     considered an integrable primitive)
  /// @param gamma The adiabatic index
  static enzo_float sound_speed_(const enzo_float prim_vals[],
				 const EnzoAdvectionFieldLUT lut,
				 const enzo_float pressure,
				 const enzo_float gamma) throw()
  { return std::sqrt(gamma * pressure / prim_vals[lut.density]); }

};

#endif /* ENZO_ENZO_RIEMANN_HPP */
