// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_GrackleFacade.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Feb 27 2023
/// @brief    [\ref Enzo] Declaration of GrackleFacade class
///
/// This class provides a facade to the Grackle library. It is responsible for
/// managing the lifetime of all associated objects

#ifndef ENZO_ENZO_GRACKLE_FACADE_HPP
#define ENZO_ENZO_GRACKLE_FACADE_HPP

enum class GracklePropertyEnum
  { cooling_time, dust_temperature, gamma, pressure, temperature };

class GrackleFacade : public PUP::able {

  /// @class    GrackleFacade
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Grackle's state and operations
  ///
  /// This class is named for the Facade design pattern (the name could be
  /// improved). Instances of this class manage Grackle's state, hide away
  /// certain details (e.g. users of the class never have to worry about the
  /// ``code_units`` or ``chemistry_data_storage``) and provides a clean
  /// interface for each of Grackle's operations, which can be categorized as:
  ///    - solving chemistry over an integration timestep
  ///    - computing chemistry-related (e.g. cooling time) or thermodynamic
  ///      quantities (e.g. pressure, temperature)
  /// This class is somewhat analogous to the python extension class provided
  /// by Grackle, but it enforces stronger invariants.
  ///
  /// The key invariant of this class is that an instance is ALWAYS properly
  /// initialized (i.e. all Grackle structs are properly configured). At the
  /// end of the description, we include a note detailing how this is
  /// accomplished with deserialization.
  ///
  /// To simplify this class's semantics and to help enforce the above
  /// invariant, the attributes of a given instance are treated as immutable
  /// after construction or deserialization (see below for more discussion).
  /// This is hard to enforce with the type-system (especially while supporting
  /// ``pup`` routines), so it is accomplished with the API.
  ///
  /// This choice of immutability has two subtleties that deserve further
  /// elaboration. Both of these relate to the other structs involved in
  /// initializing the ``chemistry_data_storage`` struct. They include:
  ///
  ///   1. Certain parameters within the ``chemistry_data`` struct (managed by
  ///      the ``GrackleChemistryData`` class) are overwritten. For example,
  ///      ``HydrogenFractionByMass`` is always overwritten when primordial
  ///      chemistry is ``0``. Likewise, ``UVbackground_redshift_off`` and
  ///      ``UVbackground_redshift_drop`` also can be overwritten. While
  ///      nothing is currently done about these parameters, if we want to
  ///      introduce handling for them, we should make sure to build-in the
  ///      adjustments to these values of ``GrackleChemistryData`` to happen
  ///      right after ``chemistry_data_storage`` is initialized.
  ///
  ///   2. As a general rule, the fields in the ``code_units`` struct shouldn't
  ///      change after initializing ``chemistry_data_storage`` (lots of data
  ///      in the storage struct are stored in terms of the specified units).
  ///      The notable exception is that the ``a_value`` field must be kept up
  ///      to date when using comoving-coordinates. To maintain our general
  ///      rule about immutability, we generate a new copy of ``code_units``
  ///      every time we need it. The overhead of this should be minimal.
  ///      However, if we want to further reduce the overhead, we can instead
  ///      mutate a copy of the stored ``code_units`` struct in the future
  ///      (since ``code_units`` is small, a ``memcpy`` will be cheap!)
  ///
  /// @note
  /// (De)serialization poses a challenge for maintaining the invariant that
  /// instances of this class are always in a valid state. In practice, we
  /// actually have the Charm++ migration constructor construct this object in
  /// an invalid state. However, this is ok because that constructor is
  /// GUARANTEED to be followed by an unpacking call to the pup method.
  ///   - This is actually the ONLY reason why this class has been declared as
  ///     PUPable
  ///   - Alternatively, we'd have to support a default constructor that could
  ///     be used at any time


public: // interface

  static bool linked_against_grackle() noexcept;

  /// public-facing constructor
  GrackleFacade(GrackleChemistryData&& my_chemistry,
                const double radiation_redshift,
                const double physics_cosmology_initial_redshift,
                const double time) noexcept;

  /// CHARM++ PUP::able declaration
  PUPable_decl(GrackleFacade);

  /// CHARM++ migration constructor for PUP::able
  GrackleFacade(CkMigrateMessage *m);

  /// CHARM++ Pack / Unpack function
  void pup(PUP::er &p);

  /// Destructor
  ~GrackleFacade() noexcept;

  /// delete default constructor
  GrackleFacade() = delete;

  // delete copy/move constructor and assignment (could be added later)
  GrackleFacade(GrackleFacade&&) = delete;
  GrackleFacade& operator=(GrackleFacade&&) = delete;
  GrackleFacade(const GrackleFacade&) = delete;
  GrackleFacade& operator=(const GrackleFacade&) = delete;

  /// returns a pointer to the stored instance of GrackleChemistryData, if the
  /// simulation is configured to actually use grackle
  ///
  /// @note
  /// We may be able to simplify this implementation. I'm pretty sure that
  /// Grackle MUST be used for an instance of GrackleFacade to exist at all
  /// (in which case, we never need to worry about returning nullptr)
  const GrackleChemistryData* try_get_chemistry() const noexcept
  {
    bool not_null = (GrackleFacade::linked_against_grackle() &&
                     (my_chemistry_.get<int>("use_grackle") == 1));
    return (not_null) ? &my_chemistry_ : nullptr;
  }

public: // low-level legacy methods - these will (probably) be removed or made
        // private in the near future

  void setup_grackle_fields(const EnzoFieldAdaptor& fadaptor,
                            grackle_field_data * grackle_fields,
                            int stale_depth = 0,
                            bool omit_cell_width = false) const noexcept;

  void setup_grackle_fields(Block * block, grackle_field_data* grackle_fields,
                            int i_hist = 0 ) const noexcept
  {
    setup_grackle_fields(EnzoFieldAdaptor(block, i_hist),
                         grackle_fields, 0, false);
  }

  void delete_grackle_fields(grackle_field_data* grackle_fields) const noexcept;

public: // wrapped grackle functions:

  /// light-weight wrapper around local_solve_chemistry function from grackle
  ///
  /// @param[in] block Holds the field data that will be passed through to
  ///     Grackle (and will be mutated).
  /// @param[in] compute_time The nominal simulation time at which this
  ///     evaluation occurs. This only matters in cosmological simulations.
  /// @param[in] dt The integration timestep in code units
  ///
  /// @note
  /// This function nominally modifies the internal energy field and the
  /// species density fields. But, based on the inputs and how Grackle is
  /// configured, other fields may be modified as well. This can happen in 1 of
  /// 2 ways:
  ///    1. Internally, Grackle may perform some transformations in place (e.g.
  ///       converting comoving values to proper values), and while it always
  ///       reverts the transformation before returning, floating point errors
  ///       could lead to slightly different field values at the end.
  ///    2. Grackle may also apply some floors to other field values
  void solve_chemistry(Block* block, double compute_time,
                       double dt) const noexcept;

  /// wrapper around the various methods for computing various grackle
  /// properties.
  ///
  /// This wraps calls to local_calculate_* family of Grackle functions.
  ///
  /// @param[in]  fadaptor Field data to be used in the calculation.
  /// @param[in]  func_choice Specifies the choice of quantity that should be
  ///     computed.
  /// @param[out] values Output buffer that is used to hold the results.
  /// @param[in]  stale_depth Specifies stale depth of the field
  /// @param[in]  grackle_fields Pointer to a precomputed grackle_field_data
  ///     struct. When this is a nullptr, a temporary is constructed.
  void compute_property(const EnzoFieldAdaptor& fadaptor,
                        GracklePropertyEnum func_choice,
                        enzo_float* values, int stale_depth = 0,
                        grackle_field_data* grackle_fields = nullptr
                        ) const noexcept;

private: // attributes
  // the attributes are treated as immutable after construction/deserialization

  /// wrapper around chemistry_data struct, which stores runtime parameters
  GrackleChemistryData my_chemistry_;

  // the following 2 attributes are primarily pointers in order to avoid
  // including ifdef statements in this header file.

  /// pointer to code_units struct
  std::unique_ptr<code_units> grackle_units_;
  /// pointer to chemistry_data_storage
  std::unique_ptr<chemistry_data_storage> grackle_rates_;

  /// the following parameter is only used in non-cosmological simulations. It
  /// specifies the redshift of the UV background. A negative value means that
  /// this is unset
  double radiation_redshift_;

};

#endif /* ENZO_ENZO_GRACKLE_FACADE_HPP */
