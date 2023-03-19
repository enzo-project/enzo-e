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

// this probably needs to be extended in the future to include
// gamma, mmw, dust_temperature
enum class GracklePropertyEnum{ cooling_time, pressure, temperature };

class GrackleFacade {

  /// @class    GrackleFacade
  /// @ingroup  Enzo
  ///
  /// An important invariant is that an instance of this object always
  /// represents a fully initialized object (it will not hold nullptrs or be
  /// partially initialized). That simplifies reasoning about the object
  /// significantly!
  ///
  /// @note
  /// Currently, this class requires that we re-initialize the code_units and
  /// grackle_data_storage structs at every simulation time. From looking at
  /// the internals of Grackle, and the way that enzo-classic uses Grackle,
  /// this isn't necessary - we should just need to initialize once. With that
  /// said, while initially introducing this class into Enzo-E, we are choosing
  /// to maintaining this existing behavior so that we are only modifying one
  /// aspect of the Grackle-wrapping code at a time. Once we are satisfied that
  /// this class works, we will look at removing this behavior (this will also
  /// let us stop tracking the simulation time at which these structs are
  /// initialized).
  ///
  /// @note
  /// It turns out that the invariant that instances are ALWAYS fully
  /// initialized, complicates serialization with the pup framework. The pup
  /// framework implicitly assumes that the object gets constructed in
  /// some default state, and then uses a pup method to update the fields of
  /// the existing object. This implies that GrackleFacade should have a
  /// default constructor (which represents a null state).
  ///
  /// @par
  /// Alternatively, we could have GrackleFacade use the PUP::able machinery.
  /// In that case, the GrackleFacade(CkMigrateMessage *m) constructor would be
  /// the only way to create an unitialized GrackleFacade instance. When that
  /// constructor is invoked, it should always be followed by a pup call.
  /// In either case, we plan to revisit this after we remove the need to
  /// reinitialize the data structures during every simulation time.

  static bool linked_against_grackle() noexcept;

  /// either pack up a GrackleFacade or unpack and allocate an instance
  ///
  /// @note
  /// This interface is a little funky to account for the class's invariant.
  /// In the overview documentation, we discuss when we may revist this choice.
  friend void pup_GrackleFacade(PUP::er &p,
                                std::unique_ptr<GrackleFacade> &ptr) noexcept;

public: // interface

  /// public-facing constructor
  GrackleFacade(GrackleChemistryData&& my_chemistry,
                const double radiation_redshift,
                const double physics_cosmology_initial_redshift,
                const double time) noexcept;

  /// Destructor
  ~GrackleFacade() noexcept;

  // delete the default constructor
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
  const GrackleChemistryData* get_chemistry() const noexcept
  {
    bool not_null = (GrackleFacade::linked_against_grackle() &&
                     (my_chemistry_.get<int>("use_grackle") == 1));
    return (not_null) ? &my_chemistry_ : nullptr;
  }

public: // low-level legacy methods - these will (probably) be removed or made
        // private in the near future

  /// sets up grackle units.
  ///
  /// @param[in]  current_time The current time. A negative value can be passed
  ///     if the current value is not known or convenient to get. In that case,
  ///     the program will abort if the current time is needed (i.e. because
  ///     this is a cosmological simulation)
  /// @param[out] grackle_units The object pointed to by this pointer is set up
  void setup_grackle_units(double current_time,
                           code_units * grackle_units) const noexcept;

  void setup_grackle_units(const EnzoFieldAdaptor& fadaptor,
                           code_units * grackle_units) const noexcept;

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

  /// this method reinitializes the object if current_time differs from the
  /// previous time when this object was initialized
  ///
  /// @note
  /// Inspection of the grackle source code strongly suggests that this method
  /// is unnecessary (we should only need to configure once). In the immediate
  /// short-term we are maintaining the current behavior (while testing this
  /// class), and then we will test removal of this functionallity in the near
  /// future.
  bool ensure_configured_for_curtime(double current_time) noexcept{
    if (time_grackle_data_initialized_ != current_time){
      init_grackle_rates_(current_time, false);
      return true;
    }
    return false;
  }

public: // wrapped grackle functions:

  /// light-weight wrapper around local_solve_chemistry function from grackle
  void solve_chemistry(Block* block, double dt) const noexcept;

  void calculate_cooling_time(const EnzoFieldAdaptor& fadaptor, enzo_float* ct,
                              int stale_depth = 0,
			      code_units* grackle_units = nullptr,
			      grackle_field_data* grackle_fields = nullptr)
    const noexcept
  {
    compute_local_property_(fadaptor, ct, stale_depth, grackle_units,
                            grackle_fields, GracklePropertyEnum::cooling_time);
  }

  void calculate_pressure(const EnzoFieldAdaptor& fadaptor,
                          enzo_float* pressure, int stale_depth = 0,
			  code_units* grackle_units = nullptr,
			  grackle_field_data* grackle_fields = nullptr)
    const noexcept
  {
    compute_local_property_(fadaptor, pressure, stale_depth, grackle_units,
                            grackle_fields, GracklePropertyEnum::pressure);
  }

  void calculate_temperature(const EnzoFieldAdaptor& fadaptor,
                             enzo_float* temperature, int stale_depth = 0,
			     code_units* grackle_units = nullptr,
			     grackle_field_data* grackle_fields = nullptr)
    const noexcept
  {
    compute_local_property_(fadaptor, temperature, stale_depth, grackle_units,
                            grackle_fields, GracklePropertyEnum::temperature);
  }



private: // helper methods


  // when grackle_units is nullptr, new values are temporarily allocated
  void compute_local_property_(const EnzoFieldAdaptor& fadaptor,
                               enzo_float* values, int stale_depth,
			       code_units* grackle_units,
			       grackle_field_data* grackle_fields,
			       GracklePropertyEnum func_choice) const noexcept;


  // this is private because its only used internally by the public constructor
  // and the pup routine during deserialization
  GrackleFacade(GrackleChemistryData&& my_chemistry,
                std::unique_ptr<code_units>&& grackle_units_,
                const double radiation_redshift,
                const double units_init_time);

  /// this (re)initializes the grackle_rates_ pointer.
  ///
  /// the following difference arise based on ``first_initialization``'s value:
  ///   - when ``true``, ``grackle_units_`` is assumed to be completely
  ///     initialized (i.e. allocated and initialized for the specified value
  ///     of ``current_time``), but ``grackle_rates_`` holds ``nullptr``.
  ///   - when ``false``, both ``grackle_units_`` and ``grackle_rates_`` are
  ///     assumed not to be ``nullptr``, but may need to be initialized based
  ///     on the value of ``current_time``.
  ///
  /// In both cases, this function modifies/(re)allocates values held by fields
  /// of ``grackle_rates_``. It also ensures that
  /// ``time_grackle_data_initialized_`` matches ``current_time``.
  void init_grackle_rates_(double current_time,
                           bool first_initialization = false) noexcept;

  void deallocate_grackle_rates_() noexcept;

private: // attributes
  /// wrapper around chemistry_data struct, which stores runtime parameters
  GrackleChemistryData my_chemistry_;

  /// pointer to code_units struct (this is only a pointer to avoid ifdef
  /// statements in this header file, it should NEVER actually hold a nullptr)
  std::unique_ptr<code_units> grackle_units_;

  /// pointer to chemistry_data_storage (this is only a pointer to avoid ifdef
  /// statements in this header file, it should NEVER actually hold a nullptr)
  std::unique_ptr<chemistry_data_storage> grackle_rates_;

  /// stores the last simulation time when grackle_rates_ was initialized.
  double time_grackle_data_initialized_;

  /// the following parameter is only used in non-cosmological simulations. It
  /// specifies the redshift of the UV background. A negative
  double radiation_redshift_;

};

#endif /* ENZO_ENZO_GRACKLE_FACADE_HPP */
