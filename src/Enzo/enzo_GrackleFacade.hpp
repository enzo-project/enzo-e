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

class GrackleFacade : public PUP::able {

  /// @class    GrackleFacade
  /// @ingroup  Enzo
  ///
  /// The key invariant of this class is that an instance is ALWAYS properly
  /// initialized (i.e. all grackle structs are properly configured).
  ///
  /// There is one minor caveat to this rule: the Charm++ migration constructor
  /// leaves it in a partially initialized state. But that's ok because the
  /// migration constructor is ALWAYS followed by unpacking.
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

private: // attributes
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
