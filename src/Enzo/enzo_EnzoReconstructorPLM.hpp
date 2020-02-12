// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoReconstructorPLM.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed May 1 2019
/// @brief    [\ref Enzo] Implementation of Enzo's Piecewise Linear
///           Reconstruction

#ifndef ENZO_ENZO_RECONSTRUCTOR_PLM_HPP
#define ENZO_ENZO_RECONSTRUCTOR_PLM_HPP

#include <type_traits>

/// @typedef limiter_function_call_signature
/// @brief   This is the expected function call signature of the function call
///          operator, `operator()`, expected for a `Limiter` functor that is
///          passed to `EnzoReconstructorPLM`.
///
/// We expect the `operator()` method of a `Limiter` functor to be declared as:
///
/// @code
///     enzo_float operator() (enzo_float vm1, enzo_float v, enzo_float vp1,
///                            enzo_float theta_limiter);
/// @endcode
///
/// The function should compute the limited slope to use in the piecewise
/// linear reconstruction of the left and right interface values (which are
/// ultimately used by the Riemann Solver) for some cell-centered variable `w`.
/// The function is passed 3 contiguous values of `w`: `w_{i-1}`, `w_i`, and
/// `w_{i+1} as the arguments `vm1`, `v` and `vp1`. The `theta_limiter`
/// argument can be used to tune the limiter (or it can be ignored). The
/// resulting limited slope, `dw` is used to reconstruct the interface
/// variables in the equations:
///     `w_{L,i+1/2} = w_i + 0.5*dw` and `w_{R,i+1/2} = w_i - 0.5*dw`

typedef enzo_float (*limiter_function_call_signature)
  (enzo_float, enzo_float, enzo_float, enzo_float);

//----------------------------------------------------------------------

/// @def      DEFINE_HAS_INSTANCE_METHOD
/// @brief    Macro for defining a struct that is used to check if a given
///           class definition has a public instance method with a specified
///           name and function signature
///
/// @param traitsName The name of the struct that will be defined
/// @param func_name  Name of the public instance method to check for
/// @param signature  The signature that the method should have (this can
///     be a function pointer).
///
/// This macro is modified slightly modified from the answer provided at
/// https://stackoverflow.com/a/16824239 . The particular answer that this is
/// adapted from allows the required instance method to be inherited. It is
/// used to help check that the the template argument of EnzoReconstructorPLM
/// in order to generate more helpful error messages
#define DEFINE_HAS_INSTANCE_METHOD(traitsName, func_name, signature)        \
  template<typename T>                                                      \
  class traitsName							    \
  {                                                                         \
    private:                                                                \
      /* primary declaration of inner class template that does the check */ \
      template<typename, typename U>                                        \
      struct inner_helper_;                                                 \
                                                                            \
      /* explicit specialization that actually performs the checking */     \
      template<typename Ret, typename... Args>                              \
      struct inner_helper_<T, Ret(*)(Args...)>                              \
      {                                                                     \
	/* attempts to call the method and compares the return type */      \
	/* against the expected return type. */                             \
	template<typename U>                                                \
	static constexpr auto check(U*) -> typename std::is_same<           \
	    decltype(std::declval<U>().func_name(std::declval<Args>()...)), \
	    Ret >::type;                                                    \
                                                                            \
	template<typename>                                                  \
	static constexpr std::false_type check(...);                        \
                                                                            \
	typedef decltype(check<T>(0)) type;                                 \
                                                                            \
      public:                                                               \
	static constexpr bool value = type::value;                          \
      };                                                                    \
    public:                                                                 \
      static constexpr bool value = inner_helper_<T, signature>::value;     \
  };

//----------------------------------------------------------------------

// define a struct to use at compile time to help make sure that the limiter
// functor has a sane value.
DEFINE_HAS_INSTANCE_METHOD(has_limiter_function_call_operator, operator(),
			   limiter_function_call_signature);

//----------------------------------------------------------------------

template <class Limiter>
class EnzoReconstructorPLM : public EnzoReconstructor
{
  /// @class    EnzoReconstructorPLM
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates piecwise linear reconstruction of
  ///           primitives at interface
  ///
  /// @tparam Limiter The functor class used to implement different types of
  ///     slope limiters. It must be default constructible. See the
  ///     documentation of `limiter_function_call_signature` for details about
  ///     the expected signature of the `operator()` method, and the
  ///     expectations about the returned value.
  ///
  /// Different specializations of EnzoReconstructorPLM implement different
  /// types of slope limiters.
  ///
  /// To implement a new type of slope limiter:
  ///   1. Define a new `Limiter` functor (e.g. `PLM_EnzoRKLimiter`).
  ///   2. It might be useful to define an alias name for the specialization of
  ///      `EnzoReconstructorPLM` that uses the new `ImplStruct`
  ///      (e.g. `using EnzoReconstructorPLMEnzoRKLim)
  ///   3. Then add the particlular specialization to enzo.CI (e.g. add the
  ///      line: `PUPable EnzoReconstructorPLM<PLM_EnzoRKLimiter>;`)
  ///   4. Update `EnzoReconstructor::construct_reconstructor` to construct the
  ///      reconstructor with the appropriate slope limiter when the correct
  ///      name is specified.
  ///   5. Update the documentation with the name of the newly available
  ///      slope limiter
  ///
  /// @note The use of an enum specifying the desired limiter coupled with a
  ///       switch statement was considered to try to swap between limiters. To
  ///       simplify the implementation, compile-time macros could optionally
  ///       be used instead of templates.
  static_assert(std::is_default_constructible<Limiter>::value,
		"The Limiter functor is not default constructable");
  static_assert(has_limiter_function_call_operator<Limiter>::value,
		"The Limiter's operator() method doesn't have the correct "
		"function signature");

public: // interface

  /// Create a new EnzoReconstructorPLM
  EnzoReconstructorPLM(std::vector<std::string> group_names,
		       enzo_float theta_limiter)
    : EnzoReconstructor(group_names),
      theta_limiter_(theta_limiter)
  { }

  /// CHARM++ PUP::able declaration
  PUPable_decl_template(EnzoReconstructorPLM<Limiter>);

  /// CHARM++ migration constructor for PUP::able
  EnzoReconstructorPLM (CkMigrateMessage *m)
    : EnzoReconstructor(m),
      theta_limiter_(0.0)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    EnzoReconstructor::pup(p);
    p | theta_limiter_;
  };

  void reconstruct_interface (Block *block, Grouping &prim_group,
			      Grouping &priml_group, Grouping &primr_group,
			      int dim, EnzoEquationOfState *eos,
			      int stale_depth);

  int total_staling_rate()
  { return 2; }

  int immediate_staling_rate()
  { return 1; }

private:
  // parameter used to adjust some slope limiters
  enzo_float theta_limiter_;
};

//----------------------------------------------------------------------

template <class Limiter>
void EnzoReconstructorPLM<Limiter>::reconstruct_interface
  (Block *block, Grouping &prim_group, Grouping &priml_group,
   Grouping &primr_group, int dim, EnzoEquationOfState *eos, int stale_depth)
{
  std::vector<std::string> group_names = this->group_names_;

  EnzoFieldArrayFactory array_factory(block, stale_depth);
  EnzoPermutedCoordinates coord(dim);

  Limiter limiter_func = Limiter();
  const enzo_float theta_limiter = theta_limiter_;

  // unecessary values are computed for inside faces of outermost ghost zone
  for (unsigned int group_ind=0; group_ind<group_names.size(); group_ind++){

    // load group name and number of fields in the group
    std::string group_name = group_names[group_ind];
    int num_fields = prim_group.size(group_name);

    // Handle possibility of having a density floor
    enzo_float prim_floor =0;
    bool use_floor = false;
    if (group_name == "density"){
      prim_floor = eos->get_density_floor();
      use_floor=true;
    } else if (group_name == "pressure"){
      prim_floor = eos->get_pressure_floor();
      use_floor=true;
    }

    // iterate over the fields in the group
    for (int field_ind=0; field_ind<num_fields; field_ind++){
      // Cast the problem as reconstructing values at:
      //   wl(k, j, i+3/2) and wr(k,j,i+1/2)

      // Prepare cell-centered arrays
      // define:   wc_left(k,j,i)   -> w(k,j,i)
      //           wc_center(k,j,i) -> w(k,j,i+1)
      //           wc_right(k,j,i)  -> w(k,j,i+2)
      EFlt3DArray wc_left, wc_center, wc_right;
      wc_left = array_factory.from_grouping(prim_group, group_name, field_ind);
      wc_center = coord.left_edge_offset(wc_left, 0, 0, 1);
      wc_right  = coord.left_edge_offset(wc_left, 0, 0, 2);

      // Prepare face-centered arrays
      // define:   wl_offset(k,j,i)-> wl(k,j,i+3/2)
      //           wr(k,j,i)       -> wr(k,j,i+1/2)
      EFlt3DArray wr, wl, wl_offset;
      wr = array_factory.reconstructed_field(primr_group, group_name, field_ind,
					     dim);
      wl = array_factory.reconstructed_field(priml_group, group_name, field_ind,
					     dim);
      wl_offset = coord.left_edge_offset(wl, 0, 0, 1);
      
      // At the interfaces between the first and second cell (second-to-
      // last and last cell), along a given axis, no need to worry about
      // initializing the left (right) interface value thanks to the adoption
      // of immediate_staling_rate

      for (int iz=0; iz<wc_right.shape(0); iz++) {
	for (int iy=0; iy<wc_right.shape(1); iy++) {
	  for (int ix=0; ix<wc_right.shape(2); ix++) {

	    // compute limited slopes
	    enzo_float val = wc_center(iz,iy,ix);
	    enzo_float dv = limiter_func(wc_left(iz,iy,ix), val,
					 wc_right(iz,iy,ix), theta_limiter);
	    enzo_float half_dv = dv*0.5;
	    enzo_float left_val, right_val;

	    if (use_floor) {
	      right_val = EnzoEquationOfState::apply_floor(val - half_dv,
							   prim_floor);
	      left_val  = EnzoEquationOfState::apply_floor(val + half_dv,
							   prim_floor);
	    } else {
	      right_val = val - half_dv;
	      left_val  = val + half_dv;
	    }

	    // face centered fields: index i corresponds to the value at i-1/2
	    wr(iz,iy,ix) = right_val;
	    wl_offset(iz,iy,ix) = left_val;
	  }
	}
      }
    }
  }
}

//----------------------------------------------------------------------

// implements sign function https://stackoverflow.com/questions/1903954
inline enzo_float sign(enzo_float val) { return (0.0 < val) - (val < 0.0); }

//----------------------------------------------------------------------

// taken from Enzo's ReconstructionRoutines.h
inline enzo_float Min(enzo_float a, enzo_float b, enzo_float c)
{
  if (a<b) {
    if (c<a)
      return c;
    else 
      return a;
  } else {
    if (c<b)
      return c;
    else 
      return b;
  }
}

//----------------------------------------------------------------------

struct PLM_EnzoRKLimiter
{
  /// @class    PLM_EnzoRKLimiter
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Functor encapsulating the same PLM slope limiter
  ///           used in the Runge-Kutta integrator from Enzo
  ///
  /// This implementation is adapted from the limiter defined in Enzo's
  /// src/enzo/hydro_rk/Rec_PLM.C file. Per the Enzo documentation,
  /// `theta_limiter` is the "flux limiter" in the minmod Van Leer formulation.
  /// The value must be between 1 (least dissipative) and 2 (most dissipative),
  /// and by default it is 1.5.
  ///
  /// To apply the limiter, the left-, right, and centered differences are
  /// computed (and selectively multiplied by `theta_limiter`:
  ///   `dw_{L,i} = (w_i     - w_{i-1}) * theta_limiter`,
  ///   `dw_{R,i} = (w_{i+1} - w_i    ) * theta_limiter`,
  ///   `dw_{C,i} = (w_{i+1} - w_{i-1}) / 2.`
  /// Then the slope is computed by applying the minmod operation on each
  /// difference: `dw = minmod(dw_{L,i},dw_{R,i},dw_{C,i})`. The minmod
  /// operation yields zero if any of it arguments have conflicting signs.
  /// Otherwise the minmod operation returns the value closest to zero.
  ///
  /// We note, `theta_limiter=1` reduces to the traditional minmod limiter:
  ///   `dw = minmod(w_i - w_{i-1}, w_{i+1} - w_i)`
  /// (in this case, `dw_{C,i}` is just the mean of `dw_{R,i}` and `dw_{L,i}`).
  /// Additionally, `theta_limiter=2` is equivalent to the traditional MC
  /// limiter (monotonized central-difference limiter). For more details,
  /// see LeVeque (2002).

  enzo_float operator()(enzo_float vm1, enzo_float v,
			enzo_float vp1, enzo_float theta_limiter)
  {
    enzo_float dv_c = 0.5*(vp1-vm1);
    enzo_float dv_l = (v-vm1) * theta_limiter;
    enzo_float dv_r = (vp1-v) * theta_limiter;

    // do NOT drop the factor of 0.5 below
    // The expression (0.5*(sign(dv_l)+sign(dv_r))) evaluates to
    //   +1   if (dv_l>0 & dv_r>0)
    //   -1   if (dv_l<0 & dv_r<0)
    //    0   if ((dv_l<0 & dv_r>0) || (dv_l<0 & dv_r>0))
    // the Min function evaluates to 0 if (dv_l == 0 || dv_r == 0)
    return (0.5*(sign(dv_l)+sign(dv_r)))*Min(fabs(dv_l),fabs(dv_r),fabs(dv_c));
  }
};

//----------------------------------------------------------------------

struct PLM_AthenaLimiter
{
  /// @class    PLM_AthenaPrimLimiter
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Functor encapsulating the same PLM slope limiter
  ///           used for the reconstruction in Athena and Athena++
  ///
  /// For this limiter, the left- and right-differences,
  /// `dw_{L,i} = (w_i-w_{i-1})` and `dw_{R,i} = (w_{i+1} - w_i)` are computed.
  /// If `SIGN(dw_{L,i}) != SIGN(dw_{R,i})`, then zero is returned. Otherwise,
  /// the returned slope is the harmonic mean of `dw_{L,i}` and `dw_{R,i}`:
  ///   `2*dw_{L,i}*dw_{R,i}/(dw_{L,i}+dw_{R,i})`.
  ///
  /// While this is directly computed for Athena++, in Athena the limited slope
  /// result is more complex. The MC slope limiter (see PLM_EnzoRKLimiter
  /// description) is compared to the above harmonic mean and the value closest
  /// to zero is used. However, one can show that the harmonic mean value is
  /// always less than or equal to the MC limited slopes when the signs of
  /// `dw_{L,i}` and `dw_{R,i}` agree. Note that the limiter described in both
  /// Stone+(2008) and Stone & Gardiner (2009) is closest to the MC slope
  /// limiter (they do not use that name). This may explain the unnecessary,
  /// but explicit comparison between limited slopes in Athena.
  ///
  /// the value of `theta_limiter` is ignored in this case.
  ///
  /// @note This is definitely used for all PLM reconstruction in Athena++.
  ///     While it's clear that in the C version of Athena that this is used
  ///     for reconstruction of primitives, it is less clear whether it is used
  ///     for reconstruction of characteristics.

  enzo_float operator()(enzo_float vm1, enzo_float v,
			enzo_float vp1, enzo_float theta_limiter)
  {
    enzo_float dv_l = (v-vm1);
    enzo_float dv_r = (vp1-v);
    enzo_float temp = dv_l*dv_r;
    if (temp <= 0.) { return 0.; }
    return 2.*temp/(dv_l + dv_r);
  }
};

//----------------------------------------------------------------------

/// @class    EnzoReconstructorPLMEnzoRKLim
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulates PLM reconstruction using the slope
///           limiter from the Runge-Kutta integrator in Enzo
using EnzoReconstructorPLMEnzoRKLim = EnzoReconstructorPLM<PLM_EnzoRKLimiter>;

//----------------------------------------------------------------------

/// @class    EnzoReconstructorPLMAthenaLim
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulates PLM reconstruction using the slope
///           limiter from Athena
using EnzoReconstructorPLMAthenaLim = EnzoReconstructorPLM<PLM_AthenaLimiter>;

#endif /* ENZO_ENZO_RECONSTRUCTOR_PLM_HPP */
