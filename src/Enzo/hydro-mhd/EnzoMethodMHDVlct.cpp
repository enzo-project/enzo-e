// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodMHDVlct.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri June 14 2019
/// @brief    [\ref Enzo] Implementation of the EnzoMethodMHDVlct class

#include <algorithm>    // std::copy

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp" // EnzoBlock,
#include "Enzo/hydro-mhd/hydro-mhd.hpp"
#include "Enzo/hydro-mhd/EnzoMHDIntegratorStageCommands.hpp"

//----------------------------------------------------------------------

static void check_field_l_(std::vector<std::string> &field_l)
{
  FieldDescr * field_descr = cello::field_descr();
  for (const std::string& field : field_l){
    ASSERT1("EnzoMethodMHDVlct", "\"%s\" must be a permanent field",
	    field.c_str(), field_descr->is_field(field));
  }
}

//----------------------------------------------------------------------

// concatenate 2 vectors of strings
static str_vec_t concat_str_vec_(const str_vec_t& vec1, const str_vec_t& vec2){
  str_vec_t out(vec2);
  out.reserve(vec1.size() + vec2.size());
  out.insert(out.begin(),vec1.begin(), vec1.end());
  return out;
}

//----------------------------------------------------------------------

static std::pair<std::string, EnzoMHDIntegratorStageArgPack*>
parse_mhdchoice_pack_pair_(ParameterGroup p)
{
  // is_defined is a function that indicates if a user defined a parameter
  auto is_defined = [&](const std::string& param_name) -> bool
  { return p.param(param_name) != nullptr; };

  // parse time-scheme & reconstruct-method (they affect backwards compat.)
  const std::string time_scheme = p.value_string("time_scheme", "vl");
  std::vector<std::string> recon_names;
  if (time_scheme == "euler") {
    recon_names = {p.value_string("reconstruct_method", "plm")};
    WARNING("parse_mhdchoice_pack_pair_",
            "\"euler\" temporal integration exists for debugging purposes. It "
            "has not been rigorously tested. USE WITH CAUTION");
  } else if (time_scheme == "vl") {
    // ALWAYS use nearest-neighbor reconstruction for partial timestep!
    recon_names = {"nn", p.value_string("reconstruct_method", "plm")};
  } else {
    ERROR1("parse_mhdchoice_pack_pair_",
           "\"%s:time_scheme\" must be \"vl\" or \"euler\"",
           p.get_group_path().c_str());
  }

  if (is_defined("half_dt_reconstruct_method") ||
      is_defined("full_dt_reconstruct_method")) {
    ERROR1("parse_mhdchoice_pack_pair_",
           "In the \"%s\" parameter-group, \"half_dt_reconstruct_method\" & "
           "\"full_dt_reconstruct_method\" have been removed. The former was "
           "only allowed to have a value of \"nn\" and the latter was "
           "replaced with \"reconstruct_method\".",
           p.get_group_path().c_str());
  }

  // allocate and initialize argpack-ptr

  if (!is_defined("mhd_choice")) {
    std::string name = p.full_name("mhd_choice");
    ERROR1("parse_mhdchoice_pack_pair_", "%s wasn't specified", name.c_str());
  }

  EnzoMHDIntegratorStageArgPack *argpack_ptr =
    new EnzoMHDIntegratorStageArgPack {p.value_string("riemann_solver","hlld"),
                                       recon_names,
                                       p.value_float("theta_limiter", 1.5),
                                       p.value_string("mhd_choice", "")};

  return {time_scheme, argpack_ptr};
}

//----------------------------------------------------------------------

EnzoMethodMHDVlct::EnzoMethodMHDVlct (ParameterGroup p,
                                      bool store_fluxes_for_corrections)
  : Method()
{

  std::pair<std::string, EnzoMHDIntegratorStageArgPack*> pair
    = parse_mhdchoice_pack_pair_(p);

  time_scheme_ = pair.first;
  integrator_arg_pack_ = pair.second;
  const double dflt_courant = (time_scheme_ == "vl") ? 0.3 : 1.0;
  this->set_courant(p.value_float("courant",dflt_courant));

  int nstages = static_cast<int>(integrator_arg_pack_->recon_names.size());

  // initialize the integrator
  integrator_ = new EnzoMHDIntegratorStageCommands(*integrator_arg_pack_);

  // Determine the lists of fields that are required to hold the integration
  // quantities and primitives and ensure that they are defined
  integration_field_list_ = integrator_->integration_quantity_keys();
  primitive_field_list_ = integrator_->primitive_quantity_keys();

  check_field_l_(integration_field_list_);
  check_field_l_(primitive_field_list_);

  // make sure "pressure" is defined (it's needed to compute the timestep)
  FieldDescr * field_descr = cello::field_descr();
  ASSERT("EnzoMethodMHDVlct", "\"pressure\" must be a permanent field",
	 field_descr->is_field("pressure"));

  bfield_method_ = integrator_->construct_bfield_method(nstages);
  if (bfield_method_ != nullptr) { bfield_method_->check_required_fields(); }

  // Check that the (cell-centered) ghost depth is large enough
  // Face-centered ghost depth could in principle be 1 smaller
  int min_gdepth_req = 0;
  for (int stage_index = 0; stage_index < nstages; stage_index++) {
    min_gdepth_req += integrator_->staling_from_stage(stage_index);
  }
  int gx,gy,gz;
  field_descr->ghost_depth(field_descr->field_id("density"), &gx, &gy, &gz);
  ASSERT1("EnzoMethodMHDVlct::compute", "ghost depth must be at least %d.",
	  min_gdepth_req, std::min(gx, std::min(gy, gz)) >= min_gdepth_req);

  // initialize other attributes
  store_fluxes_for_corrections_ = store_fluxes_for_corrections;
  if (store_fluxes_for_corrections){
    ASSERT("EnzoMethodMHDVlct::EnzoMethodMHDVlct",
           "Flux corrections are currently only supported in hydro-mode",
           integrator_->is_pure_hydro());
  }

  scratch_space_ = nullptr;

  // Finally, initialize the default Refresh object
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  // Need to refresh all fields because the fields holding passively advected
  // scalars won't necessarily be known until after all Methods have been
  // constructed and all intializers have been executed
  refresh->add_all_fields();
}

//----------------------------------------------------------------------

EnzoMethodMHDVlct::~EnzoMethodMHDVlct()
{
  delete integrator_arg_pack_;
  delete integrator_;
  if (scratch_space_ != nullptr){
    delete scratch_space_;
  }
  if (bfield_method_ != nullptr){
    delete bfield_method_;
  }
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | time_scheme_;

  if (p.isUnpacking()) {
    integrator_arg_pack_ = new EnzoMHDIntegratorStageArgPack;
  }
  p | *integrator_arg_pack_;
  if (p.isUnpacking()) {
    integrator_ = new EnzoMHDIntegratorStageCommands(*integrator_arg_pack_);
    bfield_method_ = integrator_->construct_bfield_method(2);
  }


  // skip scratch_space_. This will be freshly constructed the first time that
  // the compute method is called.

  p|integration_field_list_;
  p|primitive_field_list_;
  p|lazy_passive_list_;
  p|store_fluxes_for_corrections_;
}

//----------------------------------------------------------------------

EnzoEFltArrayMap EnzoMethodMHDVlct::get_integration_map_
(Block * block,  const str_vec_t *passive_list) const noexcept
{
  str_vec_t field_list = (passive_list == nullptr) ? integration_field_list_ :
    concat_str_vec_(integration_field_list_, *passive_list);

  Field field = block->data()->field();
  std::vector<EFlt3DArray> arrays;
  arrays.reserve(field_list.size());
  for (const std::string& field_name : field_list){
    arrays.push_back( field.view<enzo_float>(field_name) );
  }

  return EnzoEFltArrayMap("integration",field_list,arrays);
}

//----------------------------------------------------------------------

static EnzoEFltArrayMap get_accel_map_(Block* block) noexcept
{
  Field field = block->data()->field();
  if (field.field_id("acceleration_x") < 0){
    return EnzoEFltArrayMap();
  }

  str_vec_t field_list = {"acceleration_x", "acceleration_y", "acceleration_z"};
  std::vector<CelloView<enzo_float,3>> arrays
    = {field.view<enzo_float>("acceleration_x"),
       field.view<enzo_float>("acceleration_y"),
       field.view<enzo_float>("acceleration_z")};
  return EnzoEFltArrayMap("accel", field_list, arrays);
}

//----------------------------------------------------------------------

EnzoVlctScratchSpace* EnzoMethodMHDVlct::get_scratch_ptr_
(const std::array<int,3>& field_shape, const str_vec_t& passive_list) noexcept
{
  if (scratch_space_ == nullptr){
    scratch_space_ = new EnzoVlctScratchSpace
      (field_shape, integration_field_list_, primitive_field_list_,
       integrator_->dUcons_map_keys(), passive_list,
       enzo::fluid_props()->dual_energy_config().any_enabled());
  }
  return scratch_space_;
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::save_fluxes_for_corrections_
(Block * block, const EnzoEFltArrayMap &flux_map, int dim, double cell_width,
 double dt) const noexcept
{

  Field field = block->data()->field();

  // load the cell-centered shape and the ghost depth
  int density_field_id = field.field_id("density");
  int cc_mx, cc_my, cc_mz; // the values are ordered as x,y,z
  field.dimensions(density_field_id, &cc_mx, &cc_my, &cc_mz);
  int gx, gy, gz;
  field.ghost_depth(density_field_id, &gx, &gy, &gz);

  double dt_dxi = dt/cell_width;

  FluxData * flux_data = block->data()->flux_data();
  const int nf = flux_data->num_fields();

  for (int i_f = 0; i_f < nf; i_f++) {
    const int index_field = flux_data->index_field(i_f);
    const std::string field_name = field.field_name(index_field);

    // note the field_name is the same as the key
    CelloView<const enzo_float, 3> flux_arr = flux_map.at(field_name);

    FaceFluxes * left_ff = flux_data->block_fluxes(dim,0,i_f);
    FaceFluxes * right_ff = flux_data->block_fluxes(dim,1,i_f);

    int mx, my, mz;
    left_ff->get_size(&mx,&my,&mz);

    int dx_l,dy_l,dz_l,    dx_r,dy_r,dz_r;
    enzo_float* left_dest = left_ff->flux_array(&dx_l,&dy_l,&dz_l);
    enzo_float* right_dest = right_ff->flux_array(&dx_r,&dy_r,&dz_r);

    // NOTE: dx_l/dx_r, dy_l/dy_r, dz_l/dz_r have the wrong value when
    // dim is 0, 1, or 2 respectively. In each case the variables are equal to
    // 1 when they should be 0.

    if (dim == 0){
      int left_ix = gx-1;
      int right_ix = cc_mx-gx-1;

      for (int iz = 0; iz < mz; iz++){
        for (int iy = 0; iy < my; iy++){
          left_dest[dz_l*iz + dy_l*iy]
            = dt_dxi* flux_arr(gz+iz, gy+iy, left_ix);
          right_dest[dz_r*iz + dy_r*iy]
            = dt_dxi* flux_arr(gz+iz, gy+iy, right_ix);
        }
      }

    } else if (dim == 1) {
      int left_iy = gy-1;
      int right_iy = cc_my-gy-1;

      for (int iz = 0; iz < mz; iz++){
        for (int ix = 0; ix < mx; ix++){
          left_dest[dz_l*iz + dx_l*ix]
            = dt_dxi* flux_arr(gz+iz, left_iy, gx+ix);
          right_dest[dz_l*iz + dx_l*ix]
            = dt_dxi* flux_arr(gz+iz, right_iy, gx+ix);
        }
      }
    } else {
      int left_iz = gz-1;
      int right_iz = cc_mz-gz-1;

      for (int iy = 0; iy < my; iy++){
        for (int ix = 0; ix < mx; ix++){
          left_dest[dy_l*iy + dx_l*ix]
            = dt_dxi* flux_arr(left_iz, gy + iy, gx+ix);
          right_dest[dy_l*iy + dx_l*ix]
            = dt_dxi* flux_arr(right_iz, gy + iy, gx+ix);
        }
      }

    }
  }
}

//----------------------------------------------------------------------

static void allocate_FC_flux_buffer_(Block * block) throw()
{
  Field field = block->data()->field();
  // this could be better integrated with fields required by the solver
  auto field_names = field.groups()->group_list("conserved");
  const int nf = field_names.size();
  std::vector<int> field_list;
  field_list.resize(nf);
  for (int i=0; i<nf; i++) {
    field_list[i] = field.field_id(field_names[i]);
  }

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  // this needs to be allocated every cycle
  block->data()->flux_data()->allocate (nx,ny,nz,field_list,
                                        true /* = single_flux_array */ );
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::compute ( Block * block) throw()
{
  if (cello::is_initial_cycle(InitCycleKind::fresh_or_noncharm_restart)) {
    post_init_checks_();
  }

  if (store_fluxes_for_corrections_){ allocate_FC_flux_buffer_(block); }

  if (block->is_leaf()) {
    // load the list of keys for the passively advected scalars
    const str_vec_t passive_list = *(lazy_passive_list_.get_list());

    // initialize map that holds arrays wrapping the Cello Fields holding each
    // of the integration quantities. Additionally, this also includes
    // temporary arrays used to hold the specific form of the passive scalar
    //
    // by the very end of EnzoMethodMHDVlct::compute, the arrays in this map
    // will be updated with their new values
    EnzoEFltArrayMap external_integration_map = get_integration_map_
      (block, &passive_list);

    // get maps of arrays and stand-alone arrays that serve as scratch space.
    // (first, retrieve the pointer to the scratch space struct)
    const std::array<int,3> shape = {external_integration_map.array_shape(0),
                                     external_integration_map.array_shape(1),
                                     external_integration_map.array_shape(2)};
    EnzoVlctScratchSpace* const scratch = get_scratch_ptr_(shape, passive_list);

    // map used for storing integration values at the half time-step. This
    // includes key,array pairs for each entry in external_integration_map
    // (there should be no aliased arrays shared between maps)
    EnzoEFltArrayMap temp_integration_map = scratch->temp_integration_map;

    // Map of arrays used to temporarily store the cell-centered primitive
    // quantities that are subsequently reconstructed. This includes arrays for
    // storing the specific form of each of the passively advected scalars.
    EnzoEFltArrayMap primitive_map = scratch->primitive_map;

    // holds left and right reconstructed primitives (scratch-space)
    EnzoEFltArrayMap priml_map = scratch->priml_map;
    EnzoEFltArrayMap primr_map = scratch->primr_map;

    // maps used to store fluxes
    std::array<EnzoEFltArrayMap, 3> flux_maps_xyz = {scratch->xflux_map,
                                                     scratch->yflux_map,
                                                     scratch->zflux_map};

    // map of arrays used to accumulate the changes to the conserved forms of
    // the integration quantities and passively advected scalars. In other
    // words, at the start of the (partial) timestep, the fields are all set to
    // zero and are used to accumulate the flux divergence and source terms. If
    // CT is used, it won't have space to store changes in the magnetic fields.
    EnzoEFltArrayMap dUcons_map = scratch->dUcons_map;

    // initialize the map that wraps the fields holding the acceleration
    // components (these are nominally computed from gravity). This data is
    // used for the gravity source term calculation. An empty map indicates
    // that the gravity source term is not included.
    const EnzoEFltArrayMap accel_map = get_accel_map_(block);

    // array used to temporary velocity data for computing the internal
    // energy source terms (while using the dual-energy formalism)
    EFlt3DArray interface_vel_arr = scratch->interface_vel_arr;

    // allocate constrained transport object
    if (bfield_method_ != nullptr) {
      bfield_method_->register_target_block(block);
    }

    const std::array<enzo_float,3> cell_widths_xyz =
      { enzo::block(block)->CellWidth[0],
        enzo::block(block)->CellWidth[1],
        enzo::block(block)->CellWidth[2],
      };

    const double dt = block->state()->dt();

    // stale_depth indicates the number of field entries from the outermost
    // field value that the region including "stale" values (need to be
    // refreshed) extends over.
    int stale_depth = 0;

    // NOTE: it would take minimal effort to extend this code to support
    // Runge-Kutta integration (e.g. 2nd order or 3rd order). We just need to
    // adopt the 2 register method discussed in the Athena++ method paper.
    // That change boils down to 4 parts (for unigrid):
    //   1. we probably need to introduce a memory copy method to copy data
    //      between the integration map (OR we may be able to use Cello's field
    //      history capacity)
    //   2. we need to tweak the signature of
    //      EnzoMHDIntegratorStageCommands::compute_update_stage
    //   3. to EnzoIntegrationQuanUpdate needs a few minor tweaks
    //   4. Adding MHD support requires some extra steps
    unsigned short nstages;
    if (time_scheme_ == "euler") {
      nstages = 1;
    } else if (time_scheme_ == "vl") {
      nstages = 2;
    } else {
      ERROR("EnzoMethodMHDVlct::compute", "unrecognized time_scheme_");
    }

    // repeat the following loop twice (for half time-step and full time-step)
    for (unsigned short stage_index = 0; stage_index < nstages; stage_index++){

      const bool is_final_stage = (stage_index + 1) == nstages;
      const double cur_dt = (!is_final_stage) ? dt/2. : dt;
      EnzoEFltArrayMap cur_stage_integration_map =
        (stage_index == 0) ? external_integration_map : temp_integration_map;
      EnzoEFltArrayMap out_integration_map =
        (is_final_stage) ? external_integration_map : temp_integration_map;

      integrator_->compute_update_stage
        (external_integration_map,  // holds values from start of the timestep
         cur_stage_integration_map, // holds values from start of current stage
         out_integration_map,       // where to write results of current stage
         primitive_map, priml_map, primr_map,
         flux_maps_xyz, dUcons_map, accel_map, interface_vel_arr,
         passive_list, this->bfield_method_, stage_index,
         cur_dt, stale_depth, cell_widths_xyz);

      // update the stale_depth from the current stage
      stale_depth += integrator_->staling_from_stage((int)stage_index);

      if (is_final_stage && store_fluxes_for_corrections_) {
        // Dual Energy Formalism Note:
        // - the interface velocities on the edge of the blocks will be
        //   different if using SMR/AMR. This means that the internal energy
        //   source terms won't be fully self-consistent along the edges. This
        //   same effect is also present in the Ppm Solver
        for (int i = 0; i < 3; i++) {
          save_fluxes_for_corrections_(block, flux_maps_xyz[i], i,
                                       cell_widths_xyz[i], cur_dt);
        }
      }

      if (bfield_method_ != nullptr) {
        bfield_method_->increment_partial_timestep();
      }

    }
  }

  block->compute_done();
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::post_init_checks_() const noexcept
{
  if (enzo::grackle_chemistry() != nullptr){
    // we can remove this check if:
    // - EnzoMethodGrackle is modified to no longer require internal_energy to
    //   be a permanent field
    // - OR if EnzoMethodMHDVlct is modified to always track internal_energy
    //   (even if not using dual-energy formalism)
    ASSERT("EnzoMethodMHDVlct::determine_quantities_",
           ("Grackle cannot currently be used alongside this integrator "
            "unless the dual-energy formalism is in use"),
           enzo::fluid_props()->dual_energy_config().any_enabled());
  }

  ASSERT("EnzoMethodMHDVlct::post_init_checks_",
         "This solver isn't currently compatible with cosmological sims",
         enzo::cosmology() == nullptr);

  const Problem* problem = enzo::problem();

  // problems would arise relating to particle-mesh deposition (relating to
  // particle drift before deposition) and in cosmological simulations if the
  // the VL+CT method were to precede the gravity method.
  ASSERT("EnzoMethodMHDVlct::post_init_checks_",
         "when the gravity method exists, it must precede this method.",
         problem->method_precedes("gravity", "mhd_vlct") |
         (!problem->method_exists("gravity")) );

  // the following checks address some problems I've encountered in the past
  // (they probably need to be revisited when we add Bfield flux corrections)
  bool fc_exists = problem->method_exists("flux_correct");
  if (fc_exists & !problem->method_precedes("mhd_vlct", "flux_correct")) {
    ERROR("EnzoMethodMHDVlct::post_init_checks_",
          "this method can't precede the flux_correct method");
  } else if (fc_exists & !store_fluxes_for_corrections_) {
    ERROR("EnzoMethodMHDVlct::post_init_checks_",
          "the flux_correct method exists, but this method isn't saving "
          "fluxes to be used for corrections.");
  } else if (store_fluxes_for_corrections_ & !fc_exists) {
    ERROR("EnzoMethodMHDVlct::post_init_checks_",
          "this method is saving fluxes for corrections, but the flux_correct "
          "method doesn't exist");
  }
}

//----------------------------------------------------------------------

double EnzoMethodMHDVlct::timestep ( Block * block ) throw()
{
  // analogous to ppm timestep calulation, probably want to require that cfast
  // is no smaller than some tiny positive number.

  // Constructs a map containing the field data for each integration quantity
  // This includes each passively advected scalar (as densities)
  EnzoEFltArrayMap integration_map = get_integration_map_
    (block, (lazy_passive_list_.get_list()).get());

  EnzoPhysicsFluidProps* fluid_props = enzo::fluid_props();

  if (fluid_props->dual_energy_config().any_enabled()){
    // synchronize eint and etot.
    // This is only strictly necessary after problem initialization and when
    // there is an inflow boundary condition
    fluid_props->apply_floor_to_energy_and_sync(integration_map, 0);
  }

  // Compute thermal pressure
  Field field = block->data()->field();
  CelloView<enzo_float, 3> pressure = field.view<enzo_float>("pressure");
  fluid_props->pressure_from_integration(integration_map, pressure, 0);

  // widths of cells
  EnzoBlock * enzo_block = enzo::block(block);
  const double dx = enzo_block->CellWidth[0];
  const double dy = enzo_block->CellWidth[1];
  const double dz = enzo_block->CellWidth[2];

  // compute the actual timestep
  double dtBaryons = integrator_->timestep(integration_map, pressure,
                                           dx, dy, dz);

  // Multiply resulting dt by CourantSafetyNumber (for extra safety!).
  // This should be less than 0.5 for standard algorithm
  return dtBaryons * courant_;
}
