// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodFluxCorrect.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-03-19
/// @brief    Implementation of the flux_correct method

#include "problem.hpp"
#include "charm_simulation.hpp"
#include "test.hpp"

//----------------------------------------------------------------------

MethodFluxCorrect::MethodFluxCorrect
(std::string group, bool enable,
 const std::vector<std::string>& min_digits_fields,
 const std::vector<double>& min_digits_vals) throw() 
  : Method (),
    ir_pre_(-1),
    group_(group),
    enable_(enable),
    min_digits_map_(),
    field_sum_(),
    field_sum_0_(),
    scratch_()
{
  // Set up post-refresh to refresh all conserved fields in group_
  cello::simulation()->refresh_set_name(ir_post_,name());
  Grouping * groups = cello::field_groups();
  const int nf=groups->size(group_);
  for (int i_f=0; i_f<nf; i_f++) {
    cello::refresh(ir_post_)->add_field(groups->item(group_,i_f));
  }

  //  neighbor_flux causes synchronization errors; using neighbor_leaf:
  //  ir_pre_ = add_refresh_(neighbor_flux);
  ir_pre_ = add_refresh_(neighbor_leaf);
  Refresh * refresh_pre = cello::refresh(ir_pre_);
  cello::simulation()->refresh_set_name(ir_pre_,name()+"_fluxes");
  refresh_pre->set_callback(CkIndex_Block::p_method_flux_correct_refresh());
  refresh_pre->add_all_fluxes();
  // Also ensure conserved fields are themselves refreshed
  ASSERT1("MethodFluxCorrect::MethodFluxCorrect",
          "Must add field(s) to the field group \"%s\"",
          group_.c_str(),
          (nf > 0));
  for (int i_f=0; i_f<nf; i_f++) {
    refresh_pre->add_field(groups->item(group_,i_f));
  }

  // Set up min_digits_pairs_
  ASSERT("MethodFluxCorrect::MethodFluxCorrect",
         "min_digits_fields and min_digits_values must have the same length",
         min_digits_fields.size() == min_digits_vals.size());
  for (std::size_t i = 0; i < min_digits_fields.size(); i++){
    const std::string& field = min_digits_fields[i];
    ASSERT2("MethodFluxCorrect::MethodFluxCorrect",
            ("Can't check conserved digits for the \"%s\" field since it's "
             "not in the field group \"%s\"."),
            field.c_str(), group_.c_str(), groups->is_in(field,group_));
    min_digits_map_[field] = min_digits_vals[i];
  }

  // sum mass, momentum, energy

  field_sum_.resize(nf);
  field_sum_0_.resize(nf);
}

//----------------------------------------------------------------------

void MethodFluxCorrect::compute ( Block * block) throw()
{
  cello::refresh(ir_pre_)->set_active(block->is_leaf());

  block->refresh_start
    (ir_pre_, CkIndex_Block::p_method_flux_correct_refresh());

}

//----------------------------------------------------------------------

void Block::Block::p_method_flux_correct_refresh()
{
  static_cast<MethodFluxCorrect*>
    (this->method())->compute_continue_refresh(this);
}

//----------------------------------------------------------------------

void MethodFluxCorrect::compute_continue_refresh( Block * block ) throw()
{
  // accumulate local sums of conserved fields for global sum reduction

  flux_correct_ (block);

  Field field = block->data()->field();
  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth (0,&gx,&gy,&gz);

  cello_float * values;

  FluxData * flux_data = block->data()->flux_data();

  const int nf = flux_data->num_fields();
  long double * reduce = new long double [nf+1];
  std::fill_n(reduce,nf+1,0.0);
  reduce[0] = nf;

  if (block->is_leaf()) {

    cello_float * density = (cello_float *) field.values("density");

    Grouping * groups = cello::field_groups();

    for (int i_f=0; i_f<nf; i_f++) {

      const int index_field = flux_data->index_field(i_f);

      const bool scale_by_density =
        groups->is_in(field.field_name(index_field),"make_field_conservative");

      values = (cello_float *) field.values(index_field);

      if (scale_by_density) {
        for (int iz=gz; iz<mz-gz; iz++) {
          for (int iy=gy; iy<my-gy; iy++) {
            for (int ix=gx; ix<mx-gx; ix++) {
              int i=ix + mx*(iy + my*iz);
              reduce[i_f+1] += values[i]*density[i];
            }
          }
        }
      } else {
        for (int iz=gz; iz<mz-gz; iz++) {
          for (int iy=gy; iy<my-gy; iy++) {
            for (int ix=gx; ix<mx-gx; ix++) {
              int i=ix + mx*(iy + my*iz);
              reduce[i_f+1] += values[i];
            }
          }
        }
      }

      // scale by relative mesh cell volume/area

      const int level = block->level();

      const int w = 1 << level*cello::rank();
      reduce[i_f+1] /= w;
    }
  }

  CkCallback callback (CkIndex_Block::r_method_flux_correct_sum_fields(nullptr), 
                       block->proxy_array());

  block->contribute
    ((nf+1)*sizeof(long double), reduce, sum_long_double_n_type, callback);

  delete [] reduce;
}

//----------------------------------------------------------------------

void Block::Block::r_method_flux_correct_sum_fields(CkReductionMsg * msg)
{
  static_cast<MethodFluxCorrect*>
    (this->method())->compute_continue_sum_fields(this,msg);
}

//----------------------------------------------------------------------

void MethodFluxCorrect::compute_continue_sum_fields
( Block * block, CkReductionMsg * msg) throw()
{
  FluxData * flux_data = block->data()->flux_data();
  const int nf = flux_data->num_fields();
  long double * data = (long double *) msg->getData();
  for (int i_f=0; i_f<nf; i_f++) {
    field_sum_[i_f] = data[i_f+1];
  }
  delete msg;

  Field field = block->data()->field();

  // Write conserved field sums to output (root block only)
  
  if (block->index().is_root()) {

    // for each conserved field
    for (int i_f=0; i_f<nf; i_f++) {

      const int index_field = flux_data->index_field(i_f);

      // save initial sum
      if (block->cycle() == 0) {
        field_sum_0_[i_f] = field_sum_[i_f];
      }
      const int precision = field.precision (index_field);
      const double digits =
        -log10(cello::err_rel(field_sum_0_[i_f],field_sum_[i_f]));
      const std::string& field_name = field.field_name(index_field);
      cello::monitor()->print
        ("Method", "Field %s sum %20.16Le conserved to %g digits of %d",
         field_name.c_str(),
         field_sum_[i_f],
         digits,
         cello::digits_max(precision));

      auto search = min_digits_map_.find(field_name);
      if (search != min_digits_map_.end()){
        const double& min_digits = search->second;
        std::string test_name = "MethodFluxCorrect prec: " + field_name;
        unit_func(test_name.c_str());
        unit_assert (digits >= min_digits);
      }
    }
  }

  block->data()->flux_data()->deallocate();

  block->compute_done();
}

//======================================================================


static void flux_correct_helper_(cello_float * const field_array,
                                 int mx, int my, int mz,
                                 int nx, int ny, int nz,
                                 int i_f, FluxData * flux_data,
                                 const int rank,
                                 const bool (&perform_correction)[3][2])
{
  int ix,iy,iz;

  int axis,face,level_face;

  axis=0;
  // X axis
  for (face=0; face<2; face++) {
    if (perform_correction[axis][face]) {
      int dbx,dby,dbz;
      int dnx,dny,dnz;
      auto block_fluxes    = flux_data->block_fluxes(axis,face,i_f);
      auto neighbor_fluxes = flux_data->neighbor_fluxes(axis,face,i_f);
      auto block_flux_array = block_fluxes->flux_array(&dbx,&dby,&dbz);
      auto neighbor_flux_array = neighbor_fluxes->flux_array(&dnx,&dny,&dnz);

      ix = (face == 0) ? 0 : nx-1;
      for (iz=0; iz<nz; iz++) {
        for (iy=0; iy<ny; iy++) {
          int i=ix+mx*(iy+my*iz);
          int ib = iy*dby+iz*dbz;
          int in = iy*dny+iz*dnz;

          field_array[i] +=
            (2*face-1)*(block_flux_array[ib] - neighbor_flux_array[in]);
        }
      }
    }
  }

  if (rank >= 2) {
    axis=1;
    // Y axis
    for (face=0; face<2; face++) {
      if (perform_correction[axis][face]) {
        int dbx,dby,dbz;
        int dnx,dny,dnz;
        auto block_fluxes    = flux_data->block_fluxes(axis,face,i_f);
        auto neighbor_fluxes = flux_data->neighbor_fluxes(axis,face,i_f);
        auto block_flux_array =
          block_fluxes->flux_array(&dbx,&dby,&dbz);
        auto neighbor_flux_array =
          neighbor_fluxes->flux_array(&dnx,&dny,&dnz);

        iy = (face == 0) ? 0 : ny-1;
        for (iz=0; iz<nz; iz++) {
          for (ix=0; ix<nx; ix++) {
            int i=ix+mx*(iy+my*iz);
            int ib = ix*dbx+iz*dbz;
            int in = ix*dnx+iz*dnz;
            field_array[i] +=
              (2*face-1)*(block_flux_array[ib] - neighbor_flux_array[in]);
          }
        }
      }
    }

  } // rank >= 2

  if (rank >= 3) {
    axis=2;
    for (face=0; face<2; face++) {
      if (perform_correction[axis][face]) {
        int dbx,dby,dbz;
        int dnx,dny,dnz;
        auto block_fluxes    = flux_data->block_fluxes(axis,face,i_f);
        auto neighbor_fluxes = flux_data->neighbor_fluxes(axis,face,i_f);
        auto block_flux_array = block_fluxes->flux_array(&dbx,&dby,&dbz);
        auto neighbor_flux_array = neighbor_fluxes->flux_array(&dnx,&dny,&dnz);

        iz = (face == 0) ? 0 : nz-1;
        for (iy=0; iy<ny; iy++) {
          for (ix=0; ix<nx; ix++) {
            int i=ix+mx*(iy+my*iz);
            int ib = ix*dbx+iy*dby;
            int in = ix*dnx+iy*dny;
            field_array[i] +=
              (2*face-1)*(block_flux_array[ib] - neighbor_flux_array[in]);
          }
        }
      } // perform_correction[axis][face]
    } // face
  } // rank >= 3
}

/// divide each cell in input_array for which the flux correction was performed
/// and store the result in the corresponding cell of output_array 
static void divide_flux_corrected_cells_by_density_
(const cello_float * const input_array, const cello_float * const density_field,
 cello_float* const output_array, const bool (&correction_performed)[3][2],
 const int rank, int mx, int my, int mz,
 int nx, int ny, int nz)
{

  {
    const int axis = 0;
    for (int face=0; face<2; face++) {
      if (correction_performed[axis][face]){

        int ix = (face == 0) ? 0 : nx-1;
        for (int iz=0; iz<nz; iz++) {
          for (int iy=0; iy<ny; iy++) {
            int i=ix+mx*(iy+my*iz);
            output_array[i] = input_array[i] / density_field[i];
          }
        }
      }
    }

  } // rank >= 1

  if (rank >= 2) {
    const int axis = 1;
    for (int face=0; face<2; face++) {
      if (correction_performed[axis][face]){

        int iy = (face == 0) ? 0 : ny-1;
        for (int iz=0; iz<nz; iz++) {
          for (int ix=0; ix<nx; ix++) {
            int i=ix+mx*(iy+my*iz);
            output_array[i] = input_array[i] / density_field[i];
          }
        }
      }
    }

  } // rank >= 2

  if (rank >= 3) {
    const int axis=2;
    for (int face=0; face<2; face++) {
      if (correction_performed[axis][face]){

        int iz = (face == 0) ? 0 : nz-1;
        for (int iy=0; iy<ny; iy++) {
          for (int ix=0; ix<nx; ix++) {
            int i=ix+mx*(iy+my*iz);
            output_array[i] = input_array[i] / density_field[i];
          }
        }
      } // correction_performed[axis][face]
    } // face
  } // rank >= 3

}


void MethodFluxCorrect::flux_correct_(Block * block)
{

  Field field = block->data()->field();
  FluxData * flux_data = block->data()->flux_data();
  const int nf = flux_data->num_fields();

  // Perform flux-correction
  if (enable_ && block->is_leaf()) {
    if (nf  == 0){
      return;
    }

    const int level = block->level();

    const int rank = cello::rank();

    int gx,gy,gz;
    field.ghost_depth (0,&gx,&gy,&gz);

    int nx,ny,nz;
    field.size(&nx,&ny,&nz);

    int mx = nx + 2*gx;
    int my = ny + 2*gy;
    int mz = nz + 2*gz;

    // determine which faces require flux corrections
    bool perform_correction[3][2];
    for (int axis=0; axis < 3; axis++){
      for (int face = 0; face < 2; face++){
        if ((axis < rank) && (block->face_level(axis,face) > level)){
          perform_correction[axis][face] = true;
        } else {
          perform_correction[axis][face] = false;
        }
      }
    }
    int i_f_density = -1; // will be used to store i_f for density
    for (int i_f=0; i_f<nf; i_f++) {
      const int index_field = flux_data->index_field(i_f);
      const std::string field_name = field.field_name(index_field);

      if (field_name == "density"){
        i_f_density = i_f;
      }
    }

    // Allocate scratch space (if not already allocated necessary). It's
    // technically possible to do all of these operations in-place, (without
    // scratch space), but that gets complex (you need to be very careful
    // iterating over cells on the edges of the active zone multiple times).
    std::size_t num_field_elements = (std::size_t) (mx * my * mz);
    if (scratch_.size() == 0) { scratch_.resize(2 * num_field_elements); }

    std::size_t active_zone_offset = (std::size_t) (gx + mx*(gy + my*gz));
    cello_float *old_rho_array, *temp_cons_array;
    old_rho_array = scratch_.data() + active_zone_offset;
    temp_cons_array = scratch_.data() + num_field_elements + active_zone_offset;

    // load the density array
    cello_float* density_array = nullptr;
    if (field.is_field("density")){
      density_array = (cello_float*) field.unknowns("density");

      // copy the values in the density_array (we could be more selective about
      // what we copy)
      for (int iz=0; iz<nz; iz++) {
        for (int iy=0; iy<ny; iy++) {
          for (int ix=0; ix<nx; ix++) {
            int i=ix + mx*(iy + my*iz);
            old_rho_array[i] = density_array[i];
          }
        }
      }
    }

    // perform the density flux correction
    if (i_f_density > -1){
      flux_correct_helper_(density_array, mx,my,mz, nx,ny,nz,
                           i_f_density, flux_data, rank, perform_correction);
    }


    // perform the flux corrections for the other fields
    Grouping * groups = cello::field_groups();
    for (int i_f=0; i_f<nf; i_f++) {
      const int index_field = flux_data->index_field(i_f);
      const std::string field_name = field.field_name(index_field);

      if (i_f == i_f_density){ // density flux correction already happened
        continue;
      }

      cello_float* field_array = (cello_float*) field.unknowns(index_field);

      if (groups->is_in(field_name, "make_field_conservative")){
        // Handle flux corrections for fields that must be multiplied by the
        // density to be made conservative

        ASSERT1("MethodFluxCorrect::flux_correct_",
                ("The \"density\" field must exist to perform flux "
                 "corrections on \"%s\"."), field_name.c_str(),
                density_array != nullptr);

        // compute the conserved quantity
        for (int iz=0; iz<nz; iz++) {
          for (int iy=0; iy<ny; iy++) {
            for (int ix=0; ix<nx; ix++) {
              int i=ix + mx*(iy + my*iz);
              temp_cons_array[i] = old_rho_array[i]*field_array[i];
            }
          }
        }

        // now perform the actual flux correction
        flux_correct_helper_(temp_cons_array, mx,my,mz, nx,ny,nz,
                             i_f, flux_data, rank, perform_correction);

        // now, divide the updated quantities by appropriate values in
        // density_array and write the result to field_array
        divide_flux_corrected_cells_by_density_(temp_cons_array, density_array,
                                                field_array, perform_correction,
                                                rank, mx,my,mz, nx,ny,nz);
      } else {
        // the quantity is already in conserved form. We can just apply the
        // flux correction directly on the field_data
        flux_correct_helper_(field_array, mx,my,mz, nx,ny,nz,
                             i_f, flux_data, rank, perform_correction);
      }
    }
  }
}
