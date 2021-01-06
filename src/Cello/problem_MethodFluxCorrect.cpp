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
(std::string group, bool enable, double min_digits) throw() 
  : Method (),
    group_(group),
    enable_(enable),
    min_digits_(min_digits),
    field_sum_(),
    field_sum_0_(),
    ir_pre_(-1)
{
  // Set up post-refresh to refresh all conserved fields in group_
  cello::simulation()->new_refresh_set_name(ir_post_,name());
  Grouping * groups = cello::field_groups();
  const int nf=groups->size(group_);
  for (int i_f=0; i_f<nf; i_f++) {
    cello::refresh(ir_post_)->add_field(groups->item(group_,i_f));
  }
  
  //  ir_pre_ = add_new_refresh_(neighbor_flux);
  ir_pre_ = add_new_refresh_(neighbor_leaf);
  Refresh * refresh_pre = cello::refresh(ir_pre_);
  cello::simulation()->new_refresh_set_name(ir_pre_,name()+"_fluxes");
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
  // sum mass, momentum, energy

  field_sum_.resize(nf);
  field_sum_0_.resize(nf);
}

//----------------------------------------------------------------------

void MethodFluxCorrect::compute ( Block * block) throw()
{
  cello::refresh(ir_pre_)->set_active(block->is_leaf());

  block->new_refresh_start
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

    const int index_density = field.field_id("density");
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
      cello::monitor()->print
        ("Method", "Field %s sum %20.16Le conserved to %g digits of %d",
         field.field_name(index_field).c_str(),
         field_sum_[i_f],
         digits,
         cello::digits_max(precision));
      unit_func("MethodFluxCorrect precision (density)");
      if (index_field == index_density) unit_assert (digits >= min_digits_);
    }
  }

  block->data()->flux_data()->deallocate();

  block->compute_done();
}

//======================================================================

void MethodFluxCorrect::flux_correct_(Block * block)
{
    
  Field field = block->data()->field();
  FluxData * flux_data = block->data()->flux_data();
  const int nf = flux_data->num_fields();

  // Perform flux-correction
  if (enable_ && block->is_leaf()) {

    const int level = block->level();

    const int rank = cello::rank();

    int ix,iy,iz;
    int nx,ny,nz;
    field.size(&nx,&ny,&nz);

    for (int i_f=0; i_f<nf; i_f++) {

      const int index_field = flux_data->index_field(i_f);
      cello_float * field_array = (cello_float*) field.unknowns(index_field);
      int mx,my,mz;
      int gx,gy,gz;
      field.dimensions (index_field,&mx,&my,&mz);
      field.ghost_depth(index_field,&gx,&gy,&gz);

      int axis,face,level_face;

      axis=0;
      // X axis
      for (face=0; face<2; face++) {
        level_face = block->face_level(axis,face);
        if (level_face > level) {
          int dbx,dby,dbz;
          int dnx,dny,dnz;
          auto block_fluxes    = flux_data->block_fluxes(axis,face,i_f);
          auto neighbor_fluxes = flux_data->neighbor_fluxes(axis,face,i_f);
          auto block_flux_array =
            block_fluxes->flux_array(&dbx,&dby,&dbz);
          auto neighbor_flux_array =
            neighbor_fluxes->flux_array(&dnx,&dny,&dnz);
          
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
          level_face = block->face_level(axis,face);
          if (level_face > level) {
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
          level_face = block->face_level(axis,face);
          if (level_face > level) {
            int dbx,dby,dbz;
            int dnx,dny,dnz;
            auto block_fluxes    = flux_data->block_fluxes(axis,face,i_f);
            auto neighbor_fluxes = flux_data->neighbor_fluxes(axis,face,i_f);
            auto block_flux_array =
              block_fluxes->flux_array(&dbx,&dby,&dbz);
            auto neighbor_flux_array =
              neighbor_fluxes->flux_array(&dnx,&dny,&dnz);
          
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
          } // level_face > level
        } // face
      } // rank >= 3
    }
  }
}
