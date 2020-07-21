// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodFluxCorrect.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-03-19
/// @brief    Implementation of the flux_correct method

#include "problem.hpp"
#include "charm_simulation.hpp"

// #define DEBUG_FLUX_CORRECT

//----------------------------------------------------------------------

MethodFluxCorrect::MethodFluxCorrect (std::string group, int sign) throw() 
  : Method (),
    group_(group),
    sign_(sign),
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
  // BUG WORKAROUND: also add a dummy field to prevent refresh hang
  refresh_pre->add_field("density");
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

  union {
    char        * values;
    float       * values_4;
    double      * values_8;
    long double * values_16;
  };

  FluxData * flux_data = block->data()->flux_data();

  const int nf = flux_data->num_fields();
  long double * reduce = new long double [nf+1];
  std::fill_n(reduce,nf+1,0.0);
  reduce[0] = nf;

  if (block->is_leaf()) {

    for (int i_f=0; i_f<nf; i_f++) {

      const int index_field = flux_data->index_field(i_f);

      values = field.values(index_field);

      const int precision = field.precision(index_field);

      if (precision == precision_single) {
        for (int iz=gz; iz<mz-gz; iz++) {
          for (int iy=gy; iy<my-gy; iy++) {
            for (int ix=gx; ix<mx-gx; ix++) {
              int i=ix + mx*(iy + my*iz);
              reduce[i_f+1] += values_4[i];
            }    
          }    
        }
      } else if (precision == precision_double) {
        double sum = 0.0;
        for (int iz=gz; iz<mz-gz; iz++) {
          for (int iy=gy; iy<my-gy; iy++) {
            for (int ix=gx; ix<mx-gx; ix++) {
              int i=ix + mx*(iy + my*iz);
              reduce[i_f+1] += values_8[i];
              sum += values_8[i];
            }    
          }    
        }
      } else if (precision == precision_quadruple) {
        for (int iz=gz; iz<mz-gz; iz++) {
          for (int iy=gy; iy<my-gy; iy++) {
            for (int ix=gx; ix<mx-gx; ix++) {
              int i=ix + mx*(iy + my*iz);
              reduce[i_f+1] += values_16[i];
            }    
          }    
        }
      }

      // scale by relative mesh cell volume/area

      const int level = block->level();
      
      const int w = 1 << level*cello::rank();
      reduce[i_f+1] /= w;
#ifdef DEBUG_FLUX_CORRECT      
      CkPrintf ("DEBUG_FLUX_CORRECT weight = %d\n",w);

      CkPrintf ("DEBUG_FLUX_CORRECT %s %d field %d field_sum %Lg\n",
                block->name().c_str(),block->cycle(),index_field,reduce[i_f+1]);
#endif      
    }
  }
  
  CkCallback callback (CkIndex_Block::r_method_flux_correct_sum_fields(NULL), 
                       block->proxy_array());

  block->contribute
    ((nf+1)*sizeof(long double), reduce, sum_long_double_n_type, callback);
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
      int precision = field.precision (index_field);
      const int digits_max =
        (precision == precision_single) ? 7 :
        (precision == precision_double) ? 16 : 34;
      cello::monitor()->print
        ("Method", "Field %s sum %Lg conserved to %Lg digits of %d",
         field.field_name(index_field).c_str(),
         field_sum_[i_f],
         -log10(std::min(cello::err_abs(field_sum_[i_f],field_sum_0_[i_f]),
                         cello::err_rel(field_sum_[i_f],field_sum_0_[i_f]))),
         digits_max);
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
  int c3[3]={0,0,0};
  double s3[3] = {0,0,0};
  if (block->is_leaf()) {

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

      long double sum = 0.0;
      axis=0;
      // X axis
      long double sum_block,sum_neighbor;
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
          int sign = sign_*(1 - face*2);
          sum_block = sum_neighbor = 0.0;
          for (iz=0; iz<nz; iz++) {
            for (iy=0; iy<ny; iy++) {
              int i=ix+mx*(iy+my*iz);
              int ib = iy*dby+iz*dbz;
              int in = iy*dny+iz*dnz;
              sum_block += block_flux_array[ib];
              sum_neighbor += neighbor_flux_array[in];
              double update = sign*
                (block_flux_array[ib] - neighbor_flux_array[in]);
              field_array[i] += update;
              ++c3[axis];
              s3[axis]+=update*update;
            }
          }
#ifdef DEBUG_FLUX_CORRECT      
          CkPrintf ("DEBUG_FLUX_CORRECT %s %d field %d axis %d face %d flux_sum %Lg %Lg\n",
                    block->name().c_str(),block->cycle(),index_field,axis,face,sum_block,sum_neighbor);
#endif
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
            int sign = sign_*(1 - face*2);
            sum_block = sum_neighbor = 0.0;
            for (iz=0; iz<nz; iz++) {
              for (ix=0; ix<nx; ix++) {
                int i=ix+mx*(iy+my*iz);
                int ib = ix*dbx+iz*dbz;
                int in = ix*dnx+iz*dnz;
                sum_block += block_flux_array[ib];
                sum_neighbor += neighbor_flux_array[in];
                double update = sign*
                  (block_flux_array[ib] - neighbor_flux_array[in]);
                field_array[i] += update;
                ++c3[axis];
                s3[axis]+=update*update;
              }
            }
#ifdef DEBUG_FLUX_CORRECT      
            CkPrintf ("DEBUG_FLUX_CORRECT %s %d field %d axis %d face %d flux_sum %Lg %Lg\n",
                      block->name().c_str(),block->cycle(),index_field,axis,face,sum_block,sum_neighbor);
#endif
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
            int sign = sign_*(1 - face*2);
            sum_block = sum_neighbor = 0.0;
            for (iy=0; iy<ny; iy++) {
              for (ix=0; ix<nx; ix++) {
                int i=ix+mx*(iy+my*iz);
                int ib = ix*dbx+iy*dby;
                int in = ix*dnx+iy*dny;
                sum_block += block_flux_array[ib];
                sum_neighbor += neighbor_flux_array[in];
                double update = sign*
                  (block_flux_array[ib] - neighbor_flux_array[in]);
                field_array[i] += update;
                ++c3[axis];
                s3[axis]+=update*update;
              }
            }
#ifdef DEBUG_FLUX_CORRECT      
            CkPrintf ("DEBUG_FLUX_CORRECT %s %d field %d axis %d face %d flux_sum %Lg %Lg\n",
                      block->name().c_str(),block->cycle(),index_field,axis,face,sum_block,sum_neighbor);
#endif
          } // level_face > level
        } // face
      } // rank >= 3
    }
  }

}
