// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodFluxCorrect.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-03-19
/// @brief    Implementation of the flux_correct method

#include "problem.hpp"
#include "charm_simulation.hpp"

// #define DEBUG_METHOD_FLUX_CORRECT

#ifdef DEBUG_METHOD_FLUX_CORRECT
#   define TRACE_FLUX_CORRECT(BLOCK,MSG)                                \
  CkPrintf ("DEBUG_METHOD_FLUX_CORRECT level-%d %s %d %s MethodFluxCorrect::compute()\n", \
            BLOCK->level(),MSG,CkMyPe(),BLOCK->name().c_str());         \
  fflush(stdout);
#else
#   define TRACE_FLUX_CORRECT(BLOCK,MSG)  /* ... */  
#endif

//----------------------------------------------------------------------

MethodFluxCorrect::MethodFluxCorrect (std::string group) throw() 
  : Method (),
    group_(group),
    field_sum_(),
    field_sum_0_(),
    ir_pre_(-1)
{
  cello::simulation()->new_refresh_set_name(ir_post_,name());
  Refresh * refresh_post = cello::refresh(ir_post_);
  refresh_post->add_all_fields();

  ir_pre_ = add_new_refresh_(neighbor_flux);
  Refresh * refresh_pre = cello::refresh(ir_pre_);
  refresh_pre->set_callback(CkIndex_Block::p_method_flux_correct_refresh());
  refresh_pre->add_all_fields();

  // sum mass, momentum, energy
  field_sum_.resize(3);
  field_sum_0_.resize(3);
}

//----------------------------------------------------------------------

void MethodFluxCorrect::compute ( Block * block) throw()
{
  TRACE_FLUX_CORRECT(block,"1 ENTER");
  cello::refresh(ir_pre_)->set_active(block->is_leaf());
  TRACE_FLUX_CORRECT(block,"1 EXIT  ");
  
  block->new_refresh_start
    (ir_pre_, CkIndex_Block::p_method_flux_correct_refresh());
  //  block->debug_new_refresh(__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void Block::Block::p_method_flux_correct_refresh()
{
  //  this->debug_new_refresh(__FILE__,__LINE__);
  TRACE_FLUX_CORRECT(this,"2 ENTER");
  static_cast<MethodFluxCorrect*>
    (this->method())->compute_continue_refresh(this);
  TRACE_FLUX_CORRECT(this,"2 EXIT ");
}

//----------------------------------------------------------------------

void MethodFluxCorrect::compute_continue_refresh( Block * block ) throw()
{
  //   block->debug_new_refresh(__FILE__,__LINE__);
  TRACE_FLUX_CORRECT(block,"3 ENTER");

  // @@@@@ TEMPORARY @@@@@
  //  TRACE_FLUX_CORRECT(block,"3 EXIT ");
  //  block->compute_done();
  //  return;
  // @@@@@ TEMPORARY @@@@@

  // block sum mass, momentum, energy
  long double reduce[3] = {0.0, 0.0, 0.0};

  Field field = block->data()->field();
  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth (0,&gx,&gy,&gz);

  int precision = field.precision (field.field_id("density"));

  union {
    char        * d_c;
    float       * d_s;
    double      * d_d;
    long double * d_q;
  };
  
  d_c  = field.values("density");

  const int rank = cello::rank();

  reduce[0] = 1;
  
  if (precision == precision_single) {
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
        for (int ix=gx; ix<mx-gx; ix++) {
          int i=ix + mx*(iy + my*iz);
          reduce[1] += d_s[i];
        }    
      }    
    }
  } else if (precision == precision_double) {
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
        for (int ix=gx; ix<mx-gx; ix++) {
          int i=ix + mx*(iy + my*iz);
          reduce[1] += d_d[i];
        }    
      }    
    }
  } else if (precision == precision_quadruple) {
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
        for (int ix=gx; ix<mx-gx; ix++) {
          int i=ix + mx*(iy + my*iz);
          reduce[1] += d_q[i];
        }    
      }    
    }
  }

  // scale by relative mesh cell volume/area
  const int w = (1 << block->level())*rank;
  reduce[1] /= w;

  CkCallback callback (CkIndex_Block::r_method_flux_correct_sum_fields(NULL), 
                       block->proxy_array());

  block->contribute
    (2*sizeof(long double), &reduce, sum_long_double_n_type, callback);
  TRACE_FLUX_CORRECT(block,"3 EXIT ");
}

//----------------------------------------------------------------------

void Block::Block::r_method_flux_correct_sum_fields(CkReductionMsg * msg)
{
  //   this->debug_new_refresh(__FILE__,__LINE__);
  TRACE_FLUX_CORRECT(this,"4 ENTER");
  static_cast<MethodFluxCorrect*>
    (this->method())->compute_continue_sum_fields(this,msg);
  TRACE_FLUX_CORRECT(this,"4 EXIT ");
}

//----------------------------------------------------------------------

void MethodFluxCorrect::compute_continue_sum_fields
( Block * block, CkReductionMsg * msg) throw()
{
  //   block->debug_new_refresh(__FILE__,__LINE__);
  TRACE_FLUX_CORRECT(block,"5 ENTER");

  long double * data = (long double *) msg->getData();

  field_sum_[1] = data[1];

  delete msg;

  if (block->index().is_root()) {
    // save initial total mass
    if (block->cycle() == 0) {
      field_sum_0_[1] = field_sum_[1];
    }
    Field field = block->data()->field();
    int precision = field.precision (field.field_id("density"));
    const int precision_max =
      (precision == precision_single) ? 7 :
      (precision == precision_double) ? 16 : 34;
    cello::monitor()->print
      ("Method", "Mass %Lg conserved to within %Lg digits out of %d",
       field_sum_[1],
       -log10(std::min(cello::err_abs(field_sum_[1],field_sum_0_[1]),
                       cello::err_rel(field_sum_[1],field_sum_0_[1]))),
       precision_max);
  }

  const int rank = cello::rank();

  if (block->is_leaf()) {

    FluxData * flux_data = block->data()->flux_data();

    auto field_names = block->data()->field().groups()->group_list("conserved");
    
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    // Loop over neighbors in different levels

    // ... only include facet neighbors, not edges or corners

    ItNeighbor it_neighbor =
      block->it_neighbor (cello::rank() - 1,block->index(), neighbor_flux);
    
    Grouping * groups = cello::field_groups();
    int ng=groups->size(group_);
    
    // Loop over each face with a level jump
    int of3[3];
    while (it_neighbor.next(of3)) {
      Index index_neighbor = it_neighbor.index();
      for (int ig=0; ig<ng; ig++) {
        std::string field_name = groups->item(group_,ig);
      }
    }
  }  
  TRACE_FLUX_CORRECT(block,"5 EXIT ");
  block->compute_done();
}
