// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodFluxCorrect.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-03-19
/// @brief    Implementation of the flux_correct method

#include "problem.hpp"
#include "charm_simulation.hpp"

// #define DEBUG_FLUX_CORRECT

#ifdef DEBUG_FLUX_CORRECT
#   define TRACE_FLUX(BLOCK,MSG)                        \
  CkPrintf ("TRACE_FLUX %10s %d %s\n",                  \
            BLOCK->name().c_str(),CkMyPe(),MSG);        \
  fflush(stdout);
#else
#   define TRACE_FLUX(BLOCK,MSG)  /* ... */
#endif

//----------------------------------------------------------------------

MethodFluxCorrect::MethodFluxCorrect (std::string group) throw() 
  : Method (),
    group_(group),
    ir_pre_(-1)
{
  cello::simulation()->new_refresh_set_name(ir_post_,name());
  Refresh * refresh_post = cello::refresh(ir_post_);
  refresh_post->add_all_fields();

  ir_pre_ = add_new_refresh_();
  Refresh * refresh_pre = cello::refresh(ir_pre_);
  refresh_pre->set_callback(CkIndex_Block::p_method_flux_correct());
  refresh_pre->add_all_fields();
}

//----------------------------------------------------------------------

void MethodFluxCorrect::compute ( Block * block) throw()
{
  TRACE_FLUX(block,"01 compute()");
  cello::refresh(ir_pre_)->set_active(block->is_leaf());
  block->new_refresh_start
    (ir_pre_, CkIndex_Block::p_method_flux_correct());
}

//----------------------------------------------------------------------

void Block::Block::p_method_flux_correct()
{
  TRACE_FLUX(this,"02 p_method_flux_correct()");
  MethodFluxCorrect * method = static_cast<MethodFluxCorrect*> (this->method());
  method->compute_continue(this);
}

//----------------------------------------------------------------------

void MethodFluxCorrect::compute_continue( Block * block ) throw()
{
  TRACE_FLUX(block,"03 compute_continue()");
  const int rank = cello::rank();
  if (block->is_leaf()) {

    // Loop over neighbors in different levels

    // ... only include facet neighbors, not edges or corners
    const int min_face_rank = rank - 1;
    Index index = block->index();
    int neighbor_type = neighbor_flux;
    int min_level = 0;
    int root_level = 0;
    ItNeighbor it_neighbor = block->it_neighbor
      (min_face_rank,index, neighbor_type,min_level,root_level);
    
    Grouping * groups = cello::field_groups();
    int ng=groups->size(group_);
    for (int ig=0; ig<ng; ig++) {
      std::string field_name = groups->item(group_,ig);
#ifdef DEBUG_FLUX_CORRECT      
      CkPrintf ("DEBUG_FLUX MethodFluxCorrect field %s\n",
                field_name.c_str());
#endif
      
    }
  }  
  block->compute_done();
}
