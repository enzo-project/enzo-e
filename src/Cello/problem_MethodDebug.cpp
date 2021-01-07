// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodDebug.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-03-19
/// @brief    Implementation of the debug method

#include "problem.hpp"
#include "charm_simulation.hpp"
#include "test.hpp"

// #define DEBUG_METHOD_DEBUG
// #define WRITE_FILES
// #define WRITE_FILES_CYCLE_START 90
// #define WRITE_FILES_CYCLE_STOP 95
//----------------------------------------------------------------------

MethodDebug::MethodDebug (int num_fields) throw() 
  : Method (),
    field_sum_(),
    field_min_(),
    field_max_(),
    field_count_()
{
  // Set up post-refresh to refresh all conserved fields in group_

  cello::simulation()->new_refresh_set_name(ir_post_,name());
  cello::refresh(ir_post_)->add_all_fields();

  field_sum_.resize(num_fields);
  field_min_.resize(num_fields);
  field_max_.resize(num_fields);
  field_count_.resize(num_fields);
}

//----------------------------------------------------------------------

void MethodDebug::compute ( Block * block) throw()
{
  // accumulate local reductions for global 

  Field field = block->data()->field();
  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth (0,&gx,&gy,&gz);

  const int nf = field.num_permanent();
  long double * reduce = new long double [1+4*nf];
  reduce[0] = nf;
  for (int i_f=0; i_f<nf; i_f++) {
    reduce[4*i_f+1] = std::numeric_limits<long double>::max();
    reduce[4*i_f+2] = -std::numeric_limits<long double>::max();
    reduce[4*i_f+3] = 0;
    reduce[4*i_f+4] = 0;
  }

  if (block->is_leaf()) {

    const double rel_vol = cello::relative_cell_volume (block->level());
    for (int i_f=0; i_f<nf; i_f++) {

      cello_float * values = (cello_float *) field.values(i_f);

      for (int iz=gz; iz<mz-gz; iz++) {
        for (int iy=gy; iy<my-gy; iy++) {
          for (int ix=gx; ix<mx-gx; ix++) {
            int i=ix + mx*(iy + my*iz);
            reduce[4*i_f+1] = std::min(reduce[4*i_f+1],(long double)values[i]);
            reduce[4*i_f+2] = std::max(reduce[4*i_f+2],(long double)values[i]);
            reduce[4*i_f+3] += values[i];
            ++reduce[4*i_f+4];
          }    
        }    
      }
      reduce[4*i_f+3] *= rel_vol;
    }
  }
#ifdef DEBUG_METHOD_DEBUG
  CkPrintf ("DEBUG_METHOD_DEBUG nf = %d\n",nf);
  for (int i_f=0; i_f<nf; i_f++) {
    CkPrintf ("DEBUG_METHOD_DEBUG %s %s %g %g %g %g\n",block->name().c_str(),
              field.field_name(i_f).c_str(),reduce[4*i_f+1],
              reduce[4*i_f+2],reduce[4*i_f+3],reduce[4*i_f+1]);
  }    
#endif    
  
  CkCallback callback (CkIndex_Block::r_method_debug_sum_fields(NULL), 
                       block->proxy_array());

  block->contribute
    ((4*nf+1)*sizeof(long double), reduce, r_reduce_method_debug_type, callback);

  delete [] reduce;
}

//----------------------------------------------------------------------

void Block::Block::r_method_debug_sum_fields(CkReductionMsg * msg)
{
  static_cast<MethodDebug*>
    (this->method())->compute_continue_sum_fields(this,msg);
}

//----------------------------------------------------------------------

void MethodDebug::compute_continue_sum_fields
( Block * block, CkReductionMsg * msg) throw()
{

  Field field = block->data()->field();
  const int nf = field.num_permanent();
  long double * data = (long double *) msg->getData();
  int id = 1;
  int i_f = 0;
  for (int i_f=0; i_f<nf; i_f++) {
    field_min_[i_f] = data[id++];
    field_max_[i_f] = data[id++];
    field_sum_[i_f] = data[id++];
    field_count_[i_f] = data[id++];
  }
  delete msg;


  // Write conserved field sums to output (root block only)
  
#ifdef WRITE_FILES

  const int cycle = block->cycle();
  if (block->is_leaf() &&
      ( WRITE_FILES_CYCLE_START <= cycle &&
        cycle <= WRITE_FILES_CYCLE_STOP)) {
    int mx,my,mz;
    int gx,gy,gz;
    field.dimensions (0,&mx,&my,&mz);
    field.ghost_depth (0,&gx,&gy,&gz);
    for (int i_f=0; i_f<nf; i_f++) {
      char buffer[81];
      snprintf(buffer,80,"field-%03d-%02d-%s.data",block->cycle(),CkMyPe(),
               field.field_name(i_f).c_str());
      FILE * fp = fopen (buffer,"a");
      cello_float * values = (cello_float*)field.values(i_f);
      for (int iz=gz; iz<mz-gz; iz++) {
        for (int iy=gy; iy<my-gy; iy++) {
          for (int ix=gx; ix<mx-gx; ix++) {
            int i=ix + mx*(iy + my*iz);
            fprintf(fp,"%20.18f\n",values[i]);
          }
        }
      }
      fclose(fp);
    }
  }
#endif
    
  if (block->index().is_root()) {

    int nx,ny,nz;
    cello::hierarchy()->root_size(&nx,&ny,&nz);
    long int root_cells = nx*ny*nz;
    for (int i_f=0; i_f<nf; i_f++) {

      cello::monitor()->print
        ("Method", "Field %s min %20.16Le",field.field_name(i_f).c_str(),
         field_min_[i_f]);
      cello::monitor()->print
        ("Method", "Field %s max %20.16Le",field.field_name(i_f).c_str(),
         field_max_[i_f]);
      cello::monitor()->print
        ("Method", "Field %s avg %20.16Le", field.field_name(i_f).c_str(),
         field_sum_[i_f]/root_cells);
    }
  }

  block->compute_done();
}

