// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodDebug.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-03-19
/// @brief    Implementation of the debug method

#include "problem.hpp"
#include "charm_simulation.hpp"
#include "test.hpp"

// #define DEBUG_DEBUG
//----------------------------------------------------------------------

MethodDebug::MethodDebug
(int num_fields,
 int num_particles,
 bool l_print,
 bool l_coarse,
 bool l_ghost
 ) throw()
  : Method (),
    num_fields_(num_fields),
    num_particles_(num_particles),
    field_sum_(),
    field_min_(),
    field_max_(),
    field_count_(),
    particle_sum_(),
    particle_min_(),
    particle_max_(),
    particle_count_(),
    l_print_(l_print),
    l_coarse_(l_coarse),
    l_ghost_(l_ghost)
{
  // Set up post-refresh to refresh all fields

  cello::simulation()->refresh_set_name(ir_post_,name());
  cello::refresh(ir_post_)->add_all_fields();

  field_sum_.resize(num_fields);
  field_min_.resize(num_fields);
  field_max_.resize(num_fields);
  field_count_.resize(num_fields);

  for (int i=0; i<3; i++) {
    particle_sum_[i].resize(num_particles);
    particle_min_[i].resize(num_particles);
    particle_max_[i].resize(num_particles);
    particle_count_[i].resize(num_particles);
  }
}

//----------------------------------------------------------------------

void MethodDebug::compute ( Block * block) throw()
{
  const int num_reduce = 4*(num_fields_+3*num_particles_);
  long double * reduce = new long double [1+num_reduce];
  reduce[0] = num_reduce+1;
  const int kmin=0;
  const int kmax=1;
  const int ksum=2;
  const int knum=3;
  for (int k=1; k<num_reduce; k+=4) {
    reduce[k+kmin] = std::numeric_limits<long double>::max();
    reduce[k+kmax] = -std::numeric_limits<long double>::max();
    reduce[k+ksum] = 0;
    reduce[k+knum] = 0;
  }

  if (block->is_leaf()) {

    // accumulate local reductions for global

    Field field = block->data()->field();
    int mx,my,mz;
    int gx,gy,gz;
    field.dimensions (0,&mx,&my,&mz);
    field.ghost_depth (0,&gx,&gy,&gz);

    const double rel_vol = cello::relative_cell_volume (block->level());
    int k=1;
    for (int index_field=0; index_field<num_fields_; index_field++) {

      cello_float * values = (cello_float *) field.values(index_field);

      for (int iz=gz; iz<mz-gz; iz++) {
        for (int iy=gy; iy<my-gy; iy++) {
          for (int ix=gx; ix<mx-gx; ix++) {
            int i=ix + mx*(iy + my*iz);
            reduce[k+kmin] = std::min(reduce[k],(long double)(values[i]));
            reduce[k+kmax] = std::max(reduce[k+1],(long double)(values[i]));
            reduce[k+ksum] += values[i];
            reduce[k+knum] += rel_vol;
          }
        }
      }
      k += 4;
    }

    // particles
    Particle particle = block->data()->particle();
    const int mb = particle.batch_size();
    std::vector<double> position[3];
    position[0].resize(mb); position[0].clear();
    position[1].resize(mb); position[1].clear();
    position[2].resize(mb); position[2].clear();
    for (int it=0; it<num_particles_; it++) {
      const int nb = particle.num_batches(it);
      for (int ib=0; ib<nb; ib++) {
        particle.position(it,ib,position[0].data(),position[1].data(),position[2].data());
        const int np = particle.num_particles(it,ib);
        for (int i=0; i<cello::rank(); i++) {
          for (int ip=0; ip<np; ip++) {
            double value = position[i][ip];
            reduce[k+4*i+kmin] = std::min(reduce[k+4*i+0],(long double)(value));
            reduce[k+4*i+kmax] = std::max(reduce[k+4*i+1],(long double)(value));
            reduce[k+4*i+ksum] += value;
            reduce[k+4*i+knum] += 1;
          }
        }
      }
      k += 4*3;
    }

    ASSERT2("MethodDebug::compute()",
            "reduce array mismatch %d != %d",
            k,num_reduce+1,(k == num_reduce+1));
  }

#ifdef DEBUG_DEBUG  
  {
    int id = 0;
    Field field = block->data()->field();
    for (int i_f=0; i_f<num_fields_; i_f++) {
      std::string name = field.field_name(i_f).c_str();
      cello::monitor()->print
        ("Method", "Field %s %s min %Lg max %Lg sum %Lg cnt %Lg",
         name.c_str(),block->name().c_str(),
         reduce[id],reduce[id+1],reduce[id+2],reduce[id+3]);
      id+=4;
    }
    Particle particle = block->data()->particle();
    for (int it=0; it<num_particles_; it++) {
      const std::string name = particle.type_name(it).c_str();
      cello::monitor()->print
        ("Method", "Particle %s %s num_particles %d",
         name.c_str(),block->name().c_str(),particle.num_particles(it));
      for (int i=0; i<3; i++) {
        const char axis[3] = {'X','Y','Z'};
        cello::monitor()->print
          ("Method", "Particle %s %c %s min %Lg max %Lg sum %Lg cnt %Lg",
           name.c_str(),axis[i],block->name().c_str(),
           reduce[id],reduce[id+1],reduce[id+2],reduce[id+3]);
        id+=4;
      }
    }
  }
#endif
    CkCallback callback (CkIndex_Block::r_method_debug_sum_fields(NULL),
                       block->proxy_array());

  block->contribute
    ((1+num_reduce)*sizeof(long double), reduce,
     r_reduce_method_debug_type, callback);

  delete [] reduce;
}

//----------------------------------------------------------------------

void Block::r_method_debug_sum_fields(CkReductionMsg * msg)
{
  static_cast<MethodDebug*>
    (this->method())->compute_continue(this,msg);
}

//----------------------------------------------------------------------

void MethodDebug::compute_continue
( Block * block, CkReductionMsg * msg) throw()
{

  long double * data = (long double *) msg->getData();
  int id = 1;
  for (int index_field=0; index_field<num_fields_; index_field++) {
    field_min_[index_field] = data[id];
    field_max_[index_field] = data[id+1];
    field_sum_[index_field] = data[id+2];
    field_count_[index_field] = data[id+3];
    id+=4;
  }
  for (int it=0; it<num_particles_; it++) {
    for (int i=0; i<3; i++) {
      particle_min_[i][it] = data[id];
      particle_max_[i][it] = data[id+1];
      particle_sum_[i][it] = data[id+2];
      particle_count_[i][it] = data[id+3];
      id+=4;
    }
  }
  const int num_reduce = 4*(num_fields_+3*num_particles_);
  ASSERT2("MethodDebug::compute_continue()",
          "reduce array mismatch %d != %d",
          id,num_reduce+1,(id == num_reduce+1));

  delete msg;

  Field field = block->data()->field();
  Particle particle = block->data()->particle();
  if (block->index().is_root()) {
    int nx,ny,nz;
    cello::hierarchy()->root_size(&nx,&ny,&nz);
    //    long int root_cells = nx*ny*nz;
    for (int i_f=0; i_f<num_fields_; i_f++) {
      std::string name = field.field_name(i_f).c_str();
      cello::monitor()->print
        ("Method", "Field %s min %20.16Lg avg %20.16Lg max %20.16Lg",name.c_str(),
         field_min_[i_f],field_sum_[i_f]/field_count_[i_f],field_max_[i_f]);
    }
    for (int it=0; it<num_particles_; it++) {
      for (int i=0; i<3; i++) {
        const char axis[3] = {'X','Y','Z'};
        const std::string name = particle.type_name(it).c_str();
        cello::monitor()->print
          ("Method", "Particle %s %c min avg max %20.16Lg %20.16Lg %20.16Lg",
           name.c_str(),axis[i],
           particle_min_[i][it],
           particle_sum_[i][it]/particle_count_[i][it],
           particle_max_[i][it]);
      }
    }
  }

  if (block->is_leaf()) {

    for (int i_f=0; i_f<num_fields_; i_f++) {

      int mx,my,mz;
      field.dimensions (i_f,&mx,&my,&mz);
      cello_float * values = (cello_float*)field.values(i_f);

      int gx=0,gy=0,gz=0;
      if (!l_ghost_) field.ghost_depth (i_f,&gx,&gy,&gz);

      if (l_print_) {
        // Write field sums to output
        char buffer[256];
        snprintf(buffer,255,"field-%s-%s-%03d.data",
                 field.field_name(i_f).c_str(),
                 block->name().c_str(),block->state()->cycle());
        FILE * fp = fopen (buffer,"a");
        for (int iz=gz; iz<mz-gz; iz++) {
          for (int iy=gy; iy<my-gy; iy++) {
            for (int ix=gx; ix<mx-gx; ix++) {
              int i=ix + mx*(iy + my*iz);
              fprintf(fp,"%d %d %d %20.18f\n",ix,iy,iz,values[i]);
            }
          }
        }
        fclose(fp);
      }
      // write coarse field if needed
      if (l_coarse_) {

        if (l_print_) {
          // write coarse field
          char buffer[256];
          snprintf(buffer,255,"FIELD-%s-%s-%03d.data",
                   field.field_name(i_f).c_str(),
                   block->name().c_str(),block->state()->cycle());
          FILE * fp = fopen (buffer,"a");

          int mx,my,mz;
          field.coarse_dimensions (i_f,&mx,&my,&mz);
          cello_float * values = (cello_float*)field.coarse_values(i_f);
          for (int iz=gz; iz<mz-gz; iz++) {
            for (int iy=gy; iy<my-gy; iy++) {
              for (int ix=gx; ix<mx-gx; ix++) {
                int i=ix + mx*(iy + my*iz);
                fprintf(fp,"%d %d %d %20.18f\n",ix,iy,iz,values[i]);
              }
            }
          }
          fclose(fp);
        }
      }
    }
  }

  block->compute_done();
}

