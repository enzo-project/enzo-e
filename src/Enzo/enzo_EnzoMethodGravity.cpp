// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravity.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-11-07
/// @brief    Implements the EnzoMethodGravity class


#include "cello.hpp"
#include "enzo.hpp"

#include "enzo.decl.h"

// #define DEBUG_METHOD

// #define READ_ENZO_DENSITY
// #define PRINT_DENSITY_TOTAL
// #define USE_ENZO_DENSITY
// #define WRITE_ACCUM_DENSITY

// #define READ_ENZO_POTENTIAL
// #define USE_ENZO_POTENTIAL

// #define DEBUG_COPY_DENSITY
// #define DEBUG_COPY_B
// #define READ_PARTICLES

#define CK_TEMPLATES_ONLY
#include "enzo.def.h"
#undef CK_TEMPLATES_ONLY

#ifdef DEBUG_METHOD
#   define TRACE_METHOD(METHOD,BLOCK)					\
  CkPrintf ("%d %s:%d %s TRACE %s %p\n",CkMyPe(),__FILE__,__LINE__, \
	    BLOCK->name().c_str(),METHOD,this);			    \
  fflush(stdout);
#else
#   define TRACE_METHOD(METHOD,BLOCK) /*  */ 
#endif

//----------------------------------------------------------------------

EnzoMethodGravity::EnzoMethodGravity
(const FieldDescr * field_descr,
 int index_solver,
 double grav_const,
 bool accumulate)
  : Method(),
    index_solver_(index_solver),
    grav_const_(grav_const)
{
  const int id  = field_descr->field_id("density");
  const int idt = field_descr->field_id("density_total");
  const int ib  = field_descr->field_id("B");
  const int idensity = (idt != -1) ? idt : id;

  // Refresh adds density_total field faces and one layer of ghost
  // zones to "B" field


  const int ir = add_refresh(4,0,neighbor_leaf,sync_neighbor,
			     enzo_sync_id_method_gravity);

  refresh(ir)->add_field(field_descr->field_id("acceleration_x"));
  refresh(ir)->add_field(field_descr->field_id("acceleration_y"));
  refresh(ir)->add_field(field_descr->field_id("acceleration_z"));
  refresh(ir)->add_field(field_descr->field_id("density"));

  if (idt != -1 && accumulate) {
    // WARNING: accumulate cannot be used with AMR yet [170818]
    // Accumulate is used when particles are deposited into density_total.
    refresh(ir)->add_field_src_dst(idt,ib);
    refresh(ir)->set_accumulate(true);
  }
}

//----------------------------------------------------------------------

void EnzoMethodGravity::compute(Block * block) throw()
{
  TRACE_METHOD("compute()",block);

  // Initialize the linear system

  Field field = block->data()->field();

  int order=4;
  Matrix * A = new EnzoMatrixLaplace(order);
  const int ix = field.field_id ("potential");

  /// access problem-defining fields for eventual RHS and solution
  const int ib = field.field_id ("B");
#ifdef DEBUG_COPY_B  
  const int ib_copy = field.field_id ("B_copy");
#endif  
  const int id  = field.field_id("density");
  const int idt = field.field_id("density_total");
  const int idensity = (idt != -1) ? idt : id;
  
  // Solve the linear system
  int gx,gy,gz;
  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth(0,&gx,&gy,&gz);

  enzo_float * B = (enzo_float*) field.values (ib);
#ifdef DEBUG_COPY_B  
  enzo_float * B_copy = (enzo_float*) field.values (ib_copy);
#endif  
  enzo_float * D = (enzo_float*) field.values (idensity);

  TRACE_FIELD("B",B,1.0);

  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
	int i = ix + mx*(iy + my*iz);
#ifdef DEBUG_COPY_B  
	B_copy[i] = B[i];
#endif	
	D[i]  += B[i];
      }
    }
  }

#ifdef WRITE_ACCUM_DENSITY
  {
    char buffer[80];
    sprintf (buffer,"bc-enzop-%03d.data",block->cycle());
    printf ("DEBUG_ACCUM cycle=%d\n",block->cycle());
    FILE * fp = fopen(buffer,"w");
    field.ghost_depth(0,&gx,&gy,&gz);
    gx=gy=gz=1;
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
  	for (int ix=gx; ix<mx-gx; ix++) {
	  int i = ix + mx*(iy + my*iz);
	  fprintf (fp,"%d %d %d %20.16g\n",ix-gx,iy-gy,iz-gz,B_copy[i]);
	}
      }
    }
    fclose(fp);
  }
#endif

  TRACE_FIELD("density-total",D,1.0);

  EnzoPhysicsCosmology * cosmology = (EnzoPhysicsCosmology * )
    block->simulation()->problem()->physics("cosmology");

#ifdef READ_ENZO_DENSITY
  {
    enzo_float * dt = (enzo_float*) field.values("density_total");
    enzo_float * dt_enzo = (enzo_float*) field.values("density_total_enzo");
    // save Enzo-p's computed density_total first
    char buffer[80];
    sprintf (buffer,"dt-enzo-%03d.data",block->cycle());
    printf ("DEBUG_DENSITY_TOTAL cycle=%d\n",block->cycle());
    FILE * fp = fopen(buffer,"r");
    int gx,gy,gz;
    field.ghost_depth(0,&gx,&gy,&gz);
    gx=gy=gz=1;
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
  	for (int ix=gx; ix<mx-gx; ix++) {
  	  int i=ix+mx*(iy+my*iz);
	  int jx,jy,jz;
	  double value;
  	  fscanf (fp,"%d %d %d %lf\n",&jx,&jy,&jz,&value);
	  dt_enzo[i] = value;
#ifdef USE_ENZO_DENSITY_TOTAL
	  dt[i] = value;
#endif	  
  	}
      }
    }
    CkPrintf ("DEBUG_DENSITY dte[0] = %20.15g\n",dt_enzo[1+mx*(1+my*1)]);
    fclose(fp);
    fp=NULL;
  }
#endif

  if (block->is_leaf()) {
    if (cosmology) {

      enzo_float a=0.0,dadt=0.0;
      double time = block->time();
      double dt   = block->dt();
      cosmology-> compute_expansion_factor (&a,&dadt,time+0.5*dt);
#ifdef PRINT_DENSITY_TOTAL
      char buffer[80];
    sprintf (buffer,"dt-enzop-%03d.data",block->cycle());
    printf ("DEBUG_DENSITY_TOTAL cycle=%d\n",block->cycle());
    FILE * fp = fopen(buffer,"w");
#endif    
    int gx,gy,gz;
    field.ghost_depth(0,&gx,&gy,&gz);
#ifdef READ_ENZO_DENSITY    
    enzo_float * dt_diff = (enzo_float*) field.values("density_total_diff");
    enzo_float * dt_enzo = (enzo_float*) field.values("density_total_enzo");
#endif    
    gx=gy=gz=1;
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
  	for (int ix=gx; ix<mx-gx; ix++) {
	  int i = ix + mx*(iy + my*iz);
	  D[i]=-(D[i]-1.0);
#ifdef READ_ENZO_DENSITY    
	  dt_diff[i] = D[i] - dt_enzo[i];
#endif	  
	  B[i]  = D[i];
#ifdef PRINT_DENSITY_TOTAL
	  fprintf (fp,"%d %d %d %20.16g\n",ix-gx,iy-gy,iz-gz,D[i]);
#endif	    
	}
      }
    }
#ifdef PRINT_DENSITY_TOTAL
    fclose(fp);
#endif      
    } else {


      field.scale(ib, -4.0 * (cello::pi) * grav_const_, idensity);

    }

  } else {

    for (int i=0; i<mx*my*mz; i++) B[i] = 0.0;

  }
  
  //  TRACE_FIELD("density-shift",D,1.0);

  TRACE_FIELD("density-rhs",B,-1.0);

  Solver * solver = block->simulation()->problem()->solver(index_solver_);
  
  // May exit before solve is done...
  solver->set_callback (CkIndex_EnzoBlock::r_method_gravity_continue());


  solver->apply (A, ix, ib, block);

}

//----------------------------------------------------------------------

void EnzoBlock::r_method_gravity_continue()
{

  TRACE_METHOD("r_method_gravity_end()",this);

  // So do refresh with barrier synch (note barrier instead of
  // neighbor synchronization otherwise will conflict with Method
  // refresh ("Charm++ fatal error: mis-matched client callbacks in
  // reduction messages")

  Refresh refresh (4,0,neighbor_leaf, sync_barrier,
		   enzo_sync_id_method_gravity_continue);
  
  refresh.set_active(is_leaf());
  refresh.add_all_fields();

  refresh_enter(CkIndex_EnzoBlock::r_method_gravity_end(NULL),&refresh);

}

//----------------------------------------------------------------------

void EnzoBlock::r_method_gravity_end(CkReductionMsg * msg)
{
  TRACE_METHOD("r_method_gravity_end()",this);
  
  delete msg;
  
  Field field = data()->field();
  int gx,gy,gz;
  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth(0,&gx,&gy,&gz);
  enzo_float * potential = (enzo_float*) field.values ("potential");
  TRACE_FIELD("potential",potential,-1.0);
  
  EnzoPhysicsCosmology * cosmology = (EnzoPhysicsCosmology * )
    simulation()->problem()->physics("cosmology");

  if (cosmology) {

    enzo_float a=0.0,dadt=0.0,dt=0.0;
    dt   = this->timestep();
    double time = this->time();
    cosmology-> compute_expansion_factor (&a,&dadt,time+0.5*dt);
    //    cosmology-> compute_expansion_factor (&a,&dadt,time);

    for (int i=0; i<mx*my*mz; i++) potential[i] /= a;
  }
  
#ifdef READ_ENZO_POTENTIAL  
  {
    enzo_float * po = (enzo_float*) field.values("potential");
    enzo_float * po_enzo = (enzo_float*) field.values("potential_enzo");
    enzo_float * po_diff = (enzo_float*) field.values("potential_diff");
    // save Enzo-p's computed potential first
    char buffer[80];
    sprintf (buffer,"po-enzo-%03d.data",cycle_);
    printf ("DEBUG_POTENTIAL cycle=%d\n",cycle_);
    FILE * fp = fopen(buffer,"r");
    int gx=1,gy=1,gz=1;
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
  	for (int ix=gx; ix<mx-gx; ix++) {
  	  int i=ix+mx*(iy+my*iz);
	  int jx,jy,jz;
	  double value;
  	  fscanf (fp,"%d %d %d %lf\n",&jx,&jy,&jz,&value);
	  po_enzo[i] = value;
	  po_diff[i] = po[i] - po_enzo[i];
	  //	  printf ("%20.15g %20.15g %20.15g\n",po_diff[i],po[i],value);
#ifdef USE_ENZO_POTENTIAL
	  po[i] = value;
#endif	  
  	}
      }
    }
    fclose(fp);
    fp=NULL;
  }
#endif

#ifdef READ_PARTICLES
  if ((this->cycle() % 20) == 0) {
    Particle particle = data()->particle();

    int it = particle.type_index("enzo");
    int ia_x = particle.attribute_index(it,"x");
    int ia_y = particle.attribute_index(it,"y");
    int ia_z = particle.attribute_index(it,"z");
    int ia_vx = particle.attribute_index(it,"vx");
    int ia_vy = particle.attribute_index(it,"vy");
    int ia_vz = particle.attribute_index(it,"vz");
    int ia_ax = particle.attribute_index(it,"ax");
    int ia_ay = particle.attribute_index(it,"ay");
    int ia_az = particle.attribute_index(it,"az");
    int nx,ny,nz;
    field.size (&nx,&ny,&nz);
    int mb = particle.batch_size();
    char buffer[80];
    sprintf (buffer,"particles-%03d.data",this->cycle());
    CkPrintf ("file = %s\n",buffer);
    int np = nx*ny*nz;
    FILE * fp = fopen(buffer,"r");
    if (particle.num_particles(it) == 0) {
      particle.insert_particles(it,np);
    }
    for (int i=0; i<np; i++) {
      int ib,ip;
      particle.index(i,&ib,&ip);
      enzo_float * x = (enzo_float *) particle.attribute_array(it,ia_x,ib);
      enzo_float * y = (enzo_float *) particle.attribute_array(it,ia_y,ib);
      enzo_float * z = (enzo_float *) particle.attribute_array(it,ia_z,ib);
      enzo_float * vx = (enzo_float *) particle.attribute_array(it,ia_vx,ib);
      enzo_float * vy = (enzo_float *) particle.attribute_array(it,ia_vy,ib);
      enzo_float * vz = (enzo_float *) particle.attribute_array(it,ia_vz,ib);
      enzo_float * ax = (enzo_float *) particle.attribute_array(it,ia_ax,ib);
      enzo_float * ay = (enzo_float *) particle.attribute_array(it,ia_ay,ib);
      enzo_float * az = (enzo_float *) particle.attribute_array(it,ia_az,ib);
      int j=i%mb;
      fscanf (fp,"%f %f %f %f %f %f %f %f %f\n",
	      &x[j],&y[j],&z[j],&vx[j],&vy[j],&vz[j],&ax[j],&ay[j],&az[j]);
    }
    // for (int i=0; i<np; i++) {
    //   int ib,ip;
    //   particle.index(i,&ib,&ip);
    //   enzo_float * x = (enzo_float *) particle.attribute_array(it,ia_x,ib);
    //   enzo_float * y = (enzo_float *) particle.attribute_array(it,ia_y,ib);
    //   enzo_float * z = (enzo_float *) particle.attribute_array(it,ia_z,ib);
    //   CkPrintf ("B %d %g %g %g\n",i,x[i%mb],y[i%mb],z[i%mb]);
    // }
      
    fclose (fp);
  }
  
#endif

  /// compute acceleration fields from potential
  int order;
  EnzoComputeAcceleration compute_acceleration(data()->field().field_descr(),
					       rank(), order=4);

  compute_acceleration.compute(this);

  // Clear "B" and "density_total" fields for next call
  // Note density_total may not be defined

  enzo_float * B =         (enzo_float*) field.values("B");

#ifdef DEBUG_COPY_DENSITY  
  enzo_float * B_temp =    (enzo_float*) field.values("B_temp");
  for (int i=0; i<mx*my*mz; i++) {
    B_temp[i] = B[i];
  }
#endif  

  for (int i=0; i<mx*my*mz; i++) {
    B[i] = 0.0;
  }

  enzo_float * de_t =      (enzo_float*) field.values("density_total");

  if (de_t) {
#ifdef DEBUG_COPY_DENSITY  
    enzo_float * de_t_temp = (enzo_float*) field.values("density_total_temp");
    CkPrintf ("DEBUG_DENSITY dt[0] = %20.15g\n",de_t[1+mx*(1+my*1)]);
    for (int i=0; i<mx*my*mz; i++) {
      de_t_temp[i] = de_t[i];
    }
#endif  
    for (int i=0; i<mx*my*mz; i++) {
      de_t[i] = 0.0;
    }
  }
  if (de_t) {
#ifdef DEBUG_COPY_DENSITY  
    enzo_float * po_temp = (enzo_float*) field.values("potential_temp");
    for (int i=0; i<mx*my*mz; i++) {
      po_temp[i] = potential[i];
    }
#endif  
    for (int i=0; i<mx*my*mz; i++) {
      potential[i] = 0.0;
    }
  }

  // wait for all Blocks before continuing
  compute_done();
}

//----------------------------------------------------------------------

double EnzoMethodGravity::timestep (Block * block) const throw()
{
  return timestep_(block);
}

//----------------------------------------------------------------------

double EnzoMethodGravity::timestep_ (Block * block) const throw()
{
  Field field = block->data()->field();

  int nx,ny,nz;
  int mx,my,mz;
  int gx,gy,gz;
  field.size         (&nx,&ny,&nz);
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth(0,&gx,&gy,&gz);

  enzo_float * ax = (enzo_float*) field.values ("acceleration_x");
  enzo_float * ay = (enzo_float*) field.values ("acceleration_y");
  enzo_float * az = (enzo_float*) field.values ("acceleration_z");

  enzo_float dt = std::numeric_limits<enzo_float>::max();

  double hx,hy,hz;
  block->cell_width(&hx,&hy,&hz);
  
  if (ax) {
    for (int ix=gx; ix<nx+gx; ix++) {
      for (int iy=gy; iy<ny+gy; iy++) {
	for (int iz=gz; iz<nz+gz; iz++) {
	  int i=ix + mx*(iy + iz*my);
	  dt = std::min(enzo_float(dt),enzo_float(sqrt(hx/(fabs(ax[i]+1e-20)))));
	}
      }
    }
  }
  if (ay) {
    for (int ix=gx; ix<nx+gx; ix++) {
      for (int iy=gy; iy<ny+gy; iy++) {
	for (int iz=gz; iz<nz+gz; iz++) {
	  int i=ix + mx*(iy + iz*my);
	  dt = std::min(enzo_float(dt),enzo_float(sqrt(hy/(fabs(ay[i]+1e-20)))));
	}
      }
    }
  }
  if (az) {
    for (int ix=gx; ix<nx+gx; ix++) {
      for (int iy=gy; iy<ny+gy; iy++) {
	for (int iz=gz; iz<nz+gz; iz++) {
	  int i=ix + mx*(iy + iz*my);
	  dt = std::min(enzo_float(dt),enzo_float(sqrt(hz/(fabs(az[i]+1e-20)))));
	}
      }
    }
  }

  return 0.5*dt;
}
