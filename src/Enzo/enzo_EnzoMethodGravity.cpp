// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravity.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-11-07
/// @brief    Implements the EnzoMethodGravity class


#include "cello.hpp"
#include "enzo.hpp"

#include "enzo.decl.h"

// #define DEBUG_FIELD

#ifdef DEBUG_FIELD

#   define TRACE_FIELD_GM(NAME,FIELD,SCALE,gx,gy,gz,mx,my,mz,PLOT)	\
  {									\
    int gx=3,gy=3,gz=3;							\
    double sum_all=0.0;							\
    double sum_real=0.0;						\
    double sum_abs=0.0;							\
    double sum_mean=0.0;						\
    double sum_var=0.0;							\
    double sum2_all=0.0;						\
    double sum2_real=0.0;						\
    for (int iz=0; iz<mz; iz++) {					\
      for (int iy=0; iy<my; iy++) {					\
	for (int ix=0; ix<mx; ix++) {					\
	  int i = ix + mx*(iy + my*iz);					\
	  double value=SCALE*FIELD[i];					\
	  sum_all+=value;						\
	  sum2_all += value * value;					\
	}								\
      }									\
    }									\
    for (int iz=gz; iz<mz-gz; iz++) {					\
      for (int iy=gy; iy<my-gy; iy++) {					\
	for (int ix=gx; ix<mx-gx; ix++) {				\
	  int i = ix + mx*(iy + my*iz);					\
	  double value=SCALE*FIELD[i];					\
	  sum_real+=value;						\
	  sum_abs+=std::abs(value);					\
	  sum2_real+=value*value;					\
	}								\
      }									\
    }									\
    double mean=sum_real/((mx-2*gx)*(my-2*gy)*(mz-2*gz));		\
    for (int iz=gz; iz<mz-gz; iz++) {					\
      for (int iy=gy; iy<my-gy; iy++) {					\
	for (int ix=gx; ix<mx-gx; ix++) {				\
	  int i = ix + mx*(iy + my*iz);					\
	  double value=SCALE*FIELD[i];					\
	  sum_mean +=std::abs(value-mean);				\
	  sum_var  +=(value-mean)*(value-mean);				\
	}								\
      }									\
    }									\
    CkPrintf ("DEBUG_FIELD (%g) [%g] | (%g : %g %g %g) %s:%d %s\n",	\
	      sum_real,sum_all,					\
	      SCALE*(FIELD[gx+mx*(gy+my*gz)]),				\
	      SCALE*(FIELD[(gx+1)+mx*(gy+my*gz)]),			\
	      SCALE*(FIELD[gx+mx*((gy+1)+my*gz)]),			\
	      SCALE*(FIELD[gx+mx*(gy+my*(gz+1))]),			\
	      __FILE__,__LINE__,NAME);					\
    CkPrintf ("DEBUG_FIELD %s %20.18g %20.18g %20.18g\n",  NAME,sum_abs,sum_mean,sum_var); \
    fflush(stdout);							\
    char filename[80];						\
    sprintf (filename,"renzo-p-%s.png",NAME);				\
    png_array (filename,(float*)(FIELD),gx,gy,gz,mx,my,mz,__FILE__,__LINE__,2,16,16,SCALE); \
    sprintf (filename,"enzo-p-%s.png",NAME);				\
    png_array (filename,(float*)(FIELD),0,0,0,mx,my,mz,__FILE__,__LINE__,2,16,16,SCALE); \
  }
#   define TRACE_FIELD_(NAME,FIELD,SCALE) TRACE_FIELD_GM(NAME,FIELD,SCALE,gx_,gy_,gz_,mx_,my_,mz_,false)
#   define TRACE_FIELD(NAME,FIELD,SCALE)  TRACE_FIELD_GM(NAME,FIELD,SCALE,gx,gy,gz,mx,my,mz,false)
#   define PRINT_FIELD_(NAME,FIELD,SCALE) TRACE_FIELD_GM(NAME,FIELD,SCALE,gx_,gy_,gz_,mx_,my_,mz_,true)
#   define PRINT_FIELD(NAME,FIELD,SCALE)  TRACE_FIELD_GM(NAME,FIELD,SCALE,gx,gy,gz,mx,my,mz,true)

#define TRACE_PARTICLE(ATTRIBUTE,particle,TYPE_NAME,ATTR_NAME)		\
  {									\
    int it = particle.type_index(TYPE_NAME);				\
    int ia = particle.attribute_index(it,ATTR_NAME);			\
    int nb = particle.num_batches(it);					\
    double sum_abs = 0.0;						\
    for (int ib=0; ib<nb; ib++) {					\
      int np = particle.num_particles(it,ib);				\
      enzo_float * array = (enzo_float *) particle.attribute_array(it,ia,ib); \
      for (int ip=0; ip<np; ip++) {					\
	sum_abs += std::abs(array[ip]);					\
      }									\
    }									\
    CkPrintf ("DEBUG_PARTICLE %s %20.18g\n",  ATTRIBUTE,sum_abs);		\
  }
#else
#   define TRACE_FIELD_GM(NAME,FIELD,SCALE,gx,gy,gz,mx,my,mz,PLOT)	\
  /* ... */
#   define TRACE_FIELD_(NAME,FIELD,SCALE) TRACE_FIELD_GM(NAME,FIELD,SCALE,gx_,gy_,gz_,mx_,my_,mz_,false) \	/* ... */
#   define TRACE_FIELD(NAME,FIELD,SCALE)  TRACE_FIELD_GM(NAME,FIELD,SCALE,gx,gy,gz,mx,my,mz,false) \
  /* ... */
#   define PRINT_FIELD_(NAME,FIELD,SCALE) TRACE_FIELD_GM(NAME,FIELD,SCALE,gx_,gy_,gz_,mx_,my_,mz_,true) \
  /* ... */
#   define PRINT_FIELD(NAME,FIELD,SCALE)  TRACE_FIELD_GM(NAME,FIELD,SCALE,gx,gy,gz,mx,my,mz,true) \
  /* ... */
#define TRACE_PARTICLE(ATTRIBUTE,particle,TYPE_NAME,ATTR_NAME)	\
  /* ... */
#endif

// #define DEBUG_FIELD_FACE
// #define DEBUG_COPY_DENSITY
#define DEBUG_COPY_POTENTIAL
// #define DEBUG_COPY_B

// #define DEBUG_METHOD

// #define READ_ENZO_DENSITY
// #define PRINT_DENSITY_TOTAL
// #define USE_ENZO_DENSITY
// #define WRITE_ACCUM_DENSITY

// #define READ_ENZO_POTENTIAL
// #define USE_ENZO_POTENTIAL

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
(int index_solver,
 double grav_const,
 int order,
 bool accumulate)
  : Method(),
    index_solver_(index_solver),
    grav_const_(grav_const),
    order_(order)
{
  FieldDescr * field_descr = cello::field_descr();

  // Change this if fields used in this routine change
  // declare required fields
  this->required_fields_ = std::vector<std::string>
                           {"density","density_total","B","potential",
                            "acceleration_x","acceleration_y","acceleration_z"};
#ifdef DEBUG_FIELD_FACE
  this->required_fields_.insert(this->required_fields_.end(),
                                {"debug_1","debug_2"});
#endif
#ifdef DEBUG_COPY_B
  this->required_fields_.push_back("B_copy");
#endif
#ifdef DEBUG_COPY_POTENTIAL
  this->required_fields_.push_back("potential_copy");
#endif
#ifdef DEBUG_COPY_DENSITY
  this->required_fields_.push_back("density_total_copy");
#endif


  if (accumulate){
    this->required_fields_.insert(this->required_fields_.end(),
                                  {"density_particle","density_particle_accumulate"});
  }

  // now define fields if they do not exist
  this->define_fields();

  const int id  = field_descr->field_id("density");
  const int idt = field_descr->field_id("density_total");
  const int ib  = field_descr->field_id("B");
  const int iax = field_descr->field_id("acceleration_x");
  const int iay = field_descr->field_id("acceleration_y");
  const int iaz = field_descr->field_id("acceleration_z");

  // Refresh adds density_total field faces and one layer of ghost
  // zones to "B" field

#ifdef DEBUG_FIELD_FACE
  int idebug1 = field_descr->field_id("debug_1");
  int idebug2 = field_descr->field_id("debug_2");
#endif

  const int ir = add_refresh(4,0,neighbor_leaf,sync_neighbor,
			     enzo_sync_id_method_gravity);

  refresh(ir)->add_field(iax);
  refresh(ir)->add_field(iay);
  refresh(ir)->add_field(iaz);
  refresh(ir)->add_field(id);

  // Accumulate is used when particles are deposited into density_total

  if (accumulate) {

    refresh(ir)->set_accumulate(true);

    const int idp = field_descr->field_id("density_particle");
    const int idpa = field_descr->field_id("density_particle_accumulate");

    refresh(ir)->add_field_src_dst(idp,idpa);
    refresh(ir)->add_field_src_dst(idt,ib);

#ifdef DEBUG_FIELD_FACE
    refresh(ir)->add_field_src_dst(idebug1,idebug2);
#endif
  }
}

//----------------------------------------------------------------------

void EnzoMethodGravity::compute(Block * block) throw()
{
  TRACE_METHOD("compute()",block);

  // Initialize the linear system

  Field field = block->data()->field();

  /// access problem-defining fields for eventual RHS and solution
  const int ib  = field.field_id ("B");
  const int id  = field.field_id("density");
  const int idt = field.field_id("density_total");
  const int idensity = (idt != -1) ? idt : id;

  // Solve the linear system
  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth(0,&gx,&gy,&gz);

  const int m = mx*my*mz;

  enzo_float * B = (enzo_float*) field.values (ib);
#ifdef DEBUG_COPY_B
  const int ib_copy = field.field_id ("B_copy");
  enzo_float * B_copy = (enzo_float*) field.values (ib_copy);
#endif
  enzo_float * D = (enzo_float*) field.values (idensity);

  TRACE_FIELD("B",B,1.0);

  for (int i=0; i<m; i++) D[i] += B[i];

  // Add density_particle values to density_particle_accumulate ghosts

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

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

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

    fclose(fp);
    fp=NULL;
  }
#endif

  if (block->is_leaf()) {
    if (cosmology) {

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

#ifdef DEBUG_COPY_B
  for (int i=0; i<m; i++) B_copy[i] = B[i];
#endif

  Solver * solver = enzo::problem()->solver(index_solver_);

  // May exit before solve is done...
  solver->set_callback (CkIndex_EnzoBlock::r_method_gravity_continue());

  const int ix = field.field_id ("potential");

  std::shared_ptr<Matrix> A (std::make_shared<EnzoMatrixLaplace>(order_));

  solver->set_field_x(ix);
  solver->set_field_b(ib);

  solver->apply (A, block);

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
  refresh.add_field(data()->field().field_id("potential"));

  refresh_enter(CkIndex_EnzoBlock::r_method_gravity_end(),&refresh);

}

//----------------------------------------------------------------------

void EnzoBlock::r_method_gravity_end()
{
  TRACE_METHOD("r_method_gravity_end()",this);

  EnzoMethodGravity * method = static_cast<EnzoMethodGravity*> (this->method());
  method->compute_accelerations(this);
}

void EnzoMethodGravity::compute_accelerations (EnzoBlock * enzo_block) throw()
{

  Field field = enzo_block->data()->field();
  int gx,gy,gz;
  int mx,my,mz;
  field.ghost_depth(0,&gx,&gy,&gz);
  field.dimensions (0,&mx,&my,&mz);
  const int m = mx*my*mz;
  enzo_float * potential = (enzo_float*) field.values ("potential");
  TRACE_FIELD("potential",potential,-1.0);

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  if (cosmology) {

    enzo_float cosmo_a = 1.0;
    enzo_float cosmo_dadt = 0.0;
    double dt   = enzo_block->timestep();
    double time = enzo_block->time();
    cosmology-> compute_expansion_factor (&cosmo_a,&cosmo_dadt,time+0.5*dt);
    //    cosmology-> compute_expansion_factor (&a,&dadt,time);

    for (int i=0; i<m; i++) potential[i] /= cosmo_a;
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
  if ((enzo_block->cycle() % 20) == 0) {
    Particle particle = enzo_block->data()->particle();

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
    sprintf (buffer,"particles-%03d.data",enzo_block->cycle());
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

  EnzoComputeAcceleration compute_acceleration(cello::rank(), order_);

  compute_acceleration.compute(enzo_block);

  // Clear "B" and "density_total" fields for next call
  // Note density_total may not be defined

  enzo_float * B =         (enzo_float*) field.values("B");

#ifdef DEBUG_COPY_DENSITY
  enzo_float * B_copy =    (enzo_float*) field.values("B_copy");
  for (int i=0; i<m; i++) B_copy[i] = B[i];
#endif

  for (int i=0; i<m; i++) B[i] = 0.0;

  enzo_float * de_t =      (enzo_float*) field.values("density_total");

  if (de_t != NULL) {
#ifdef DEBUG_COPY_DENSITY
    enzo_float * de_t_copy = (enzo_float*) field.values("density_total_copy");
    for (int i=0; i<m; i++) de_t_copy[i] = de_t[i];
#endif
    for (int i=0; i<m; i++) de_t[i] = 0.0;
  }
  if (potential) {
#ifdef DEBUG_COPY_POTENTIAL
    enzo_float * po_copy = (enzo_float*) field.values("potential_copy");
    if (po_copy != NULL) {
      for (int i=0; i<m; i++) po_copy[i] = potential[i];
    }
#endif
    for (int i=0; i<m; i++) potential[i] = 0.0;
  }

  // wait for all Blocks before continuing
  enzo_block->compute_done();
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

  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth(0,&gx,&gy,&gz);

  enzo_float * ax = (enzo_float*) field.values ("acceleration_x");
  enzo_float * ay = (enzo_float*) field.values ("acceleration_y");
  enzo_float * az = (enzo_float*) field.values ("acceleration_z");

  enzo_float dt = std::numeric_limits<enzo_float>::max();

  double hx,hy,hz;
  block->cell_width(&hx,&hy,&hz);

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  if (cosmology) {
    const int rank = cello::rank();
    enzo_float cosmo_a = 1.0;
    enzo_float cosmo_dadt = 0.0;
    double dt   = block->dt();
    double time = block->time();
    cosmology-> compute_expansion_factor (&cosmo_a,&cosmo_dadt,time+0.5*dt);
    if (rank >= 1) hx*=cosmo_a;
    if (rank >= 2) hy*=cosmo_a;
    if (rank >= 3) hz*=cosmo_a;
  }

  if (ax) {
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
	for (int ix=gx; ix<mx-gx; ix++) {
	  int i=ix + mx*(iy + iz*my);
	  dt = std::min(enzo_float(dt),enzo_float(sqrt(hx/(fabs(ax[i]+1e-20)))));
	}
      }
    }
  }
  if (ay) {
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
	for (int ix=gx; ix<mx-gx; ix++) {
	  int i=ix + mx*(iy + iz*my);
	  dt = std::min(enzo_float(dt),enzo_float(sqrt(hy/(fabs(ay[i]+1e-20)))));
	}
      }
    }
  }
  if (az) {
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
	for (int ix=gx; ix<mx-gx; ix++) {
	  int i=ix + mx*(iy + iz*my);
	  dt = std::min(enzo_float(dt),enzo_float(sqrt(hz/(fabs(az[i]+1e-20)))));
	}
      }
    }
  }

  return 0.5*dt;
}
