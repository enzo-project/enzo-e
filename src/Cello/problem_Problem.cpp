// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Problem.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-03
/// @brief    Implementation of the Problem container class

#include "problem.hpp"

//----------------------------------------------------------------------

Problem::Problem() throw()
  : is_periodic_(true),
    stopping_(NULL),
    prolong_(NULL),
    restrict_(NULL),
    index_refine_(0),
    index_output_(0),
    index_boundary_(0)
{
  
}

//----------------------------------------------------------------------

Problem::~Problem() throw()
{
  deallocate_();
}

//----------------------------------------------------------------------

void Problem::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  PUP::able::pup(p);

  bool pk = p.isPacking();
  bool up = p.isUnpacking();

  int n;

  if (pk) n=boundary_list_.size();
  p | n;
  if (up) boundary_list_.resize(n);
  for (int i=0; i<n; i++) {
    p | boundary_list_[i]; // PUP::able
  }

  p | is_periodic_;

  if (pk) n=initial_list_.size();
  p | n;
  if (up) initial_list_.resize(n);
  for (int i=0; i<n; i++) {
    p | initial_list_[i]; // PUP::able
  }

  if (pk) n=refine_list_.size();
  p | n;
  if (up) refine_list_.resize(n);
  for (int i=0; i<n; i++) {
    p | refine_list_[i]; // PUP::able
  }

  p | stopping_;

  // if (pk) n=timestep_list_.size();
  // p | n;
  // if (up) timestep_list_.resize(n);
  // for (int i=0; i<n; i++) {
  //   p | timestep_list_[i]; // PUP::able
  // }

  if (pk) n=method_list_.size();
  p | n;
  if (up) method_list_.resize(n);
  for (int i=0; i<n; i++) {
    p | method_list_[i]; // PUP::able
  }

  if (pk) n=solver_list_.size();
  p | n;
  if (up) solver_list_.resize(n);
  for (int i=0; i<n; i++) {
    p | solver_list_[i]; // PUP::able
  }

  if (pk) n=output_list_.size();
  p | n;
  if (up) output_list_.resize(n);
  for (int i=0; i<n; i++) {
    p | output_list_[i]; // PUP::able
  }

  p | prolong_; // PUP::able
  p | restrict_; // PUP::able

  p | index_refine_;
  p | index_output_;
  p | index_boundary_;
}

//----------------------------------------------------------------------

void Problem::initialize_boundary(Config * config, 
				  Parameters * parameters) throw()
{
  for (int index=0; index < config->num_boundary; index++) {

    std::string type = config->boundary_type[index];

    Boundary * boundary = create_boundary_(type,index,config,parameters);

    ASSERT1("Problem::initialize_boundary",
	  "Boundary type %s not recognized",
	  type.c_str(),  boundary != NULL);

    boundary_list_.push_back(boundary);
  }

}

//----------------------------------------------------------------------

void Problem::initialize_initial(Config * config,
				 Parameters * parameters,
				 const FieldDescr * field_descr) throw()
{

  for (int index=0; index < config->num_initial; index++) {

    std::string type = config->initial_list[index];

    Initial * initial = create_initial_
      (type,index,config,parameters,field_descr);

    ASSERT1("Problem::initialize_initial",
	    "Initial type %s not recognized",
	    config->initial_list[index].c_str(),
	    initial != NULL);

    initial_list_.push_back( initial );
  }
}

//----------------------------------------------------------------------

void Problem::initialize_refine(Config * config,
				Parameters * parameters,
				const FieldDescr * field_descr) throw()
{
  for (int i=0; i<config->num_adapt; i++) {

    std::string name = config->adapt_type[i];

    Refine * refine = create_refine_ 
      (name,config,parameters,field_descr,i);

    if (refine) {
      refine_list_.push_back( refine );
    } else {
      ERROR1("Problem::initialize_refine",
	     "Cannot create Refine type %s",name.c_str());
    }
  }
}

//----------------------------------------------------------------------

void Problem::initialize_stopping(Config * config) throw()
{
  stopping_ = create_stopping_("default",config);

  ASSERT("Problem::initialize_stopping",
	  "Stopping object not successfully created",
	  stopping_ != NULL);
}

//----------------------------------------------------------------------

void Problem::initialize_prolong(Config * config) throw()
{
  prolong_ = create_prolong_(config->field_prolong,config);

  ASSERT1("Problem::initialize_prolong",
	  "Prolong type %s not recognized",
	  config->field_prolong.c_str(),
	  prolong_ != NULL);
}

//----------------------------------------------------------------------

void Problem::initialize_restrict(Config * config) throw()
{
  restrict_ = create_restrict_(config->field_restrict,config);

  ASSERT1("Problem::initialize_restrict",
	  "Restrict type %s not recognized",
	  config->field_restrict.c_str(),
	  restrict_ != NULL);
}

//----------------------------------------------------------------------

void Problem::initialize_output
(Config * config,
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr,
 const Factory * factory) throw()
{

  for (int index=0; index < config->num_output; index++) {

    std::string type       = config->output_type[index];

    Output * output = create_output_ 
      (type,index, config,field_descr,particle_descr,factory);

    if (output == NULL) {
      ERROR2("Problem::initialize_output",
	     "Unknown parameter type Output:%s:type = %s",
	     config->output_list[index].c_str(),type.c_str());
    }

    if (config->output_name[index].size() > 0) {
      std::string file_name = config->output_name[index][0];

      std::vector<std::string> file_args;

      for (size_t i=1; i<config->output_name[index].size(); i++) {
	file_args.push_back(config->output_name[index][i]);
      }

      output->set_filename (file_name,file_args);
    }

    if (config->output_dir[index].size() > 0) {
      std::string dir_name = config->output_dir[index][0];

      std::vector<std::string> dir_args;

      for (size_t i=1; i<config->output_dir[index].size(); i++) {
	dir_args.push_back(config->output_dir[index][i]);
      }

      output->set_dir (dir_name,dir_args);
    }

    //--------------------------------------------------
    // field_list
    //--------------------------------------------------

    const int num_fields = config->output_field_list[index].size();

    if (num_fields == 0) {
      // field_list not initialized: assuming none
      ItIndexList * it_field = new ItIndexList;
      output->set_it_field_index(it_field);
    } else {
      // if any fields are "*", default is all fields
      bool all_fields = false;
      for (int i=0; i<num_fields; i++) {
	if (config->output_field_list[index][i] == "*") {
	  all_fields = true;
	  int field_count = field_descr->field_count();
	  ItIndexRange * it_field = new ItIndexRange(field_count);
	  output->set_it_field_index(it_field);
	  break;
	}
      } 
      // if no fields are "*", create field index iterator
      if (! all_fields) {
	ItIndexList * it_field = new ItIndexList;
	for (int i=0; i<num_fields; i++) {
	  std::string field_name = config->output_field_list[index][i];
	  int index_field = field_descr->field_id(field_name);
	  it_field->append(index_field);
	}
	output->set_it_field_index(it_field);
      }
    }

    //--------------------------------------------------
    // particle_list
    //--------------------------------------------------

    const int num_particles = config->output_particle_list[index].size();

    if (num_particles == 0) {
      // particle_list not initialized: assuming none
      ItIndexList * it_particle = new ItIndexList;
      output->set_it_particle_index(it_particle);
    } else {
      // if any particles are "*", default is all particles
      bool all_particles = false;
      for (int i=0; i<num_particles; i++) {
	if (config->output_particle_list[index][i] == "*") {
	  all_particles = true;
	  int particle_count = particle_descr->num_types();
	  ItIndexRange * it_particle = new ItIndexRange(particle_count);
	  output->set_it_particle_index(it_particle);
	  break;
	}
      } 
      // if no particles are "*", create particle index iterator
      if (! all_particles) {
	ItIndexList * it_particle = new ItIndexList;
	for (int i=0; i<num_particles; i++) {
	  std::string particle_name = config->output_particle_list[index][i];
	  int index_particle = particle_descr->type_index(particle_name);
	  it_particle->append(index_particle);
	}
	output->set_it_particle_index(it_particle);
      }
    }


    //--------------------------------------------------
    // Scheduling parameters
    //--------------------------------------------------

    int index_schedule = config->output_schedule_index[index];

    output->set_schedule
      (Schedule::create( config->schedule_var[index_schedule],
			 config->schedule_type[index_schedule],
			 config->schedule_start[index_schedule],
			 config->schedule_stop[index_schedule],
			 config->schedule_step[index_schedule],
			 config->schedule_list[index_schedule]));


    //--------------------------------------------------
    // Image parameters
    //--------------------------------------------------

    OutputImage * output_image = dynamic_cast<OutputImage *> (output);

    if (output_image != NULL) {

      // COLORMAP

      int n = config->output_colormap[index].size() / 3;

      if (n > 0) {
	double * r = new double [n];
	double * g = new double [n];
	double * b = new double [n];

	for (int i=0; i<n; i++) {

	  int ir=3*i+0;
	  int ig=3*i+1;
	  int ib=3*i+2;

	  r[i] = config->output_colormap[index][ir];
	  g[i] = config->output_colormap[index][ig];
	  b[i] = config->output_colormap[index][ib];

	}

	output_image->set_colormap(n,r,g,b);

	delete [] r;
	delete [] g;
	delete [] b;

      }

    }

    // Add the initialized Output object to the Simulation's list of
    // output objects

    output_list_.push_back(output); 

  } // (for index)

}

//----------------------------------------------------------------------

void Problem::initialize_method
(
 Config * config,
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr
 ) throw()
{
  const size_t num_method = config->method_list.size();

  Method::courant_global = config->method_courant_global;
  
  for (size_t index_method=0; index_method < num_method ; index_method++) {

    std::string name = config->method_list[index_method];

    Method * method = create_method_(name, config, index_method, 
				     field_descr, particle_descr);

    if (method) {

      method_list_.push_back(method); 

      int index_schedule = config->method_schedule_index[index_method];

      if (index_schedule != -1) {
	method->set_schedule
	  (Schedule::create( config->schedule_var[index_schedule],
			     config->schedule_type[index_schedule],
			     config->schedule_start[index_schedule],
			     config->schedule_stop[index_schedule],
			     config->schedule_step[index_schedule],
			     config->schedule_list[index_schedule]));
      }

    } else {
      ERROR1("Problem::initialize_method",
	     "Unknown Method %s",name.c_str());
    }
  }
}

//----------------------------------------------------------------------

void Problem::initialize_solver
(
 Config * config,
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr
 ) throw()
{
  const size_t num_solver = config->solver_list.size();

  for (size_t index_solver=0; index_solver < num_solver ; index_solver++) {

    std::string type = config->solver_type[index_solver];

    Solver * solver = create_solver_(type, config, index_solver, 
				     field_descr, particle_descr);

    if (solver) {

      solver_list_.push_back(solver); 

    } else {
      ERROR1("Problem::initialize_method",
	     "Unknown Method %s",type.c_str());
    }

  }
}

//----------------------------------------------------------------------

Solver * Problem::solver(size_t i) const throw()
{
  ASSERT2("Problem::solver",
	  "Solver id %d out of range [0,%d)",
	  i < solver_list_.size(),
	  i,solver_list_.size());
  return solver_list_.at(i);
}

//======================================================================

void Problem::deallocate_() throw()
{
  for (size_t i=0; i<boundary_list_.size(); i++) {
    delete boundary_list_[i];    boundary_list_[i] = 0;
  }

  for (size_t i=0; i<initial_list_.size(); i++) {
    delete initial_list_[i];    initial_list_[i] = 0;
  }
  for (size_t i=0; i<refine_list_.size(); i++) {
    delete refine_list_[i];    refine_list_[i] = 0;
  }
  delete stopping_;      stopping_ = 0;
  // for (size_t i=0; i<timestep_list_.size(); i++) {
  //   delete timestep_list_[i];     timestep_list_[i] = 0;
  // }
  for (size_t i=0; i<output_list_.size(); i++) {
    delete output_list_[i];    output_list_[i] = 0;
  }
  for (size_t i=0; i<solver_list_.size(); i++) {
    delete solver_list_[i];    solver_list_[i] = 0;
  }
  for (size_t i=0; i<method_list_.size(); i++) {
    delete method_list_[i];    method_list_[i] = 0;
  }
}

//----------------------------------------------------------------------

Boundary * Problem::create_boundary_
(
 std::string type,
 int index,
 Config * config,
 Parameters * parameters
 ) throw ()
{
  if (type != "periodic") is_periodic_ = false;

  if (type == "inflow") {

    std::string param_str = 
      "Boundary:" + config->boundary_list[index] + ":value";

    int param_type = parameters->type(param_str);

    if (! (param_type == parameter_list ||
	   param_type == parameter_float ||
	   param_type == parameter_float_expr)) {
      ERROR2("Problem::create_boundary_()",
	     "Parameter %s is of incorrect type %d",
	     param_str.c_str(),param_type);
    }

    Value * value = new Value (parameters, param_str);

    axis_enum axis = (axis_enum) config->boundary_axis[index];
    face_enum face = (face_enum) config->boundary_face[index];

    return new BoundaryValue (axis,face,value,
			      config->boundary_field_list[index]);

  } else if (type == "periodic") {

    axis_enum axis = (axis_enum) config->boundary_axis[index];
    face_enum face = (face_enum) config->boundary_face[index];

    return new BoundaryPeriodic(axis,face);

  }
  return NULL;
}

//----------------------------------------------------------------------

Initial * Problem::create_initial_
(
 std::string  type,
 int index,
 Config * config,
 Parameters * parameters,
 const FieldDescr * field_descr
 ) throw ()
{ 
  //--------------------------------------------------
  // parameter: Initial : cycle
  // parameter: Initial : time
  //--------------------------------------------------

  Initial * initial = NULL;

  if (type == "file") {
    initial = new InitialFile (parameters,
			       config->initial_cycle,
			       config->initial_time);;
  } else if (type == "value") {
    initial = new InitialValue(parameters,field_descr,
			       config->initial_cycle,
			       config->initial_time);
  } else if (type == "trace") {
    initial = new InitialTrace (config->initial_trace_field,
				config->initial_trace_mpp,
				config->initial_trace_dx,
                                config->initial_trace_dy,
                                config->initial_trace_dz);
  }

  return initial;
}

//----------------------------------------------------------------------

Refine * Problem::create_refine_
(
 std::string        type,
 Config *           config,
 Parameters *       parameters,
 const FieldDescr * field_descr,
 int                index
 ) throw ()
{ 

  if (type == "density") {

    return new RefineDensity
      (config->adapt_min_refine[index],
       config->adapt_max_coarsen[index],
       config->adapt_max_level[index],
       config->adapt_include_ghosts[index],
       config->adapt_output[index]);

  } else if (type == "slope") {

    return new RefineSlope 
      (field_descr,
       config->adapt_min_refine[index],
       config->adapt_max_coarsen[index],
       config->adapt_field_list[index],
       config->adapt_max_level[index],
       config->adapt_include_ghosts[index],
       config->adapt_output[index]);

  } else if (type == "shear") {

    return new RefineShear 
      (config->adapt_min_refine[index],
       config->adapt_max_coarsen[index],
       config->adapt_max_level[index],
       config->adapt_include_ghosts[index],
       config->adapt_output[index]);

  } else if (type == "mask") {

    std::string param_str = "Adapt:" + config->adapt_list[index] + ":value";

    return new RefineMask 
      (parameters,
       param_str,
       config->adapt_max_level[index],
       config->adapt_include_ghosts[index],
       config->adapt_output[index]);

  } else if (type == "mass") {

    double root_cell_volume = 1.0;
    for (int i=0; i<config->mesh_root_rank; i++) {
      double upper = config->domain_upper[i] ;
      double lower = config->domain_lower[i];
      int     root = config->mesh_root_size[i];

      root_cell_volume *= (upper - lower) / (root);

    }

    return new RefineMass 
      (config->adapt_min_refine[index],
       config->adapt_max_coarsen[index],
       config->adapt_level_exponent[index],
       root_cell_volume,
       config->adapt_max_level[index],
       config->adapt_include_ghosts[index],
       config->adapt_output[index]);

  } else if (type == "particle_mass") {

    return new RefineParticleMass
      (config->adapt_min_refine[index],
       config->adapt_max_coarsen[index],
       config->adapt_max_level[index],
       config->adapt_include_ghosts[index],
       config->adapt_output[index],
       config->adapt_level_exponent[index] );

  } else if (type == "particle_count") {

    return new RefineParticleCount
      (config->adapt_min_refine[index],
       config->adapt_max_coarsen[index],
       config->adapt_max_level[index],
       config->adapt_include_ghosts[index],
       config->adapt_output[index],
       config->adapt_level_exponent[index] );
  }
  return NULL;
}

//----------------------------------------------------------------------

Stopping * Problem::create_stopping_ 
(
 std::string  type,
 Config * config
 ) throw ()
/// @param type   Type of the stopping criterion to create (ignored)
/// @param config  Configuration parameter class
{
  return new Stopping(config->stopping_cycle,
		      config->stopping_time,
		      config->stopping_seconds);
}

//----------------------------------------------------------------------

Solver * Problem::create_solver_ 
( std::string  name,
  Config * config,
  int index_solver,
  const FieldDescr * field_descr,
  const ParticleDescr * particle_descr) throw ()
{
  TRACE1("Problem::create_solver %s",name.c_str());

  // No default solver
  Solver * solver = NULL;

  return solver;
}

//----------------------------------------------------------------------

Method * Problem::create_method_ 
( std::string  name,
  Config * config,
  int index_method,
  const FieldDescr * field_descr,
  const ParticleDescr * particle_descr) throw ()
{
  TRACE1("Problem::create_method %s",name.c_str());

  // No default method
  Method * method = NULL;

  if (name == "trace") {
    method = new MethodTrace(field_descr, particle_descr,
			     config->method_courant[index_method],
			     config->method_timestep[index_method]);
  }
  return method;
}

//----------------------------------------------------------------------

Output * Problem::create_output_ 
(
 std::string    name,
 int index,
 Config *  config,
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr,
 const Factory * factory
 ) throw ()
/// @param name           Name of Output object to create
/// @param index          Index of output object in Object list
/// @param config         Configuration parameter object
/// @param field_descr    Field descriptor 
/// @param particle_descr Particle descriptor 
/// @param factory        Factory object for creating Io objects of correct type
{ 

  Output * output = NULL;

  if (name == "image") {

    //--------------------------------------------------
    // parameter: Mesh : root_size
    // parameter: Mesh : root_blocks
    //--------------------------------------------------

    int nx = config->mesh_root_size[0];
    int ny = config->mesh_root_size[1];
    int nz = config->mesh_root_size[2];

    int nbx = config->mesh_root_blocks[0];
    int nby = config->mesh_root_blocks[1];
    int nbz = config->mesh_root_blocks[2];

    // NOTE: assumes cube for non-z axis images

    std::string image_type       = config->output_image_type[index];
    int         image_size_x     = config->output_image_size[index][0];
    int         image_size_y     = config->output_image_size[index][1];
    int         image_block_size = config->output_image_block_size[index];
    bool        image_ghost      = config->output_image_ghost[index];
    bool        image_log        = config->output_image_log[index];
    bool        image_abs        = config->output_image_abs[index];
    int         image_face_rank  = config->output_image_face_rank[index];
    int         min_level        = config->output_min_level[index];
    int         max_level        = std::min(config->output_max_level[index],
					    config->mesh_max_level);
    bool        leaf_only        = config->output_leaf_only[index];
    std::string image_reduce_type = config->output_image_reduce_type[index];
    std::string image_mesh_color  = config->output_image_mesh_color[index];
    std::string image_color_particle_attribute =
      config->output_image_color_particle_attribute[index];
    double      image_min = config->output_image_min[index];
    double      image_max = config->output_image_max[index];

    double image_lower[3] = { config->output_image_lower[index][0],
			      config->output_image_lower[index][1],
			      config->output_image_lower[index][2] };
    double image_upper[3] = { config->output_image_upper[index][0],
			      config->output_image_upper[index][1],
			      config->output_image_upper[index][2] };
    // AXIS

    int image_axis = config->output_axis[index][0] - 'x';

    output = new OutputImage (index,factory,
			      field_descr,
			      particle_descr,
			      CkNumPes(),
			      nx,ny,nz, 
			      nbx,nby,nbz,
			      min_level,
			      max_level,
			      leaf_only,
			      image_type,
			      image_size_x,image_size_y,
			      image_reduce_type,
			      image_mesh_color,
			      image_color_particle_attribute,
			      image_block_size,
			      image_lower, image_upper,
			      image_face_rank,
			      image_axis,
			      image_log,
			      image_abs,
			      image_ghost,
			      image_min, image_max);

  } else if (name == "data") {

    output = new OutputData (index,factory,
			     field_descr,
			     particle_descr,
			     config);

  } else if (name == "checkpoint") {

    output = new OutputCheckpoint (index,factory,
				   field_descr,
				   particle_descr,
				   config,CkNumPes());

  }

  return output;

}

//----------------------------------------------------------------------

Prolong * Problem::create_prolong_ ( std::string  name ,
				     Config * config) throw ()
{
  Prolong * prolong = 0;

  if (name == "linear") {

    prolong = new ProlongLinear;

  } else if (name == "inject") {

    prolong = new ProlongInject;

  } else {
    
    ERROR1("Problem::create_prolong_",
	  "Unrecognized Field:prolong parameter %s",name.c_str());

  }

  return prolong;
  
}

//----------------------------------------------------------------------

Restrict * Problem::create_restrict_ ( std::string  name ,
				     Config * config) throw ()
{
  Restrict * restrict = 0;

  if (name == "linear") {

    restrict = new RestrictLinear;

  } else {
    
    ERROR1("Problem::create_restrict_",
	  "Unrecognized Field:restrict parameter %s",name.c_str());

  }

  return restrict;
  
}

//----------------------------------------------------------------------
