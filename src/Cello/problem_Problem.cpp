// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Problem.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-03
/// @brief    Implementation of the Problem container class

#include "problem.hpp"

//----------------------------------------------------------------------

Problem::Problem() throw()
  : boundary_list_(),
    is_periodic_(true),
    initial_list_(),
    physics_list_(),
    refine_list_(),
    stopping_(nullptr),
    solver_list_(),
    method_list_(),
    output_list_(),
    prolong_list_(),
    restrict_list_(),
    units_(nullptr),
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

  if (pk) n=physics_list_.size();
  p | n;
  if (up) physics_list_.resize(n);
  for (int i=0; i<n; i++) {
    p | physics_list_[i]; // PUP::able
  }

  if (pk) n=refine_list_.size();
  p | n;
  if (up) refine_list_.resize(n);
  for (int i=0; i<n; i++) {
    p | refine_list_[i]; // PUP::able
  }

  p | stopping_;

  p | units_;

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

  if (pk) n=prolong_list_.size();
  p | n;
  if (up) prolong_list_.resize(n);
  for (int i=0; i<n; i++) {
    p | prolong_list_[i]; // PUP::able
  }

  if (pk) n=restrict_list_.size();
  p | n;
  if (up) restrict_list_.resize(n);
  for (int i=0; i<n; i++) {
    p | restrict_list_[i]; // PUP::able
  }
  
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
	  type.c_str(),  boundary != nullptr);

    boundary_list_.push_back(boundary);
  }

}

//----------------------------------------------------------------------

void Problem::initialize_initial(Config * config,
				 Parameters * parameters) throw()
{

  for (int index=0; index < config->num_initial; index++) {

    std::string type = config->initial_list[index];

    Initial * initial = create_initial_ (type,index,config,parameters);

    ASSERT1("Problem::initialize_initial",
	    "Initial type %s not recognized",
	    config->initial_list[index].c_str(),
	    (initial != nullptr) );

    initial_list_.push_back( initial );
  }
}

//----------------------------------------------------------------------

void Problem::initialize_physics(Config * config,
				 Parameters * parameters) throw()
{

  for (int index=0; index < config->num_physics; index++) {

    std::string type = config->physics_list[index];

    Physics * physics = create_physics_
      (type,index,config,parameters);

    ASSERT1("Problem::initialize_physics",
	    "Physics type %s not recognized",
	    config->physics_list[index].c_str(),
	    (physics != nullptr) );

    physics_list_.push_back( physics );
  }

}

//----------------------------------------------------------------------

void Problem::initialize_refine(Config * config,
				Parameters * parameters) throw()
{
  for (int index=0; index<config->num_adapt; index++) {

    std::string name = config->adapt_type[index];

    Refine * refine = create_refine_ 
      (name,index,config,parameters);

    if (refine) {
      refine_list_.push_back( refine );
      int index_schedule = config->adapt_schedule_index[index];

      if (index_schedule >= 0) {
	refine->set_schedule
	  (Schedule::create( config->schedule_var[index_schedule],
			     config->schedule_type[index_schedule],
			     config->schedule_start[index_schedule],
			     config->schedule_stop[index_schedule],
			     config->schedule_step[index_schedule],
			     config->schedule_list[index_schedule]));
      }
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
	  stopping_ != nullptr);
}

//----------------------------------------------------------------------

void Problem::initialize_prolong(Config * config) throw()
{
  // default prolongation
  ASSERT ("Problem::initialize_prolong()",
          "Initial default prolongation must be added to Problem::prolong_list_ first",
          (prolong_list_.size() == 0));
  
  Prolong * prolong = create_prolong_(config->field_prolong,config);

  ASSERT1("Problem::initialize_prolong",
	  "Prolong type %s not recognized",
	  config->field_prolong.c_str(),
	  prolong != nullptr);

  prolong_list_.push_back(prolong);

}  

//----------------------------------------------------------------------

void Problem::initialize_restrict(Config * config) throw()
{
  // default restriction
  ASSERT ("Problem::initialize_restrict()",
          "Initial default restriction must be added to Problem::restrict_list_ first",
          (restrict_list_.size() == 0));

  Restrict * restrict = create_restrict_(config->field_restrict,config);

  ASSERT1("Problem::initialize_restrict",
	  "Restrict type %s not recognized",
	  config->field_restrict.c_str(),
	  restrict != nullptr);

  restrict_list_.push_back(restrict);

}

//----------------------------------------------------------------------

void Problem::initialize_output
(Config * config, const Factory * factory) throw()
{
  FieldDescr * field_descr = cello::field_descr();
  
  for (int index=0; index < config->num_output; index++) {

    std::string type       = config->output_type[index];

    Output * output = create_output_ (type,index, config,factory);

    if (output == nullptr) {
      ERROR2("Problem::initialize_output",
	     "Unknown parameter type Output:%s:type = %s",
	     config->output_list[index].c_str(),type.c_str());
    } else {

      if (config->output_name[index].size() > 0) {
        std::string file_name = config->output_name[index][0];

        std::vector<std::string> file_args;

        for (size_t i=1; i<config->output_name[index].size(); i++) {
          file_args.push_back(config->output_name[index][i]);
        }

        output->set_filename (file_name,file_args);
      }

      if (config->output_dir[index].size() > 0 ||
          config->output_dir_global != ".") {

        std::string dir_name;
        std::vector<std::string> dir_args;
      
        dir_name = config->output_dir_global;

        if (config->output_dir[index].size() > 0) {
          dir_name = dir_name + "/" + config->output_dir[index][0];

          for (size_t i=1; i<config->output_dir[index].size(); i++) {
            dir_args.push_back(config->output_dir[index][i]);
          }
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

      ParticleDescr * particle_descr = cello::particle_descr();

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

      if (output_image != nullptr) {

        // COLORMAP

        int n = config->output_colormap[index].size() / 3;

        if (n > 0) {
          std::vector<float> colormap[3];
          for (int i=0; i<3; i++) {
            colormap[i].resize(n,0.0);
          }

          for (int i=0; i<n; i++) {

            for (int rgb=0; rgb<3; rgb++) {
              int index_colormap = 3*i+rgb;
              colormap[rgb][i] =
                config->output_colormap[index][index_colormap];
            }

          }

          output_image->set_colormap(colormap);
        }

      }

    }

    // Add the initialized Output object to the Simulation's list of
    // output objects

    output_list_.push_back(output); 

  } // (for index)

}

//----------------------------------------------------------------------

void Problem::initialize_method
( Config * config, const Factory * factory ) throw()
{
  const size_t num_method = config->method_list.size();

  Method::courant_global = config->method_courant_global;

  method_list_.push_back(new MethodNull(config->method_null_dt)); 
  
  for (size_t index_method=0; index_method < num_method ; index_method++) {

    std::string name = config->method_type[index_method];

    Method * method = create_method_(name, index_method, config, factory);

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

void Problem::initialize_solver( Config * config ) throw()
{
  const size_t num_solver = config->solver_list.size();

  for (size_t index_solver=0; index_solver < num_solver ; index_solver++) {

    std::string type = config->solver_type[index_solver];

    Solver * solver = create_solver_(type, index_solver, config);

    if (solver) {

      solver_list_.push_back(solver); 

    } else {
      ERROR1("Problem::initialize_solver",
	     "Unknown Solver %s",type.c_str());
    }
  }
}

//----------------------------------------------------------------------

void Problem::initialize_units(Config * config) throw()
{
  units_ = create_units_(config);

  ASSERT("Problem::initialize_units",
	  "Units object not successfully created",
	  units_ != nullptr);
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
  delete units_;         units_ = nullptr;
  for (size_t i=0; i<physics_list_.size(); i++) {
    delete physics_list_[i];   physics_list_[i] = nullptr;
  }
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

    return new BoundaryPeriodic(axis);

  }
  return nullptr;
}

//----------------------------------------------------------------------

Initial * Problem::create_initial_
(
 std::string  type,
 int index,
 Config * config,
 Parameters * parameters
 ) throw ()
{ 
  //--------------------------------------------------
  // parameter: Initial : cycle
  // parameter: Initial : time
  //--------------------------------------------------

  Initial * initial = nullptr;

  if (type == "file") {
    initial = new InitialFile (parameters,
			       config->initial_cycle,
			       config->initial_time);;
  } else if (type == "value") {
    initial = new InitialValue(parameters,
			       config->initial_cycle,
			       config->initial_time);
  } else if (type == "trace") {
    initial = new InitialTrace (config->initial_trace_name,
				config->initial_trace_field,
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
 int                index,
 Config *           config,
 Parameters *       parameters
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
      (config->adapt_min_refine[index],
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

  } else if (type == "particle_count") {

    return new RefineParticleCount
      (config->adapt_min_refine[index],
       config->adapt_max_coarsen[index],
       config->adapt_max_level[index],
       config->adapt_include_ghosts[index],
       config->adapt_output[index],
       config->adapt_level_exponent[index] );
  }
  return nullptr;
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

Units * Problem::create_units_ 
(
 Config * config
 ) throw ()
/// @param type   Type of the units criterion to create (ignored)
/// @param config  Configuration parameter class
{
  Units * units = new Units;
  
  if (config->units_mass == 1.0) {

    units->set_using_density (config->units_length,
			       config->units_density,
			       config->units_time);
    
  } else if (config->units_density == 1.0) {

    units->set_using_mass (config->units_length,
			    config->units_mass,
			    config->units_time);
  } else {
    
    ERROR("Problem::create_units_",
	  "Cannot set both Units:density and Units:time parameters");
  }

  return units;
}

//----------------------------------------------------------------------

Solver * Problem::create_solver_ 
( std::string  type,
  int index_solver,
  Config * config
  ) throw ()
{
  TRACE1("Problem::create_solver %s",type.c_str());

  Solver * solver = nullptr;

  if (type == "null") {

    Prolong * prolong = create_prolong_
      (config->solver_prolong[index_solver],config);
    Restrict * restrict = create_restrict_
      (config->solver_restrict[index_solver],config);

    const int index_prolong = prolong_list_.size();
    const int index_restrict = restrict_list_.size();
    prolong_list_.push_back(prolong);
    restrict_list_.push_back(restrict);
    
    solver = new SolverNull
      (config->solver_list         [index_solver],
       config->solver_field_x      [index_solver],
       config->solver_field_b      [index_solver],
       0, // restart cycle
       solve_leaf,
       index_prolong,
       index_restrict,
       config->solver_min_level[index_solver],
       config->solver_max_level[index_solver]);
  }

  if (solver) solver->set_index(index_solver);
  
  return solver;
}

//----------------------------------------------------------------------

Physics * Problem::create_physics_ 
( std::string type,
  int index,
  Config * config,
  Parameters * parameters) throw ()
{
  TRACE1("Problem::create_physics %s",type.c_str());

  // No default physics
  Physics * physics = nullptr;

  return physics;
}

//----------------------------------------------------------------------

Physics * Problem::physics (std::string type) const throw()
{
  for (size_t i=0; i<physics_list_.size(); i++) {
    if (physics_list_[i]->type() == type) return physics_list_[i];
  }
  return nullptr;
}

//----------------------------------------------------------------------

Method * Problem::method (std::string name) const throw()
{
  for (size_t i=0; i<method_list_.size(); i++) {
    if (method_list_[i]->name() == name) return method_list_[i];
  }
  return nullptr;
}

//----------------------------------------------------------------------

bool Problem::method_exists(const std::string& name) const throw() {
  for (size_t i = 0; i < method_list_.size(); i++) {
    if (method_list_[i]->name() == name) return true;
  }
  return false;
}

//----------------------------------------------------------------------

bool Problem::method_precedes(const std::string& name1,
                              const std::string& name2) const throw() {
  size_t ind_1 = 0;
  size_t ind_2 = 0;

  bool method_1_found = false;
  bool method_2_found = false;

  bool no_repeats = true;
  for (size_t i = 0; i < method_list_.size(); i++) {
    if (method_list_[i]->name() == name1) {
      if (method_1_found) {
        no_repeats = false;
        break;
      } else {
        method_1_found = true;
        ind_1 = i;
      }
    }
    if (method_list_[i]->name() == name2) {
      if (method_2_found) {
        no_repeats = false;
        break;
      } else {
        method_2_found = true;
        ind_2 = i;
      }
    }
  }

  return no_repeats && method_1_found && method_2_found && (ind_1 < ind_2);
}

//----------------------------------------------------------------------

Compute * Problem::create_compute
  ( std::string name,
    Config * config ) throw ()
{
  TRACE1("Problem::create_compute %s", name.c_str());

  Compute * compute = nullptr;

  return compute;
}

//----------------------------------------------------------------------

Method * Problem::create_method_ 
( std::string  name,
  int index_method,
  Config * config,
  const Factory * factory
  ) throw ()
{
  TRACE1("Problem::create_method %s",name.c_str());

  // No default method
  Method * method = nullptr;

  if (name == "trace") {
    method = new MethodTrace(config->method_courant[index_method],
			     config->method_timestep[index_method],
			     config->method_trace_name[index_method]);
  } else if (name == "null") {

    method = new MethodNull
      (config->method_null_dt);

  } else if (name == "flux_correct") {

    method = new MethodFluxCorrect
      (config->method_flux_correct_group[index_method],
       config->method_flux_correct_enable[index_method],
       config->method_flux_correct_min_digits_fields[index_method],
       config->method_flux_correct_min_digits_values[index_method]);

  } else if (name == "checkpoint") {

    method = new MethodCheckpoint
      ( config->method_path_name[index_method]);

  } else if (name == "output") {

    ASSERT("Problem::create_method_()",
           "MethodOutput must have 'file_name' parameter set",
           config->method_file_name[index_method].size() > 0);
    ASSERT("Problem::create_method_()",
           "MethodOutput must have 'path_name' parameter set",
           config->method_path_name[index_method].size() > 0);

    method = new MethodOutput
      ( factory,
        config->method_file_name[index_method],
        config->method_path_name[index_method],
        config->method_field_list[index_method],
        config->method_particle_list[index_method],
        config->method_ghost_depth[index_method],
        config->method_min_face_rank[index_method],
        config->method_all_fields[index_method],
        config->method_all_particles[index_method],
        config->method_output_all_blocks[index_method],
        config->method_output_blocking[0][index_method],
        config->method_output_blocking[1][index_method],
        config->method_output_blocking[2][index_method]);

  } else if (name == "refresh") {

    method = new MethodRefresh
      (config->method_field_list[index_method],
       config->method_particle_list[index_method],
       config->method_ghost_depth[index_method],
       config->method_min_face_rank[index_method],
       config->method_all_fields[index_method],
       config->method_all_particles[index_method]);

  } else if (name == "debug") {

    method = new MethodDebug
      (config->num_fields,
       config->method_debug_print[index_method],
       config->method_debug_coarse[index_method],
       config->method_debug_ghost[index_method]);

  } else if (name == "close_files") {

    method = new MethodCloseFiles
      (config->method_close_files_seconds_stagger[index_method],
       config->method_close_files_seconds_delay[index_method],
       config->method_close_files_group_size[index_method]);
      
  }
  return method;
}

//----------------------------------------------------------------------

Output * Problem::create_output_ 
(
 std::string    name,
 int index,
 Config *  config,
 const Factory * factory
 ) throw ()
/// @param name           Name of Output object to create
/// @param index          Index of output object in Object list
/// @param config         Configuration parameter object
/// @param factory        Factory object for creating Io objects of correct type
{ 

  Output * output = nullptr;

  if (name == "image") {

    // NOTE: assumes cube for non-z axis images

    std::string image_type       = config->output_image_type[index];
    bool        image_ghost      = config->output_image_ghost[index];
    bool        image_log        = config->output_image_log[index];
    bool        image_abs        = config->output_image_abs[index];
    int         image_face_rank  = config->output_image_face_rank[index];
    int         min_level        = config->output_min_level[index];
    int         max_level        = std::min(config->output_max_level[index],
					    config->mesh_max_level);
    bool        leaf_only        = config->output_leaf_only[index];
    int         image_size[2] = { config->output_image_size[index][0],
                                  config->output_image_size[index][1] };
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

    bool use_min_max =
      (image_min != std::numeric_limits<double>::max()) &&
      (image_max != -std::numeric_limits<double>::max());
    
    output = new OutputImage (index,factory,
			      CkNumPes(),
			      config->mesh_root_size,
			      config->mesh_root_blocks,
			      min_level,
			      max_level,
			      leaf_only,
			      image_type,
			      image_size,
			      image_reduce_type,
			      image_mesh_color,
			      image_color_particle_attribute,
			      image_lower, image_upper,
			      image_face_rank,
			      image_axis,
			      image_log,
			      image_abs,
			      image_ghost,
                              use_min_max,
			      image_min, image_max);

  } else if (name == "data") {

    output = new OutputData (index,factory,config);

  } else if (name == "checkpoint") {

    output = new OutputCheckpoint (index,factory,
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
