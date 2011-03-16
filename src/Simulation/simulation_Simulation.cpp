// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      simulation_Simulation.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2010-11-10
/// @brief     Implementation of the Simulation class

#include "cello.hpp"

#include "simulation.hpp"

/// Initialize the Simulation object
Simulation::Simulation
(
 Factory *      factory,
 Parameters *   parameters,
 GroupProcess * group_process
)
  : factory_(factory),
    dimension_(0),
    group_process_(group_process),
    parameters_(parameters),
    mesh_(0),
    field_descr_(0),
    stopping_(0),
    timestep_(0),
    initial_(0),
    boundary_(0),
    method_list_()
{
  ASSERT("Simulation::Simulation","Parameters is NULL", parameters != NULL);

  lower_[0] = 0.0;
  lower_[1] = 0.0;
  lower_[2] = 0.0;

  upper_[0] = 1.0;
  upper_[1] = 1.0;
  upper_[2] = 1.0;

}

//----------------------------------------------------------------------

Simulation::~Simulation() throw()
{
  deallocate_();
}

//----------------------------------------------------------------------

Simulation::Simulation(const Simulation & simulation) throw()
{
  INCOMPLETE("Simulation::Simulation(simulation)");
}

//----------------------------------------------------------------------


Simulation & Simulation::operator= (const Simulation & simulation) throw()
{
  INCOMPLETE("Simulation::operatior = (simulation)");
  return (*this);
}

//----------------------------------------------------------------------

void Simulation::initialize() throw()
{

  // Initialize parameters

  initialize_simulation_();

  // Initialize simulation components

  initialize_data_();
  initialize_mesh_();

  initialize_stopping_();
  initialize_timestep_();
  initialize_initial_();
  initialize_boundary_();
  initialize_method_();

}

//----------------------------------------------------------------------

void Simulation::finalize() throw()
{
  deallocate_();
}

//----------------------------------------------------------------------

int Simulation::dimension() const throw()
{
  return dimension_;
}

//----------------------------------------------------------------------

void Simulation::lower (double * xm, double *ym, double *zm) const throw()
{
  if (xm) *xm = lower_[0];
  if (ym) *ym = lower_[1];
  if (zm) *zm = lower_[2];
}

//----------------------------------------------------------------------

void Simulation::upper (double * xp, double *yp, double *zp) const throw()
{
  if (xp) *xp = upper_[0];
  if (yp) *yp = upper_[1];
  if (zp) *zp = upper_[2];
}

//----------------------------------------------------------------------

Mesh * Simulation::mesh() const throw()
{
  return mesh_;
}
  
//----------------------------------------------------------------------

FieldDescr * Simulation::field_descr() const throw()
{
  return field_descr_;
}

//----------------------------------------------------------------------

Stopping * Simulation::stopping() const throw() 
{ return stopping_; }

//----------------------------------------------------------------------

Timestep * Simulation::timestep() const throw() 
{ return timestep_; }

//----------------------------------------------------------------------

Initial * Simulation::initial() const throw() 
{ return initial_; }

//----------------------------------------------------------------------

Boundary * Simulation::boundary() const throw() 
{ return boundary_; }

//----------------------------------------------------------------------

int Simulation::num_method() const throw()
{ return method_list_.size(); }

//----------------------------------------------------------------------

Method * Simulation::method(int i) const throw()
{ return method_list_[i]; }

//----------------------------------------------------------------------

Factory * Simulation::factory () const throw()
{ return factory_; }

//======================================================================

void Simulation::initialize_simulation_() throw()
{

  //--------------------------------------------------
  parameters_->set_current_group("Physics");
  //--------------------------------------------------

  // parameter: Physics::dimensions

  dimension_ = parameters_->value_integer("dimensions",0);


  //--------------------------------------------------
  parameters_->set_current_group("Domain");
  //--------------------------------------------------

  // parameter: Domain::dimensions
  // parameter: Domain::extent

  // get extent_length and check consistency

  int extent_length = parameters_->list_length("extent");

  bool valid_extent_length = ((extent_length == 2) ||
			      (extent_length == 4) ||
			      (extent_length == 6));
  ASSERT ("Simulation::initialize_simulation_",
	  "Parameter Domain:extent list must have length 2, 4, or 6",
	  valid_extent_length);

  ASSERT ("Simulation::initialize_simulation_",
	  "Parameter mismatch between Physics:dimensions and Domain:extent",
	  (dimension_ == extent_length / 2));

  for (int i=0; i<extent_length/2; i++) {

    lower_[i] = parameters_->list_value_scalar(2*i,   "extent", 0.0);
    upper_[i] = parameters_->list_value_scalar(2*i+1, "extent", 1.0);

    if (lower_[i] > upper_[i]) {
      char error_message[ERROR_LENGTH];
      sprintf (error_message,
	       "Parameter Domain::extent[%d]==%g not lower than "
	       "Domain::extent[%d]==%g",
	       2*i,lower_[i],2*i+1,upper_[i]);
      ERROR ("Simulation::Simulation", error_message);
    }
  }

}

//----------------------------------------------------------------------

void Simulation::initialize_data_() throw()
{

  /* MEMORY LEAK */
  field_descr_ = new FieldDescr;

  //--------------------------------------------------
  parameters_->set_current_group("Field");
  //--------------------------------------------------

  // parameter: Field::fields

  // Add data fields

  int i;
  for (i=0; i<parameters_->list_length("fields"); i++) {
    field_descr_->insert_field
      (parameters_->list_value_string(i, "fields"));
  }

  // Define default ghost zone depth for all fields, default value of 1

  // @@@ WARNING: REPEATED CODE: SEE enzo_EnzoBlock.cpp

  // parameter: Field::ghosts

  int gx = 1;
  int gy = 1;
  int gz = 1;

  if (parameters_->type("ghosts") == parameter_integer) {
    gx = gy = gz = parameters_->value_integer("ghosts",1);
  } else if (parameters_->type("ghosts") == parameter_list) {
    gx = parameters_->list_value_integer(0,"ghosts",1);
    gy = parameters_->list_value_integer(1,"ghosts",1);
    gz = parameters_->list_value_integer(2,"ghosts",1);
  }

  if (dimension_ < 2) gy = 0;
  if (dimension_ < 3) gz = 0;

  printf ("field count = %d\n",field_descr_->field_count());
  for (i=0; i<field_descr_->field_count(); i++) {
    field_descr_->set_ghosts(i,gx,gy,gz);
    int gx2,gy2,gz2;
    field_descr_->ghosts(i,&gx2,&gy2,&gz2);
    printf ("%s:%d %d %d %d\n",__FILE__,__LINE__,gx2,gy2,gz2);
  }

  // Set precision

  // parameter: Field::precision

  std::string precision_str = parameters_->value_string("precision","default");

  precision_enum precision = precision_default;

  if (precision_str == "single")
    precision = precision_single;
  else if (precision_str == "double")
    precision = precision_double;
  else if (precision_str == "quadruple")
    precision = precision_quadruple;

  for (i=0; i<field_descr_->field_count(); i++) {
    field_descr_->set_precision(i,precision);
  }
}

//----------------------------------------------------------------------

void Simulation::initialize_mesh_() throw()
{

  ASSERT("Simulation::initialize_mesh_",
	 "data must be initialized before mesh",
	 field_descr_ != NULL);

  //--------------------------------------------------
  parameters_->set_current_group ("Mesh");
  //--------------------------------------------------

  int root_size[3];

  // parameter: Mesh::root_size

  root_size[0] = parameters_->list_value_integer(0,"root_size",1);
  root_size[1] = parameters_->list_value_integer(1,"root_size",1);
  root_size[2] = parameters_->list_value_integer(2,"root_size",1);

  int root_blocks[3];

  // parameter: Mesh::root_blocks

  root_blocks[0] = parameters_->list_value_integer(0,"root_blocks",1);
  root_blocks[1] = parameters_->list_value_integer(1,"root_blocks",1);
  root_blocks[2] = parameters_->list_value_integer(2,"root_blocks",1);

  mesh_ = factory()->create_mesh 
    (group_process_,
     root_size[0],root_size[1],root_size[2],
     root_blocks[0],root_blocks[1],root_blocks[2]);

  // Allocate and insert the root patch, using all processes

  Patch * root_patch = factory()->create_patch
    (group_process_,
     root_size[0],root_size[1],root_size[2],
     root_blocks[0],root_blocks[1],root_blocks[2]);

  mesh_->insert_patch(root_patch);

  // Domain extents

  //--------------------------------------------------
  parameters_->set_current_group("Domain");
  //--------------------------------------------------

  double domain_lower[3];

  // parameter: Domain::extent

  domain_lower[0] = parameters_->list_value_scalar(0,"extent",0.0);
  domain_lower[1] = parameters_->list_value_scalar(2,"extent",0.0);
  domain_lower[2] = parameters_->list_value_scalar(4,"extent",0.0);

  mesh_->set_lower(domain_lower[0],
		   domain_lower[1],
		   domain_lower[2]);

  double domain_upper[3];

  domain_upper[0] = parameters_->list_value_scalar(1,"extent",1.0);
  domain_upper[1] = parameters_->list_value_scalar(3,"extent",1.0);
  domain_upper[2] = parameters_->list_value_scalar(5,"extent",1.0);

  mesh_->set_upper(domain_upper[0],
		   domain_upper[1],
		   domain_upper[2]);

  // Create the root patch

  root_patch->set_lower(domain_lower[0],
			domain_lower[1],
			domain_lower[2]);

  root_patch->set_upper(domain_upper[0],
			domain_upper[1],
			domain_upper[2]);

  root_patch->allocate_blocks(field_descr_);

  // Parallel layout of the root patch
  
  Layout * layout = root_patch->layout();

  //--------------------------------------------------
  parameters_->set_current_group("Mesh");
  //--------------------------------------------------

  // parameter: Mesh::root_process_first
  // parameter: Mesh::root_process_count

  int process_first = parameters_->value_integer("root_process_first",0);
  int process_count = parameters_->value_integer("root_process_count",1);

  layout->set_process_range(process_first, process_count);

  // Mesh AMR settings

  // parameter: Mesh::refine
  // parameter: Mesh::max_level
  // parameter: Mesh::balanced
  // parameter: Mesh::backfill
  // parameter: Mesh::coalesce

  mesh_->set_refine_factor (parameters_->value_integer("refine",    2));
  mesh_->set_max_level     (parameters_->value_integer("max_level", 0));
  mesh_->set_balanced      (parameters_->value_logical("balanced",  true));
  mesh_->set_backfill      (parameters_->value_logical("backfill",  true));
  mesh_->set_coalesce      (parameters_->value_logical("coalesce",  true));

}

//----------------------------------------------------------------------

void Simulation::initialize_stopping_() throw()
{
  stopping_ = create_stopping_("ignored");
}

//----------------------------------------------------------------------

void Simulation::initialize_timestep_() throw()
{
  timestep_ = create_timestep_("ignored");
}

//----------------------------------------------------------------------

void Simulation::initialize_initial_() throw()
{
  //--------------------------------------------------
  parameters_->set_current_group("Initial");
  //--------------------------------------------------

  // parameter: Initial::problem

  std::string name = parameters_->value_string("problem","default");

  initial_ = create_initial_(name);

  if (initial_ == 0) {
    initial_ = new InitialDefault(parameters_);
  }
}

//----------------------------------------------------------------------

void Simulation::initialize_boundary_() throw()
{
  //--------------------------------------------------
  parameters_->set_current_group("Boundary");
  //--------------------------------------------------

  // parameter: Boundary::name

  std::string name = parameters_->value_string("name","");
  boundary_ = create_boundary_(name);
}


//----------------------------------------------------------------------

void Simulation::initialize_method_() throw()
{
  //--------------------------------------------------
  parameters_->set_current_group("Method");
  //--------------------------------------------------

  int method_count = parameters_->list_length("sequence");

  if (method_count == 0) {
    ERROR ("Simulation::initialize_method_",
		   "List parameter 'Method sequence' must have length greater than zero");
  }

  for (int i=0; i<method_count; i++) {

    // parameter: Method::sequence

    std::string method_name = parameters_->list_value_string(i,"sequence");

    Method * method = create_method_(method_name);
    if (method) {

      method_list_.push_back(method); 

    } else {
      char error_message[ERROR_LENGTH];
      sprintf (error_message,"Unknown Method %s",method_name.c_str());
      ERROR ("Simulation::initialize_method_",
		     error_message);
    }
  }


  INCOMPLETE("Simulation::initialize_method_");
}

//----------------------------------------------------------------------

void Simulation::deallocate_() throw()
{
  delete mesh_;
  mesh_ = 0;
  delete field_descr_;
  field_descr_ = 0;
  delete stopping_;
  stopping_ = 0;
  delete timestep_;
  timestep_ = 0;
  delete initial_;
  initial_ = 0;
  delete boundary_;
  boundary_ = 0;
  delete factory_;
  factory_ = 0;
  for (size_t i=0; i<method_list_.size(); i++) {
    delete method_list_[i];
    method_list_[i] = 0;
  }
}
