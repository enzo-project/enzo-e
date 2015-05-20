// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravityMg.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-21 17:25:09
/// @brief    Implements the EnzoMethodGravityMg class
///
/// Multigrid method on an adaptive mesh.  We use the MLAT (Multilevel
/// Adaptive Technique) of Brandt '77, based on the FAS formulation of
/// Multigrid.
///
///======================================================================
///
///  "Coarse" view of Multigrid solver
///
///   @code
///    MG(A,X,B)
///
///    enter_solver()          // initialize solver
///    begin_cycle()           // initialize a V-cycle
///    send_faces()            // send face data to neighbors
///    p_mg_receive_face()     // receive face daat from neighbors
///    pre_smooth()            // apply the smoother
///    p_mg_restrict()         // receive restricted data from parent
///    evaluate_b()            // compute the right-hand side B
///    coarse_solve()          // solve the coarse-grid equation
///    p_mg_end_coarse_solve() // callback after the coarse solver
///    p_mg_prolong()          // receive prolonged data from children
///    update_x()              // update the solution X
///    post_smooth()           // apply the smoother
///    end_cycle()             // finalize V-cycle
///    exit_solver()           // finalize solver
///  @endcode
///
///
/// ======================================================================
///
///  "Fine" view of Multigrid algorithm
///
/// @code
///  enter_solver()
///
///    initialize iter_ = 0
///    call begin_cycle()
/// @endcode
///
/// @code
///  begin_cycle()
///
///   if converged {
///      call exit_solver()
///   } else {
///      if leaf block {
///         call send_faces()
///      }
///   }
/// @endcode
///
/// @code
///  send_faces()
///
///    <have updated X>
///
///   initialize call counter
///
///   if leaf block, then
///      for each neighbor
///          if not coarse neighbor
///             pack face data
///          if neighbor level same
///             remote call p_receive_face() on neighbor
///          if neighbor level coarser
///             remote call p_receive_face() on neighbor
///             remote call p_receive_face() on neighbor parent
///   else
///      for each face
///         pack face data
///         remote call p_receive_face() on neighbor
///
///   counter = determine_count ()
///
///   local call receive_face(counter)
/// @endcode
///    
/// @code
///  determine_count()
///
///   if leaf block, then
///      for each face
///         if non-repeated coarse face
///            increment counter
///   else
///      count number neighbors in same level 
///
///
///   if not leaf block, then
///      count number of children
///
///   count self
/// @endcode
///	    
/// @code
///  p_receive_face()
///
///    receive_face(counter)
///
///    unpack ghost data
///    if (number of calls >= counter), then
///       call compute_correction()
/// @endcode
///
///  @code
///  compute_correction()
///
///    if coarsest level, then
///        call coarse_solve()
///    else 
///        if not leaf block, then
///            call evaluate_b()
///        call pre_smooth()
/// @endcode
///
/// @code
///  pre_smooth()
///
///     apply the smoother to A X = B
///     compute Y = A*X  [ need refresh? ]
///     pack B, X, Y fields
///     remote call p_restrict() on parent
/// @endcode
///
/// @code
///  p_restrict(B,X,Y)
///
///     unpack vectors B,X,Y
///     call send_faces()
///     if (number of calls >= counter), then
///        call compute_correction()
/// @endcode
///
/// @code
///  evaluate_b()
///
///     B = B - Y
///     compute Y = A*X [ need refresh? ]
///     B = B + Y
/// @endcode
///
/// @code
///  coarse_solve()
///     
///    [need refreshed X?]
///    initiate the coarse solver with p_coarse_solved() callback
///
///    [ NOTE: only coarse blocks call this ]
/// @endcode
///
/// @code
///  p_coarse_solved()
///
///    if leaf, then
///       call end_cycle()
///    else
///       for each child
///           pack correction X
///           remote call p_prolong() on child
/// @endcode
///   
/// @code
///  p_prolong()
///
///    unpack correction Y
///    apply correction to solution X:
///       X = X + Y
///    call post_smooth()  [ need refresh? ]
/// @endcode
///
/// @code
///  post_smooth()
///
///    apply smoother to A X = B
///    if leaf node, then
///       call end_cycle()
///    else
///       call prolong()
/// @endcode
///
/// @code
///  end_cycle()
///
///    increment iter_
///    call begin_cycle()
/// @endcode
///
/// @code
///  exit_solver()
///
///    deallocate temporaries
///    end_compute() with quiescence
/// @endcode
///
///======================================================================
///
/// Required Fields
///
/// - B                          linear system right-hand side
/// - R                          residual R = B - A*X
/// - X                          current solution to A*X = B
/// - Y                          temporary vector for A*X on fine
/// - potential                  computed gravitational potential
/// - density                    density field
/// - acceleration_x             acceleration along X-axis
/// - acceleration_y (rank >= 2) acceleration along Y-axis
/// - acceleration_z (rank >= 3) acceleration along Z-axis

#include "cello.hpp"
#include "enzo.hpp"
#include "enzo.decl.h"

#define CK_TEMPLATES_ONLY
#include "enzo.def.h"
#undef CK_TEMPLATES_ONLY

#define MG_VERBOSE(FUNCTION)				\
  Monitor * monitor = enzo_block->simulation()->monitor();	\
  monitor->verbose(stdout,"Method", "%s "FUNCTION,enzo_block->name().c_str());

#define VERBOSE(FUNCTION)				\
  Monitor * monitor = simulation()->monitor();	\
  monitor->verbose(stdout,"Method", "%s "FUNCTION,name().c_str());

//----------------------------------------------------------------------

extern CkReduction::reducerType r_method_gravity_mg_type;

//----------------------------------------------------------------------

EnzoMethodGravityMg::EnzoMethodGravityMg 
(FieldDescr * field_descr, int rank,
 double grav_const, int iter_max, double res_tol, int monitor_iter,
 bool is_singular,
 Compute * smooth,
 Restrict * restrict,
 Prolong * prolong,
 int level_min,
 int level_max) 
  : Method(), 
    A_(new EnzoMatrixLaplace),
    smooth_(smooth),
    restrict_(restrict),
    prolong_(prolong),
    is_singular_(is_singular),
    rank_(rank),
    grav_const_(grav_const),
    iter_max_(iter_max), 
    res_tol_(res_tol),
    monitor_iter_(monitor_iter),
    rr_(0),rr0_(0), rr_min_(0),rr_max_(0),
    idensity_(0),  ipotential_(0),
    ib_(0), ir_(0), ix_(0), iy_(0),
    iter_(0),
    level_min_(level_min),
    level_max_(level_max)
{
  idensity_   = field_descr->field_id("density");
  ipotential_ = field_descr->field_id("potential");

  ib_ = field_descr->field_id("B");
  ir_ = field_descr->field_id("R");
  ix_ = field_descr->field_id("X");
  iy_ = field_descr->field_id("Y");
}

//----------------------------------------------------------------------

EnzoMethodGravityMg::~EnzoMethodGravityMg () throw()
{
  delete smooth_;
  delete prolong_;
  delete restrict_;

  smooth_ = NULL;
  prolong_ = NULL;
  restrict_ = NULL;
}


//----------------------------------------------------------------------

void EnzoMethodGravityMg::compute ( Block * block) throw()
{
  Field field = block->data()->field();

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  MG_VERBOSE("compute()");

  precision_ = field.precision(idensity_);

  if      (precision_ == precision_single)    
    enter_solver_<float>      (enzo_block);
  else if (precision_ == precision_double)    
    enter_solver_<double>     (enzo_block);
  else if (precision_ == precision_quadruple) 
    enter_solver_<long double>(enzo_block);
  else 
    ERROR1("EnzoMethodGravityMg()", "precision %d not recognized", precision_);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg::enter_solver_ (EnzoBlock * enzo_block) throw()
//     initialize iter_
//     call begin_cycle()
{

  MG_VERBOSE("enter_solver_()");
  
  iter_ = 0;

  Data * data = enzo_block->data();
  Field field = data->field();

  T * density = (T*) field.values(idensity_);
    
  T * B = (T*) field.values(ib_);
  T * X = (T*) field.values(ix_);
  T * R = (T*) field.values(ir_);
  T * Y = (T*) field.values(iy_);

  // Initialize B, X, R, Y

  // X = 0
  // B = -h^2 * 4 * PI * G * density
  // R = B ( residual with X = 0 )
  // Y = 0

  int mx,my,mz;
  field.dimensions(idensity_,&mx,&my,&mz);

  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {
	int i = ix + mx*(iy + my*iz);
	B[i] = - 4.0 * (cello::pi) * grav_const_ * density[i];
	X[i] = 0.0;
	R[i] = B[i];
	Y[i] = 0.0;
      }
    }
  }


  begin_cycle_<T>(enzo_block);

}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg::begin_cycle_(EnzoBlock * enzo_block) throw()
//    if converged, then
//       call exit_solver()
//    else 
//       if leaf block, then
//          call send_faces()
{
  MG_VERBOSE("begin_cycle_()");

  //    if converged, then return
  if (iter_ >= iter_max_) {

    exit_solver_<T>(enzo_block, return_error_max_iter_reached);

  } else {

    //  else call send_faces()
    // (only leaf blocks have data initially)

    if (enzo_block->is_leaf()) {

      send_faces_(enzo_block);
    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodGravityMg::send_faces_(EnzoBlock * enzo_block) throw()
//    <have updated X>
//
//    initialize call counter
//
//    if leaf block, then
//       for each neighbor
//           if not coarse neighbor
//              pack face data
//           if neighbor level same
//              remote call p_receive_face() on neighbor
//           if neighbor level finer
//              remote call p_receive_face() on neighbor
//              remote call p_receive_face() on neighbor parent
//    else
//       for each face
//          pack face data
//          remote call p_receive_face() on neighbor
//
//    counter = determine_count ()
//
//    local call receive_face(counter)
{
  MG_VERBOSE("send_faces_()");

  const int min_face_rank = 0;
  const int level   = enzo_block->level();
  const Index index = enzo_block->index();

  if (enzo_block->is_leaf()) {

    // face data
    int n; 
    char * array;
    bool lghost[3] = {false,false,false};
    std::vector<int> field_list;
    field_list.push_back(ix_);
    int ic3[3];
    int if3[3];

    // for each neighbor
    ItNeighbor it_neighbor = enzo_block->it_neighbor(min_face_rank,index);

    while (it_neighbor.next()) {
      Index index_neighbor = it_neighbor.index();
      printf ("ItNeighbor %s : %s\n",
	      index.bit_string(4,2).c_str(),
	      index_neighbor.bit_string(4,2).c_str());
      int level_neighbor = index_neighbor.level();

      it_neighbor.child(ic3);
      it_neighbor.face (if3);
      int of3[3] = {if3[0],if3[1],if3[2]};
      FieldFace * field_face;
      int type_op_array = op_array_unknown;

      //  if neighbor level not coarser
      if (! (level_neighbor > level)) {
	// pack face data

	if (level_neighbor == level)   type_op_array = op_array_copy;
	if (level_neighbor == level+1) type_op_array = op_array_prolong;
	if (level_neighbor == level-1) type_op_array = op_array_restrict;

	field_face = enzo_block->load_face (&n, &array,
				 of3, ic3, lghost,
				 type_op_array,
				 field_list);
      } 
      //  if neighbor level same
      if (level_neighbor == level) {

	CProxy_EnzoBlock enzo_block_proxy = 
	  (CProxy_EnzoBlock) enzo_block->proxy_array();
	// remote call p_receive_face() on neighbor
	enzo_block_proxy[index_neighbor].p_mg_receive_face
	  (n,array, refresh_same, of3, ic3);
	
      } 
      //  if neighbor level finer
      else if (level_neighbor > level) {
	CProxy_EnzoBlock enzo_block_proxy = 
	  (CProxy_EnzoBlock) enzo_block->proxy_array();
	// remote call p_receive_face() on neighbor
	enzo_block_proxy[index_neighbor].p_mg_receive_face
	  (n,array, refresh_fine, of3, ic3);
	// remote call p_receive_face() on neighbor parent
	Index index_parent = index_neighbor.index_parent();
	enzo_block_proxy[index_neighbor].p_mg_receive_face
	  (n,array, refresh_same, of3, ic3);
      }

    } 
  } else { // not a leaf
    ItFace it_face = enzo_block->it_face(min_face_rank,index);
    // for each face
    int if3[3];
    while (it_face.next(if3)) {
      Index index_face = it_face.index();
      int level_face = index_face.level();
      printf ("ItFace %s : %s\n",
	      index.bit_string(4,2).c_str(),
	      index_face.bit_string(4,2).c_str());
      // pack face data
      // remote call p_receive_face() on neighbor
    }
  }
  
  //    counter = determine_count ()
  //
  //    local call receive_face(counter)
}

//----------------------------------------------------------------------

int EnzoMethodGravityMg::determine_count_(EnzoBlock * enzo_block) throw()
//    if leaf block, then
//       for each face
//          if non-repeated coarse face
//             increment counter
//    else
//       count number neighbors in same level 
//
//
//    if not leaf block, then
//       count number of children
//
//    count self
{
  int count = 0;

  return count;
}

//----------------------------------------------------------------------

void EnzoBlock::p_mg_receive_face
(int n, char buffer[],  int type_refresh, int if3[3], int ic3[3])
///
///    receive_face(counter)
///
///    unpack ghost data
///    if (number of calls >= counter), then
///       call compute_correction()
{
  printf ("%s p_mg_receive_face %s  face %d %d %d  child %d %d %d\n",
	  name().c_str(),type_refresh == refresh_same ? "same" :
	  (type_refresh == refresh_fine ? "fine" : "coarse"),
	   if3[0],if3[1],if3[2],ic3[0],ic3[1],ic3[2]);
  //  VERBOSE("p_mg_receive_face()");
}

//----------------------------------------------------------------------

void EnzoMethodGravityMg::compute_correction_(EnzoBlock * enzo_block) throw()
//     if coarsest level, then
//         call coarse_solve()
//     else 
//         if not leaf block, then
//             call evaluate_b()
//         call pre_smooth()
{
  MG_VERBOSE("compute_correction_()");
}

//----------------------------------------------------------------------

void EnzoMethodGravityMg::pre_smooth_(EnzoBlock * enzo_block) throw()
{
  MG_VERBOSE("pre_smooth_()");
//      apply the smoother to A X = B
//      compute Y = A*X  [ need refresh? ]
//      pack B, X, Y fields
//      remote call p_restrict() on parent
  smooth_->compute(enzo_block);

  A_->matvec (iy_, ix_, enzo_block);

}

//----------------------------------------------------------------------

void EnzoMethodGravityMg::p_restrict_() throw()
//      unpack vectors B,X,Y
//      call send_faces()
//      if (number of calls >= counter), then
//         call compute_correction()
{
  //  MG_VERBOSE("p_restrict_()");
}

//----------------------------------------------------------------------

void EnzoMethodGravityMg::evaluate_b_(EnzoBlock * enzo_block) throw()
//      B = B - Y
//      Y = A*X [ need refresh? ]
//      B = B + Y
{
  MG_VERBOSE("evaluate_b_()");
}

//----------------------------------------------------------------------

void EnzoMethodGravityMg::solve_coarse_(EnzoBlock * enzo_block) throw()
//     [ need refreshed X?]
//     initiate the coarse solver with p_coarse_solved() callback
//
//     [ NOTE: only coarse blocks call this ]
{
  MG_VERBOSE("solve_coarse_()");
}

//----------------------------------------------------------------------

void EnzoMethodGravityMg::p_coarse_solved_() throw()
//     if leaf, then
//        call end_cycle()
//     else
//        for each child
//            pack correction X
//            remote call p_prolong() on child
{
  //  MG_VERBOSE("p_coarse_solved_()");
}

//----------------------------------------------------------------------

void EnzoMethodGravityMg::p_prolong_() throw()
//     unpack correction Y
//     apply correction to solution X:
//        X = X + Y
//     call post_smooth()  [ need refresh? ]
{
  //  MG_VERBOSE("p_prolong_()");
}

//----------------------------------------------------------------------

void EnzoMethodGravityMg::post_smooth_(EnzoBlock * enzo_block) throw()
//     apply smoother to A X = B
//     if leaf node, then
//        call end_cycle()
//     else
//        call prolong()
{
  MG_VERBOSE("post_smooth_()");
}

//----------------------------------------------------------------------

void EnzoMethodGravityMg::end_cycle_(EnzoBlock * enzo_block) throw()
//     increment iter_
//     call begin_cycle()
{
  MG_VERBOSE("end_cycle_()");
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodGravityMg::exit_solver_ 
( EnzoBlock * enzo_block, int retval ) throw ()
//     deallocate temporaries
//     end_compute() with quiescence
{
  MG_VERBOSE("exit_solver_()");

  enzo_block->clear_refresh();

  Data * data = enzo_block->data();
  Field field = data->field();

  T * X         = (T*) field.values(ix_);
  T * potential = (T*) field.values(ipotential_);

  int mx,my,mz;
  field.dimensions(idensity_,&mx,&my,&mz);

  copy_(potential,X,mx,my,mz,enzo_block->is_leaf());

  FieldDescr * field_descr = field.field_descr();

  EnzoComputeAcceleration compute_acceleration (field_descr,rank_, true,2);

  compute_acceleration.compute(enzo_block);

  monitor_output_ (enzo_block);

  enzo_block->compute_done();
}

//======================================================================

void EnzoMethodGravityMg::monitor_output_(EnzoBlock * enzo_block) throw()
{
  if (enzo_block->index().is_root()) {

  Monitor * monitor = enzo_block->simulation()->monitor();

  monitor->print ("Enzo", "MG iter %04d  rr %g [%g %g]",
		  iter_,
		  (double)(rr_    / rr0_),
		  (double)(rr_min_/ rr0_),
		  (double)(rr_max_/ rr0_));
  }
}

//----------------------------------------------------------------------

