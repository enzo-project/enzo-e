//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Interface routines to HYPRE

/**
 * 
 * @file      hypre.cpp
 * @brief     Implementation of the Hypre class
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
 */

#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <math.h>
#include <string>
#include <vector>
#include <map>

#include "HYPRE_sstruct_ls.h"
#include "HYPRE_krylov.h"

#include "hdf5.h"

#include "newgrav-hypre-solve.h"

//----------------------------------------------------------------------

// Constants

const bool debug = false;
const bool trace = false;

//----------------------------------------------------------------------

// Typedefs

typedef int int3[3];

//----------------------------------------------------------------------

#include "newgrav-mpi.h"
#include "newgrav-scalar.h"
#include "newgrav-error.h"
#include "newgrav-constants.h"
#include "newgrav-point.h"
#include "newgrav-faces.h"
#include "newgrav-domain.h"
#include "newgrav-grid.h"
#include "newgrav-level.h"
#include "newgrav-hierarchy.h"
#include "newgrav-parameters.h"
#include "newgrav-problem.h"
#include "newgrav-hypre.h"
#include "newgrav-error.h"

//======================================================================

// Coefficient for Poisson problem

inline Scalar acoef(Scalar x, Scalar y, Scalar z)
{
  return 1.0;
}

//======================================================================
// PUBLIC MEMBER FUNCTIONS
//======================================================================

/// Hypre constructor

Hypre::Hypre (Hierarchy  & hierarchy,
	      Parameters & parameters)
  : grid_(0),
    graph_(0),
    stencil_(0),
    A_(0),
    B_(0),
    X_(0),
    solver_(0),
    parameters_(&parameters),
    hierarchy_(&hierarchy),
    resid_(-1.0),
    iter_(-1),
    r_factor_(2),
    matrix_scale_(1.0)
{
} // Hypre::Hypre

//----------------------------------------------------------------------

/// Hypre destructor

Hypre::~Hypre ()
{
}

//----------------------------------------------------------------------

/// Initialize the Grid Hierarchy

/** Creates a hypre grid, with one part per level and one box per Grid
    patch object, for an AMR problem.  Sets grid box extents, grid
    part variables, and periodicity. */

void Hypre::init_hierarchy (Mpi & mpi)
{

  int dim       = hierarchy_->dimension();
  int num_parts = hierarchy_->num_levels();

  // Create the hypre grid
  
  HYPRE_SStructGridCreate (MPI_COMM_WORLD, dim, num_parts, &grid_);

  ItHierarchyLevels itl (*hierarchy_);

  while (Level * level = itl++) {

    int part = level->index();

    ItLevelGridsLocal itgl (*level);

    // Set extents for boxes that comprise the hypre grid

    while (Grid * grid = itgl++) {

      int lower[3] = {grid->index_lower(0),
		      grid->index_lower(1),
		      grid->index_lower(2)};
      int upper[3] = {grid->index_upper(0),
		      grid->index_upper(1),
		      grid->index_upper(2)};

      HYPRE_SStructGridSetExtents(grid_, part, lower, upper);
      
    } // while grid = itgl++

    // Create a single cell-centered variable for each grid part (level)

    HYPRE_SStructVariable variable_types[] = { HYPRE_SSTRUCT_VARIABLE_CELL };
    const int numvars = 1;

    HYPRE_SStructGridSetVariables(grid_, part, numvars, variable_types);

    // Set periodicity of the grid part

    int period[3] = { hierarchy_->period_index(0,part),
		      hierarchy_->period_index(1,part),
		      hierarchy_->period_index(2,part) };

    HYPRE_SStructGridSetPeriodic (grid_, part, period);

  } // while level = itl++

  // When finished, assemble the hypre grid

  HYPRE_SStructGridAssemble (grid_);
  
} // Hypre::init_hierarchy()

//----------------------------------------------------------------------

/// Initialize the discretization stencils.  

/** Creates and initializes a hypre stencil object. */

void Hypre::init_stencil ()

{

  int dim = hierarchy_->dimension();

  HYPRE_SStructStencilCreate (dim,dim*2+1,&stencil_);

  int entries[][3] = { {  0, 0, 0 },     // center
		       {  1, 0, 0 },     // X+
		       { -1, 0, 0 },     // X-
		       {  0, 1, 0 },     // Y+
		       {  0,-1, 0 },     // Y-
		       {  0, 0, 1 },     // Z+
		       {  0, 0,-1 } };   // Z-

  if (dim >= 1) HYPRE_SStructStencilSetEntry (stencil_, 0, entries[0], 0);
  if (dim >= 1) HYPRE_SStructStencilSetEntry (stencil_, 1, entries[1], 0);
  if (dim >= 1) HYPRE_SStructStencilSetEntry (stencil_, 2, entries[2], 0);
  if (dim >= 2) HYPRE_SStructStencilSetEntry (stencil_, 3, entries[3], 0);
  if (dim >= 2) HYPRE_SStructStencilSetEntry (stencil_, 4, entries[4], 0);
  if (dim >= 3) HYPRE_SStructStencilSetEntry (stencil_, 5, entries[5], 0);
  if (dim >= 3) HYPRE_SStructStencilSetEntry (stencil_, 6, entries[6], 0);

} // Hypre::init_stencil()

//----------------------------------------------------------------------

/// Initialize the graph.

/** Creates a graph containing the matrix nonzero structure.  Graph
    edges include both those for nonzeros from the stencil within each
    part (level), and nonzeros for graph entries connecting linked
    parts.  The matrix nonzero structure is generally nonsymmetric.
    Only the stencil step is required for unigrid problems.
*/

void Hypre::init_graph ()

{
  // Create the hypre graph object

  HYPRE_SStructGraphCreate (MPI_COMM_WORLD, grid_, &graph_);
  
  HYPRE_SStructGraphSetObjectType (graph_, HYPRE_SSTRUCT);

  ItHierarchyLevels itl (*hierarchy_);

  while (Level * level = itl++) {

    int part = level->index();

    // Define stencil connections within each level

    HYPRE_SStructGraphSetStencil (graph_, part, 0, stencil_);

    // Define graph connections between levels

    if (part > 0) {

      ItLevelGridsAll itag (*level);

      while (Grid * grid = itag++) {

	init_graph_nonstencil_(*grid);

      } // while grid = itag++

    } // if part > 0

    ItLevelGridsAll itag (*level);

    while (Grid * grid = itag++) {

      // Initialize face counters for subsequent matrix inter-level entries

      int dim = hierarchy_->dimension();
      grid->init_counter(dim*2+1);

    } // while grid = itag++
  } // while level = itl++

  // Assemble the hypre graph

  HYPRE_SStructGraphAssemble (graph_);

} // Hypre::init_graph()


//----------------------------------------------------------------------

/// Initialize the matrix A and right-hand-side vector b

/** Creates a matrix with a given nonzero structure, and sets nonzero
    values.
*/

void Hypre::init_elements (std::vector<Point *>  points)

{
  // Create the hypre matrix A_, solution X_, and right-hand side B_ objects

  HYPRE_SStructMatrixCreate (MPI_COMM_WORLD, graph_, &A_);
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid_,  &X_);
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid_,  &B_);

  // Set the object types

  if (parameters_->value("solver") == "bicgstab-boomer") {
    HYPRE_SStructMatrixSetObjectType (A_,HYPRE_PARCSR);
    HYPRE_SStructVectorSetObjectType (X_,HYPRE_PARCSR);
    HYPRE_SStructVectorSetObjectType (B_,HYPRE_PARCSR);
  } else {
    HYPRE_SStructMatrixSetObjectType (A_,HYPRE_SSTRUCT);
    HYPRE_SStructVectorSetObjectType (X_,HYPRE_SSTRUCT);
    HYPRE_SStructVectorSetObjectType (B_,HYPRE_SSTRUCT);
  }

  // Initialize the hypre matrix and vector objects

  HYPRE_SStructMatrixInitialize (A_);
  HYPRE_SStructVectorInitialize (X_);
  HYPRE_SStructVectorInitialize (B_);

  //--------------------------------------------------
  // Initialize the matrix A_ elements
  //--------------------------------------------------

  init_elements_matrix_ ();

  //--------------------------------------------------
  // Initialize B_ elements 
  //--------------------------------------------------

  init_elements_rhs_ (points);

  // Assemble the matrix and vectors

  HYPRE_SStructMatrixAssemble (A_);
  HYPRE_SStructVectorAssemble (B_);
  HYPRE_SStructVectorAssemble (X_);

  // Optionally write the matrix and vectors to a file for debugging

  if (parameters_->value("dump_hypre") == "true") {
    HYPRE_SStructMatrixPrint ("A",A_,0);
    HYPRE_SStructVectorPrint ("B",B_,0);
    HYPRE_SStructVectorPrint ("X0",X_,0);
  }

} // Hypre::init_elements()

//----------------------------------------------------------------------

/// Initialize and solve the linear solver

void Hypre::solve ()

{
  std::string solver = parameters_->value("solver");
  int         levels = hierarchy_->num_levels();

  int    itmax  = 0;
  double restol = 0.0;

  // Check solver parameters
  std::string sitmax  = parameters_->value("solver_itmax");
  std::string srestol = parameters_->value("solver_restol");

  // If not defined, then define them
  if (sitmax == "")  parameters_->add_parameter ("solver_itmax","200");
  if (srestol == "") parameters_->add_parameter ("solver_restol","1e-6");

  // Set local variables
  itmax  = atoi(sitmax.c_str());
  restol = atof(srestol.c_str());

  if        (solver == "pfmg" && levels == 1) {

    solve_pfmg_(itmax,restol);

  } else if (solver == "fac"  && levels > 1) {

    solve_fac_(itmax,restol);

  } else if (solver == "bicgstab") {

    solve_bicgstab_(itmax,restol);

  } else if (solver == "bicgstab-boomer") {

    solve_bicgstab_boomer_(itmax,restol);

  } else if (solver == "gmres") {

    solve_gmres_(itmax,restol);

  } else {

    char error_message[100];
    sprintf (error_message, "Hypre::solve called with illegal combination of "
	     "solver %s on %d levels", solver.c_str(),levels);
    ERROR(error_message);
  }
  
  if (parameters_->value("dump_hypre") == "true") {
    HYPRE_SStructVectorPrint ("X",X_,1);
  }

} // Hypre::solve()

//----------------------------------------------------------------------

/// Evaluate the success of the solve

void Hypre::evaluate ()

{

  char filename[80];

  ItHierarchyGridsLocal itg(*hierarchy_);
  while (Grid * grid = itg++) {

    int level = grid->level();

    int lower[3],upper[3];
    grid->get_limits(lower,upper);

    grid->allocate();

    HYPRE_SStructVectorGetBoxValues (X_,level,lower,upper,0,grid->values());  

    sprintf (filename,"X.%d",grid->id());
    grid->write(filename);
    
    HYPRE_SStructVectorGetBoxValues (B_,level,lower,upper,0,grid->values());  

    sprintf (filename,"B.%d",grid->id());
    grid->write(filename);

  } // grid = itg++


  // --------------------------------------------------
  // Evaluate ||X|| and sum(X) / ||X|| and display result
  // --------------------------------------------------

  Scalar xsum_local  = 0.0;
  Scalar x2sum_local = 0.0;

  while (Grid * grid = itg++) {

    int level = grid->level();

    int lower[3],upper[3];

    grid->get_limits(lower,upper);

    grid->allocate();

    HYPRE_SStructVectorGetBoxValues (X_,level,lower,upper,0,grid->values());  

    Scalar * x = grid->values();

    int n3[3];

    grid->get_size(n3);

    for (int i2=0; i2<n3[2]; i2++) {
      for (int i1=0; i1<n3[1]; i1++) {
	for (int i0=0; i0<n3[0]; i0++) {
	  Scalar xval = x[grid->index(i0,i1,i2,n3[0],n3[1],n3[2])];
	  xsum_local  += xval;
	  x2sum_local += xval*xval;
	}
      }
    }
  }

  Scalar xsum = 0.0;
  MPI_Allreduce (&xsum_local, &xsum, 1, 
		 MPI_SCALAR, MPI_SUM, MPI_COMM_WORLD);

  Scalar x2sum = 0.0;
  MPI_Allreduce (&x2sum_local, &x2sum, 1, 
		 MPI_SCALAR, MPI_SUM, MPI_COMM_WORLD);


  if (pmpi->is_root()) { 

    bool success = true;

    // Residual too high

    double restol = atof(parameters_->value("solver_restol").c_str());

    if (resid_ > restol) {
      printf ("Diverged: %g > %g\n", resid_,restol);
      success = false;
    }

    // Iterations reached limit

    int itmax     = atoi(parameters_->value("solver_itmax").c_str());
    if (iterations() >= itmax) {
      printf ("Stalled: %d >= %d\n", iterations(),itmax);
      success = false;
    }

    // Appears to have completed successfully

    if (success) {
      printf ("Success!\n"); fflush(stdout); 
    }
  }
  if (pmpi->is_root()) printf ("norm(X)        = %g\n",sqrt(x2sum));
  if (pmpi->is_root()) printf ("sum(X)/norm(X) = %g\n",xsum/sqrt(x2sum));

} // Hypre::evaluate()

//======================================================================
// PRIVATE MEMBER FUNCTIONS
//======================================================================

/// init_nonstencil_() is called twice: first by
/// init_graph_nonstencil_() with phase == phase_graph to set nonstencil
/// graph entries, and again by init_matrix_nonstencil_() with phase
/// == phase_matrix to set nonstencil matrix entries.

void Hypre::init_nonstencil_ (Grid & grid, phase_enum phase)
{

  // Input parameter check

  if ( ! (phase == phase_graph || phase == phase_matrix) ) {
    char error_message[80];
    sprintf (error_message,"init_matrix_nonstencil_ called with phase = %d",
	     int(phase));
    ERROR(error_message);
  } // if phase unexpected

  // global grid index limits

  int index_global[3][2];
  grid.indices(index_global);
  
  // Determine the discretization method (constant, linear, quadratic)
  // Currently only constant is implemented

  enum discret_type_enum {discret_type_unknown, discret_type_const};
  discret_type_enum discret_type;

  if (parameters_->value("discret") == "constant") 
    discret_type = discret_type_const;

  // Loop over each face zone in the grid, adding non-stencil graph
  // entries wherever a zone is adjacent to a coarse zone.  Both
  // fine-to-coarse and coarse-to-fine entries are added.

  int level_fine   = grid.level();
  int level_coarse = grid.level() - 1;
  assert (level_coarse >= 0);

  for (int axis=0; axis<3; axis++) {

    // axis0:       axis normal to face
    // axis1,axis2: axes within face

    int axis0 = axis;
    int axis1 = (axis+1)%3;
    int axis2 = (axis+2)%3;

    // n0     grid size normal to face
    // n1,n2: grid size within face

    int n0 = index_global[axis0][1] - index_global[axis0][0];
    int n1 = index_global[axis1][1] - index_global[axis1][0];
    int n2 = index_global[axis2][1] - index_global[axis2][0];

    // index_global[][] should be divisible by r_factor_**level.  Just
    // test r_factor_ here.

    bool l0 = (index_global[axis1][0]/r_factor_)*r_factor_ == index_global[axis1][0];
    bool l1 = (index_global[axis1][1]/r_factor_)*r_factor_ == index_global[axis1][1];

    if (!l0) printf ("index_global[%d][0] = %d\n",axis1,index_global[axis1][0]);
    assert (l0);
    if (!l1) printf ("index_global[%d][1] = %d\n",axis1,index_global[axis1][1]);
    assert (l1);

    for (int face=0; face<2; face++) {

      // Loop over face zones that are aligned with coarse zones (hence "+= r")

      for (int index1=0; index1<n1; index1 += r_factor_) {
	for (int index2=0; index2<n2; index2 += r_factor_) {

	  Grid * adjacent   = grid.faces().adjacent(axis,face,index1,index2);

	  Faces::Label & fz = grid.faces().label(axis,face,index1,index2);

	  // Add graph entries iff grid or adjacent grid is local, and if
	  // adjacent grid (if it exists) is in the next-coarser level

	  bool is_local =
	    (adjacent != NULL) &&  (adjacent->is_local() || grid.is_local());

	  bool is_coarse = fz == Faces::_coarse_;

	  if (is_local && is_coarse) {

	    // (fine) grid global indices

	    int index_fine[3]; 

	    index_fine[axis0] = index_global[axis0][0] + face*(n0 - r_factor_);
	    index_fine[axis1] = index_global[axis1][0] + index1;
	    index_fine[axis2] = index_global[axis2][0] + index2;

	    // (coarse) adjacent global indices

	    int index_coarse[3]; 

	    index_coarse[axis0] = (index_fine[axis0]) / r_factor_  + (face*r_factor_-1);
	    index_coarse[axis1] = (index_fine[axis1]) / r_factor_;
	    index_coarse[axis2] = (index_fine[axis2]) / r_factor_;

	    // adjust for periodicity

	    if (hierarchy_->is_periodic(axis0)) {
	      int period = hierarchy_->period_index(axis0,level_coarse);
 	      index_coarse[axis0] = (index_coarse[axis0] + period) % period;
 	    }

	    //--------------------------------------------------
	    // GRAPH ENTRY: FINE-TO-COARSE 
	    //--------------------------------------------------

	    if (discret_type == discret_type_const) {

	      update_fine_coarse_const_(face,grid,axis0,phase,
					level_fine,level_coarse,
					index_fine,index_coarse);

	      //--------------------------------------------------
	      // GRAPH ENTRY: COARSE-TO-FINE
	      //--------------------------------------------------
	      
	      if (adjacent->is_local()) {

		update_coarse_fine_const_(face,*adjacent,axis0,phase,
					  level_fine,level_coarse,
					  index_fine,index_coarse);

	      }

	    } else {
	      char error_message[80];
	      strcpy (error_message,"Unknown parameter discret = ");
	      strcat (error_message,parameters_->value("discret").c_str());
	      ERROR(error_message);
	    } // if discret unexpected
	  } // if is_local && fz == Faces::_coarse_
	} // for index2
      } // for index1
    } // for face
  } // for axis
} // Hypre::init_nonstencil_()

//------------------------------------------------------------------------

/// Initialize matrix stencil and graph entries

void Hypre::init_elements_matrix_ ()

{

  ItHierarchyLevels itl (*hierarchy_);

  while (Level * level = itl++) {

    int part = level->index();

    ItLevelGridsLocal itlg (*level);
    
    // 1. Set stencil values within level

    while (Grid * grid = itlg++) {

      init_matrix_stencil_(*grid);

    } // while grid = itlg++

    if (part > 0) {

      // *** WARNING: POSSIBLE SCALING ISSUE.  Below we loop over all
      // *** grids; however, we only need to loop over "parent-child
      // *** pairs such that either child or parent is local to this
      // *** MPI process."
 
      // Set matrix values between levels

      ItLevelGridsAll itag (*level);

      while (Grid * grid = itag++) {

 	init_matrix_nonstencil_(*grid);

      } // while grid = itag++

    } // while level > 0

  } // while level = itl++

  for (int part = 1; part < hierarchy_->num_levels(); part++) {

    // Clean up stencil connections between levels

    init_matrix_clear_(part);

  } // for part

} // init_elements_matrix_()

//------------------------------------------------------------------------

/// Set right-hand-side elements

void Hypre::init_elements_rhs_ (std::vector<Point *>  & points)
{

  Scalar local_shift_b_sum = 0.0;
  long long shift_b_count  = 0.0;

  bool enzo_density = parameters_->value("enzo_density") == "true";

  if (enzo_density) {

    // Either set density according to Enzo grid files...

    bool        enzo_packed = parameters_->value("enzo_packed") == "true";
    std::string enzo_prefix = parameters_->value("enzo_prefix");
    local_shift_b_sum += init_vector_density_ (enzo_prefix,enzo_packed);

  } else {

    // ...or set density according to point masses

    local_shift_b_sum += init_vector_points_  (points);

  }


  if ( parameters_->value("boundary") == "periodic" ) {

    // Clear under overlapped areas

    int parts[hierarchy_->num_levels()];
    int3 *refinements = new int3 [hierarchy_->num_levels()];
    for (int part = 0; part < hierarchy_->num_levels(); part++) {
      parts[part] = part;
      refinements[part][0] = r_factor_;
      refinements[part][1] = r_factor_;
      refinements[part][2] = r_factor_;
    }
    HYPRE_SStructFACZeroAMRVectorData (B_, parts, refinements);

    // Accumulate local sums

    local_shift_b_sum  = 0.0;

    // Compute total number of variables (excluding overlap)

    shift_b_count = 0; 
    const int r_factor3 = r_factor_*r_factor_*r_factor_;

    ItHierarchyLevels itl (*hierarchy_);
    while (Level * level = itl++) {
      
      int part = level->index();
      ItLevelGridsAll itgl (*level);
      while (Grid * grid = itgl++) {

	// Adjust shift_b_count: add grid; subtract overlap

	shift_b_count += grid->num_unknowns();
	if (part > 0) shift_b_count -= grid->num_unknowns() / r_factor3;

	if (grid->is_local()) {

 	  // Get grid info

 	  int lower[3],upper[3];
	  grid->get_limits(lower, upper);
	  int n = grid->num_unknowns();

	  // Create space for the patch

 	  double * values = new double [n];
 	  for (int i=0; i<n; i++) values[i] = 0.0;

 	  // Copy vector values to the array
 	  HYPRE_SStructVectorGetBoxValues (B_,part,lower,upper,0,values);  

	  // Accumulate the sum
	  for (int i=0; i<n; i++) local_shift_b_sum += values[i];

	  // Delete the zeroed values
 	  delete [] values;
	}
      }
    }

    // Get global sum from local sums

    Scalar shift_b_sum = 0.0;

    MPI_Allreduce (&local_shift_b_sum, &shift_b_sum, 1, 
		   MPI_SCALAR, MPI_SUM, MPI_COMM_WORLD);

    // Compute the shift given the global sum and global count

    Scalar shift_b_amount = - shift_b_sum / shift_b_count;

    // Perform the shift of B 
    
    while (Level * level = itl++) {

      int part = level->index();

      ItLevelGridsLocal itg (*level);
      while (Grid * grid = itg++) {

	int lower[3],upper[3];
	grid->get_limits(lower, upper);
	Scalar * values = new Scalar[grid->num_unknowns()];

	for (int i=0; i<grid->num_unknowns(); i++) values[i] = shift_b_amount;

	HYPRE_SStructVectorAddToBoxValues (B_,part,lower,upper,0,values);

	delete [] values;

      } // while grid = itg++

    } // while level = itl++

    // Re-clear under overlapped area that got shifted

    HYPRE_SStructFACZeroAMRVectorData (B_, parts, refinements);

  } // if periodic

} // init_elements_rhs_()

//------------------------------------------------------------------------

/// Set matrix stencil values for the grid interior

void Hypre::init_matrix_stencil_ (Grid & grid)

{
  int n          = grid.num_unknowns();
  int entries[7] = { 0,1,2,3,4,5,6 };
  double h3[3]   = {grid.h(0),grid.h(1),grid.h(2)};
  int    n3[3]   = {grid.n(0),grid.n(1),grid.n(2)};

  double h120 = h3[1]*h3[2] / h3[0];
  double h201 = h3[2]*h3[0] / h3[1];
  double h012 = h3[0]*h3[1] / h3[2];

  double * v0;         // Diagonal elements
  double * v1[3][2];   // Off-diagonal elements

  // Allocate storage

  v0 = new double [n];
  for (int axis=0; axis<3; axis++) {
    for (int face=0; face<2; face++) {
      v1[axis][face] = new double [n];
    } // for face
  } // for axis

  //-----------------------------------------------------------
  // Set stencil for all unknowns, ignoring boundary conditions
  //-----------------------------------------------------------

  double lower[3];
  grid.x_lower(lower[0],lower[1],lower[2]);

  int i0,i1,i2,i;
  for (i2 = 0; i2 < n3[2]; i2++) {
    for (i1 = 0; i1 < n3[1]; i1++) {
      for (i0 = 0; i0 < n3[0]; i0++) {

	// DIFFUSION COEFFICIENTS HERE
	// (NOTE: inefficient due to recomputing)

	// Compute position of the zone center

	double x0 = lower[0]+(i0 + 0.5)*h3[0];
	double x1 = lower[1]+(i1 + 0.5)*h3[1];
	double x2 = lower[2]+(i2 + 0.5)*h3[2];

	double axm = acoef(x0-0.5*h3[0],x1,x2);
	double axp = acoef(x0+0.5*h3[0],x1,x2);
	double aym = acoef(x0,x1-0.5*h3[1],x2);
	double ayp = acoef(x0,x1+0.5*h3[1],x2);
	double azm = acoef(x0,x1,x2-0.5*h3[2]);
	double azp = acoef(x0,x1,x2+0.5*h3[2]);

	i = Grid::index(i0,i1,i2,n3[0],n3[1],n3[2]);

	v1[0][0][i] = matrix_scale_ * h120 * axm;
	v1[0][1][i] = matrix_scale_ * h120 * axp;
	v1[1][0][i] = matrix_scale_ * h201 * aym;
	v1[1][1][i] = matrix_scale_ * h201 * ayp;
	v1[2][0][i] = matrix_scale_ * h012 * azm;
	v1[2][1][i] = matrix_scale_ * h012 * azp;

	v0[i] = -( v1[0][0][i] + v1[0][1][i] +
		   v1[1][0][i] + v1[1][1][i] + 
		   v1[2][0][i] + v1[2][1][i] );

      } // for i0
    } // for i1
  } // for i2

  //-----------------------------------------------------------
  // Adjust stencil at grid boundaries
  //-----------------------------------------------------------

  //   Faces & faces = grid.faces();

  int level = grid.level();

  int index_lower[3] = { grid.index_lower(0), 
			 grid.index_lower(1), 
			 grid.index_lower(2) };
  int index_upper[3] = { grid.index_upper(0), 
			 grid.index_upper(1), 
			 grid.index_upper(2) };

  // Update matrix stencil values

  HYPRE_SStructMatrixSetBoxValues 
    (A_,level,index_lower,index_upper,0,1,&entries[0],v0);
  HYPRE_SStructMatrixSetBoxValues 
    (A_,level,index_lower,index_upper,0,1,&entries[1],v1[0][1]);
  HYPRE_SStructMatrixSetBoxValues 
    (A_,level,index_lower,index_upper,0,1,&entries[2],v1[0][0]);
  HYPRE_SStructMatrixSetBoxValues 
    (A_,level,index_lower,index_upper,0,1,&entries[3],v1[1][1]);
  HYPRE_SStructMatrixSetBoxValues 
    (A_,level,index_lower,index_upper,0,1,&entries[4],v1[1][0]);
  HYPRE_SStructMatrixSetBoxValues 
    (A_,level,index_lower,index_upper,0,1,&entries[5],v1[2][1]);
  HYPRE_SStructMatrixSetBoxValues 
    (A_,level,index_lower,index_upper,0,1,&entries[6],v1[2][0]);

  // Deallocate arrays

  delete [] v0;
  for (int axis=0; axis<3; axis++) {
    for (int face=0; face<2; face++) {
      delete [] v1[axis][face];
    } // for face
  } // for axis

} // Hypre::init_matrix_stencil_()

//------------------------------------------------------------------------

/// Clean up stencil connections between parts for FAC solver

void Hypre::init_matrix_clear_ (int part)
{

  if (part > 0) {

    int r_factors[3] = {r_factor_,r_factor_,r_factor_}; 

    HYPRE_SStructFACZeroAMRMatrixData (A_, part-1, r_factors);

  }

} // Hypre::init_matrix_clear_()

//------------------------------------------------------------------------

/// Add contributions from point sources to right-hand side B

Scalar Hypre::init_vector_points_ (std::vector<Point *> & points)

{

  const Scalar scaling0 = -4.0*Constants::G()*Constants::pi();

  Scalar shift_b_sum = 0.0;

  int i;
  for (i=0; i<int(points.size()); i++) {

    Point & point      = *points[i];
    Grid & grid        = hierarchy_->return_grid(point.igrid());
    if (grid.is_local()) {

      Scalar cell_volume = grid.h(0) * grid.h(1) * grid.h(2);
      Scalar density     = point.mass() / cell_volume;
      Scalar value       = scaling0 * density;

      // Add contribution of the point to the right-hand side vector

      int index[3];
      Scalar lower[3],upper[3];
      grid.x_lower(lower[0],lower[1],lower[2]);
      grid.x_upper(upper[0],upper[1],upper[2]);
      for (int k=0; k<3; k++) {
	Scalar ap = point.x(k)      - lower[k];
	Scalar ag = upper[k] - lower[k];
	int    ig = grid.num_unknowns(k);
	int    i0 = grid.index_lower(k);
	index[k] = int (ap/ag*ig) + i0;
      } // for k=0,2
      if (index[0] < grid.index_lower(0) || grid.index_upper(0) < index[0] ||
	  index[1] < grid.index_lower(1) || grid.index_upper(1) < index[1] ||
	  index[2] < grid.index_lower(2) || grid.index_upper(2) < index[2]) {
	printf ("WARNING: Point apparently not in grid: \n");
	printf ("WARNING:    Point: (%g,%g,%g)\n",point.x(0),point.x(1),point.x(2));
	printf ("WARNING:    Grid:  (%g,%g,%g) - (%g,%g,%g)\n",
		lower[0],lower[1],lower[2],upper[0],upper[1],upper[2]);
      } // if index not in grid
      if (debug) {
	point.print();
	grid.print();
	printf ("Point index  = %d %d %d)\n",index[0],index[1],index[2]);
	printf ("Cell size    = %g %g %g\n",grid.h(0),grid.h(1),grid.h(2));
	printf ("Cell volume  = %g\n",cell_volume);
	printf ("Cell density = %g\n",density);
	printf ("RHS contribution = %g\n",value);
      } // if debug
    
      shift_b_sum += value;

      HYPRE_SStructVectorAddToValues (B_, grid.level(), index, 0, &value);

    } // if grid.is_local()
  } // for i=0 to # points

  return shift_b_sum;
} // Hypre::init_vector_points_()

//------------------------------------------------------------------------

/// Add contributions from Density in enzo HDF5 files to right-hand side B

Scalar Hypre::init_vector_density_ (std::string             file_prefix,
				    bool                    enzo_packed)

{
  ItHierarchyGridsLocal itg (*hierarchy_);
  char error_message[80];

  herr_t status;
  hid_t  file_id;
  hid_t  dataset_id;

  while (Grid * grid = itg++) {

    // Open the HDF5 grid file

    char grid_num_str[10];
    sprintf (grid_num_str,"%04d",grid->id() + 1);
    std::string file_name = file_prefix + ".grid" + grid_num_str;
    file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
      strcpy (error_message,"H5Fopen cannot open file ");
      strcat (error_message,file_name.c_str());
      ERROR(error_message);
    } else {
      printf ("DEBUG %s:%d %s opened successfully\n",
	      __FILE__,__LINE__,file_name.c_str());
    }

    // Open the dataset Density

    dataset_id = H5Dopen(file_id, "Density");
    if (dataset_id < 0) {
      strcpy (error_message,"H5Dopen cannot open dataset Density");
      ERROR(error_message);
    } else {
      printf ("DEBUG %s:%d Density opened successfully\n",
	      __FILE__,__LINE__);
    }

    // Read the dataset

    double       *values = new double [grid->n()];

    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
		     H5P_DEFAULT, values);

    if (status < 0) {
      strcpy (error_message,"H5Dread exited with status ");
      char status_str[10];
      sprintf (status_str,"%d",status);
      strcat (error_message,status_str);
      ERROR(error_message);
    } else {
      printf ("DEBUG %s:%d Density read successfully\n",
	      __FILE__,__LINE__);
      printf ("%d %d %d  %d  [%g %g]\n",
	      grid->n(0),grid->n(1),grid->n(2),grid->n(),
	      values[0],values[grid->n()-1]);
    }

    // Copy the values to the hypre vector

    int part = grid->level();
    int lower[3] = { grid->index_lower(0), 
		     grid->index_lower(1), 
		     grid->index_lower(2) };
    int upper[3] = { grid->index_upper(0), 
		     grid->index_upper(1), 
		     grid->index_upper(2) };
    HYPRE_SStructVectorAddToBoxValues (B_,part,lower,upper,0,values);

    delete [] values;

    // Close the HDF5 grid file

    status = H5Fclose(file_id);

    if (status < 0) {
      strcpy (error_message,"H5Fclose exited with status ");
      char status_str[10];
      sprintf (status_str,"%d",status);
      strcat (error_message,status_str);
      ERROR(error_message);
    }
    
  }

  return 0.0;
}

//------------------------------------------------------------------------

/// Initialize the PFMG hypre solver

void Hypre::solve_pfmg_ (int itmax, double restol)

{

  // Create and initialize the solver

  HYPRE_SStructSysPFMGCreate    (MPI_COMM_WORLD, &solver_);

  // stopping criteria

  if (itmax != 0 )   HYPRE_SStructSysPFMGSetMaxIter(solver_,itmax);
  if (restol != 0.0) HYPRE_SStructSysPFMGSetTol(solver_,    restol);

  HYPRE_SStructSysPFMGSetLogging(solver_, 1);
  HYPRE_SStructSysPFMGSetup     (solver_,A_,B_,X_);

  // Solve the linear system

  HYPRE_SStructSysPFMGSolve     (solver_,A_,B_,X_);

  // Write out some diagnostic info about the solve

  HYPRE_SStructSysPFMGGetNumIterations (solver_,&iter_);
  HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm (solver_,&resid_);

  printf ("HYPRE_SStructSysPFMGSolve num iterations: %d\n",iter_);
  printf ("HYPRE_SStructSysPFMGSolve final relative residual norm: %g\n",resid_);

  // Delete the solver

  HYPRE_SStructSysPFMGDestroy (solver_);

} // Hypre::solve_pfmg_()

//------------------------------------------------------------------------

/// Initialize the FAC hypre solver

void Hypre::solve_fac_ (int itmax, double restol)

{
  int i;

  // Create the solver

  HYPRE_SStructFACCreate(MPI_COMM_WORLD, &solver_);

  // Initialize parts

  int num_parts = hierarchy_->num_levels();

  HYPRE_SStructFACSetMaxLevels(solver_,  num_parts);

  int *parts  = new int [num_parts];
  for (i=0; i<num_parts; i++) parts[i] = i;

  HYPRE_SStructFACSetPLevels(solver_, num_parts, parts);

  // Initialize refinement factors

  int3 *refinements = new int3 [num_parts];
  
  for (i=0; i<num_parts; i++) {
    refinements[i][0] = r_factor_;
    refinements[i][1] = r_factor_;
    refinements[i][2] = r_factor_;
  }

  HYPRE_SStructFACSetPRefinements(solver_, num_parts, refinements);

  // solver parameters

  int npre   = 2;
  int npost  = 2;
  int csolve = 1;
  int relax  = 2;

  HYPRE_SStructFACSetNumPreRelax(solver_,      npre);
  HYPRE_SStructFACSetNumPostRelax(solver_,     npost);
  HYPRE_SStructFACSetCoarseSolverType(solver_, csolve);
  HYPRE_SStructFACSetRelaxType(solver_,        relax);

  // stopping criteria

  if (itmax != 0 ) {
    HYPRE_SStructFACSetMaxIter(solver_,itmax);
  }
  if (restol != 0.0) {
    HYPRE_SStructFACSetTol(solver_,    restol);
  }

  // output amount

  HYPRE_SStructFACSetLogging(solver_, 1);

  // prepare for solve

  HYPRE_SStructFACSetup2(solver_, A_, B_, X_);

  // Solve the linear system

  HYPRE_SStructFACSolve3(solver_, A_, B_, X_);

  // Write out some diagnostic info about the solve

  HYPRE_SStructFACGetNumIterations(solver_, &iter_);
  HYPRE_SStructFACGetFinalRelativeResidualNorm(solver_, &resid_);

  printf ("HYPRE_SStructFACSolve num iterations: %d\n",iter_);
  printf ("HYPRE_SStructFACSolve final relative residual norm: %g\n",resid_);

  // Delete the solver

  HYPRE_SStructFACDestroy2(solver_);

  // Delete local dynamic storage

  delete [] parts;
  delete [] refinements;

} // Hypre::solve_fac_()

//------------------------------------------------------------------------

/// Initialize the BICGSTAB hypre solver

void Hypre::solve_bicgstab_ (int itmax, double restol)

{

  // Create the solver

  HYPRE_SStructBiCGSTABCreate(MPI_COMM_WORLD, &solver_);

  // stopping criteria

  if (itmax != 0 )   HYPRE_SStructBiCGSTABSetMaxIter(solver_,itmax);
  if (restol != 0.0) HYPRE_SStructBiCGSTABSetTol(solver_,    restol);

  // output amount

  HYPRE_SStructBiCGSTABSetLogging(solver_, 1);

  // Initialize the solver

  HYPRE_SStructBiCGSTABSetup(solver_, A_, B_, X_);

  // Solve the linear system

  HYPRE_SStructBiCGSTABSolve(solver_, A_, B_, X_);

  // Write out some diagnostic info about the solve

  HYPRE_SStructBiCGSTABGetNumIterations(solver_, &iter_);
  HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(solver_, &resid_);

  if (debug) printf ("HYPRE_SStructBiCGSTABSolve num iterations: %d\n",iter_);
  if (debug) printf ("HYPRE_SStructBiCGSTABSolve final relative residual norm: %g\n",resid_);


  // Delete the solver

  HYPRE_SStructBiCGSTABDestroy(solver_);

} // Hypre::solve_bicgstab_()

//------------------------------------------------------------------------

/// Initialize the BiCGSTAB hypre solver with BoomerAMG preconditioning

void Hypre::solve_bicgstab_boomer_ (int itmax, double restol)

{

  HYPRE_Solver solver;

  // Create the solver

  HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &solver);

  // stopping criteria

  if (itmax != 0 )   HYPRE_BiCGSTABSetMaxIter(solver,itmax);
  if (restol != 0.0) HYPRE_BiCGSTABSetTol(solver,    restol);

  // output amount

  HYPRE_BiCGSTABSetLogging(solver, 1);

  // Set BoomerAMG preconditioner

  HYPRE_Solver par_precond;

  HYPRE_BoomerAMGCreate(&par_precond);
  HYPRE_BoomerAMGSetCoarsenType(par_precond, 6);
  HYPRE_BoomerAMGSetStrongThreshold(par_precond, 0.25);
  HYPRE_BoomerAMGSetTol(par_precond, 0.0);
  HYPRE_BoomerAMGSetPrintLevel(par_precond, 1);
  HYPRE_BoomerAMGSetPrintFileName(par_precond, "sstruct.out.log");
  HYPRE_BoomerAMGSetMaxIter(par_precond, 1);

  HYPRE_BiCGSTABSetPrecond (solver,
			    (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
			    (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
			    par_precond);

  // Initialize the solver

  HYPRE_BiCGSTABSetup(solver, 
		      (HYPRE_Matrix) A_, 
		      (HYPRE_Vector) B_, 
		      (HYPRE_Vector) X_);

  // Solve the linear system

  HYPRE_BiCGSTABSolve(solver, (HYPRE_Matrix) A_, (HYPRE_Vector)B_,(HYPRE_Vector) X_);

  // Write out some diagnostic info about the solve

  HYPRE_BiCGSTABGetNumIterations(solver, &iter_);
  HYPRE_BiCGSTABGetFinalRelativeResidualNorm(solver, &resid_);

  if (debug) printf ("HYPRE_BiCGSTABSolve num iterations: %d\n",iter_);
  if (debug) printf ("HYPRE_BiCGSTABSolve final relative residual norm: %g\n",resid_);


  // Delete the solver

  HYPRE_BiCGSTABDestroy(solver);

} // Hypre::solve_bicgstab_boomer_()

//------------------------------------------------------------------------

/// Initialize the GMRES hypre solver

void Hypre::solve_gmres_ (int itmax, double restol)

{

  // Create the solver

  HYPRE_SStructGMRESCreate(MPI_COMM_WORLD, &solver_);

  // stopping criteria

  if (itmax != 0 )   HYPRE_SStructGMRESSetMaxIter(solver_,itmax);
  if (restol != 0.0) HYPRE_SStructGMRESSetTol(solver_,    restol);

  // output amount

  HYPRE_SStructGMRESSetLogging(solver_, 1);

  // Initialize the solver

  HYPRE_SStructGMRESSetup(solver_, A_, B_, X_);

  // Solve the linear system

  HYPRE_SStructGMRESSolve(solver_, A_, B_, X_);

  // Write out some diagnostic info about the solve

  HYPRE_SStructGMRESGetNumIterations(solver_, &iter_);
  HYPRE_SStructGMRESGetFinalRelativeResidualNorm(solver_, &resid_);

  if (debug) printf ("HYPRE_SStructGMRESSolve num iterations: %d\n",iter_);
  if (debug) printf ("HYPRE_SStructGMRESSolve final relative residual norm: %g\n",resid_);


  // Delete the solver

  HYPRE_SStructGMRESDestroy(solver_);

} // Hypre::solve_gmres_()

//------------------------------------------------------------------------

/// Update the matrix at a coarse face zone adjacent to fine face
/// zones using a simple piecewise-constant finite volume
/// discretization.  Called in two phases, with phase == phase_graph
/// (via init_graph_nonstencil_()) for the nonzero structure, and with
/// phase == phase_matrix (via init_matrix_nonstencil_()) for the
/// matrix nonzeros.

void Hypre::update_fine_coarse_const_
(   int        face, 
    Grid &     grid, 
    int        axis0, 
    phase_enum phase,
    int        level_fine, 
    int        level_coarse,
    int        index_fine[3], 
    int        index_coarse[3])
{

  int axis1 = (axis0+1)%3;
  int axis2 = (axis0+2)%3;

  //--------------------------------------------------
  // (*) CONSTANT
  //     Scale        = 2/3
  //     Coefficients = 1
  //--------------------------------------------------

  Scalar val_h_fine   = grid.h(axis1) * grid.h(axis2) / grid.h(axis0);

  int index_increment[][3] = {{face*(r_factor_-1),0,0},
		    {0,1,0},
		    {0,0,1},
		    {0,-1,0},
		    {-face*(r_factor_-1),0,-1}};

  if (grid.is_local()) {

    if (phase == phase_graph) {

      int k = 0;

      index_fine[axis0] += index_increment[k][0];
      index_fine[axis1] += index_increment[k][1];
      index_fine[axis2] += index_increment[k][2];

      for (k=1; k<5; k++) {

	HYPRE_SStructGraphAddEntries 
	  (graph_, level_fine, index_fine, 0, level_coarse, index_coarse, 0);

	index_fine[axis0] += index_increment[k][0];
	index_fine[axis1] += index_increment[k][1];
	index_fine[axis2] += index_increment[k][2];

      } // for k = 1:4

    } else if (phase == phase_matrix) {

      // fine->coarse off-diagonal

      double val_s = 2. / 3.;

      int entry;
      double val_a;
      double val,value;

      int k=0;

      index_fine[axis0] += index_increment[k][0];
      index_fine[axis1] += index_increment[k][1];
      index_fine[axis2] += index_increment[k][2];

      for (int k=1; k<5; k++) {

	// (x,y,z) = ???

	val_a = 1.0; // DIFFUSION COEFFICIENT GOES HERE

	// Update off-diagonal

	entry = grid.counter(index_fine)++;
	val   = matrix_scale_ * val_h_fine * val_s * val_a;

	value = val;
	HYPRE_SStructMatrixAddToValues 
	  (A_, level_fine, index_fine, 0, 1, &entry, &value);

	// Update diagonal

	entry = 0;
	value = -val;

	HYPRE_SStructMatrixAddToValues 
	  (A_, level_fine, index_fine, 0, 1, &entry, &value);

	// Clear coarse-fine stencil values

	val   = matrix_scale_ * val_h_fine * val_a;

	entry = 2*axis0 + 1 + (1-face); // stencil xp=1,xm,yp,ym,zp,zm=6
	value = - val;

	HYPRE_SStructMatrixAddToValues 
	  (A_, level_fine, index_fine, 0, 1, &entry, &value);

	entry = 0; // diagonal
	value = val;

	HYPRE_SStructMatrixAddToValues 
	  (A_, level_fine, index_fine, 0, 1, &entry, &value);

	// Update indices

	index_fine[axis0] += index_increment[k][0];
	index_fine[axis1] += index_increment[k][1];
	index_fine[axis2] += index_increment[k][2];

      } // for k = 1,4



    } // if phase == phase_matrix
  } // if grid.is_local()
} // if discret_type_const

//------------------------------------------------------------------------

/// Update the matrix at a coarse face zone adjacent to fine face
/// zones using a simple piecewise-constant finite volume
/// discretization.  Called in two phases, with phase == phase_graph
/// (via init_graph_nonstencil_()) for the nonzero structure, and with
/// phase == phase_matrix (via init_matrix_nonstencil_()) for the
/// matrix nonzeros.

void Hypre::update_coarse_fine_const_
(   int        face, 
    Grid &     grid_coarse, 
    int        axis0,
    phase_enum phase,
    int        level_fine, 
    int        level_coarse,
    int        index_fine[3], 
    int        index_coarse[3])
{

    int axis1 = (axis0+1)%3;
    int axis2 = (axis0+2)%3;

  Scalar val_h_coarse  = grid_coarse.h(axis1) * grid_coarse.h(axis2) / grid_coarse.h(axis0);

  if (phase == phase_graph) {

    int index_increment[][3] = {{1,0,0},
		      {0,1,0},
		      {-1,0,0},
		      {0,0,1},
		      {1,0,0},
		      {0,-1,0},
		      {-1,0,0},
		      {0,0,-1}};

    for (int k=0; k<8; k++) {

      HYPRE_SStructGraphAddEntries 
	(graph_, level_coarse, index_coarse, 0, level_fine, index_fine, 0);

      index_fine[0] += index_increment[k][0];
      index_fine[1] += index_increment[k][1];
      index_fine[2] += index_increment[k][2];

    } // for k=0,7

  } else if (phase == phase_matrix) {

    double val_s = 1./8.;
    double val_a = 1.0; // DIFFUSION COEFFICIENT GOES HERE
    double val   = matrix_scale_ * val_h_coarse * val_s * val_a;
    int    entry;

    double value;

    // Adjust coarse-fine nonstencil values

    for (int i=0; i<8; i++) {

      entry = grid_coarse.counter(index_coarse)++;
      value = val;

      // Set new nonstencil coarse-fine entry

      HYPRE_SStructMatrixAddToValues 
	(A_, level_coarse,index_coarse, 0, 1, &entry, &value);

      // Adjust stencil diagonal

      entry = 0;
      value = - val;
      HYPRE_SStructMatrixAddToValues 
	(A_, level_coarse, index_coarse, 0, 1, &entry, &value);

    } // for i=0,7

    val   = matrix_scale_ * val_h_coarse * val_a;

    // Clear coarse-fine stencil values

    // (note: "face" is for fine grid, but we want coarse)

    entry = 2*axis0 + 1 + face;     //stencil xp=1,xm,yp,ym,zp,zm=6
    value = - val;

    HYPRE_SStructMatrixAddToValues 
      (A_, level_coarse, index_coarse, 0, 1, &entry, &value);

    // Adjust stencil values

    entry = 0; // diagonal
    value = val;

    HYPRE_SStructMatrixAddToValues 
          (A_, level_coarse, index_coarse, 0, 1, &entry, &value);


  } // if phase == phase_matrix
}
