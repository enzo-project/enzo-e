//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Interface routines to HYPRE

/**
 * 
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 *
 */

#include <mpi.h>
#include <assert.h>
#include <string>
#include <vector>
#include <map>

#include "HYPRE_sstruct_ls.h"

#include "hypre-solve.hpp"

#include "mpi.hpp"
#include "scalar.hpp"
#include "error.hpp"
#include "constants.hpp"
#include "point.hpp"
#include "faces.hpp"
#include "domain.hpp"
#include "grid.hpp"
#include "level.hpp"
#include "hierarchy.hpp"
#include "sphere.hpp"
#include "parameters.hpp"
#include "problem.hpp"
#include "hypre.hpp"
#include "error.hpp"

const int debug  = 0;
const int trace  = 0;
const int trace_hypre  = 0;
const int trace_graph  = 0;

const Scalar matrix_scale = 1.0;  // 1.0:  1 1 1 -6 1 1 1
const Scalar matrix_diag = 0e0;  // + cu 

//======================================================================
// PUBLIC MEMBER FUNCTIONS
//======================================================================

/// Hypre constructor

Hypre::Hypre (Parameters & parameters)
  : grid_(0),
    stencil_(0),
    graph_(0),
    A_(0),
    B_(0),
    X_(0),
    solver_(0),
    parameters_(parameters)
{
  
}

//----------------------------------------------------------------------

/// Initialize the Grid Hierarchy

/** Creates a hypre grid, with one part per level and one box per Grid
    patch object, for an AMR problem.  Sets grid box extents, grid
    part variables, and periodicity of the root-level grid part. */

void Hypre::init_hierarchy (Parameters & parameters,
			    Hierarchy  & hierarchy, 
			    Mpi        & mpi)
{

  int dim       = hierarchy.dimension();
  int num_parts = hierarchy.num_levels();

  // Create the hypre grid
  
  HYPRE_SStructGridCreate (MPI_COMM_WORLD, dim, num_parts, &grid_);

  ItHierarchyLevels itl (hierarchy);

  int part = 0;

  while (Level * level = itl++) {

    ItLevelGridsLocal itg (*level);

    while (Grid * grid = itg++) {

      int lower[3] = {grid->i_lower(0),grid->i_lower(1),grid->i_lower(2)};
      int upper[3] = {grid->i_upper(0),grid->i_upper(1),grid->i_upper(2)};

      // Set extents for boxes that comprise the hypre grid

      HYPRE_SStructGridSetExtents(grid_, part, lower, upper);
      
    }

    // Create a single cell-centered variable for each grid part (level)

    HYPRE_SStructVariable variable_types[] = { HYPRE_SSTRUCT_VARIABLE_CELL };
    const int numvars = 1;

    HYPRE_SStructGridSetVariables(grid_, part, numvars, variable_types);

    // Set grid part to be periodic, with periodicity determined by the root
    // level size, current level, and refinement factor (ASSUMED TO BE 2)

    const int r = 2;
    int periodicity[3];

    // Determine periodicity of Level

    for (int i=0; i<3; i++) {
      periodicity[i] = hierarchy.level(0).zones(i);
      for (int k=0; k < part; k++) periodicity[i] *= r;
    }

    if (parameters.value("boundary") == "dirichlet") {
      periodicity[0] = 0;
      periodicity[1] = 0;
      periodicity[2] = 0;
    }

    if (debug) printf ("%s:%d  Level = %d Periodicity = (%d,%d,%d)\n",
		       __FILE__,__LINE__, part,
		       periodicity[0],periodicity[1],periodicity[2]);

    HYPRE_SStructGridSetPeriodic (grid_, part, periodicity);

    ++ part;
  }

  // When finished, assemble the hypre grid

  HYPRE_SStructGridAssemble (grid_);
  
}

//----------------------------------------------------------------------

/// Initialize the discretization stencils.  

/** Creates and initializes a stencil object.  Supports 1, 2, or 3
    dimensional stencils. */

void Hypre::init_stencil (Hierarchy & hierarchy)

{

  int dim = hierarchy.dimension();

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

}

//----------------------------------------------------------------------

/// Initialize the graph.

/** Creates a graph containing the matrix non-zero structure.  Graph
    edges include both those for non-zeros from the stencil within
    each part (level), and non-zeros for graph entries connecting
    linked parts.

    Setting up the matrix nonzero structure is performed using the
    following steps:

     1. Define stencil connections within each level

     2. Define connections for unknowns adjacent to coarse unknowns

     3. Define connections for unknowns adjacent to fine unknowns 

    The matrix nonzero structure is generally nonsymmetric.  Only step
    1 is required for unigrid problems.

*/

void Hypre::init_graph (Hierarchy & hierarchy)

{
  // Create the hypre graph object

  HYPRE_SStructGraphCreate (MPI_COMM_WORLD, grid_, &graph_);

  HYPRE_SStructGraphSetObjectType (graph_, HYPRE_SSTRUCT);

  ItHierarchyLevels itl (hierarchy);

  int part = 0;

  while (Level * level = itl++) {

    // 1. Define stencil connections within each level

    HYPRE_SStructGraphSetStencil (graph_, part, 0, stencil_);

    ++ part;

    // 2. Define matrix nonzero structure connecting grids in level
    // with next-coarser.

    if (level->index() > 0) {
      ItLevelGridsAll itag (*level);

      while (Grid * grid = itag++) {

	// Define nonstencil entries for the grid

	init_graph_nonstencil_(*grid);

      }
    }

    ItLevelGridsAll itag (*level);

    while (Grid * grid = itag++) {

      // Clear the nonstencil entry counter for subsequent matrix
      // nonstencil entries

      int dim = hierarchy.dimension();
      grid->init_counter(dim*2+1);
    }
  }

  // Assemble the graph

  HYPRE_SStructGraphAssemble (graph_);

}


//----------------------------------------------------------------------

/// Initialize the matrix A and right-hand-side vector b

/** Creates a matrix with a given non-zero structure, and sets nonzero
    values.

    Setting up the matrix elements is done with the following
    steps:

     1. Set stencil values

     2. Clean up stencil connections between parts and under overlapped grids

     3. Set values for unknowns between parent and children

    The matrix is generally nonsymmetric.  Only step 1 is required for
    unigrid problems.

*/

void Hypre::init_linear (Parameters          & parameters,
			 Hierarchy           & hierarchy,
			 std::vector<Point *>  points,
			 std::vector<Sphere *> spheres)

{
  // Create the hypre matrix A_, solution X_, and right-hand side B_ objects

  HYPRE_SStructMatrixCreate (MPI_COMM_WORLD, graph_, &A_);
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid_,  &X_);
  HYPRE_SStructVectorCreate (MPI_COMM_WORLD, grid_,  &B_);

  // Set the object types

  HYPRE_SStructMatrixSetObjectType (A_,HYPRE_SSTRUCT);
  HYPRE_SStructVectorSetObjectType (X_,HYPRE_SSTRUCT);
  HYPRE_SStructVectorSetObjectType (B_,HYPRE_SSTRUCT);

  // Initialize the hypre matrix and vector objects

  HYPRE_SStructMatrixInitialize (A_);
  HYPRE_SStructVectorInitialize (X_);
  HYPRE_SStructVectorInitialize (B_);

  //--------------------------------------------------
  // Initialize the matrix A_
  //--------------------------------------------------

  ItHierarchyLevels itl (hierarchy);

  while (Level * level = itl++) {

    ItLevelGridsLocal itlg (*level);
    
    // 1. Set stencil values within level

    while (Grid * grid = itlg++) {

      init_matrix_stencil_(*grid);
    }

    if (level->index() > 0) {

      // 2. Clean up stencil connections between levels

      init_matrix_clear_(*level);

      // 3. Set matrix values between levels

      // WARNING: POSSIBLE SCALING ISSUE.  Below we loop over all grids;
      // however, we only need too loop over parent-child pairs such
      // that either child or parent is local to this MPI process.
 
      ItLevelGridsAll itag (*level);

      while (Grid * grid = itag++) {

 	init_matrix_nonstencil_(*grid);

      }
    }
  }

  //--------------------------------------------------
  // Initialize B_ according to density
  //--------------------------------------------------

  Scalar local_shift_b_sum = 0.0;

  local_shift_b_sum += init_vector_points_  (hierarchy,points);
  local_shift_b_sum += init_vector_spheres_ (hierarchy,spheres);

  Scalar shift_b_sum = 0.0;

  MPI_Allreduce (&local_shift_b_sum, &shift_b_sum, 1, 
		 MPI_SCALAR, MPI_SUM, MPI_COMM_WORLD);

  // Shift B to zero out the null space if problem is periodic

  if ( parameters.value("boundary") == "periodic" ) {

    // Compute the shift

    int part = 0;
    long long shift_b_count = 0;
    ItHierarchyLevels itl (hierarchy);
    while (Level * level = itl++) {
      ItLevelGridsAll itg (*level);
      while (Grid * grid = itg++) {
	shift_b_count += grid->num_unknowns();
      }
    }

    Scalar shift_b_amount = - shift_b_sum / shift_b_count;
    if (debug) printf ("Periodic shift = %g\n",shift_b_amount);

    // Perform the shift
    
    part = 0;
    while (Level * level = itl++) {
      ItLevelGridsLocal itg (*level);
      while (Grid * grid = itg++) {
	int lower[3] = { grid->i_lower(0), grid->i_lower(1), grid->i_lower(2) };
	int upper[3] = { grid->i_upper(0), grid->i_upper(1), grid->i_upper(2) };
	Scalar * values = new Scalar[grid->num_unknowns()];
	for (int i=0; i<grid->num_unknowns(); i++) values[i] = shift_b_amount;
	HYPRE_SStructVectorAddToBoxValues (B_,part,lower,upper,0,values);

	delete [] values;
      }


      ++part;
    }
    if (debug) printf ("%s:%d shift (count,sum,amount) = (%lld,%g,%g)\n",
		       __FILE__,__LINE__,shift_b_count,shift_b_sum,
		       shift_b_amount);
  }

  // Assemble the matrix and vectors

  HYPRE_SStructMatrixAssemble (A_);
  HYPRE_SStructVectorAssemble (B_);
  HYPRE_SStructVectorAssemble (X_);

  // Write the vector to a file for debugging

  if (parameters.value("dump_a") == "true") HYPRE_SStructMatrixPrint ("A",A_,1);
  if (parameters.value("dump_b") == "true") HYPRE_SStructVectorPrint ("B",B_,1);
  if (parameters.value("dump_x") == "true") HYPRE_SStructVectorPrint ("X0",X_,1);

}

//----------------------------------------------------------------------

/// Initialize and solve the linear solver

void Hypre::solve (Parameters & parameters,
		   Hierarchy & hierarchy)

{
  std::string solver = parameters.value("solver");
  int         levels = hierarchy.num_levels();

  int    itmax  = 0;
  double restol = 0.0;
  
  std::string sitmax  = parameters.value("solver_itmax");
  std::string srestol = parameters.value("solver_restol");
  if (sitmax != "")  itmax  = atoi(sitmax.c_str());
  if (srestol != "") restol = atof(srestol.c_str());

  if        (solver == "fac"  && levels > 1) {
    solve_fac_(hierarchy,itmax,restol);
  } else if (solver == "bicgstab") {
    solve_bicgstab_(hierarchy,itmax,restol);
  } else if (solver == "pfmg" && levels == 1) {
    solve_pfmg_(hierarchy,itmax,restol);
  } else {
    char error_message[100];
    sprintf (error_message, "Hypre::solve called with illegal combination of "
	     "solver %s on %d levels", solver.c_str(),levels);
    ERROR(error_message);
  }
  
  if (parameters.value("dump_x") == "true") HYPRE_SStructVectorPrint ("X",X_,1);

}

//----------------------------------------------------------------------

/// Evaluate the success of the solve

void Hypre::evaluate (Hierarchy & hierarchy)

{
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  NOT_IMPLEMENTED("Hypre::evaluate()");
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
}

//======================================================================
// PRIVATE MEMBER FUNCTIONS
//======================================================================

/// init_nonstencil_() is called twice: first by
/// init_graph_nonstencil_(), with phase == "graph", to set nonstencil
/// graph entries, and again by init_matrix_nonstencil_(), with phase
/// == "matrix", to set nonstencil matrix entries.

void Hypre::init_nonstencil_ (Grid & grid, std::string phase)
{

  if (phase != "graph" && phase != "matrix") {
    char error_message[80];
    strcpy (error_message,"init_matrix_nonstencil_ called with phase = ");
    strcat (error_message,phase.c_str());
    ERROR(error_message);
  }

  const int r = 2; // WARNING: hard-coded refinement factor r = 2

  int face,axis,ig1,ig2;

  // global grid index limits

  int ig3[3][2];
  grid.indices(ig3);
  
  // Loop over each face zone in the grid, adding non-stencil graph
  // entries wherever a zone is adjacent to a coarse zone.  Both
  // fine-to-coarse and coarse-to-fine entries are added.

  int level_fine   = grid.level();
  int level_coarse = grid.level() - 1;

  assert (level_coarse >= 0);

  for (axis=0; axis<3; axis++) {

    // j0:    axis normal to face
    // j1,j2: axes within face

    int j0 = axis;
    int j1 = (axis+1)%3;
    int j2 = (axis+2)%3;

    // n0     grid size normal to face
    // n1,n2: grid size within face

    int n0 = ig3[j0][1] - ig3[j0][0];
    int n1 = ig3[j1][1] - ig3[j1][0];
    int n2 = ig3[j2][1] - ig3[j2][0];

    double h0 = grid.h(j0);
    double h1 = grid.h(j1);
    double h2 = grid.h(j2);

    double H0 = r*grid.h(j0);
    double H1 = r*grid.h(j1);
    double H2 = r*grid.h(j2);

    // ig3[][] should be divisible by r**level.  Just test r here.

    if ((ig3[j1][0]/r)*r != ig3[j1][0]) printf ("ig3[%d][0] = %d)\n",j1,ig3[%d][0]);
    assert ((ig3[j1][0]/r)*r == ig3[j1][0]);
    if ((ig3[j1][1]/r)*r != ig3[j1][1]) printf ("ig3[%d][1] = %d)\n",j1,ig3[%d][1]);
    assert ((ig3[j1][1]/r)*r == ig3[j1][1]);

    for (face=0; face<2; face++) {

      // Loop over face zones that are aligned with coarse zones (hence "+= r")

      for (ig1=0; ig1<n1; ig1 += r) {
	for (ig2=0; ig2<n2; ig2 += r) {

	  Grid * adjacent   = grid.faces().adjacent(axis,face,ig1,ig2);

	  Faces::Label & fz = grid.faces().label(axis,face,ig1,ig2);

	  // Add graph entries iff grid or adjacent grid is local, and if
	  // adjacent grid (if it exists) is in the next-coarser level

	  bool is_local = 
	    (adjacent != NULL) &&  (adjacent->is_local() || grid.is_local());

	  if (is_local && fz == Faces::_coarse_) {

	    // (fine) grid global indices

	    int igg3[3]; 

	    igg3[j0] = ig3[j0][0] + face*(n0 - r);
	    igg3[j1] = ig3[j1][0] + ig1;
	    igg3[j2] = ig3[j2][0] + ig2;

	    // (coarse) adjacent global indices

	    int ign3[3]; 

	    ign3[j0] = (igg3[j0]) / r  + (face*r-1);
	    ign3[j1] = (igg3[j1]) / r;
	    ign3[j2] = (igg3[j2]) / r;

	    //--------------------------------------------------
	    // GRAPH ENTRY: FINE-TO-COARSE 
	    //--------------------------------------------------

	    if (parameters_.value("discret") == "constant") {

	      //--------------------------------------------------
	      // (*) CONSTANT
	      //     Scale        = 2/3
	      //     Coefficients = 1
	      //--------------------------------------------------

	      int diggs[][3] = {{face*(r-1),0,0},
				{0,1,0},
				{0,0,1},
				{0,-1,0},
				{-face*(r-1),0,-1}};

	      if (phase == "graph") {
		int k=0;

		igg3[j0]+=diggs[k][0];
		igg3[j1]+=diggs[k][1];
		igg3[j2]+=diggs[k][2];
		for (int k=1; k<5; k++) {
		  if (trace_graph) {
		    printf ("TRACE GRAPH (%d %d %d; %d) : (%d %d %d; %d)\n",
			    igg3[0],igg3[1],igg3[2], level_fine,  
			    ign3[0],ign3[1],ign3[2], level_coarse);
		  }
		  HYPRE_SStructGraphAddEntries 
		    (graph_, level_fine, igg3, 0, level_coarse, ign3, 0);
		  _TRACE_;
		  igg3[j0]+=diggs[k][0];
		  igg3[j1]+=diggs[k][1];
		  igg3[j2]+=diggs[k][2];
		}

	      } else if (phase == "matrix") {

		// fine->coarse off-diagonal

		double val_h = h1*h2/h0;
		double val_s = 2. / 3.;

		int entry;
		double val_a;
		double val;

		int k=0;

		igg3[j0]+=diggs[k][0];
		igg3[j1]+=diggs[k][1];
		igg3[j2]+=diggs[k][2];

		for (int k=1; k<5; k++) {

		  val_a = 1.0; // DIFFUSION COEFFICIENT GOES HERE

		  // Update off-diagonal

		  entry = grid.counter(igg3)++;
		  val   = matrix_scale * val_h * val_s * val_a;
		  HYPRE_SStructMatrixAddToValues 
		    (A_, level_fine, igg3, 0, 1, &entry, &val);

		  // Update diagonal

		  entry = 0;
		  val = -val;
		  HYPRE_SStructMatrixAddToValues 
		    (A_, level_fine, igg3, 0, 1, &entry, &val);

		  igg3[j0]+=diggs[k][0];
		  igg3[j1]+=diggs[k][1];
		  igg3[j2]+=diggs[k][2];
		}
	      }

	    } else {
	      char error_message[80];
	      strcpy (error_message,"Unknown parameter discret = ");
	      strcat (error_message,parameters_.value("discret").c_str());
	      ERROR(error_message);
	    }

	    //--------------------------------------------------
	    // GRAPH ENTRY: COARSE-TO-FINE
	    //--------------------------------------------------

	    if (phase == "graph") {

	      int diggs[][3] = {{1,0,0},
				{0,1,0},
				{-1,0,0},
				{0,0,1},
				{1,0,0},
				{0,-1,0},
				{-1,0,0},
				{0,0,-1}};

	      for (int k=0; k<8; k++) {
		HYPRE_SStructGraphAddEntries 
		  (graph_, level_coarse, ign3, 0, level_fine, igg3, 0);
		  if (trace_graph) {
		    printf ("TRACE GRAPH (%d %d %d; %d) : (%d %d %d; %d)\n",
			    ign3[0],ign3[1],ign3[2], level_coarse,  
			    igg3[0],igg3[1],igg3[2], level_fine);
		  }
		igg3[0] += diggs[k][0];
		igg3[1] += diggs[k][1];
		igg3[2] += diggs[k][2];
	      }

 	    } else if (phase == "matrix") {

 	      double val_h = H1*H2/H0;
	      double val_s = 1.;
 	      double val_a = 1.0; // DIFFUSION COEFFICIENT GOES HERE
 	      double val   = matrix_scale * val_h * val_s * val_a;
	      int    entry;
 	      double value;

	      // coarse->fine off-diagonal

 	      for (int i=0; i<8; i++) {
 		entry = adjacent->counter(ign3)++;
 		value = (1./8.) * val;
		HYPRE_SStructMatrixAddToValues 
		  (A_, level_coarse,ign3, 0, 1, &entry, &value);
 	      }

	      // coarse->coarse diagonal

	      double val_diag = -val;
	      entry = 0;

	      //	      _TEMPORARY_;
	      //	      val_diag*=0.00;
	      // HYPRE_SStructMatrixAddToValues 
	      // (A_, level_coarse, ign3, 0, 1, &entry, &val_diag);

	    }
	  }
	}
      }
    }
  }
}

//------------------------------------------------------------------------

/// Set matrix stencil values for the grid interior

void Hypre::init_matrix_stencil_ (Grid & grid)

{
  int lower[3]    = { grid.i_lower(0), grid.i_lower(1), grid.i_lower(2) };
  int upper[3]    = { grid.i_upper(0), grid.i_upper(1), grid.i_upper(2) };
  int n           = grid.num_unknowns();
  int entries[7]  = { 0,1,2,3,4,5,6 };
  double h3[3]    = {grid.h(0),grid.h(1),grid.h(2)};
  int    n3[3]    = {grid.n(0),grid.n(1),grid.n(2)};

  double h120 = h3[1]*h3[2] / h3[0];
  double h201 = h3[2]*h3[0] / h3[1];
  double h012 = h3[0]*h3[1] / h3[2];

  double hhh = h3[0]*h3[1]*h3[2];

  double * v0  = new double [n];
  double * vxp  = new double [n];
  double * vxm  = new double [n];
  double * vyp  = new double [n];
  double * vym  = new double [n];
  double * vzp  = new double [n];
  double * vzm  = new double [n];

  int i=0;
  int level = grid.level();

  for (int i2 = 0; i2 < n3[2]; i2++) {

    // DIFFUSION COEFFICIENTS HERE

    double azp = (level==0) || (i2 < n3[2]-1) ? 1.0 : 0.0;
    double azm = (level==0) || (i2 >       0) ? 1.0 : 0.0;

    for (int i1 = 0; i1 < n3[1]; i1++) {

      // DIFFUSION COEFFICIENTS HERE

      double ayp = (level==0) || (i1 < n3[1]-1) ? 1.0 : 0.0;
      double aym = (level==0) || (i1 >       0) ? 1.0 : 0.0;

      for (int i0 = 0; i0 < n3[0]; i0++) {

	// DIFFUSION COEFFICIENTS HERE

	double axp = (level==0) || (i0 < n3[0]-1) ? 1.0 : 0.0;
	double axm = (level==0) || (i0 >       0) ? 1.0 : 0.0;

	int i = i0 + n3[0]*(i1 + n3[1]*i2);

	vxp[i] = matrix_scale * h120 * axp;
	vxm[i] = matrix_scale * h120 * axm;
	vyp[i] = matrix_scale * h201 * ayp;
	vym[i] = matrix_scale * h201 * aym;
	vzp[i] = matrix_scale * h012 * azp;
	vzm[i] = matrix_scale * h012 * azm;


	v0[i] = -( vxp[i] + vxm[i] + vyp[i] + vym[i] + vzp[i] + vzm[i] 
		   + hhh*matrix_diag );

      }
    }
  }

  HYPRE_SStructMatrixSetBoxValues (A_,level,lower,upper,0,1,&entries[0],v0);
  HYPRE_SStructMatrixSetBoxValues (A_,level,lower,upper,0,1,&entries[1],vxp);
  HYPRE_SStructMatrixSetBoxValues (A_,level,lower,upper,0,1,&entries[2],vxm);
  HYPRE_SStructMatrixSetBoxValues (A_,level,lower,upper,0,1,&entries[3],vyp);
  HYPRE_SStructMatrixSetBoxValues (A_,level,lower,upper,0,1,&entries[4],vym);
  HYPRE_SStructMatrixSetBoxValues (A_,level,lower,upper,0,1,&entries[5],vzp);
  HYPRE_SStructMatrixSetBoxValues (A_,level,lower,upper,0,1,&entries[6],vzm);

  delete [] vzm;
  delete [] vzp;
  delete [] vym;
  delete [] vyp;
  delete [] vxm;
  delete [] vxp;
  delete [] v0;

}

//------------------------------------------------------------------------

/// Clean up stencil connections between parts for FAC solver

void Hypre::init_matrix_clear_ (Level & level)
{
  // WARNING: hard-coding refinement factor of 2
  int r_factors[3] = {2,2,2}; 
  int part = level.index();
  if (part > 0) {

    // Clear stencil values from coarse to fine part

    HYPRE_SStructFACZeroCFSten (A_,grid_, part, r_factors);

    // Clear stencil values from fine to coarse part

    HYPRE_SStructFACZeroFCSten (A_,grid_, part);

    // Set overlapped areas of part with identity

    HYPRE_SStructFACZeroAMRMatrixData (A_, part-1, r_factors);

    // Need to clear under rhs also
    //   HYPRE_SStructFACZeroAMRVectorData(B_, plevels, prefinements);
  }
}

//------------------------------------------------------------------------

/// Add contributions from point sources to right-hand side B

Scalar Hypre::init_vector_points_ (Hierarchy            & hierarchy,
				   std::vector<Point *> & points)

{

  const Scalar scaling0 = -4.0*Constants::G()*Constants::pi();

  Scalar shift_b_sum = 0.0;

  int i;
  for (i=0; i<int(points.size()); i++) {
    Point & point      = *points[i];
    Grid & grid        = hierarchy.grid(point.igrid());
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
	int    i0 = grid.i_lower(k);
	index[k] = int (ap/ag*ig) + i0;
     
      }
      if (debug) {
	point.print();
	grid.print();
	printf ("Point index  = %d %d %d)\n",index[0],index[1],index[2]);
	printf ("Cell size    = %g %g %g\n",grid.h(0),grid.h(1),grid.h(2));
	printf ("Cell volume  = %g\n",cell_volume);
	printf ("Cell density = %g\n",density);
	printf ("RHS contribution = %g\n",value);
      }
    
      shift_b_sum += value;

      HYPRE_SStructVectorAddToValues (B_, grid.level(), index, 0, &value);
    }
  }
  return shift_b_sum;
}

//------------------------------------------------------------------------

/// Add contributions from sphere sources to right-hand side B

Scalar Hypre::init_vector_spheres_ (Hierarchy             & hierarchy,
				    std::vector<Sphere *> & spheres)

{

  if (spheres.size() > 0) {  
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    NOT_IMPLEMENTED("Contribution of sphere mass to right-hand side");
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  }

  return 0.0;
}

//------------------------------------------------------------------------

/// Initialize the PFMG hypre solver

void Hypre::solve_pfmg_ (Hierarchy & hierarchy, int itmax, double restol)

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

  int num_iterations;
  HYPRE_SStructSysPFMGGetNumIterations (solver_,&num_iterations);
  if (debug) printf ("HYPRE_SStructSysPFMGSolve num iterations: %d\n",num_iterations);

  double residual;
  HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm (solver_,&residual);
  if (debug) printf ("HYPRE_SStructSysPFMGSolve final relative residual norm: %g\n",residual);

  // Delete the solver

  HYPRE_SStructSysPFMGDestroy (solver_);
  solver_ = 0;

  
}

//------------------------------------------------------------------------

/// Initialize the FAC hypre solver

void Hypre::solve_fac_ (Hierarchy & hierarchy, int itmax, double restol)

{
  int i;

  const int r = 2; // WARNING: hard-coded refinement factor r = 2

  // Create the solver

  HYPRE_SStructFACCreate(MPI_COMM_WORLD, &solver_);

  // Initialize parts

  int num_parts = hierarchy.num_levels();
  int *parts  = new int [num_parts];
  for (i=0; i<num_parts; i++) parts[i] = i;
  HYPRE_SStructFACSetMaxLevels(solver_,  num_parts);
  HYPRE_SStructFACSetPLevels(solver_, num_parts, parts);

  // Initialize refinement factors

  typedef int int3[3];
  int3 *refinements = new int3 [num_parts];
  
  for (i=0; i<num_parts; i++) {
    refinements[i][0] = r;
    refinements[i][1] = r;
    refinements[i][2] = r;
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

  if (itmax != 0 )   HYPRE_SStructFACSetMaxIter(solver_,itmax);
  if (restol != 0.0) HYPRE_SStructFACSetTol(solver_,    restol);

  // output amount

  HYPRE_SStructFACSetLogging(solver_, 1);

  // prepare for solve

  HYPRE_SStructFACSetup2(solver_, A_, B_, X_);

  // Solve the linear system

  HYPRE_SStructFACSolve3(solver_, A_, B_, X_);

  // Write out some diagnostic info about the solve

  int num_iterations;
  HYPRE_SStructFACGetNumIterations(solver_, &num_iterations);
  double residual;
  HYPRE_SStructFACGetFinalRelativeResidualNorm(solver_, &residual);

  if (debug) printf ("HYPRE_SStructFACSolve3 num iterations: %d\n",num_iterations);

  if (debug) printf ("HYPRE_SStructFACSolve3 final relative residual norm: %g\n",residual);


  // Delete the solver

  HYPRE_SStructFACDestroy2(solver_);
  solver_ = 0;

  // Delete local dynamic storage
  delete [] parts;
  delete [] refinements;
}

//------------------------------------------------------------------------

/// Initialize the BICGSTAB hypre solver

void Hypre::solve_bicgstab_ (Hierarchy & hierarchy, int itmax, double restol)

{
  _TRACE_;

  // Create the solver

  HYPRE_SStructBiCGSTABCreate(MPI_COMM_WORLD, &solver_);

  _TRACE_;

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

  int num_iterations;
  HYPRE_SStructBiCGSTABGetNumIterations(solver_, &num_iterations);
  double residual;
  HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(solver_, &residual);

  if (debug) printf ("HYPRE_SStructBiCGSTABSolve3 num iterations: %d\n",num_iterations);

  if (debug) printf ("HYPRE_SStructBiCGSTABSolve3 final relative residual norm: %g\n",residual);


  // Delete the solver

  HYPRE_SStructBiCGSTABDestroy(solver_);
  solver_ = 0;

}

