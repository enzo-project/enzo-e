
/// Hypre class include file

/**
 * 
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 *
 */

class Hypre {

public:

  Hypre ();

  void init_hierarchy (Hierarchy & hierarchy, Mpi & mpi);
  void init_stencil   (Hierarchy & hierarchy);
  void init_graph     (Hierarchy & hierarchy);
  void init_matrix    (Hierarchy & hierarchy);
  void init_rhs       (Hierarchy & hierarchy);
  void init_solver    (Hierarchy & hierarchy);
  void solve          (Hierarchy & hierarchy);
  void evaluate       (Hierarchy & hierarchy);


private:

  HYPRE_SStructStencil stencil_;

private:

  void init_hierarchy_create_grid_        (Grid & grid, int dim, int levels);
  void init_hierarchy_set_grid_extents_   (Grid & grid);
  void init_hierarchy_set_grid_variables_ (Grid & grid);
  void init_hierarchy_set_grid_neighbors_ (Grid & grid);
  void init_hierarchy_assemble_grids_     (Grid & grid);
  
  void init_graph_stencil                 (Grid & grid);
  void init_graph_zero_covered            (Grid & grid);
  void init_graph_neighbors               (Grid & grid);
  void init_graph_boundary                (Grid & grid);
  void init_graph_children                (Grid & grid);
  void init_graph_parent                  (Grid & grid);
  void init_graph_neighbors_children      (Grid & grid);
  void init_graph_parents_neighbor        (Grid & grid);
  void init_graph_check                   (Grid & grid);

};
