// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_NodeTrace.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-09
/// @brief    [\ref Mesh] Declaration of the NodeTrace class
///
/// The NodeTrace class is used for storing information about
/// the location of a Node in a Tree, and provides functions
/// for generating and accessing a NodeTrace.  Traces are limited
/// to Tree's with at most 64 levels.

#ifndef MESH_NODE_TRACE_HPP
#define MESH_NODE_TRACE_HPP

class NodeTrace {

  /// @class    NodeTrace
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  NodeTrace(int d, int r) throw();

  /// Trace one level down the Tree according to the provided x,y,z indices from
  /// the given Node; return true when we've reached a leaf node
  Node * trace (Node * node, int ix, int iy=0, int iz=0);

  /// Reset the trace variables for the next trace
  void reset_trace() 
  { trace_[0]=0;
    trace_[1]=0; 
    trace_[2]=0; };

  /// Reverse the trace_[] bit values after they're accumulated so
  /// that they can be used
  void reverse_trace();

  unsigned long long trace(axis_enum axis) const
  { return trace_[axis]; }

private: // functions

private: // attributes

  /// Dimensionality of the Tree the Node is contained in
  int d_;

  /// Refinement amout of the Tree the Node is contained in
  int r_;

  /// Level of the node in the Tree, with Tree's root being level 0
  int level_;

  /// Number of bits used per level (1 for octree, 2 for 3^2-tree or
  /// 4^2-tree, etc.)
  int bits_;

  /// bit trace from Tree root, enough for 3D Tree's of at least 64 levels
  unsigned long long trace_[3]; 

};

#endif /* MESH_NODE_TRACE_HPP */

