// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_NodeTrace.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-09
/// @brief    [\ref Mesh] Declaration of the NodeTrace class
///
/// The NodeTrace class is used for storing information about the
/// location of a Node in a Tree with respect to the root node.  It
/// functions like a stack, with operations for "pushing" down to
/// child nodes, and "popping" back up to the parent.  Child indexing
/// is the same as how Nodes are stored (Morton- or Z-ordering).
/// Other orderings such as Hilbert or Moore space-filling curves are
/// implemented at a higher level in the ItNode class hierarchy.

#ifndef MESH_NODE_TRACE_HPP
#define MESH_NODE_TRACE_HPP

class Tree;
class Node;

class NodeTrace {

  /// @class    NodeTrace
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  NodeTrace(Node * root) throw();

  /// Copy constructor
  NodeTrace(const NodeTrace &) throw();

  /// assignment operator
  NodeTrace & operator=(const NodeTrace &) throw();

  /// Destructor
  ~NodeTrace () throw();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    p |  level_;
    p |  index_;
    WARNING("NodeTrace::pup","not pup'ing node_ vector of pointers");
    //    p |  node_;
  }
#endif
  
  /// Move one level "up" the Tree to the specified child (Z-ordering)
  inline Node * push (int index) 
  {
    Node * child = node_[level_]->child(index);

    if (child != NULL) {
      ++level_;
      index_.push_back(index);
      node_.push_back(child);
    }

    return child;
  };

  /// Move one level "down" the Tree to the parent
  inline Node * pop ()
  {
    --level_;
    node_.pop_back();
    index_.pop_back();
    return node_[level_];
  };


  /// Reset the node trace to the root
  void reset()
  { level_ = 0; 
    node_.resize(1);
    index_.resize(1);
  };

  /// Return the level of the traced node
  int level() const
  { return level_; }

  /// Return the last node in the trace
  Node * node() const
  { return node_[level_]; }

  /// Return the node in the node trace at the specified level
  Node * node_level(int level) const;

  /// Return the parent's child index of the last node in the trace
  int index() const;

  /// Return the parent's child index of the node in the node trace at the
  /// specified level.  Root parent's child index is undefined.
  int index_level(int level) const;

private: // functions

  /// Copy node_trace into this NodeTrace; called by constructor and
  /// assignment
  void copy_ (const NodeTrace & node_trace);

private: // attributes

  /// Level of the current Node in the Tree, with Tree's root being level 0
  int level_;

  /// Trace of child node indices
  std::vector<int> index_;

  /// Trace of node pointers
  std::vector<Node *> node_;

};

#endif /* MESH_NODE_TRACE_HPP */

