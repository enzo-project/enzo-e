// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_ParamNode.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon May 10 12:43:27 PDT 2010
/// @brief [\ref Parameters] Node for representing parameters in a
///           tree of groups / subgroups /parameters /values

#ifndef PARAMETERS_PARAM_NODE_HPP
#define PARAMETERS_PARAM_NODE_HPP

class ParamNode {

  /// @class    ParamNode
  /// @ingroup  Parameters
  /// @brief    [\ref Parameters] Node representing a subtree of parameters

public: // interface

  /// Constructor
  ParamNode(std::string name) throw()
    : name_(name),
      subnodes_()
  {};

  //----------------------------------------------------------------------
  // Big Five
  //----------------------------------------------------------------------

  /// Destructor
  ~ParamNode() throw()
  {
    for (std::pair<const std::string, ParamNode*>& key_val : subnodes_) {
      delete key_val.second;
    }
  };

  // No copy constructor or copy assignment:
  ParamNode(const ParamNode&) = delete;
  ParamNode & operator= (const ParamNode &) = delete;

  // force the compiler to generate move constructor and move assignment
  ParamNode(ParamNode&&) = default;
  ParamNode & operator= (ParamNode &&) = default;

public: // interface

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change

    p | name_;

    // pup std::map<std::string,ParamNode*> subnodes_
    
    int n = this->size();
    p | n;
    if (!p.isUnpacking()) {
      for (std::pair<const std::string, ParamNode*>& key_val : subnodes_) {
        std::string name = key_val.first;
        p | name;
        ParamNode* subnode = key_val.second;
        p | *subnode;
      }
    } else {
      for (int i=0; i<n; i++) {
        std::string name;
        p | name;
	ParamNode * subnode = new ParamNode(name);
        p | *subnode;
        subnodes_[name] = subnode;
      }
    }
  }

  /// Return the node name
  std::string name() const {return name_;};

  /// Return the number of subgroups
  int size() const { return subnodes_.size(); }

  /// Return the given subnode, returning 0 if it doesn't exist
  const ParamNode * subnode(std::string subgroup) const
  {
    auto search = subnodes_.find(subgroup);
    return (search != subnodes_.end()) ? search->second : nullptr;
  }

  /// Return the given subgroup, creating a new one if it doesn't exist
  ParamNode * new_subnode(std::string subgroup)
  {
    if (subnodes_[subgroup] == 0) {
      subnodes_[subgroup] = new ParamNode(subgroup);
    }

    return subnodes_[subgroup];
  }

private: // attributes

  /// Subnodes of the tree
  std::string name_;
  std::map<std::string, ParamNode *> subnodes_;

};

#endif /* PARAMETERS_PARAM_NODE_HPP */
