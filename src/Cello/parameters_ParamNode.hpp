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
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~ParamNode() throw()
  {
    std::map<std::string,ParamNode *>::iterator it_param;
    for (it_param =  subnodes_.begin();
	 it_param != subnodes_.end();
	 ++it_param) {
      delete it_param->second;
    }
  };

private: // No copy or assign

  /// Copy constructor
  ParamNode(const ParamNode & param_node) throw()
  {  }

  /// Assignment operator
  ParamNode & operator= (const ParamNode & param_node) throw()
  {
    return *this;
  }

public: // interface

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    p | name_;
    WARNING("ParamNode::pup","skipping subnodes_");
    //    p | subnodes_;
  }
#endif

  /// Return the node name
  std::string name() const {return name_;};

  /// Return the number of subgroups
  int size()
  {return subnodes_.size(); }

  /// Return the ith subgroup
  std::string subgroup (int group_index)
  {
    if (0 <= group_index && group_index < size()) {
      std::map<std::string,ParamNode *>::iterator it_param;
      int i;
      for (i=0,it_param =  subnodes_.begin();
	   it_param != subnodes_.end();
	   ++it_param,i++) {
	if (group_index == i) {
	  return it_param->first;
	}
      }
    }
    return "";
  };

  /// Return the given subnode, returning 0 if it doesn't exist
  ParamNode * subnode(std::string subgroup)
  {
    return subnodes_[subgroup];
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
