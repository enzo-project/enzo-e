// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MsgOrder.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2023-07-22
/// @brief    [\ref Charm] Declaration of the MsgOrder Charm++ message

#include "charm.hpp"

//----------------------------------------------------------------------

long MsgOrder::counter[CONFIG_NODE_SIZE] = {0};

