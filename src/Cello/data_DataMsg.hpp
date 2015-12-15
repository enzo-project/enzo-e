// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsg.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-07-06
/// @brief    [\ref Data] Declaration of the DataMsg Charm++ Message

#ifndef DATA_DATA_MSG_HPP
#define DATA_DATA_MSG_HPP

class DataMsg : public CMessage_DataMsg {

public:
  /// Array length
  int n;

  /// Array data
  char * a;

  /// Child indices
  int ic3[3];
};

#endif /* DATA_DATA_MSG_HPP */

