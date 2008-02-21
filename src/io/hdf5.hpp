#ifndef CLASS_HDF5
#define CLASS_HDF5
// $Id$
/**
 * @file
 * @brief Header file for the Hdf5 class
 * @author James Bordner 
 * @version 1.0
 *
 * Attributes
 *
 * ( ) 
 *
 * Operations
 *
 *     public:
 *
 * ( ) Hdf5();
 * ( ) void file_open ()
 * ( ) void file_read (Array & array);
 * ( ) void file_write (Array & array);
 * ( ) void file_close ()
 * ( ) void group_open ()
 * ( ) void group_close ()
 * ( ) void dataset_open ()
 * ( ) void dataset_close ()
 * ( ) void append (Array & array);
 *
 * ( ) Hdf5()
 *
 *     private:
 *
 */
// $Log$

/// Brief description

/**
Detailed description
*/

 
class Hdf5 {

  // PRIVATE ATTRIBUTES

private:

  // PUBLIC OPERATIONS

  Hdf5();
  void read   (Array &, std::string);
  void write  (Array &, std::string);
  void append (Array &, std::string);

public:

  // PRIVATE OPERATIONS

private:

};

#endif
