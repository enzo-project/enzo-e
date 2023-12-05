// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_IoEnzoSimulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-05-12
/// @brief    [\ref Enzo] Declaration of the IoEnzoSimulation class

#ifndef IO_IO_ENZO_SIMULATION_HPP
#define IO_IO_ENZO_SIMULATION_HPP

class EnzoSimulation;

class IoEnzoSimulation : public IoSimulation {

  /// @class    IoEnzoSimulation
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Class for interfacing between EnzoSimulation and Output objects

public: // interface

  /// Constructor
  IoEnzoSimulation(const EnzoSimulation * simulation) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(IoEnzoSimulation);

  /// CHARM++ migration constructor
  IoEnzoSimulation(CkMigrateMessage *m) : IoSimulation(m) {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    IoSimulation::pup(p);
    TRACEPUP;

    p | index_enzo_;

    p | enzo_turbou_real_state_;
    p | enzo_turbou_int_state_;
    p | num_real_;
    p | num_int_;

  }

  /// Return the ith metadata item associated with the EnzoSimulation object
  void meta_value
  (int index,
   void ** buffer, std::string * name, int * type,
   int * nxd=0, int * nyd=0, int * nzd=0) throw();

  /// Return the ith data item associated with the EnzoSimulation object
  void data_value
  (int index,
   void ** buffer, std::string * name, int * type,
   int * nxd=0, int * nyd=0, int * nzd=0) throw();

  /// PACKING / UNPACKING

  /// Return the number of bytes required to serialize the data object
  virtual int data_size () const;

  /// Serialize the object into the provided empty memory buffer.
  virtual char * save_data (char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  virtual char * load_data (char * buffer);

  /// Copy the values to the object
  virtual void save_to (void *); 

private:

  int index_enzo_;

  std::vector<double> enzo_turbou_real_state_;
  std::vector<int>    enzo_turbou_int_state_;
  int num_real_;
  int num_int_;
};

#endif /* IO_IO_ENZO_SIMULATION_HPP */

