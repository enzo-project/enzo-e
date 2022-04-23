/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodFire2Feedback.hpp
/// @author     Andrew Emerick (aemerick11@gmail.com)
/// @date
/// @brief  Implements the FIRE2 model for stellar Feedback
///         as implemented into Enzo from Azton Wells
///


#ifdef NOTDEFINED

#ifndef ENZO_ENZO_METHOD_FIRE2_FEEDBACK
#define ENZO_ENZO_METHOD_FIRE2_FEEDBACK

class EnzoMethodFire2Feedback : public Method {

  /// @class   EnzoMethodFire2Feedback
  /// @ingroup Enzo
  /// @btief   [\ref Enzo] Encapsulate Feedback Routines

public:

  EnzoMethodFire2Feedback();

  /// Destructor
  virtual ~EnzoMethodFire2Feedback() throw() {};

  /// Charm++ Pup::able declarations
  PUPable_decl(EnzoMethodFire2Feedback);

  /// Charm++ Pup::able migration Constructor
  EnzoMethodFire2Feedback (CkMigrateMessage *m)
   : Method (m)
   {  }

   /// Charm++ Pack / Unpack function
   void pup(PUP::er &p);

   /// Apply the method
   virtual void compute (Block * block) throw();

   void compute_ (Block * block);

   /// name
   virtual std::string name() throw()
   { return "feedback"; }

   // Compute the maximum timestep for this method
   virtual double timestep (Block * block) throw();

protected:

  int sf_minimum_level_;
  int single_sn_;


}
#endif
