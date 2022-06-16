/// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodFeedbackSTARSS.hpp
/// @author     Will Hicks (whicks@ucsd.edu)
/// @date
/// @brief  Implements the STARSS model for stellar Feedback
///         as implemented into Enzo from Azton Wells
///


//#ifdef NOTDEFINED //uncomment when ready to compile

#ifndef ENZO_ENZO_METHOD_FEEDBACK_STARSS
#define ENZO_ENZO_METHOD_FEEDBACK_STARSS

class EnzoMethodFeedbackSTARSS : public Method {

  /// @class   EnzoMethodFeedbackSTARSS 
  /// @ingroup Enzo
  /// @btief   [\ref Enzo] Encapsulate Feedback Routines

public:

  EnzoMethodFeedbackSTARSS();

  /// Destructor
  virtual ~EnzoMethodFeedbackSTARSS() throw() {};

  /// Charm++ Pup::able declarations
  PUPable_decl(EnzoMethodFeedbackSTARSS);

  /// Charm++ Pup::able migration Constructor
  EnzoMethodFeedbackSTARSS (CkMigrateMessage *m)
   : Method (m)
   , ir_feedback_(-1)
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

   int determineSN (double age_Myr, int * nSNII, int * nSNIA,
                    double mass_Msun, double tunit, float dt);
   
   int determineWinds(double age_Myr, double * eWinds, double * mWinds, double * zWinds,
                      double mass_Msun, double metallicity_Zsun, double tunit, double dt); 

   // this can raise errors -- remove const throw() ???
   void deposit_feedback (Block * block, 
                          double ejectaEnergy, double ejectaMass, double ejectaMetals,
                          const enzo_float up, const enzo_float vp, const enzo_float wp,
                          const enzo_float xp, const enzo_float yp, const enzo_float zp,
                          const int ix, const int iy, const int iz,
                          const int winds, const int nSNII, const int nSNIa,
                          enzo_float starZ) const throw();

   void transformComovingWithStar(enzo_float * density,
                                  enzo_float * velocity_x, enzo_float * velocity_y, enzo_float * velocity_z,
                                  const enzo_float up, const enzo_float vp, const enzo_float wp,
                                  const int mx, const int my, const int mz, int direction) const throw();

   void add_accumulate_fields(EnzoBlock * enzo_block) throw();

   // window function for calculating CiC fractions. May make more sense to put this in Cello somewhere
   double Window(double xd, double yd, double zd, double width) const throw();

protected:

  int sf_minimum_level_;
  int single_sn_;
  int NEvents;
  double tiny_number; 

  // Refresh ID
  int ir_feedback_;

};

#endif
