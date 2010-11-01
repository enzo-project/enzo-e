#ifndef BASIC_BLOCK_HPP
#define BASIC_BLOCK_HPP

// enum face_type {
//   face_xm,
//   face_xp,
//   face_ym,
//   face_yp,
//   face_zm,
//   face_zp
// };

// class Block {
// public:
//   Block(int n,
// 	double xm, double xp,
// 	double ym, double yp,
// 	double zm, double zp ) 
//     : n_(n),
//       cycle_values_(0)
//   {
//     PARALLEL_PRINTF ("Creating Block(%d  %g %g  %g %g  %g %g)\n",n,
// 		     xm,xp,ym,yp,zm,zp);
//     values_ = new double [n*n*n];
//     for (int k=0; k<n*n*n; k++) values_[k] = 0.0;
//     for (int i=0; i<6; i++) {
//       cycle_ghosts_[i] = 0;
//       ghosts_[i]       = new double [n*n];
//       for (int k=0; k<n*n; k++) ghosts_[i][k] = 0.0;
//     }
//     lower_[0] = xm;
//     lower_[1] = ym;
//     lower_[2] = zm;
//     upper_[0] = xp;
//     upper_[1] = yp;
//     upper_[2] = zp;
//   };

//   ~Block()
//   {
//     PARALLEL_PRINTF ("Deleting Block(%d)\n",n_);
//     delete [] values_;
//     for (int i=0; i<6; i++) {
//       delete ghosts_[i];
//     }
//   };

//   Block (const Block & block)
//   {
//     PARALLEL_PRINTF ("Creating Block(Block(%d))\n",block.n_);
//     int n3 = n_*n_*n_; 
//     n_ = block.n_;
//     values_ = new double [n_*n_*n_];
//     for (int k=0; k<n3; k++) values_[k] = block.values_[k];
//     cycle_values_ = block.cycle_values_;
//     for (int i=0; i<6; i++) {
//       ghosts_[i] = new double [n_*n_];
//       for (int k=0; k<n3; k++) ghosts_[i][k] = block.ghosts_[i][k];
//     }
//     lower_[0] = block.lower_[0];
//     lower_[1] = block.lower_[1];
//     lower_[2] = block.lower_[2];
//     upper_[0] = block.upper_[0];
//     upper_[1] = block.upper_[1];
//     upper_[2] = block.upper_[2];
//   };


//   Block & operator= (const Block & block)
//   {
//     PARALLEL_PRINTF ("Assigning Block = Block(%d))\n",block.n_);
//     int n3 = n_*n_*n_; 
//     n_ = block.n_;
//     values_ = new double [n_*n_*n_];
//     for (int k=0; k<n3; k++) values_[k] = block.values_[k];
//     cycle_values_ = block.cycle_values_;
//     for (int i=0; i<6; i++) {
//       ghosts_[i] = new double [n_*n_];
//       for (int k=0; k<n3; k++) ghosts_[i][k] = block.ghosts_[i][k];
//     }
//     lower_[0] = block.lower_[0];
//     lower_[1] = block.lower_[1];
//     lower_[2] = block.lower_[2];
//     upper_[0] = block.upper_[0];
//     upper_[1] = block.upper_[1];
//     upper_[2] = block.upper_[2];
//     return *this;
//   };

// private:
//   int n_;
//   double * values_;      // values
//   double * ghosts_[6];   // ghosts
//   int cycle_values_;     // cycle number of values
//   int cycle_ghosts_[6];  // cycle number of ghosts
//   double lower_[3];      // lower corner
//   double upper_[3];      // upper corner
// };

#endif
