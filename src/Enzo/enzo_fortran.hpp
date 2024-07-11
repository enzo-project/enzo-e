#ifndef ENZO_FORTRAN_HPP
#define ENZO_FORTRAN_HPP

extern "C" void FORTRAN_NAME(interp3d)
  (enzo_float * parent, enzo_float * work, 
   int *dim1, int *dim2, int *dim3, 
   int *start1, int *start2, int *start3,
   int *end1, int *end2, int *end3, 
   int *refine1, int *refine2, int *refine3, 
   enzo_float *grid,
   int *gdim1, int *gdim2, int *gdim3, 
   int *gstart1, int *gstart2, int *gstart3,
   int *wdim1, int *wdim2, int *wdim3,
   int *ierror);

#endif /* ENZO_FORTRAN_HPP */
