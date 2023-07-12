#ifndef ENZO_FORTRAN_HPP
#define ENZO_FORTRAN_HPP

extern "C" void FORTRAN_NAME(interpolate)
  (int *rank, enzo_float *pfield, int pdim[],
   int pis[], int pie[], int r[],
   enzo_float *field, int dim[], int is[], enzo_float *work,
   int *imethod, int *posflag,
   int *ierror);

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
 
extern "C" void FORTRAN_NAME(expand_terms)(
   int *rank, int *isize, int *idual, enzo_float *coef,
   int *imethod, enzo_float *gamma,
   enzo_float *p,  enzo_float *d, enzo_float *e, enzo_float *ge,
   enzo_float *u, enzo_float *v, enzo_float *w,
   enzo_float *dold, enzo_float *eold, enzo_float *geold,
   enzo_float *uold, enzo_float *vold, enzo_float *wold,
   int *icr, enzo_float *ecr, enzo_float *ecrold);

extern "C" void FORTRAN_NAME(turboinit)
  (int *rank, int *nbox,
   enzo_float *u, enzo_float *v, enzo_float *w,
   int *in, int *jn, int *kn,
   int *ig, int *jg, int *kg);

extern "C" void FORTRAN_NAME(turboinit2d)
  (int *rank, int *nbox,
   enzo_float *u, enzo_float *v,
   int *in, int *jn,
   int *ig, int *jg);

#endif /* ENZO_FORTRAN_HPP */
