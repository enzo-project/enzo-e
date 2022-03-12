// See LICENSE_CELLO file for license and copyright information

/// @file     problem_ProlongLinear.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    Implentation of default linear prolongation

// #define CELLO_TRACE

// #define DEBUG_ARRAY
#include "problem.hpp"

// #define TRACE_SUMS

//----------------------------------------------------------------------

ProlongLinear::ProlongLinear() throw()
  : Prolong ()
{
  //  CkPrintf ("TRACE_PROLONG ProlongLinear()\n");
}

//----------------------------------------------------------------------

void ProlongLinear::apply 
( precision_type precision,
  void *       values_f, int mf3[3], int of3[3], int nf3[3],
  const void * values_c, int mc3[3], int oc3[3], int nc3[3],
  bool accumulate)
{
  TRACE6("ProlongLinear fine   %d:%d %d:%d %d:%d",
         of3[0],nf3[0]+of3[0],
         of3[1],nf3[1]+of3[1],
         of3[2],nf3[2]+of3[2]);

  TRACE6("ProlongLinear coarse %d:%d %d:%d %d:%d",
         oc3[0],nc3[0]+oc3[0],
         oc3[1],nc3[1]+oc3[1],
         oc3[2],nc3[2]+oc3[2]);

#ifdef TRACE_PROLONG
  CkPrintf("TRACE_PROLONG accum %d\n",accumulate?1:0);
  CkPrintf("TRACE_PROLONG nf %d %d %d\n", nf3[0],nf3[1],nf3[2]);
  CkPrintf("TRACE_PROLONG mf %d %d %d\n", mf3[0],mf3[1],mf3[2]);
  CkPrintf("TRACE_PROLONG of %d %d %d\n", of3[0],of3[1],of3[2]);
  CkPrintf("TRACE_PROLONG nc %d %d %d\n", nc3[0],nc3[1],nc3[2]);
  CkPrintf("TRACE_PROLONG mc %d %d %d\n", mc3[0],mc3[1],mc3[2]);
  CkPrintf("TRACE_PROLONG oc %d %d %d\n", oc3[0],oc3[1],oc3[2]);
  CkPrintf ("\n");
#endif

  switch (precision)  {

  case precision_single:

    apply_((float *)       values_f, mf3, of3, nf3,
           (const float *) values_c, mc3, oc3, nc3,
           accumulate);

    break;

  case precision_double:

    apply_((double *)       values_f, mf3, of3, nf3,
           (const double *) values_c, mc3, oc3, nc3,
           accumulate);

    break;

  default:

    ERROR1 ("ProlongLinear::apply()",
            "Unknown precision %d",
            precision);
  }
}

//----------------------------------------------------------------------

template <class T>
void ProlongLinear::apply_
(  T * values_f, int mf3[3], int of3[3], int nf3[3],
   const T * values_c, int mc3[3], int oc3[3], int nc3[3],
   bool accumulate)
{
  const int dcx = 1;
  const int dcy = mc3[0];
  const int dcz = mc3[0]*mc3[1];

  int rank = (mf3[1] == 1) ? 1 : ( (mf3[2] == 1) ? 2 : 3 );

  for (int i=0; i<rank; i++) {
    const char * xyz = "xyz";
    ASSERT3 ("ProlongLinear::apply_",
             "fine array %c-axis %d must be 2 times coarse axis %d",
             xyz[i],nf3[i],nc3[i],
             nf3[i]==2*nc3[i] || nf3[i]==2*(nc3[i]-2));

    // ASSERT2 ("ProlongLinear::apply_",
    //          "fine grid %c-axis %d must be divisible by 2",
    //          xyz[i],nf3[i],nf3[i] % 2 == 0);
    // ASSERT2 ("ProlongLinear::apply_",
    //          "fine grid %c-axis %d must be at least 4",
    //          xyz[i],nf3[i],nf3[i] >= 4);
  }

  // adjustment if coarse ghost cells available
  // NOTE:1 if ghosts not available , 0 if ghosts available
  
  int gcx = (nf3[0]==2*nc3[0]) ? 1 : 0;
  int gcy = (nf3[1]==2*nc3[1]) ? 1 : 0;
  int gcz = (nf3[2]==2*nc3[2]) ? 1 : 0;

#ifdef TRACE_SUMS
  T cmin=1e30,cmax=-1e30,cavg=0.0,ccount=0;
  T fmin=1e30,fmax=-1e30,favg=0.0,fcount=0;
#endif

  const int mf = mf3[0]*mf3[1]*mf3[2];
  
  T * temp_f = (accumulate) ? (new T [mf]) : values_f;

  if (rank == 1) {

    const int ofx = of3[0];
    const int ocx = oc3[0];
    const int nfx = nf3[0];

    for (int ifx = 0; ifx<nfx; ifx++) {

      int icx = ((ifx+1) >> 1) - gcx;

      // Default weighting factor
      int wx[2] = { 1, 3 };

      // Update weights if no ghosts and on edges
      if (ifx==0)     { icx += gcx; }
      if (ifx==nfx-1) { icx -= gcx; }
      if (ifx==0 || ifx==nfx-1) {
        wx[0] += 4*gcx;
        wx[1] -= 4*gcx;
      }

      T wx0 = 0.25*wx[ ifx&1];
      T wx1 = 0.25*wx[~ifx&1];

      int i_c = (ocx+icx) ;
      int i_f = (ofx+ifx) ;

      temp_f[i_f] = wx0*values_c[i_c]
        +             wx1*values_c[i_c + dcx ];
    }

  } else if (rank == 2) {

    const int ofy = of3[1];
    const int ocy = oc3[1];
    const int nfy = nf3[1];
    
    for (int ify = 0; ify<nfy; ify++) {

      int icy = ((ify+1) >> 1) - gcy;

      // Default weighting factor
      int wy[2] = { 1, 3 };

      // Update weights if no ghosts and on edges
      if (ify==0)     { icy += gcy; }
      if (ify==nfy-1) { icy -= gcy; }
      if (ify==0 || ify==nfy-1) {
        wy[0] += 4*gcy;
        wy[1] -= 4*gcy;
      }

      T wy0 = 0.25*wy[ ify&1];
      T wy1 = 0.25*wy[~ify&1];

      const int ofx = of3[0];
      const int ocx = oc3[0];
      const int nfx = nf3[0];

      const int mcx = mc3[0];
      const int mfx = mf3[0];

      for (int ifx = 0; ifx<nfx; ifx++) {

        int icx = ((ifx+1) >> 1) - gcx;

        // Default weighting factor
        int wx[2] = { 1, 3 };

        // Update weights if no ghosts and on edges
        if (ifx==0)     { icx += gcx; }
        if (ifx==nfx-1) { icx -= gcx; }
        if (ifx==0 || ifx==nfx-1) {
          wx[0] += 4*gcx;
          wx[1] -= 4*gcx;
        }

        T wx0 = 0.25*wx[ ifx&1];
        T wx1 = 0.25*wx[~ifx&1];

        int i_c = (ocx+icx) + mcx * ( (ocy+icy) );
        int i_f = (ofx+ifx) + mfx * ( (ofy+ify) );

        temp_f[i_f] = wx0*wy0*values_c[i_c]
          +             wx1*wy0*values_c[i_c + dcx ]
          +             wx0*wy1*values_c[i_c       + dcy ]
          +             wx1*wy1*values_c[i_c + dcx + dcy ];
      }
    }


  } else { // rank == 3

    const int ofz = of3[2];
    const int ocz = oc3[2];
    const int nfz = nf3[2];

    for (int ifz = 0; ifz<nfz; ifz++) {

      int icz = ((ifz+1) >> 1) - gcz;

      // Default weighting factor
      int wz[2] = { 1, 3 };

      // Update weights if no ghosts and on edges
      if (ifz==0)     { icz += gcz; }
      if (ifz==nfz-1) { icz -= gcz; }
      if (ifz==0 || ifz==nfz-1) {
        wz[0] += 4*gcz;
        wz[1] -= 4*gcz;
      }

      T wz0 = 0.25*wz[ ifz&1];
      T wz1 = 0.25*wz[~ifz&1];

      const int ofy = of3[1];
      const int ocy = oc3[1];
      const int nfy = nf3[1];

      for (int ify = 0; ify<nfy; ify++) {

        int icy = ((ify+1) >> 1) - gcy;

        // Default weighting factor
        int wy[2] = { 1, 3 };

        // Update weights if no ghosts and on edges
        if (ify==0)     { icy += gcy; }
        if (ify==nfy-1) { icy -= gcy; }
        if (ify==0 || ify==nfy-1) {
          wy[0] += 4*gcy;
          wy[1] -= 4*gcy;
        }

        T wy0 = 0.25*wy[ ify&1];
        T wy1 = 0.25*wy[~ify&1];

        const int ofx = of3[0];
        const int ocx = oc3[0];
        const int nfx = nf3[0];

        const int mcx = mc3[0];
        const int mfx = mf3[0];
        const int mcy = mc3[1];
        const int myf = mf3[1];
    
        for (int ifx = 0; ifx<nfx; ifx++) {

          int icx = ((ifx+1) >> 1) - gcx;

          // Default weighting factor
          int wx[2] = { 1, 3 };

          // Update weights if no ghosts and on edges
          if (ifx==0)     { icx += gcx; }
          if (ifx==nfx-1) { icx -= gcx; }
          if (ifx==0 || ifx==nfx-1) {
            wx[0] += 4*gcx;
            wx[1] -= 4*gcx;
          }

          T wx0 = 0.25*wx[ ifx&1];
          T wx1 = 0.25*wx[~ifx&1];

          int i_c = (ocx+icx) + mcx*( (ocy+icy) + mcy*(ocz+icz) );
          int i_f = (ofx+ifx) + mfx*( (ofy+ify) + myf*(ofz+ifz) );
	  
          temp_f[i_f] = wx0*wy0*wz0*values_c[i_c]
            +             wx1*wy0*wz0*values_c[i_c + dcx ]
            +             wx0*wy1*wz0*values_c[i_c       + dcy ]
            +             wx1*wy1*wz0*values_c[i_c + dcx + dcy ]
            +             wx0*wy0*wz1*values_c[i_c             + dcz ]
            +             wx1*wy0*wz1*values_c[i_c + dcx       + dcz ]
            +             wx0*wy1*wz1*values_c[i_c       + dcy + dcz ]
            +             wx1*wy1*wz1*values_c[i_c + dcx + dcy + dcz ];
        }
      }
    }
  }

#ifdef TRACE_SUMS
  if (accumulate) {
    DEBUG_PRINT_ARRAY0("values_c",values_c,mc3,nc3,oc3);
    DEBUG_PRINT_ARRAY0("values_f",values_f,mf3,nf3,of3);
  }
  //  if (accumulate) {
  for (int icz=0; icz<nc3[2]; icz++) {
    int iz=icz+oc3[2];
    for (int icy=0; icy<nc3[1]; icy++) {
      int iy=icy+oc3[1];
      for (int icx=0; icx<nc3[0]; icx++) {
        ccount++;
        int ix=icx+oc3[0];
        int i_c=ix + mc3[0]*(iy + mc3[1]*iz);
        cavg += values_c[i_c];
        cmin = std::min(cmin,values_c[i_c]);
        cmax = std::max(cmax,values_c[i_c]);
      }
    }
  }
  for (int ifz=0; ifz<nf3[2]; ifz++) {
    int iz=ifz+of3[2];
    for (int ify=0; ify<nf3[1]; ify++) {
      int iy=ify+of3[1];
      for (int ifx=0; ifx<nf3[0]; ifx++) {
        fcount++;
        int ix=ifx+of3[0];
        int i_f=ix + mf3[0]*(iy + mf3[1]*iz);
        favg += values_f[i_f];
        fmin = std::min(fmin,values_f[i_f]);
        fmax = std::max(fmax,values_f[i_f]);
      }
    }
  }
  
  favg/=fcount;
  cavg/=ccount;
  if (accumulate == 1) {
    CkPrintf ("TRACE_SUMS accum %d nc %d %d %d mn/avg/max %g %g %g\n",
              accumulate,
              nc3[0],nc3[1],nc3[2],
              //            oc3[0],oc3[1],oc3[2],
              //            mc3[0],mc3[1],mc3[2],
              cmin,cavg,cmax);
    CkPrintf ("TRACE_SUMS accum %d nf %d %d %d mn/avg/max %g %g %g\n",
              accumulate,
              nf3[0],nf3[1],nf3[2],
              //            of3[0],of3[1],of3[2],
              //            mf3[0],mf3[1],mf3[2],
              fmin,favg,fmax);
  }
#endif        

  if (accumulate) {
    int i0 = of3[0] + mf3[0]*(of3[1] + mf3[1]*of3[2]);
    for (int iz=0; iz<nf3[2]; iz++) {
      for (int iy=0; iy<nf3[1]; iy++) {
        for (int ix=0; ix<nf3[0]; ix++) {
          int i = i0 + ix + mf3[0]*(iy + mf3[1]*iz);
          values_f[i] += temp_f[i];
        }
      }
    }
    delete [] temp_f;
  }
}

//======================================================================

