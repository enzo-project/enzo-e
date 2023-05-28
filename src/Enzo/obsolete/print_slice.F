#include "fortran.h"
#define FORTRAN      
      
      subroutine print_slice(iphase,islice,jn,kn,gx,ncolor,
     &     dls,drs,flatten,
     &     pbar,    pls,
     &     prs, pslice,   ubar,
     &     uls,    urs,    vls,
     &     vrs,    wls,    wrs,
     &     diffcoef,     df,     ef,
     &     uf,     vf,     wf,
     &     gef,   gels,   gers,
     &     ges,
     &     dslice, eslice,grslice,
     &     uslice, vslice, wslice,
     &     geslice,
     &     colslice,  colf,
     &     colls, colrs)
      implicit none
#include "fortran_types.h"
      integer jn,kn,gx,iphase,ncolor
      R_PREC dls(jn,kn),    drs(jn,kn),flatten(jn,kn),
     &     pbar(jn,kn),    pls(jn,kn),
     &     prs(jn,kn), pslice(jn,kn),   ubar(jn,kn),
     &     uls(jn,kn),    urs(jn,kn),    vls(jn,kn),
     &     vrs(jn,kn),    wls(jn,kn),    wrs(jn,kn),
     &     diffcoef(jn,kn),     df(jn,kn),     ef(jn,kn),
     &     uf(jn,kn),     vf(jn,kn),     wf(jn,kn),
     &     gef(jn,kn),   gels(jn,kn),   gers(jn,kn),
     &     ges(jn,kn)
      R_PREC dslice(jn,kn), eslice(jn,kn),grslice(jn,kn),
     &     uslice(jn,kn), vslice(jn,kn), wslice(jn,kn),
     &     geslice(jn,kn)
      R_PREC colslice(jn,kn,ncolor),  colf(jn,kn,ncolor),
     &     colls(jn,kn,ncolor), colrs(jn,kn,ncolor)
      double precision sum(40)
      integer i,j,k,n,ic,islice

      do j=1,40
         sum(j) = 0.0D0
      end do

      n=0
C      print*, 'DEBUG_SWEEP jn,kn=',jn,kn
C      print*, 'DEBUG_SWEEP gx,ncolor,islice=',gx,ncolor,islice

      do j=1+gx,jn-gx
         do k=1+gx,kn-gx

            i=1
            sum(i) = sum(i) + abs(dls(j,k))
            i=2
            sum(i) = sum(i) + abs(dls(j,k))
            i=3
            sum(i) = sum(i) + abs(drs(j,k))
            i=4
            sum(i) = sum(i) + abs(flatten(j,k))
            
            i=5
            sum(i) = sum(i) + abs(pbar(j,k))
            i=6
            sum(i) = sum(i) + abs(pls(j,k))
            
            i=7
            sum(i) = sum(i) + abs(prs(j,k))
            i=8
            sum(i) = sum(i) + abs(pslice(j,k))
            i=9
            sum(i) = sum(i) + abs(ubar(j,k))
            
            i=10
            sum(i) = sum(i) + abs(uls(j,k))
            i=11
            sum(i) = sum(i) + abs(urs(j,k))
            i=12
            sum(i) = sum(i) + abs(vls(j,k))
            
            i=13
            sum(i) = sum(i) + abs(vrs(j,k))
            i=14
            sum(i) = sum(i) + abs(wls(j,k))
            i=15
            sum(i) = sum(i) + abs(wrs(j,k))
            
            i=16
            sum(i) = sum(i) + abs(diffcoef(j,k))
            i=17
            sum(i) = sum(i) + abs(df(j,k))
            i=18
            sum(i) = sum(i) + abs(ef(j,k))
            
            i=19
            sum(i) = sum(i) + abs(uf(j,k))
            i=20
            sum(i) = sum(i) + abs(vf(j,k))
            i=21
            sum(i) = sum(i) + abs(wf(j,k))
            
            i=22
            sum(i) = sum(i) + abs(gef(j,k))
            i=23
            sum(i) = sum(i) + abs(gels(j,k))
            i=24
            sum(i) = sum(i) + abs(gers(j,k))
            
            i=25
            sum(i) = sum(i) + abs(ges(j,k))
            i=26
            sum(i) = sum(i) + abs(dslice(j,k))
            i=27
            sum(i) = sum(i) + abs(eslice(j,k))
            i=28
            sum(i) = sum(i) + abs(grslice(j,k))
            
            i=29
            sum(i) = sum(i) + abs(uslice(j,k))
            i=30
            sum(i) = sum(i) + abs(vslice(j,k))
            i=31
            sum(i) = sum(i) + abs(wslice(j,k))
            
            i=32
            sum(i) = sum(i) + abs(geslice(j,k))

            do ic=1,ncolor
               i=33
               sum(i) = sum(i) + abs(colslice(j,k,ic))
               i=34
               sum(i) = sum(i) + abs(colf(j,k,ic))
               i=35
               sum(i) = sum(i) + abs(colls(j,k,ic))
               i=36
               sum(i) = sum(i) + abs(colrs(j,k,ic))
            end do
            n=i
            
         end do
      end do
C      print*, 'DEBUG_PPM SLICE n = ',n
      if (gx .le. islice-1 .and. islice-1 .lt.jn-gx) then
         do i=1,n
C            print*, 'DEBUG_PPM SLICE ',iphase,islice-gx-1,i,sum(i)
         end do
      end if
      
      return
      end


      subroutine print_a_slice(iphase,islice,jn,gx,slice)
      implicit none
#include "fortran_types.h"
      integer jn,kn,gx,iphase,ncolor
      R_PREC slice(jn)
      double precision sum
      integer i,j,k,n,ic,islice

      sum = 0.0D0

      do j=1+gx,jn-gx
         sum = sum + abs(slice(j))
      end do
C      print*, 'DEBUG_PPM A SLICE ',iphase,islice-gx-1,sum
      return
      end
      
