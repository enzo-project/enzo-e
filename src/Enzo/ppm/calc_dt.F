c     $Id$
c     See LICENSE_ENZO file for license and copyright information

#include "fortran.h"
c=======================================================================
c////////////////////////  SUBROUTINE CALC_DT  \\\\\\\\\\\\\\\\\\\\\\\\\
c
c
      subroutine calc_dt(rank, idim, jdim, kdim, 
     &                   i1, i2, j1, j2, k1, k2, 
     &                   dx, dy, dz, gamma, ipfree, aye,
     &                   d, p, u, v, w, dt)
c
c  COMPUTES TIME STEP FOR NEXT CYCLE
c
c     written by: Greg Bryan
c     date:       February, 1996
c     modified1:  Alexei Kritsuk, Jan. 2001; changed the expression for
c                 hydro time step for PPM_DE. Now it follows the linear 
c                 stability condition for multidimensional Godunov scheme 
c                 (Godunov 1959). It is now safe to use CourantSafetyNumber =
c                 0.8 for PPM_DE runs.
c
c  PURPOSE:  Computes the new timestep using the Courant condition.
c            (For rank < 3, the unused fields and cell widths may be
c             null)
c
c  INPUTS:
c    rank    - rank of fields
c    i,j,dim - declared dimensions of fields
c    i,j,k1  - start index of active region in fields (0 based)
c    i,j,k2  - end index of active region in fields (0 based)
c    dx,y,z  - cell widths along each dimension
c    gamma   - ratio of specific heats
c    ipfree  - pressure free flag (1 = on, 0 = off)
c    aye     - expansion factor (or 1 if not using comoving coordinates)
c    d,p     - density and pressure fields
c    u,v,w   - velocity fields (x,y,z)
c
c  OUTPUTS:
c    dt      - minimum allowed dt (without Courant safety factor)
c
c  LOCALS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c     Arguments
c
      integer idim, jdim, kdim, i1, i2, j1, j2, k1, k2, rank, ipfree
      ENZO_REAL dx, dy, dz
      ENZO_REAL dt, gamma, aye
      ENZO_REAL d(idim,jdim,kdim), p(idim,jdim,kdim), u(idim,jdim,kdim),
     &        v(idim,jdim,kdim), w(idim,jdim,kdim)
c
c     Locals
c
      integer i,j,k
      ENZO_REAL cs, dt1
      ENZO_REAL dmax,dmin,pmax,pmin,umax,umin,vmax,vmin
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c     Set initial timestep to a large number
c
C      print*, rank, idim, jdim, kdim
C      print*, i1, i2, j1, j2, k1, k2
C      print*, 'dx = ',dx, 'dy = ',dy
C      print*, gamma, ipfree, aye
C      print*, 'd(i1+1,j1+1,1) = ', d(i1+1,j1+1,1)
C      print*, 'd(i2+1,j2+1,1) = ', d(i2+1,j2+1,1)
C      print*, 'p(i1+1,j1+1,1) = ', p(i1+1,j1+1,1)
C      print*, 'p(i2+1,j2+1,1) = ', p(i2+1,j2+1,1)
C      print*, 'u(i1+1,j1+1,1) = ', u(i1+1,j1+1,1)
C      print*, 'u(i2+1,j2+1,1) = ', u(i2+1,j2+1,1)
C      print*, 'v(i1+1,j1+1,1) = ', v(i1+1,j1+1,1)
C      print*, 'v(i2+1,j2+1,1) = ', v(i2+1,j2+1,1)
C      print*, dt

      dmax=-huge
      pmax=-huge
      umax=-huge
      vmax=-huge
      dmin=huge
      pmin=huge
      umin=huge
      vmin=huge
      do j = j1+1, j2+1
         do i = i1+1, i2+1
            dmin = min(dmin,d(i,j,1))
            pmin = min(pmin,p(i,j,1))
            umin = min(umin,u(i,j,1))
            vmin = min(vmin,v(i,j,1))
            dmax = max(dmax,d(i,j,1))
            pmax = max(pmax,p(i,j,1))
            umax = max(umax,u(i,j,1))
            vmax = max(vmax,v(i,j,1))
         end do
      end do

C      print*, 'dmin = ',dmin
C      print*, 'pmin = ',pmin
C      print*, 'umin = ',umin
C      print*, 'vmin = ',vmin
C      print*, 'dmax = ',dmax
C      print*, 'pmax = ',pmax
C      print*, 'umax = ',umax
C      print*, 'vmax = ',vmax

      dt = huge
c
c     one-dimensional version
c
      if (rank .eq. 1) then
         do i = i1+1, i2+1
            cs = max(sqrt(gamma*p(i,1,1)/d(i,1,1)), tiny)
            if (ipfree .eq. 1) cs = tiny
            dt = min(dt, real(dx*aye/(cs + abs(u(i,1,1)))))
         enddo
      endif
c
c     two-dimensional version
c
      if (rank .eq. 2) then
         do j = j1+1, j2+1
            do i = i1+1, i2+1
               cs = max(sqrt(gamma*p(i,j,1)/d(i,j,1)), tiny)
               if (ipfree .eq. 1) cs = tiny

c
c 	Here is Godunov's formula (to make sure ppm works with 0.8)
c
               dt = min(dt, real(aye/((cs + abs(u(i,j,1)))/dx+
     &                       (cs + abs(v(i,j,1)))/dy) ))
C               print*, i,j,p(i,j,1),d(i,j,1),
C     &              u(i,j,1),v(i,j,1),dx,dy
c
c       The way it was originally in enzo
c
c               dt = min(dt, real(dx*aye/(cs + abs(u(i,j,1)))),
c     &                      real(dy*aye/(cs + abs(v(i,j,1)))))
            enddo
c
         enddo
      endif
c
c     three-dimensional version
c
      if (rank .eq. 3) then
         do k = k1+1, k2+1
            do j = j1+1, j2+1
               do i = i1+1, i2+1
                  if (d(i,j,k) .ne. d(i,j,k) .or.
     &                p(i,j,k) .ne. p(i,j,k))
     &               write(6,*) 'calc_dt',d(i,j,k),p(i,j,k),i,j,k
                  cs = max(sqrt(gamma*p(i,j,k)/d(i,j,k)), tiny)
                  if (ipfree .eq. 1) cs = tiny

c     
c     Godunov's formula.
c     
                  dt1 = aye/((cs + abs(u(i,j,k)))/dx +
     &                 (cs + abs(v(i,j,k)))/dy +
     &                 (cs + abs(w(i,j,k)))/dz)
c     
c     The way it was originally in enzo
c     
c     dt1 = min(dx*aye/(cs + abs(u(i,j,k))),
c     &                    dy*aye/(cs + abs(v(i,j,k))),
c     &                    dz*aye/(cs + abs(w(i,j,k))))

                  dt = min(dt, dt1)
!                  if (dt1 .lt. 1.0e-5) write(6,1000) dt1,d(i,j,k),
!    &                  p(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k)
! 1000             format('calc_dt (small dt): dt,d,p,uvw=',1p,6e12.3)
               enddo
c
            enddo
         enddo
      endif
c
      return
      end
