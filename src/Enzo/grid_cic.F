#include "fortran.h"
c=======================================================================
c///////////////////////  SUBROUTINE DEP_GRID_CIC  \\\\\\\\\\\\\\\\\\\\\
c
      subroutine dep_grid_cic
     &     (source, dest, temp, velx, vely, velz, 
     &     dt, rfield, ndim, delx, dely, delz,
     &     sdim1, sdim2, sdim3, 
     &     sstart1, sstart2, sstart3, 
     &     send1, send2, send3,
     &     offset1, offset2, offset3,
     &     ddim1, ddim2, ddim3,
     &     refine1, refine2, refine3)
      use, intrinsic :: iso_fortran_env
c     
c  DEPOSIT SOURCE GRID INTO DEST GRID USING CIC INTERPOLATION
c
c  written by: Greg Bryan
c  date:       March, 1999
c  modified1:
c
c  PURPOSE:
c
c  INPUTS:
c     source       - source field
c     rfield       - source-like field indicating if cell is refined 
c                       (1=no, 0=yes)
c     sdim1-3      - source dimension
c     ddim1-3      - destination dimension
c     ndim         - rank of fields
c     refine1-3    - refinement factors
c     sstart1-3    - source start index
c     send1-3      - source end index
c     offset1-3     - offset from this grid edge to dest grid edge
c                    (>= 0, in dest cell units)
c     velx,y,z     - velocities
c     dt           - time step
c     delx         - cell size of source grid
c     temp         - temporary field, 4*size of dest
c
c  OUTPUT ARGUMENTS: 
c     dest         - prolonged field
c
c  EXTERNALS: 
c
c  LOCALS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer ddim1, ddim2, ddim3, sdim1, sdim2, sdim3, ndim,
     &        refine1, refine2, refine3, sstart1, sstart2, sstart3,
     &        send1, send2, send3
      ENZO_REAL
     &     source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3),
     &        rfield(sdim1, sdim2, sdim3),
     &        velx(sdim1, sdim2, sdim3), vely(sdim1, sdim2, sdim3),
     &        velz(sdim1, sdim2, sdim3), dt, delx, dely, delz,
     &        offset1, offset2, offset3,
     &        temp(ddim1, ddim2, ddim3, 4)
c
c  locals
c
      integer i, j, k, i1, j1, k1, n
      ENZO_REAL   fact1, fact2, fact3, x, y, z, dx, dy, dz,
     &     weight, mass,
     &        coef1, coef2, coef3, shift1, shift2, shift3,
     &     start1, start2, start3, half, edge1, edge2, edge3, temp1,
     &     ox,oy,oz
      parameter (half = 0.500e0)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////

C      print*, 'DEBUG GRID_CIC'
C      print*, source(1,1,1), dest(1,1,1), temp(1,1,1,1)
C      print*, velx(1,1,1), vely(1,1,1), velz(1,1,1)
C      print*, dt, rfield(1,1,1), ndim, delx, dely, delz
C      print*, sdim1, sdim2, sdim3
C      print*, sstart1, sstart2, sstart3
C      print*, send1, send2, send3
C      print*, offset1, offset2, offset3
C      print*, ddim1, ddim2, ddim3
C      print*, refine1, refine2, refine3

      do k=1,ddim3
         do j=1,ddim2
            do i=1,ddim1
               dest(i,j,k) = 0.e0
            enddo
            do n=1,ndim+1
               do i=1,ddim1
                  temp(i,j,k,n) = 0.e0
               enddo
            enddo
         enddo
      enddo
c
c     Precompute some things
c
      fact1 = 1.e0/REAL(refine1,4)
      fact2 = 1.e0/REAL(refine2,4)
      fact3 = 1.e0/REAL(refine3,4)
c
                       coef1 = dt/delx*fact1
      if (ndim .gt. 1) coef2 = dt/dely*fact2
      if (ndim .gt. 2) coef3 = dt/delz*fact3
c
C      print*, 'DEBUG_CIC coef = ',coef1,coef2,coef3
C      print*, 'DEBUG_CIC dt = ',dt
      start1 = -sstart1 - 0.5e0 + offset1*REAL(refine1,4)
      start2 = -sstart2 - 0.5e0 + offset2*REAL(refine2,4)
      start3 = -sstart3 - 0.5e0 + offset3*REAL(refine3,4)
c
      edge1 = REAL(ddim1,4) - half
      edge2 = REAL(ddim2,4) - half
      edge3 = REAL(ddim3,4) - half

!     write(0,'("grid_cic: ",6f12.4)') start1,start2,start3,
!    &                                 edge1,edge2,edge3

c
c     a) 1D
c
      if (ndim .eq. 1) then
         weight = fact1
c
c        compute density and mass-weighted velocity field
c
         do i=sstart1+1, send1+1
            x = min(max((start1 + i)*fact1, half), edge1)
            i1 = int(x + 0.5e0)
            dx = REAL(i1,4) + 0.5e0 - x
            mass = (1.e0 - rfield(i,1,1))*weight*source(i,1,1)
c
            ox = (1.e0-dx)
            temp(i1  ,1  ,1, 1) = temp(i1  ,1  ,1, 1) + 
     &                                 mass*dx
            temp(i1+1,1  ,1, 1) = temp(i1+1,1  ,1, 1) + 
     &                                 mass*ox
c
            temp1 = velx(i,1,1)*mass
            temp(i1  ,1  ,1, 2) = temp(i1  ,1  ,1, 2) + 
     &                     temp1*dx
            temp(i1+1,1  ,1, 2) = temp(i1+1,1  ,1, 2) + 
     &                     temp1*ox
         enddo
c
c        Use velocity and mass field to generate mass field advanced by dt
c
         do i=1,ddim1
            shift1 = temp(i,1,1,2)/max(temp(i,1,1,1), tiny) * coef1
            x = min(max((i - 0.5e0 + shift1), half), edge1)
            i1 = int(x + 0.5e0)
            dx = REAL(i1,4) + 0.5e0 - x
            ox = (1.e0-dx)
            mass = temp(i,1,1,1)
            dest(i1  ,1  ,1) = dest(i1  ,1  ,1) + mass*dx
            dest(i1+1,1  ,1) = dest(i1+1,1  ,1) + mass*ox
         enddo
      endif
c
c     b) 2D
c
      if (ndim .eq. 2) then
         weight = fact1*fact2
c
c        compute density and mass-weighted velocity field
c
         do j=sstart2+1, send2+1
            y = min(max((start2 + j)*fact2, half), edge2)
            j1 = int(y + 0.5e0)
            dy = REAL(j1,4) + 0.5e0 - y
            do i=sstart1+1, send1+1
               x = min(max((start1 + i)*fact1, half), edge1)
               i1 = int(x + 0.5e0)
               dx = REAL(i1,4) + 0.5e0 - x
               ox = (1.e0-dx)
               oy = (1.e0-dy)
               mass = (1.e0 - rfield(i,j,1))*weight*source(i,j,1)
c
               temp(i1  ,j1  ,1,1) = temp(i1  ,j1  ,1,1) +
     &                             mass*dx*dy
               temp(i1+1,j1  ,1,1) = temp(i1+1,j1  ,1,1) +
     &                             mass*ox*dy
               temp(i1  ,j1+1,1,1) = temp(i1  ,j1+1,1,1) +
     &                             mass*dx*oy
               temp(i1+1,j1+1,1,1) = temp(i1+1,j1+1,1,1) +
     &                             mass*(1. 0-dx)*oy
c
               temp1 = velx(i,j,1)*mass
               temp(i1  ,j1  ,1,2) = temp(i1  ,j1  ,1,2) +
     &                 temp1*dx*dy
               temp(i1+1,j1  ,1,2) = temp(i1+1,j1  ,1,2) +
     &                 temp1*ox*dy
               temp(i1  ,j1+1,1,2) = temp(i1  ,j1+1,1,2) +
     &                 temp1*dx*oy
               temp(i1+1,j1+1,1,2) = temp(i1+1,j1+1,1,2) +
     &                 temp1*ox*oy
c
               temp1 = vely(i,j,1)*mass
               temp(i1  ,j1  ,1,3) = temp(i1  ,j1  ,1,3) +
     &                 temp1*dx*dy
               temp(i1+1,j1  ,1,3) = temp(i1+1,j1  ,1,3) +
     &                 temp1*ox*dy
               temp(i1  ,j1+1,1,3) = temp(i1  ,j1+1,1,3) +
     &                 temp1*dx*oy
               temp(i1+1,j1+1,1,3) = temp(i1+1,j1+1,1,3) +
     &                 temp1*ox*oy
c
            enddo
         enddo
c
c        Use velocity and mass field to generate mass field advanced by dt
c
         do j=1, ddim2
            do i=1, ddim1
               shift2 = temp(i,j,1,3)/
     &                  max(temp(i,j,1,1),tiny)*coef2
               y = min(max((j - 0.5e0 + shift2), half), edge2)
               j1 = int(y + 0.5e0)
               dy = REAL(j1,4) + 0.5e0 - y
               shift1 = temp(i,j,1,2)/max(temp(i,j,1,1),tiny)*coef1
               x = min(max((i - 0.5e0 + shift1), half), edge1)
               i1 = int(x + 0.5e0)
               dx = REAL(i1,4) + 0.5e0 - x
               ox = (1.e0-dx)
               oy = (1.e0-dy)
               mass = temp(i,j,1,1)
               dest(i1  ,j1  ,1) = dest(i1  ,j1  ,1) +
     &                             mass*dx*dy
               dest(i1+1,j1  ,1) = dest(i1+1,j1  ,1) +
     &                             mass*ox*dy
               dest(i1  ,j1+1,1) = dest(i1  ,j1+1,1) +
     &                             mass*dx*oy
               dest(i1+1,j1+1,1) = dest(i1+1,j1+1,1) +
     &                             mass*ox*oy
            enddo
         enddo
      endif
c
c     c) 3D
c
      if (ndim .eq. 3) then
         weight = fact1*fact2*fact3
c
c        compute density and mass-weighted velocity field
c
         do k=sstart3+1, send3+1
            z = min(max((start3 + k)*fact3, half), edge3)
            k1 = int(z + 0.5e0)
            dz = REAL(k1) + 0.5e0 - z
            do j=sstart2+1, send2+1
               y = min(max((start2 + j)*fact2, half), edge2)
               j1 = int(y + 0.5e0)
               dy = REAL(j1) + 0.5e0 - y
               do i=sstart1+1, send1+1
                  x = min(max((start1 + i)*fact1, half), edge1)
                  i1 = int(x + 0.5e0)
                  dx = REAL(i1) + 0.5e0 - x
                  mass = (1.e0 - rfield(i,j,k))*weight*source(i,j,k)
                  ox = (1.e0-dx)
                  oy = (1.e0-dy)
                  oz = (1.e0-dz)
c
C                  if (i.eq.sstart1+4.and.
C     &                 j.eq.sstart2+4.and.
C     &                 k.eq.sstart3+4) then
C                     print*, 'DEBUG_CIC dxyz=',dx,dy,dz
C                     print*, 'DEBUG_CIC vel123=',
C     &                    velx(i,j,k),vely(i,j,k),velz(i,j,k)
C                     print*, 'DEBUG_CIC mass=',mass
C                     print*, 'DEBUG_CIC xyz=',x,y,z
C                     print*, 'DEBUG_CIC ijk1=',i1,j1,k1
C                     print*, 'DEBUG_CIC fact123=',fact1,fact2,fact3
C                     print*, 'DEBUG_CIC start123=',start1,start2,start3
C                     print*, 'DEBUG_CIC edge123=',edge1,edge2,edge3
C                  endif
                  temp(i1  ,j1  ,k1  ,1) = temp(i1  ,j1  ,k1  ,1) +
     &                       mass*dx*dy*dz
                  temp(i1+1,j1  ,k1  ,1) = temp(i1+1,j1  ,k1  ,1) +
     &                       mass*ox*dy*dz
                  temp(i1  ,j1+1,k1  ,1) = temp(i1  ,j1+1,k1  ,1) +
     &                       mass*dx*oy*dz
                  temp(i1+1,j1+1,k1  ,1) = temp(i1+1,j1+1,k1  ,1) +
     &                       mass*ox*oy*dz
                  temp(i1  ,j1  ,k1+1,1) = temp(i1  ,j1  ,k1+1,1) +
     &                       mass*dx*dy*oz
                  temp(i1+1,j1  ,k1+1,1) = temp(i1+1,j1  ,k1+1,1) +
     &                       mass*ox*dy*oz
                  temp(i1  ,j1+1,k1+1,1) = temp(i1  ,j1+1,k1+1,1) +
     &                       mass*dx*oy*oz
                  temp(i1+1,j1+1,k1+1,1) = temp(i1+1,j1+1,k1+1,1) +
     &                 mass*ox*oy*oz
c
                  temp1 = velx(i,j,k)*mass
                  temp(i1  ,j1  ,k1  ,2) = temp(i1  ,j1  ,k1  ,2) +
     &                       temp1*dx*dy*dz
                  temp(i1+1,j1  ,k1  ,2) = temp(i1+1,j1  ,k1  ,2) +
     &                       temp1*ox*dy*dz
                  temp(i1  ,j1+1,k1  ,2) = temp(i1  ,j1+1,k1  ,2) +
     &                       temp1*dx*oy*dz
                  temp(i1+1,j1+1,k1  ,2) = temp(i1+1,j1+1,k1  ,2) +
     &                       temp1*ox*oy*dz
                  temp(i1  ,j1  ,k1+1,2) = temp(i1  ,j1  ,k1+1,2) +
     &                       temp1*dx*dy*oz
                  temp(i1+1,j1  ,k1+1,2) = temp(i1+1,j1  ,k1+1,2) +
     &                       temp1*ox*dy*oz
                  temp(i1  ,j1+1,k1+1,2) = temp(i1  ,j1+1,k1+1,2) +
     &                       temp1*dx*oy*oz
                  temp(i1+1,j1+1,k1+1,2) = temp(i1+1,j1+1,k1+1,2) +
     &                 temp1*ox*oy*oz
c
                  temp1 = vely(i,j,k)*mass
                  temp(i1  ,j1  ,k1  ,3) = temp(i1  ,j1  ,k1  ,3) +
     &                       temp1*dx*dy*dz
                  temp(i1+1,j1  ,k1  ,3) = temp(i1+1,j1  ,k1  ,3) +
     &                       temp1*ox*dy*dz
                  temp(i1  ,j1+1,k1  ,3) = temp(i1  ,j1+1,k1  ,3) +
     &                       temp1*dx*oy*dz
                  temp(i1+1,j1+1,k1  ,3) = temp(i1+1,j1+1,k1  ,3) +
     &                       temp1*ox*oy*dz
                  temp(i1  ,j1  ,k1+1,3) = temp(i1  ,j1  ,k1+1,3) +
     &                       temp1*dx*dy*oz
                  temp(i1+1,j1  ,k1+1,3) = temp(i1+1,j1  ,k1+1,3) +
     &                       temp1*ox*dy*oz
                  temp(i1  ,j1+1,k1+1,3) = temp(i1  ,j1+1,k1+1,3) +
     &                       temp1*dx*oy*oz
                  temp(i1+1,j1+1,k1+1,3) = temp(i1+1,j1+1,k1+1,3) +
     &                 temp1*ox*oy*oz
c
                  temp1 = velz(i,j,k)*mass
                  temp(i1  ,j1  ,k1  ,4) = temp(i1  ,j1  ,k1  ,4) +
     &                       temp1*dx*dy*dz
                  temp(i1+1,j1  ,k1  ,4) = temp(i1+1,j1  ,k1  ,4) +
     &                       temp1*ox*dy*dz
                  temp(i1  ,j1+1,k1  ,4) = temp(i1  ,j1+1,k1  ,4) +
     &                       temp1*dx*oy*dz
                  temp(i1+1,j1+1,k1  ,4) = temp(i1+1,j1+1,k1  ,4) +
     &                       temp1*ox*oy*dz
                  temp(i1  ,j1  ,k1+1,4) = temp(i1  ,j1  ,k1+1,4) +
     &                       temp1*dx*dy*oz
                  temp(i1+1,j1  ,k1+1,4) = temp(i1+1,j1  ,k1+1,4) +
     &                       temp1*ox*dy*oz
                  temp(i1  ,j1+1,k1+1,4) = temp(i1  ,j1+1,k1+1,4) +
     &                       temp1*dx*oy*oz
                  temp(i1+1,j1+1,k1+1,4) = temp(i1+1,j1+1,k1+1,4) +
     &                 temp1*ox*oy*oz
c
               enddo
            enddo
         enddo
c
c        Use velocity and mass field to generate mass field advanced by dt
c
         do k=1, ddim3
            do j=1, ddim2
               do i=1, ddim1
                 shift3 = temp(i,j,k,4)/max(temp(i,j,k,1),tiny)*coef3
                  z = min(max((k - 0.5e0 + shift3), half), edge3)
                  k1 = int(z + 0.5e0)
                  dz = REAL(k1) + 0.5e0 - z
c
                 shift2 = temp(i,j,k,3)/max(temp(i,j,k,1),tiny)*coef2
                  y = min(max((j - 0.5e0 + shift2), half), edge2)
                  j1 = int(y + 0.5e0)
                  dy = REAL(j1) + 0.5e0 - y
c
                 shift1 = temp(i,j,k,2)/max(temp(i,j,k,1),tiny)*coef1
                  x = min(max((i - 0.5e0 + shift1), half), edge1)
                  i1 = int(x + 0.5e0)
                  dx = REAL(i1) + 0.5e0 - x
c
c                  if (i1 .lt. 1 .or. i1 .ge. ddim1 .or.
c     &                j1 .lt. 1 .or. j1 .ge. ddim2 .or.
c     &                k1 .lt. 1 .or. k1 .ge. ddim3    )
c     &                write(6,*) i1,j1,k1,ddim1,ddim2,ddim3
c
                  mass = temp(i,j,k,1)
                  ox = (1.e0-dx)
                  oy = (1.e0-dy)
                  oz = (1.e0-dz)
                  dest(i1  ,j1  ,k1  ) = dest(i1  ,j1  ,k1  ) +
     &                       mass*dx*dy*dz
                  dest(i1+1,j1  ,k1  ) = dest(i1+1,j1  ,k1  ) +
     &                       mass*ox*dy*dz
                  dest(i1  ,j1+1,k1  ) = dest(i1  ,j1+1,k1  ) +
     &                       mass*dx*oy*dz
                  dest(i1+1,j1+1,k1  ) = dest(i1+1,j1+1,k1  ) +
     &                       mass*ox*oy*dz
                  dest(i1  ,j1  ,k1+1) = dest(i1  ,j1  ,k1+1) +
     &                       mass*dx*dy*oz
                  dest(i1+1,j1  ,k1+1) = dest(i1+1,j1  ,k1+1) +
     &                       mass*ox*dy*oz
                  dest(i1  ,j1+1,k1+1) = dest(i1  ,j1+1,k1+1) +
     &                       mass*dx*oy*oz
                  dest(i1+1,j1+1,k1+1) = dest(i1+1,j1+1,k1+1) +
     &                       mass*ox*oy*oz
               enddo
            enddo
         enddo
c
      endif
c
      return
      end
