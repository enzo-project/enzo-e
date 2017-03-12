#include "fortran.h"
c=======================================================================
c//////////////////////  SUBROUTINE CIC_DEPOSIT  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine cic_deposit(posx, posy, posz, ndim, npositions, 
     &                      mass, field, leftedge, 
     &                      dim1, dim2, dim3, cellsize, cloudsize)

c
c  PERFORMS 1/2/3D CLOUD-IN-CELL INTERPOLATION FROM FIELD TO SUMFIELD
c
c  written by: Greg Bryan
c  date:       January, 1998
c  modified1:
c
c  PURPOSE: This routine performs a three-dimension, second-order
c           interpolation from field to sumfield (without clearing sumfield
c           first) at the positions specified.
c           The CIC cloud size is set by cloudsize, which should be
c           equal to or less than cellsize.  If cloudsize is equal to
c           cellsize, then this is the typical CIC deposit, but if it
c           is less than cellsize, then the effective cloud size is
c           smaller and the overlap calculations are adjusted.
c
c  INPUTS:
c     ndim       - dimensionality
c     cellsize   - the cell size of field
c     cloudsize  - size of the particle "cloud" (must be <= cellsize)
c     dim1,2,3   - real dimensions of field
c     leftedge   - the left edge(s) of field
c     npositions - number of particles
c     posx,y,z   - particle positions
c     sumfield   - 1D field (length npositions) of masses
c
c  OUTPUT ARGUMENTS: 
c     field      - field to be deposited to
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
      integer dim1, dim2, dim3, npositions, ndim
      ENZO_REAL posx(npositions), posy(npositions), posz(npositions),
     &        leftedge(3)
      ENZO_REAL    mass(npositions), field(dim1, dim2, dim3), cellsize,
     &          cloudsize
c
c  locals
c
      integer i1, j1, k1, i2,j2,k2,n
      ENZO_REAL    xpos, ypos, zpos, dx1, dy1, dz1, dx2,dy2,dz2
      ENZO_REAL edge1, edge2, edge3, fact, shift, half, refine
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c

!     write(0,*) npositions, leftedge, dim1, dim2, dim3, cellsize

      fact = 1.0/cellsize
      refine = cellsize/cloudsize
      half = 0.5001/refine
      shift = 0.5/refine
      edge1 = REAL(dim1) - half
      edge2 = REAL(dim2) - half
      edge3 = REAL(dim3) - half
c
      if (cloudsize .gt. cellsize) then
         write(6,*) "cloudsize > cellsize in cic_deposit!"
         ERROR_MESSAGE
      endif
c
c     1D
c
      if (ndim .eq. 1) then
c
         do n=1, npositions
c
c           Compute the position of the central cell
c
            xpos = min(max((posx(n) - leftedge(1))*fact, half), edge1)
c
c           Convert this into an integer index
c
            i1  = int(xpos - shift) + 1
c
c           Compute the weights
c
            dx1 = min((REAL(i1)+shift-xpos)*refine, 1.0)
c
c           Interpolate from field into sumfield
c
            field(i1  ,1,1) = field(i1  ,1,1) + mass(n)*dx1
            field(i1+1,1,1) = field(i1+1,1,1) + mass(n)*(1.0-dx1)
c
         enddo
c
      endif
c
c     2D
c
      if (ndim .eq. 2) then
c
         do n=1, npositions
c
c           Compute the position of the central cell
c
            xpos = min(max((posx(n) - leftedge(1))*fact, half), edge1)
            ypos = min(max((posy(n) - leftedge(2))*fact, half), edge2)
c
c           Convert this into an integer index
c
            i1 = int(xpos - shift) + 1
            j1 = int(ypos - shift) + 1
            i2 = i1 + 1
            j2 = j1 + 1
c
c           Compute the weights
c
            dx1 = min((REAL(i1)+shift-xpos)*refine, 1.0)
            dy1 = min((REAL(j1)+shift-ypos)*refine, 1.0)
            dx2 = (1.0-dx1)
            dy2 = (1.0-dy1)
c
c           Interpolate from field into sumfield
c
            field(i1,j1,1) = field(i1,j1,1) + mass(n)*dx1*dy1
            field(i2,j1,1) = field(i2,j1,1) + mass(n)*dx2*dy1
            field(i1,j2,1) = field(i1,j2,1) + mass(n)*dx1*dy2
            field(i2,j2,1) = field(i2,j2,1) + mass(n)*dx2*dy2
c
         enddo
c
      endif
c
c     3D
c
      if (ndim .eq. 3) then
c     
         do n=1, npositions
c     
c     Compute the position of the central cell
c     
            xpos = min(max((posx(n) - leftedge(1))*fact, half),edge1)
            ypos = min(max((posy(n) - leftedge(2))*fact, half),edge2)
            zpos = min(max((posz(n) - leftedge(3))*fact, half),edge3)
c     
c     Convert this into an integer index
c     
            i1 = int(xpos - shift) + 1
            j1 = int(ypos - shift) + 1
            k1 = int(zpos - shift) + 1
            i2 = i1 + 1
            j2 = j1 + 1
            k2 = k1 + 1
c     
c     Compute the weights
c     
            dx1 = min((REAL(i1)+shift-xpos)*refine, 1.0)
            dy1 = min((REAL(j1)+shift-ypos)*refine, 1.0)
            dz1 = min((REAL(k1)+shift-zpos)*refine, 1.0)
            dx2 = (1.0 - dx1)
            dy2 = (1.0 - dy1)
            dz2 = (1.0 - dz1)
c     
c     Interpolate from field into sumfield
c     

            field(i1,j1,k1) = field(i1,j1,k1) + mass(n)*dx1*dy1*dz1
            field(i2,j1,k1) = field(i2,j1,k1) + mass(n)*dx2*dy1*dz1
            field(i1,j2,k1) = field(i1,j2,k1) + mass(n)*dx1*dy2*dz1
            field(i2,j2,k1) = field(i2,j2,k1) + mass(n)*dx2*dy2*dz1
            field(i1,j1,k2) = field(i1,j1,k2) + mass(n)*dx1*dy1*dz2
            field(i2,j1,k2) = field(i2,j1,k2) + mass(n)*dx2*dy1*dz2
            field(i1,j2,k2) = field(i1,j2,k2) + mass(n)*dx1*dy2*dz2
            field(i2,j2,k2) = field(i2,j2,k2) + mass(n)*dx2*dy2*dz2
c     
         enddo
      
      endif
c
      return

      end
