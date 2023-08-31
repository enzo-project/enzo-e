c     See LICENSE_PPML file for license and copyright information

#include "fortran.h"

      Subroutine Primitiv_ig(gamma,nx,ny,nz,qu1,qu2,qu3,qu4,qu5,qu6,
     &                    qu7,qu8,dn,vx,vy,vz,bx,by,bz,pn)
      Implicit NONE

      Integer nx,ny,nz,i,j,k
	  
      ENZO_REAL gamma
      ENZO_REAL dn(nx,ny,nz)  	 
      ENZO_REAL vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)  	 
      ENZO_REAL bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)  	 
      ENZO_REAL pn(nx,ny,nz)

      ENZO_REAL qu1(nx,ny,nz)  	 
      ENZO_REAL qu2(nx,ny,nz),qu3(nx,ny,nz),qu4(nx,ny,nz)  	 
      ENZO_REAL qu5(nx,ny,nz),qu6(nx,ny,nz),qu7(nx,ny,nz)  	 
      ENZO_REAL qu8(nx,ny,nz)
	  
      Do k=1,nz
         Do j=1,ny
            Do i=1,nx
               dn(i,j,k)=QU1(i,j,k)
               vx(i,j,k)=QU2(i,j,k)/QU1(i,j,k)
               vy(i,j,k)=QU3(i,j,k)/QU1(i,j,k)
               vz(i,j,k)=QU4(i,j,k)/QU1(i,j,k)
               bx(i,j,k)=QU5(i,j,k)
               by(i,j,k)=QU6(i,j,k)
               bz(i,j,k)=QU7(i,j,k)
               pn(i,j,k)=(gamma-1.)*(qu8(i,j,k)
     &         -(qu2(i,j,k)**2+qu3(i,j,k)**2+qu4(i,j,k)**2)/qu1(i,j,k)/2.-
     &         -(qu5(i,j,k)**2+qu6(i,j,k)**2+qu7(i,j,k)**2)/2.)
		  
            Enddo
         Enddo	
      Enddo
	   
      RETURN
      END
