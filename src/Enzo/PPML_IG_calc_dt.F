c     See LICENSE_PPML file for license and copyright information

#include "fortran.h"

      Subroutine calc_dt_ppml_ig(gamma,nx,ny,nz,
     1     i1,i2,j1,j2,k1,k2,
     1     dx,dy,dz,
     1     dn,vx,vy,vz,bx,by,bz,pn,
     1     b0,dt)
      
      Implicit NONE

      Integer nx,ny,nz,i,j,k,i1,j1,k1,i2,j2,k2
      
      REAL dn(nx,ny,nz)  	 
      REAL vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)  	 
      REAL bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)  	 
      REAL pn(nx,ny,nz)  	 
      
      REAL dtmx,dtmy,dtmz,rr0,bkb,vah,vax,vay,vaz,cg
      REAL cga,dsc,cfsx,cgs,taux,cfsy,tauy,cfsz,tauz,c0,dt	  
      Real dx,dy,dz,gamma,bx0,by0,bz0
      
      Real B0(3)

      DTMX=1.E+10
      DTMY=1.E+10
      DTMZ=1.E+10
      
      c0=0.8d0

      DO K=k1,k2
         DO J=j1,j2
            DO I=i1,i2

               Bx0=bx(I,J,K)+B0(1)
               By0=by(I,J,K)+B0(2)
               Bz0=bz(I,J,K)+B0(3)
               
               RR0=dn(I,J,K)
               BKB=Bx0**2+By0**2+Bz0**2
               VAH=BKB/RR0
               VAX=Bx0**2/RR0
               VAY=By0**2/RR0
               VAZ=Bz0**2/RR0
               CG=GAMMA*pn(I,J,K)/RR0
               CGA=CG+VAH
               DSC=CGA**2-4.d0*VAX*CG
               IF(DSC.LT.0.) DSC=0.
               CFSX=dsqrt(DSC)
               CGS=dsqrt((CGA+CFSX)/2.d0)
               TAUX=DX/(DABS(vx(I,J,K))+CGS)
               DSC=CGA**2-4.d0*VAY*CG
               IF(DSC.LT.0.) DSC=0.
               CFSY=dsqrt(DSC)
               CGS=dsqrt((CGA+CFSY)/2.d0)
               TAUY=DY/(DABS(vy(I,J,K))+CGS)
               DSC=CGA**2-4.d0*VAZ*CG
               IF(DSC.LT.0.) DSC=0.
               CFSZ=dsqrt(DSC)
               CGS=dsqrt((CGA+CFSZ)/2.d0)
               TAUZ=DZ/(DABS(vz(I,J,K))+CGS)

               DTMX=MIN(DTMX,TAUX)
               DTMY=MIN(DTMY,TAUY)
               DTMZ=MIN(DTMZ,TAUZ)

            ENDDO
         ENDDO
      ENDDO

      DT=C0/(1.D0/DTMX+1.D0/DTMY+1.D0/DTMZ)
      
      Return
      End
