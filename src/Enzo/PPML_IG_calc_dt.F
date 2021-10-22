c     See LICENSE_PPML file for license and copyright information

C#include "fortran.h"
C
C      Subroutine calc_dt_ppml_ig(gamma,nx,ny,nz,
C     1     i1,i2,j1,j2,k1,k2,
C     1     dx,dy,dz,
C     1     dn,vx,vy,vz,bx,by,bz,pn,
C     1     b0,dt)
C      
C      Implicit NONE
C
C      Integer nx,ny,nz,i,j,k,i1,j1,k1,i2,j2,k2
C      
C      REAL dn(nx,ny,nz)  	 
C      REAL vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)  	 
C      REAL bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)  	 
C      REAL pn(nx,ny,nz)  	 
C      
C      REAL dtmx,dtmy,dtmz,rr0,bkb,vah,vax,vay,vaz,cg
C      REAL cga,dsc,cfsx,cgs,taux,cfsy,tauy,cfsz,tauz,c0,dt	  
C      Real dx,dy,dz,gamma,bx0,by0,bz0
C      
C      Real B0(3)
C
C      DTMX=1.E+10
C      DTMY=1.E+10
C      DTMZ=1.E+10
C      
C      c0=0.8d0
C
C      DO K=k1,k2
C         DO J=j1,j2
C            DO I=i1,i2
C
C               Bx0=bx(I,J,K)+B0(1)
C               By0=by(I,J,K)+B0(2)
C               Bz0=bz(I,J,K)+B0(3)
C               
C               RR0=dn(I,J,K)
C               BKB=Bx0**2+By0**2+Bz0**2
C               VAH=BKB/RR0
C               VAX=Bx0**2/RR0
C               VAY=By0**2/RR0
C               VAZ=Bz0**2/RR0
C               CG=GAMMA*pn(I,J,K)/RR0
C               CGA=CG+VAH
C               DSC=CGA**2-4.d0*VAX*CG
C               IF(DSC.LT.0.) DSC=0.
C               CFSX=dsqrt(DSC)
C               CGS=dsqrt((CGA+CFSX)/2.d0)
C               TAUX=DX/(DABS(vx(I,J,K))+CGS)
C               DSC=CGA**2-4.d0*VAY*CG
C               IF(DSC.LT.0.) DSC=0.
C               CFSY=dsqrt(DSC)
C               CGS=dsqrt((CGA+CFSY)/2.d0)
C               TAUY=DY/(DABS(vy(I,J,K))+CGS)
C               DSC=CGA**2-4.d0*VAZ*CG
C               IF(DSC.LT.0.) DSC=0.
C               CFSZ=dsqrt(DSC)
C               CGS=dsqrt((CGA+CFSZ)/2.d0)
C               TAUZ=DZ/(DABS(vz(I,J,K))+CGS)
C
C               DTMX=MIN(DTMX,TAUX)
C               DTMY=MIN(DTMY,TAUY)
C               DTMZ=MIN(DTMZ,TAUZ)
C
C            ENDDO
C         ENDDO
C      ENDDO
C
C      DT=C0/(1.D0/DTMX+1.D0/DTMY+1.D0/DTMZ)
C      
C      Return
C      End
