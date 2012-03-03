c     See LICENSE_PPML file for license and copyright information

#include "fortran.h"

      Subroutine calc_dt_ppml
     &     (nx,ny,nz,
     &     i1,i2,j1,j2,k1,k2,
     &     dx,dy,dz,
     &     dn,vx,vy,vz,bx,by,bz,
     &     dt)
      
      Implicit NONE

      Integer nx,ny,nz,i,j,k,i1,j1,k1,i2,j2,k2

      Integer ii1,ii2,jj1,jj2,kk1,kk2

      ENZO_REAL dn(nx,ny,nz)           
      ENZO_REAL vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)           
      ENZO_REAL bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)           
      
      ENZO_REAL dtmx,dtmy,dtmz,rr0,bkb,vah,vax,vay,vaz,cg,a2
      ENZO_REAL cga,dsc,cfsx,cgs,taux,cfsy,tauy,cfsz,tauz,c0,dt          
!      ENZO_REAL dx(nx),dy(ny),dz(nz)
      ENZO_REAL dx,dy,dz
      ENZO_REAL dmin,dmax
      ENZO_REAL vxmin,vxmax,vymin,vymax,vzmin,vzmax
      ENZO_REAL bxmin,bxmax,bymin,bymax,bzmin,bzmax
!      print*, 'nx,ny,nz=',nx,ny,nz
!      print*, 'i1,i2=',i1,i2
!      print*, 'j1,j2=',j1,j2
!      print*, 'k1,k2=',k1,k2
!!      print*, 'dx(i1)=',dx(1)
!!      print*, 'dy(j1)=',dy(1)
!!      print*, 'dz(k1)=',dz(1)
!      print*, 'dx=',dx
!      print*, 'dy=',dy
!      print*, 'dz=',dz
!      print*, 'dn(10,10,1)=',dn(10,10,1)
!      print*, 'vx(10,10,1)=',vx(10,10,1)
!      print*, 'vy(10,10,1)=',vy(10,10,1)
!      print*, 'vz(10,10,1)=',vz(10,10,1)
!      print*, 'bx(10,10,1)=',bx(10,10,1)
!      print*, 'by(10,10,1)=',by(10,10,1)
!      print*, 'bz(10,10,1)=',bz(10,10,1)
!      print*, 'dt=',dt

      DTMX=1.E+10
      DTMY=1.E+10
      DTMZ=1.E+10
      
      c0=0.8d0
      a2=1.0d0

      dmin = 1.e36
      dmax = -1e36
      vxmin = 1.e36
      vxmax = -1e36
      vymin = 1.e36
      vymax = -1e36
      vzmin = 1.e36
      vzmax = -1e36
      bxmin = 1.e36
      bxmax = -1e36
      bymin = 1.e36
      bymax = -1e36
      bzmin = 1.e36
      bzmax = -1e36

!     Adjust array bounds for C / C++
!     copies since original variables are call by reference

      ii1 = i1 + 1
      ii2 = i2 + 1
      jj1 = j1 + 1
      jj2 = j2 + 1
      kk1 = k1 + 1
      kk2 = k2 + 1

      DO K=kk1,kk2
         DO J=jj1,jj2
            DO I=ii1,ii2

               RR0=dn(I,J,K)
               BKB=bx(I,J,K)**2+by(I,J,K)**2+bz(I,J,K)**2
               VAH=BKB/RR0
               VAX=bx(I,J,K)**2/RR0
               VAY=by(I,J,K)**2/RR0
               VAZ=bz(I,J,K)**2/RR0
               CG=A2
               CGA=CG+VAH
               DSC=CGA**2-4.d0*VAX*CG
               IF(DSC.LT.0.) DSC=0.
               CFSX=sqrt(DSC)
               CGS=sqrt((CGA+CFSX)/2.d0)
!               TAUX=DX(ii1)/(abs(vx(I,J,K))+CGS)
               TAUX=DX/(abs(vx(I,J,K))+CGS)
               DSC=CGA**2-4.d0*VAY*CG
               IF(DSC.LT.0.) DSC=0.
               CFSY=sqrt(DSC)
               CGS=sqrt((CGA+CFSY)/2.d0)
!               TAUY=DY(jj1)/(abs(vy(I,J,K))+CGS)
               TAUY=DY/(abs(vy(I,J,K))+CGS)
               DSC=CGA**2-4.d0*VAZ*CG
               IF(DSC.LT.0.) DSC=0.
               CFSZ=sqrt(DSC)
               CGS=sqrt((CGA+CFSZ)/2.d0)
!               TAUZ=DZ(kk1)/(abs(vz(I,J,K))+CGS)
               TAUZ=DZ/(abs(vz(I,J,K))+CGS)

               DTMX=MIN(DTMX,TAUX)
               DTMY=MIN(DTMY,TAUY)
               DTMZ=MIN(DTMZ,TAUZ)
            ENDDO
         ENDDO
      ENDDO

      DT=C0/(1.D0/DTMX+1.D0/DTMY+1.D0/DTMZ)
!         
      Return
      End
