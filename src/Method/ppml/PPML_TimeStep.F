      Subroutine calc_dt_ppml(nx,ny,nz,
     1                        i1,i2,j1,j2,k1,k2,
     1                        dx,dy,dz,
     1                        dn,vx,vy,vz,bx,by,bz,
     1                        dt)
	 
 	  Implicit NONE

      Integer nx,ny,nz,i,j,k,i1,j1,k1,i2,j2,k2
  
	  REAL dn(nx,ny,nz)  	 
      REAL vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)  	 
      REAL bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)  	 
	  
      REAL dtmx,dtmy,dtmz,rr0,bkb,vah,vax,vay,vaz,cg,a2
      REAL cga,dsc,cfsx,cgs,taux,cfsy,tauy,cfsz,tauz,c0,dt	  
	  Real dx,dy,dz

       DTMX=1.E+10
       DTMY=1.E+10
       DTMZ=1.E+10
	   
	   c0=0.8d0
       a2=1.0d0

      DO K=k1,k2
      DO J=j1,j2
      DO I=i1,i2

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
        TAUX=DX/(abs(vx(I,J,K))+CGS)
	    DSC=CGA**2-4.d0*VAY*CG
	    IF(DSC.LT.0.) DSC=0.
	    CFSY=sqrt(DSC)
	    CGS=sqrt((CGA+CFSY)/2.d0)
        TAUY=DY/(abs(vy(I,J,K))+CGS)
	    DSC=CGA**2-4.d0*VAZ*CG
	    IF(DSC.LT.0.) DSC=0.
	    CFSZ=sqrt(DSC)
	    CGS=sqrt((CGA+CFSZ)/2.d0)
        TAUZ=DZ/(abs(vz(I,J,K))+CGS)

        DTMX=MIN(DTMX,TAUX)
        DTMY=MIN(DTMY,TAUY)
        DTMZ=MIN(DTMZ,TAUZ)

      ENDDO
      ENDDO
      ENDDO

		DT=C0/(1.D0/DTMX+1.D0/DTMY+1.D0/DTMZ)
		
      Return
	  End