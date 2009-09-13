	Subroutine QDD6(Nx,Qvr,Qvl,Qrx,Qlx,Qp)
	Implicit NONE
	
	integer i,m,nx
	REAL Qvr(Nx,7),Qvl(Nx,7),Qp(Nx,7) 
	REAL Qrx(Nx,7),Qlx(Nx,7) 
	REAL DQ,Q6

	DO I=1,Nx	
	   
	   DO M=1,7

	      IF((QVR(I,M)-QP(I,M))*(QP(I,M)-QVL(I,M)).LE.0.)THEN
		 QRX(I,M)=QP(I,M)
		 QLX(I,M)=QP(I,M)
	      ELSE

		 DQ=QVR(I,M)-QVL(I,M)
		 Q6=6.0D0*(QP(I,M)-(QVR(I,M)+QVL(I,M))/2.D0)

		 QLX(I,M)=QVL(I,M)
		 QRX(I,M)=QVR(I,M)
		 
		 IF(DQ*Q6.GT.DQ**2) QLX(I,M)=3.*QP(I,M)-2.*QVR(I,M)
		 IF(DQ*Q6.LT.-DQ**2) QRX(I,M)=3.*QP(I,M)-2.*QVL(I,M)
		 
	      ENDIF  
	      
	   ENDDO 
	ENDDO

	Return
	END	
