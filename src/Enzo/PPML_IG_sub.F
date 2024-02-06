c     See LICENSE_PPML file for license and copyright information

#include "fortran.h"

      SUBROUTINE VECTLRx_ig(gamma,QU,BS,QL,QR,U)
      Implicit NONE

      ENZO_REAL QU(8),QL(8,8),QR(8,8),U(8),BS(3)  	
      ENZO_REAL DLT,X22,SK1,BVX,BVY,BVZ,VAX,VAY,VAZ
      ENZO_REAL VAH,CKV,CSR,CSV,CWH,DSCV,VF2
      ENZO_REAL VS2,CF,CS,CA,BXX,BYZ,SBB,BMG,GPP,GPB
      ENZO_REAL ASS,AS,AF,BY,BZ,BSGN,CD2,GAMMA,QU5,QU6,QU7  
      
      DLT=1.E-12
      X22=1./SQRT(2.)
	  
      QU5=QU(5)+BS(1)  
      QU6=QU(6)+BS(2)  
      QU7=QU(7)+BS(3)  
 	
      SK1=SQRT(QU(1))
      BVX=QU5/SK1
      BVY=QU6/SK1
      BVZ=QU7/SK1
      VAX=BVX**2
      VAY=BVY**2
      VAZ=BVZ**2
      VAH=VAX+VAY+VAZ
      CKV=GAMMA*QU(8)/QU(1)
      CSR=SQRT(CKV)
      CSV=CKV+VAH
      CWH=CKV-VAH
      DSCV=CSV**2-4.*CKV*VAX
      IF(DSCV.LT.0.) DSCV=0.
      DSCV=SQRT(DSCV)
      VF2=(CSV+DSCV)/2.
      VS2=(CSV-DSCV)/2.
      IF(VS2.LT.0.) VS2=0.
      CF=SQRT(VF2)
      CS=SQRT(VS2)
      CA=ABS(BVX)

      BXX=QU5**2
      BYZ=QU6**2+QU7**2
      SBB=SQRT(BYZ)
      BMG=BXX+BYZ
      GPP=GAMMA*QU(8)
      GPB=ABS(GPP-BXX)

      IF(BYZ.GT.DLT*BMG) THEN

         ASS=SQRT(CWH**2+4.*CKV*(VAY+VAZ))
         AS=SQRT(ABS((1.-CWH/ASS)/2.))
         AF=SQRT(ABS((1.+CWH/ASS)/2.))

         BY=QU6/SBB
         BZ=QU7/SBB

      ELSE

         IF(GPB.GT.DLT*GPP) THEN
            AS=SQRT(ABS((1.-SIGN(1.,CWH))/2.))
            AF=SQRT(ABS((1.+SIGN(1.,CWH))/2.))
            CS=SQRT(ABS((CKV+VAX-SIGN(1.,CWH)*CWH)/2.))
            CF=SQRT(ABS((CKV+VAX+SIGN(1.,CWH)*CWH)/2.))
         ELSE
            AS=X22
            AF=X22
            CS=CSR
            CF=CSR  
         ENDIF
         BY=X22
         BZ=X22

      ENDIF

	
      BSGN=SIGN(1.,QU5)
      CD2=2.*CKV
	
      QL(1,1)=0.
      QL(1,2)=-QU(1)*AF*CF/CD2
      QL(1,3)=QU(1)*AS*CS*BY*BSGN/CD2
      QL(1,4)=QU(1)*AS*CS*BZ*BSGN/CD2
      QL(1,5)=0.
      QL(1,6)=SK1*CSR*AS*BY/CD2
      QL(1,7)=SK1*CSR*AS*BZ/CD2
      QL(1,8)=AF/CD2

      QL(2,1)=0.
      QL(2,2)=0.
      QL(2,3)=-QU(1)*BZ*BSGN/2.
      QL(2,4)=QU(1)*BY*BSGN/2.
      QL(2,5)=0.
      QL(2,6)=-SK1*BZ/2.
      QL(2,7)=SK1*BY/2.
      QL(2,8)=0.

      QL(3,1)=0.
      QL(3,2)=-QU(1)*AS*CS/CD2
      QL(3,3)=-QU(1)*AF*CF*BY*BSGN/CD2
      QL(3,4)=-QU(1)*AF*CF*BZ*BSGN/CD2
      QL(3,5)=0.
      QL(3,6)=-SK1*CSR*AF*BY/CD2
      QL(3,7)=-SK1*CSR*AF*BZ/CD2
      QL(3,8)=AS/CD2

      QL(4,1)=1.
      QL(4,2)=0.
      QL(4,3)=0.
      QL(4,4)=0.
      QL(4,5)=0.
      QL(4,6)=0.
      QL(4,7)=0.
      QL(4,8)=-1./CKV

      QL(5,1)=0.
      QL(5,2)=QU(1)*AS*CS/CD2
      QL(5,3)=QU(1)*AF*CF*BY*BSGN/CD2
      QL(5,4)=QU(1)*AF*CF*BZ*BSGN/CD2
      QL(5,5)=0.
      QL(5,6)=-SK1*CSR*AF*BY/CD2
      QL(5,7)=-SK1*CSR*AF*BZ/CD2
      QL(5,8)=AS/CD2

      QL(6,1)=0.
      QL(6,2)=0.
      QL(6,3)=-QU(1)*BZ*BSGN/2.
      QL(6,4)=QU(1)*BY*BSGN/2.
      QL(6,5)=0.
      QL(6,6)=SK1*BZ/2.
      QL(6,7)=-SK1*BY/2.
      QL(6,8)=0.

      QL(7,1)=0.
      QL(7,2)=QU(1)*AF*CF/CD2
      QL(7,3)=-QU(1)*AS*CS*BY*BSGN/CD2
      QL(7,4)=-QU(1)*AS*CS*BZ*BSGN/CD2
      QL(7,5)=0.
      QL(7,6)=SK1*CSR*AS*BY/CD2
      QL(7,7)=SK1*CSR*AS*BZ/CD2
      QL(7,8)=AF/CD2
	
      QL(8,1)=0.
      QL(8,2)=0.
      QL(8,3)=0.
      QL(8,4)=0.
      QL(8,5)=1.
      QL(8,6)=0.
      QL(8,7)=0.
      QL(8,8)=0.

      QR(1,1)=AF
      QR(1,2)=0.
      QR(1,3)=AS
      QR(1,4)=1.
      QR(1,5)=AS
      QR(1,6)=0.
      QR(1,7)=AF
      QR(1,8)=0.

      QR(2,1)=-AF*CF/QU(1)
      QR(2,2)=0.
      QR(2,3)=-AS*CS/QU(1)
      QR(2,4)=0.
      QR(2,5)=AS*CS/QU(1)
      QR(2,6)=0.
      QR(2,7)=AF*CF/QU(1)
      QR(2,8)=0.

      QR(3,1)=AS*CS*BY*BSGN/QU(1)
      QR(3,2)=-BZ*BSGN/QU(1)
      QR(3,3)=-AF*CF*BY*BSGN/QU(1)
      QR(3,4)=0.
      QR(3,5)=AF*CF*BY*BSGN/QU(1)
      QR(3,6)=-BZ*BSGN/QU(1)
      QR(3,7)=-AS*CS*BY*BSGN/QU(1)
      QR(3,8)=0.

      QR(4,1)=AS*CS*BZ*BSGN/QU(1)
      QR(4,2)=BY*BSGN/QU(1)
      QR(4,3)=-AF*CF*BZ*BSGN/QU(1)
      QR(4,4)=0.
      QR(4,5)=AF*CF*BZ*BSGN/QU(1)
      QR(4,6)=BY*BSGN/QU(1)
      QR(4,7)=-AS*CS*BZ*BSGN/QU(1)
      QR(4,8)=0.

      QR(5,1)=0.
      QR(5,2)=0.
      QR(5,3)=0.
      QR(5,4)=0.
      QR(5,5)=0.
      QR(5,6)=0.
      QR(5,7)=0.
      QR(5,8)=1.

      QR(6,1)=AS*BY*CSR/SK1
      QR(6,2)=-BZ/SK1
      QR(6,3)=-AF*BY*CSR/SK1
      QR(6,4)=0.
      QR(6,5)=-AF*BY*CSR/SK1
      QR(6,6)=BZ/SK1
      QR(6,7)=AS*BY*CSR/SK1
      QR(6,8)=0.

      QR(7,1)=AS*BZ*CSR/SK1
      QR(7,2)=BY/SK1
      QR(7,3)=-AF*BZ*CSR/SK1
      QR(7,4)=0.
      QR(7,5)=-AF*BZ*CSR/SK1
      QR(7,6)=-BY/SK1
      QR(7,7)=AS*BZ*CSR/SK1
      QR(7,8)=0.

      QR(8,1)=AF*CKV
      QR(8,2)=0.
      QR(8,3)=AS*CKV
      QR(8,4)=0.
      QR(8,5)=AS*CKV
      QR(8,6)=0.
      QR(8,7)=AF*CKV
      QR(8,8)=0.

      U(1)=QU(2)-CF 
      U(2)=QU(2)-CA
      U(3)=QU(2)-CS
      U(4)=QU(2) 
      U(5)=QU(2)+CS
      U(6)=QU(2)+CA
      U(7)=QU(2)+CF 
      U(8)=QU(2)

      RETURN
      END

      SUBROUTINE VECTLRy_ig(gamma,QU,BS,QL,QR,U)
      Implicit NONE

      ENZO_REAL QU(8),QL(8,8),QR(8,8),U(8),BS(3)  	
      ENZO_REAL DLT,X22,SK1,BVX,BVY,BVZ,VAX,VAY,VAZ
      ENZO_REAL VAH,CKV,CSR,CSV,CWH,DSCV,VF2
      ENZO_REAL VS2,CF,CS,CA,BYY,BXZ,SBB,BMG,GPP,GPB
      ENZO_REAL ASS,AS,AF,BX,BZ,BSGN,CD2,GAMMA,QU5,QU6,QU7
      
      DLT=1.E-12
      X22=1./SQRT(2.)

      QU5=QU(5)+BS(1)  
      QU6=QU(6)+BS(2)  
      QU7=QU(7)+BS(3)  
 	
      SK1=SQRT(QU(1))
      BVX=QU5/SK1
      BVY=QU6/SK1
      BVZ=QU7/SK1
      VAX=BVX**2
      VAY=BVY**2
      VAZ=BVZ**2
      VAH=VAX+VAY+VAZ
      CKV=GAMMA*QU(8)/QU(1)
      CSR=SQRT(CKV)
      CSV=CKV+VAH
      CWH=CKV-VAH
      DSCV=CSV**2-4.*CKV*VAY
      IF(DSCV.LT.0.) DSCV=0.
      DSCV=SQRT(DSCV)
      VF2=(CSV+DSCV)/2.
      VS2=(CSV-DSCV)/2.
      IF(VS2.LT.0.) VS2=0.
      CF=SQRT(VF2)
      CS=SQRT(VS2)
      CA=ABS(BVY)

      BYY=QU6**2
      BXZ=QU5**2+QU7**2
      SBB=SQRT(BXZ)
      BMG=BYY+BXZ
      GPP=GAMMA*QU(8)
      GPB=ABS(GPP-BYY)

      IF(BXZ.GT.DLT*BMG) THEN

         ASS=SQRT(CWH**2+4.*CKV*(VAX+VAZ))
         AS=SQRT(ABS((1.-CWH/ASS)/2.))
         AF=SQRT(ABS((1.+CWH/ASS)/2.))

         BX=QU5/SBB
         BZ=QU7/SBB
          
      ELSE

         IF(GPB.GT.DLT*GPP) THEN
            AS=SQRT(ABS((1.-SIGN(1.,CWH))/2.))
            AF=SQRT(ABS((1.+SIGN(1.,CWH))/2.))
            CS=SQRT(ABS((CKV+VAY-SIGN(1.,CWH)*CWH)/2.))
            CF=SQRT(ABS((CKV+VAY+SIGN(1.,CWH)*CWH)/2.))
         ELSE
            AS=X22
            AF=X22
            CS=CSR
            CF=CSR  
         ENDIF

         BX=X22
         BZ=X22

      ENDIF

      BSGN=SIGN(1.,QU6)
      CD2=2.*CKV
	
      QL(1,1)=0.
      QL(1,2)=QU(1)*AS*CS*BX*BSGN/CD2
      QL(1,3)=-QU(1)*AF*CF/CD2
      QL(1,4)=QU(1)*AS*CS*BZ*BSGN/CD2
      QL(1,5)=SK1*CSR*AS*BX/CD2
      QL(1,6)=0.
      QL(1,7)=SK1*CSR*AS*BZ/CD2
      QL(1,8)=AF/CD2
	
      QL(2,1)=0.
      QL(2,2)=-QU(1)*BZ*BSGN/2.
      QL(2,3)=0.
      QL(2,4)=QU(1)*BX*BSGN/2.
      QL(2,5)=-SK1*BZ/2.
      QL(2,6)=0.
      QL(2,7)=SK1*BX/2.
      QL(2,8)=0.

      QL(3,1)=0.
      QL(3,2)=-QU(1)*AF*CF*BX*BSGN/CD2
      QL(3,3)=-QU(1)*AS*CS/CD2
      QL(3,4)=-QU(1)*AF*CF*BZ*BSGN/CD2
      QL(3,5)=-SK1*CSR*AF*BX/CD2
      QL(3,6)=0.
      QL(3,7)=-SK1*CSR*AF*BZ/CD2
      QL(3,8)=AS/CD2

      QL(4,1)=1.
      QL(4,2)=0.
      QL(4,3)=0.
      QL(4,4)=0.
      QL(4,5)=0.
      QL(4,6)=0.
      QL(4,7)=0.
      QL(4,8)=-1./CKV

      QL(5,1)=0.
      QL(5,2)=QU(1)*AF*CF*BX*BSGN/CD2
      QL(5,3)=QU(1)*AS*CS/CD2
      QL(5,4)=QU(1)*AF*CF*BZ*BSGN/CD2
      QL(5,5)=-SK1*CSR*AF*BX/CD2
      QL(5,6)=0.
      QL(5,7)=-SK1*CSR*AF*BZ/CD2
      QL(5,8)=AS/CD2

      QL(6,1)=0.
      QL(6,2)=-QU(1)*BZ*BSGN/2.
      QL(6,3)=0.
      QL(6,4)=QU(1)*BX*BSGN/2.
      QL(6,5)=SK1*BZ/2.
      QL(6,6)=0.
      QL(6,7)=-SK1*BX/2.
      QL(6,8)=0.

      QL(7,1)=0.
      QL(7,2)=-QU(1)*AS*CS*BX*BSGN/CD2
      QL(7,3)=QU(1)*AF*CF/CD2
      QL(7,4)=-QU(1)*AS*CS*BZ*BSGN/CD2
      QL(7,5)=SK1*CSR*AS*BX/CD2
      QL(7,6)=0.
      QL(7,7)=SK1*CSR*AS*BZ/CD2
      QL(7,8)=AF/CD2
	
      QL(8,1)=0.
      QL(8,2)=0.
      QL(8,3)=0.
      QL(8,4)=0.
      QL(8,5)=0.
      QL(8,6)=1.
      QL(8,7)=0.
      QL(8,8)=0.

      QR(1,1)=AF
      QR(1,2)=0.
      QR(1,3)=AS
      QR(1,4)=1.
      QR(1,5)=AS
      QR(1,6)=0.
      QR(1,7)=AF
      QR(1,8)=0.

      QR(2,1)=AS*CS*BX*BSGN/QU(1)
      QR(2,2)=-BZ*BSGN/QU(1)
      QR(2,3)=-AF*CF*BX*BSGN/QU(1)
      QR(2,4)=0.
      QR(2,5)=AF*CF*BX*BSGN/QU(1)
      QR(2,6)=-BZ*BSGN/QU(1)
      QR(2,7)=-AS*CS*BX*BSGN/QU(1)
      QR(2,8)=0.

      QR(3,1)=-AF*CF/QU(1)
      QR(3,2)=0.
      QR(3,3)=-AS*CS/QU(1)
      QR(3,4)=0.
      QR(3,5)=AS*CS/QU(1)
      QR(3,6)=0.
      QR(3,7)=AF*CF/QU(1)
      QR(3,8)=0.

      QR(4,1)=AS*CS*BZ*BSGN/QU(1)
      QR(4,2)=BX*BSGN/QU(1)
      QR(4,3)=-AF*CF*BZ*BSGN/QU(1)
      QR(4,4)=0.
      QR(4,5)=AF*CF*BZ*BSGN/QU(1)
      QR(4,6)=BX*BSGN/QU(1)
      QR(4,7)=-AS*CS*BZ*BSGN/QU(1)
      QR(4,8)=0.

      QR(5,1)=AS*BX*CSR/SK1
      QR(5,2)=-BZ/SK1
      QR(5,3)=-AF*BX*CSR/SK1
      QR(5,4)=0.
      QR(5,5)=-AF*BX*CSR/SK1
      QR(5,6)=BZ/SK1
      QR(5,7)=AS*BX*CSR/SK1
      QR(5,8)=0.

      QR(6,1)=0.
      QR(6,2)=0.
      QR(6,3)=0.
      QR(6,4)=0.
      QR(6,5)=0.
      QR(6,6)=0.
      QR(6,7)=0.
      QR(6,8)=1.

      QR(7,1)=AS*BZ*CSR/SK1
      QR(7,2)=BX/SK1
      QR(7,3)=-AF*BZ*CSR/SK1
      QR(7,4)=0.
      QR(7,5)=-AF*BZ*CSR/SK1
      QR(7,6)=-BX/SK1
      QR(7,7)=AS*BZ*CSR/SK1
      QR(7,8)=0.

      QR(8,1)=AF*CKV
      QR(8,2)=0.
      QR(8,3)=AS*CKV
      QR(8,4)=0.
      QR(8,5)=AS*CKV
      QR(8,6)=0.
      QR(8,7)=AF*CKV
      QR(8,8)=0.

      U(1)=QU(3)-CF 
      U(2)=QU(3)-CA
      U(3)=QU(3)-CS
      U(4)=QU(3) 
      U(5)=QU(3)+CS
      U(6)=QU(3)+CA
      U(7)=QU(3)+CF 
      U(8)=QU(3)

      RETURN
      END

      SUBROUTINE VECTLRz_ig(gamma,QU,BS,QL,QR,U)
      Implicit NONE

      ENZO_REAL QU(8),QL(8,8),QR(8,8),U(8),BS(3)  	
      ENZO_REAL DLT,X22,SK1,BVX,BVY,BVZ,VAX,VAY,VAZ
      ENZO_REAL VAH,CKV,CSR,CSV,CWH,DSCV,VF2
      ENZO_REAL VS2,CF,CS,CA,BYY,BXZ,SBB,BMG,GPP,GPB
      ENZO_REAL ASS,AS,AF,BX,BY,BSGN,CD2,GAMMA
      ENZO_REAL QU5,QU6,QU7
      
      DLT=1.E-12
      X22=1./SQRT(2.)

      QU5=QU(5)+BS(1)  
      QU6=QU(6)+BS(2)  
      QU7=QU(7)+BS(3)  
 	
      SK1=SQRT(QU(1))
      BVX=QU5/SK1
      BVY=QU6/SK1
      BVZ=QU7/SK1
      VAX=BVX**2
      VAY=BVY**2
      VAZ=BVZ**2
      VAH=VAX+VAY+VAZ
      CKV=GAMMA*QU(8)/QU(1)
      CSR=SQRT(CKV)
      CSV=CKV+VAH
      CWH=CKV-VAH
      DSCV=CSV**2-4.*CKV*VAZ
      IF(DSCV.LT.0.) DSCV=0.
      DSCV=SQRT(DSCV)
      VF2=(CSV+DSCV)/2.
      VS2=(CSV-DSCV)/2.
      IF(VS2.LT.0.) VS2=0.
      CF=SQRT(VF2)
      CS=SQRT(VS2)
      CA=ABS(BVZ)

      BYY=QU7**2
      BXZ=QU5**2+QU6**2
      SBB=SQRT(BXZ)
      BMG=BYY+BXZ
      GPP=GAMMA*QU(8)
      GPB=ABS(GPP-BYY)

      IF(BXZ.GT.DLT*BMG) THEN

         ASS=SQRT(CWH**2+4.*CKV*(VAX+VAY))
         AS=SQRT(ABS((1.-CWH/ASS)/2.))
         AF=SQRT(ABS((1.+CWH/ASS)/2.))

         BX=QU5/SBB
         BY=QU6/SBB

      ELSE

         IF(GPB.GT.DLT*GPP) THEN
            AS=SQRT(ABS((1.-SIGN(1.,CWH))/2.))
            AF=SQRT(ABS((1.+SIGN(1.,CWH))/2.))
            CS=SQRT(ABS((CKV+VAZ-SIGN(1.,CWH)*CWH)/2.))
            CF=SQRT(ABS((CKV+VAZ+SIGN(1.,CWH)*CWH)/2.))
         ELSE
            AS=X22
            AF=X22
            CS=CSR
            CF=CSR  
         ENDIF

         BX=X22
         BY=X22

      ENDIF

      BSGN=SIGN(1.,QU7)
      CD2=2.*CKV
	
      QL(1,1)=0.
      QL(1,2)=QU(1)*AS*CS*BX*BSGN/CD2
      QL(1,3)=QU(1)*AS*CS*BY*BSGN/CD2
      QL(1,4)=-QU(1)*AF*CF/CD2
      QL(1,5)=SK1*CSR*AS*BX/CD2
      QL(1,6)=SK1*CSR*AS*BY/CD2
      QL(1,7)=0.
      QL(1,8)=AF/CD2
	
      QL(2,1)=0.
      QL(2,2)=-QU(1)*BY*BSGN/2.
      QL(2,3)=QU(1)*BX*BSGN/2.
      QL(2,4)=0.
      QL(2,5)=-SK1*BY/2.
      QL(2,6)=SK1*BX/2.
      QL(2,7)=0.
      QL(2,8)=0.

      QL(3,1)=0.
      QL(3,2)=-QU(1)*AF*CF*BX*BSGN/CD2
      QL(3,3)=-QU(1)*AF*CF*BY*BSGN/CD2
      QL(3,4)=-QU(1)*AS*CS/CD2
      QL(3,5)=-SK1*CSR*AF*BX/CD2
      QL(3,6)=-SK1*CSR*AF*BY/CD2
      QL(3,7)=0.
      QL(3,8)=AS/CD2

      QL(4,1)=1.
      QL(4,2)=0.
      QL(4,3)=0.
      QL(4,4)=0.
      QL(4,5)=0.
      QL(4,6)=0.
      QL(4,7)=0.
      QL(4,8)=-1./CKV

      QL(5,1)=0.
      QL(5,2)=QU(1)*AF*CF*BX*BSGN/CD2
      QL(5,3)=QU(1)*AF*CF*BY*BSGN/CD2
      QL(5,4)=QU(1)*AS*CS/CD2
      QL(5,5)=-SK1*CSR*AF*BX/CD2
      QL(5,6)=-SK1*CSR*AF*BY/CD2
      QL(5,7)=0.
      QL(5,8)=AS/CD2

      QL(6,1)=0.
      QL(6,2)=-QU(1)*BY*BSGN/2.
      QL(6,3)=QU(1)*BX*BSGN/2.
      QL(6,4)=0.
      QL(6,5)=SK1*BY/2.
      QL(6,6)=-SK1*BX/2.
      QL(6,7)=0.
      QL(6,8)=0.

      QL(7,1)=0.
      QL(7,2)=-QU(1)*AS*CS*BX*BSGN/CD2
      QL(7,3)=-QU(1)*AS*CS*BY*BSGN/CD2
      QL(7,4)=QU(1)*AF*CF/CD2
      QL(7,5)=SK1*CSR*AS*BX/CD2
      QL(7,6)=SK1*CSR*AS*BY/CD2
      QL(7,7)=0.
      QL(7,8)=AF/CD2
	
      QL(8,1)=0.
      QL(8,2)=0.
      QL(8,3)=0.
      QL(8,4)=0.
      QL(8,5)=0.
      QL(8,6)=0.
      QL(8,7)=1.
      QL(8,8)=0.

      QR(1,1)=AF
      QR(1,2)=0.
      QR(1,3)=AS
      QR(1,4)=1.
      QR(1,5)=AS
      QR(1,6)=0.
      QR(1,7)=AF
      QR(1,8)=0.

      QR(2,1)=AS*CS*BX*BSGN/QU(1)
      QR(2,2)=-BY*BSGN/QU(1)
      QR(2,3)=-AF*CF*BX*BSGN/QU(1)
      QR(2,4)=0.
      QR(2,5)=AF*CF*BX*BSGN/QU(1)
      QR(2,6)=-BY*BSGN/QU(1)
      QR(2,7)=-AS*CS*BX*BSGN/QU(1)
      QR(2,8)=0.

      QR(3,1)=AS*CS*BY*BSGN/QU(1)
      QR(3,2)=BX*BSGN/QU(1)
      QR(3,3)=-AF*CF*BY*BSGN/QU(1)
      QR(3,4)=0.
      QR(3,5)=AF*CF*BY*BSGN/QU(1)
      QR(3,6)=BX*BSGN/QU(1)
      QR(3,7)=-AS*CS*BY*BSGN/QU(1)
      QR(3,8)=0.

      QR(4,1)=-AF*CF/QU(1)
      QR(4,2)=0.
      QR(4,3)=-AS*CS/QU(1)
      QR(4,4)=0.
      QR(4,5)=AS*CS/QU(1)
      QR(4,6)=0.
      QR(4,7)=AF*CF/QU(1)
      QR(4,8)=0.

      QR(5,1)=AS*BX*CSR/SK1
      QR(5,2)=-BY/SK1
      QR(5,3)=-AF*BX*CSR/SK1
      QR(5,4)=0.
      QR(5,5)=-AF*BX*CSR/SK1
      QR(5,6)=BY/SK1
      QR(5,7)=AS*BX*CSR/SK1
      QR(5,8)=0.

      QR(6,1)=AS*BY*CSR/SK1
      QR(6,2)=BX/SK1
      QR(6,3)=-AF*BY*CSR/SK1
      QR(6,4)=0.
      QR(6,5)=-AF*BY*CSR/SK1
      QR(6,6)=-BX/SK1
      QR(6,7)=AS*BY*CSR/SK1
      QR(6,8)=0.

      QR(7,1)=0.
      QR(7,2)=0.
      QR(7,3)=0.
      QR(7,4)=0.
      QR(7,5)=0.
      QR(7,6)=0.
      QR(7,7)=0.
      QR(7,8)=1.

      QR(8,1)=AF*CKV
      QR(8,2)=0.
      QR(8,3)=AS*CKV
      QR(8,4)=0.
      QR(8,5)=AS*CKV
      QR(8,6)=0.
      QR(8,7)=AF*CKV
      QR(8,8)=0.

      U(1)=QU(4)-CF 
      U(2)=QU(4)-CA
      U(3)=QU(4)-CS
      U(4)=QU(4) 
      U(5)=QU(4)+CS
      U(6)=QU(4)+CA
      U(7)=QU(4)+CF 
      U(8)=QU(4)

      RETURN
      END


      SUBROUTINE VECTEGx_ig(gamma,QU,BS,U)
      Implicit NONE

      ENZO_REAL QU(8),U(8),BS(3)  	
      ENZO_REAL DLT,SK1,BVX,BVY,BVZ,VAX,VAY,VAZ
      ENZO_REAL VAH,CKV,A1,A2,CSR,CSV,CWH,DSCV,VF2
      ENZO_REAL VS2,CF,CS,VA,BXX,BYZ,SBB,BMG,GPP,GPB
      ENZO_REAL ASS,AS,AF,GAMMA,QU5,QU6,QU7

      DLT=1.E-12

      QU5=QU(5)+BS(1)  
      QU6=QU(6)+BS(2)  
      QU7=QU(7)+BS(3)  
	  
      SK1=SQRT(QU(1))
      BVX=QU5/SK1
      BVY=QU6/SK1
      BVZ=QU7/SK1
      VAX=BVX**2
      VAY=BVY**2
      VAZ=BVZ**2
      VAH=VAX+VAY+VAZ
      CKV=GAMMA*QU(8)/QU(1)
      CSR=SQRT(CKV)
      CSV=CKV+VAH
      CWH=CKV-VAH
      DSCV=CSV**2-4.*CKV*VAX
      IF(DSCV.LT.0.) DSCV=0.
      DSCV=SQRT(DSCV)
      VF2=(CSV+DSCV)/2.
      VS2=(CSV-DSCV)/2.
      IF(VS2.LT.0.) VS2=0.
      CF=SQRT(VF2)
      CS=SQRT(VS2)
      VA=ABS(BVX)

      BXX=QU5**2
      BYZ=QU6**2+QU7**2
      SBB=SQRT(BYZ)
      BMG=BXX+BYZ
      GPP=GAMMA*QU(8)
      GPB=ABS(GPP-BXX)

      IF(BYZ.LE.DLT*BMG) THEN

         IF(GPB.GT.DLT*GPP) THEN
            CS=SQRT(ABS((CKV+VAX-SIGN(1.,CWH)*CWH)/2.))
            CF=SQRT(ABS((CKV+VAX+SIGN(1.,CWH)*CWH)/2.))
         ELSE
            CS=CSR
            CF=CSR
         ENDIF

      ENDIF

      U(1)=QU(2)-CF 
      U(2)=QU(2)-VA
      U(3)=QU(2)-CS
      U(4)=QU(2) 
      U(5)=QU(2)+CS
      U(6)=QU(2)+VA
      U(7)=QU(2)+CF 
      U(8)=QU(2)

      RETURN
      END

      SUBROUTINE VECTEGy_ig(gamma,QU,BS,U)
      Implicit NONE

      ENZO_REAL QU(8),U(8),BS(3)  	
      ENZO_REAL DLT,SK1,BVX,BVY,BVZ,VAX,VAY,VAZ
      ENZO_REAL VAH,CKV,A1,A2,CSR,CSV,CWH,DSCV,VF2
      ENZO_REAL VS2,CF,CS,VA,BYY,BXZ,SBB,BMG,GPP,GPB
      ENZO_REAL ASS,AS,AF,GAMMA ,QU5,QU6,QU7
	
      DLT=1.E-12

      QU5=QU(5)+BS(1)  
      QU6=QU(6)+BS(2)  
      QU7=QU(7)+BS(3)  
	  
      SK1=SQRT(QU(1))
      BVX=QU5/SK1
      BVY=QU6/SK1
      BVZ=QU7/SK1
      VAX=BVX**2
      VAY=BVY**2
      VAZ=BVZ**2
      VAH=VAX+VAY+VAZ
      CKV=GAMMA*QU(8)/QU(1)
      CSR=SQRT(CKV)
      CSV=CKV+VAH
      CWH=CKV-VAH
      DSCV=CSV**2-4.*CKV*VAY
      IF(DSCV.LT.0.) DSCV=0.
      DSCV=SQRT(DSCV)
      VF2=(CSV+DSCV)/2.
      VS2=(CSV-DSCV)/2.
      IF(VS2.LT.0.) VS2=0.
      CF=SQRT(VF2)
      CS=SQRT(VS2)
      VA=ABS(BVY)

      BYY=QU6**2
      BXZ=QU5**2+QU7**2
      SBB=SQRT(BXZ)
      BMG=BYY+BXZ
      GPP=GAMMA*QU(8)
      GPB=ABS(GPP-BYY)

      IF(BXZ.LE.DLT*BMG) THEN

         IF(GPB.GT.DLT*GPP) THEN
            CS=SQRT(ABS((CKV+VAY-SIGN(1.,CWH)*CWH)/2.))
            CF=SQRT(ABS((CKV+VAY+SIGN(1.,CWH)*CWH)/2.))
         ELSE
            CS=CSR
            CF=CSR
         ENDIF

      ENDIF

      U(1)=QU(3)-CF 
      U(2)=QU(3)-VA
      U(3)=QU(3)-CS
      U(4)=QU(3) 
      U(5)=QU(3)+CS
      U(6)=QU(3)+VA
      U(7)=QU(3)+CF 
      U(8)=QU(3)

      RETURN
      END

      SUBROUTINE VECTEGz_ig(gamma,QU,BS,U)
      Implicit NONE

      ENZO_REAL QU(8),U(8),BS(3)  	
      ENZO_REAL DLT,SK1,BVX,BVY,BVZ,VAX,VAY,VAZ
      ENZO_REAL VAH,CKV,A1,A2,CSR,CSV,CWH,DSCV,VF2
      ENZO_REAL VS2,CF,CS,VA,BXY,BZZ,SBB,BMG,GPP,GPB
      ENZO_REAL ASS,AS,AF,GAMMA,QU5,QU6,QU7 

      DLT=1.E-12

      QU5=QU(5)+BS(1)  
      QU6=QU(6)+BS(2)  
      QU7=QU(7)+BS(3)  
	  
      SK1=SQRT(QU(1))
      BVX=QU5/SK1
      BVY=QU6/SK1
      BVZ=QU7/SK1
      VAX=BVX**2
      VAY=BVY**2
      VAZ=BVZ**2
      VAH=VAX+VAY+VAZ
      CKV=GAMMA*QU(8)/QU(1)
      CSR=SQRT(CKV)
      CSV=CKV+VAH
      CWH=CKV-VAH
      DSCV=CSV**2-4.*CKV*VAZ
      IF(DSCV.LT.0.) DSCV=0.
      DSCV=SQRT(DSCV)
      VF2=(CSV+DSCV)/2.
      VS2=(CSV-DSCV)/2.
      IF(VS2.LT.0.) VS2=0.
      CF=SQRT(VF2)
      CS=SQRT(VS2)
      VA=ABS(BVZ)

      BZZ=QU7**2
      BXY=QU5**2+QU6**2
      SBB=SQRT(BXY)
      BMG=BXY+BZZ
      GPP=GAMMA*QU(8)
      GPB=ABS(GPP-BZZ)

      IF(BXY.LE.DLT*BMG) THEN

         IF(GPB.GT.DLT*GPP) THEN
            CS=SQRT(ABS((CKV+VAZ-SIGN(1.,CWH)*CWH)/2.))
            CF=SQRT(ABS((CKV+VAZ+SIGN(1.,CWH)*CWH)/2.))
         ELSE
            CS=CSR
            CF=CSR
         ENDIF

      ENDIF

      U(1)=QU(4)-CF 
      U(2)=QU(4)-VA
      U(3)=QU(4)-CS
      U(4)=QU(4) 
      U(5)=QU(4)+CS
      U(6)=QU(4)+VA
      U(7)=QU(4)+CF 
      U(8)=QU(4)

      RETURN
      END

      SUBROUTINE FLUXx_ig(gamma,QP,BS,F)
      Implicit NONE
	
      ENZO_REAL QP(8),Q(8),BS(3),F(8)
      ENZO_REAL gamma,qb5,qb6,qb7,bp,bp0,ux,uy,uz,vk,bm,pg,vb,pt

      CALL CONS_ig(gamma,QP,Q)
       
      qb5 = q(5)+bs(1)
      qb6 = q(6)+bs(2)
      qb7 = q(7)+bs(3)

      bp  = qb5**2 + qb6**2 + qb7**2
      bp0 = bs(1)**2 + bs(2)**2 + bs(3)**2
	 
      ux = q(2)/q(1)
      uy = q(3)/q(1)
      uz = q(4)/q(1)

      vk = q(2)**2 + q(3)**2 + q(4)**2
      bm = q(5)**2 + q(6)**2 + q(7)**2
      pg = (gamma-1.)*(q(8) - vk/q(1)/2. - bm/2.)

      vb = ux*q(5) + uy*q(6) + uz*q(7)
      pt = pg + bp/2.-bp0/2.


      f(1) = q(2)
      f(2) = q(2)*ux + pt - qb5**2 + bs(1)**2
      f(3) = q(3)*ux - qb5*qb6 + bs(1)*bs(2)
      f(4) = q(4)*ux - qb5*qb7 + bs(1)*bs(3)
      f(5) = 0.
      f(6) = qb6*ux - qb5*uy
      f(7) = qb7*ux - qb5*uz
      f(8) = (q(8) + pt)*ux - qb5*vb

      RETURN
      END
      
      SUBROUTINE FLUXy_ig(gamma,QP,BS,F)
      Implicit NONE
	
      ENZO_REAL QP(8),Q(8),BS(3),F(8)
      ENZO_REAL gamma,qb5,qb6,qb7,bp,bp0,ux,uy,uz,vk,bm,pg,vb,pt

      CALL CONS_ig(gamma,QP,Q)

      qb5 = q(5)+bs(1)
      qb6 = q(6)+bs(2)
      qb7 = q(7)+bs(3)

      bp  = qb5**2 + qb6**2 + qb7**2
      bp0 = bs(1)**2 + bs(2)**2 + bs(3)**2
	 
      ux = q(2)/q(1)
      uy = q(3)/q(1)
      uz = q(4)/q(1)

      vk = q(2)**2 + q(3)**2 + q(4)**2
      bm = q(5)**2 + q(6)**2 + q(7)**2
      pg = (gamma-1.)*(q(8) - vk/q(1)/2. - bm/2.)

      vb = ux*q(5) + uy*q(6) + uz*q(7)
      pt = pg + bp/2.-bp0/2.

      f(1) = q(3)
      f(2) = q(2)*uy - qb5*qb6 + bs(1)*bs(2)
      f(3) = q(3)*uy + pt - qb6**2 + bs(2)**2
      f(4) = q(4)*uy - qb6*qb7 + bs(2)*bs(3)
      f(5) = qb5*uy - qb6*ux
      f(6) = 0.
      f(7) = qb7*uy - qb6*uz
      f(8) = (q(8) + pt)*uy - qb6*vb

      RETURN
      END

      SUBROUTINE FLUXz_ig(gamma,QP,BS,F)
      Implicit NONE
      
      ENZO_REAL QP(8),Q(8),BS(3),F(8)
      ENZO_REAL gamma,qb5,qb6,qb7,bp,bp0,ux,uy,uz,vk,bm,pg,vb,pt

      CALL CONS_ig(gamma,QP,Q)

      qb5 = q(5)+bs(1)
      qb6 = q(6)+bs(2)
      qb7 = q(7)+bs(3)

      bp  = qb5**2 + qb6**2 + qb7**2
      bp0 = bs(1)**2 + bs(2)**2 + bs(3)**2
	 
      ux = q(2)/q(1)
      uy = q(3)/q(1)
      uz = q(4)/q(1)

      vk = q(2)**2 + q(3)**2 + q(4)**2
      bm = q(5)**2 + q(6)**2 + q(7)**2
      pg = (gamma-1.)*(q(8) - vk/q(1)/2. - bm/2.)

      vb = ux*q(5) + uy*q(6) + uz*q(7)
      pt = pg + bp/2.-bp0/2.

      f(1) = q(4)
      f(2) = q(2)*uz - qb5*qb7 + bs(1)*bs(3)
      f(3) = q(3)*uz - qb6*qb7 + bs(2)*bs(3)
      f(4) = q(4)*uz + pt - qb7**2 + bs(3)**2
      f(5) = qb5*uz - qb7*ux
      f(6) = qb6*uz - qb7*uy
      f(7) = 0.
      f(8) = (q(8) + pt)*uz - qb7*vb

      RETURN
      END

      SUBROUTINE FPMR_ig(Y,QL,QR,QP,QC) 
      Implicit NONE
      
      ENZO_REAL QL(8),QR(8),QP(8),QC(8)
      ENZO_REAL DQ,Q6,Y
      integer M
      
      DO M=1,8
         DQ=QR(M)-QL(M)
	 Q6=6.0*(QP(M)-(QR(M)+QL(M))/2.)
	 QC(M)=QL(M)+Y/2.*(DQ+(1.-2./3.*Y)*Q6)
      ENDDO

      RETURN
      END
	 
      SUBROUTINE FPML_ig(Y,QL,QR,QP,QC) 
      Implicit NONE

      integer m
      ENZO_REAL QL(8),QR(8),QP(8),QC(8)
      ENZO_REAL DQ,Q6,y
      
      DO M=1,8
         DQ=QR(M)-QL(M)
	 Q6=6.0*(QP(M)-(QR(M)+QL(M))/2.)
	 QC(M)=QR(M)-Y/2.*(DQ-(1.-2./3.*Y)*Q6)
      ENDDO

      RETURN
      END
	
      SUBROUTINE FQML_ig(Y,QL,QR,QP,QC) 
      Implicit NONE
      
      Integer M
      ENZO_REAL QL(8),QR(8),QP(8),QC(8)
      ENZO_REAL DQ,Q6,Y

      DO M=1,8
         DQ=QR(M)-QL(M)
	 Q6=6.0*(QP(M)-(QR(M)+QL(M))/2.)
	 QC(M)=QL(M)+Y*(DQ+(1.-Y)*Q6)
      ENDDO

      RETURN
      END


      Subroutine QDD6_ig(Nx,Qvr,Qvl,Qrx,Qlx,Qp)
      Implicit NONE

      Integer m,nx,i
      ENZO_REAL Qvr(Nx,8),Qvl(Nx,8),Qp(Nx,8) 
      ENZO_REAL Qrx(Nx,8),Qlx(Nx,8) 
      ENZO_REAL DQ,Q6

      Do I=1,Nx	
	
         DO M=1,8

            IF((QVR(I,M)-QP(I,M))*(QP(I,M)-QVL(I,M)).LE.0.)THEN
               QRX(I,M)=QP(I,M)
               QLX(I,M)=QP(I,M)
            ELSE

               DQ=QVR(I,M)-QVL(I,M)
               Q6=6.0*(QP(I,M)-(QVR(I,M)+QVL(I,M))/2.)

               QLX(I,M)=QVL(I,M)
               QRX(I,M)=QVR(I,M)
		 
               IF(DQ*Q6.GT.DQ**2) QLX(I,M)=3.*QP(I,M)-2.*QVR(I,M)
               IF(DQ*Q6.LT.-DQ**2) QRX(I,M)=3.*QP(I,M)-2.*QVL(I,M)
         
            ENDIF  
	 
         ENDDO 
      ENDDO

      Return
      End	
	 
      SUBROUTINE MATR_AX_ig(gamma,QU,BS,QL)
      Implicit NONE

      ENZO_REAL QU(8),QL(8,8),BS(3)
      ENZO_REAL GAMMA,QU5,QU6,QU7

      QU5=QU(5)+BS(1)  
      QU6=QU(6)+BS(2)  
      QU7=QU(7)+BS(3)  
	  
      QL(1,1)=QU(2)
      QL(1,2)=QU(1)
      QL(1,3)=0.
      QL(1,4)=0.
      QL(1,5)=0.
      QL(1,6)=0.
      QL(1,7)=0.
      QL(1,8)=0.
	
      QL(2,1)=0.
      QL(2,2)=QU(2)
      QL(2,3)=0.
      QL(2,4)=0.
      QL(2,5)=0.
      QL(2,6)=QU6/QU(1)
      QL(2,7)=QU7/QU(1)
      QL(2,8)=1./QU(1)

      QL(3,1)=0.
      QL(3,2)=0.
      QL(3,3)=QU(2)
      QL(3,4)=0.
      QL(3,5)=0.
      QL(3,6)=-QU5/QU(1)
      QL(3,7)=0.
      QL(3,8)=0.

      QL(4,1)=0.
      QL(4,2)=0.
      QL(4,3)=0.
      QL(4,4)=QU(2)
      QL(4,5)=0.
      QL(4,6)=0.
      QL(4,7)=-QU5/QU(1)
      QL(4,8)=0.

      QL(5,1)=0.
      QL(5,2)=0.
      QL(5,3)=0.
      QL(5,4)=0.
      QL(5,5)=QU(2)
      QL(5,6)=0.
      QL(5,7)=0.
      QL(5,8)=0.

      QL(6,1)=0.
      QL(6,2)=QU6
      QL(6,3)=-QU5
      QL(6,4)=0.
      QL(6,5)=0.
      QL(6,6)=QU(2)
      QL(6,7)=0.
      QL(6,8)=0.

      QL(7,1)=0.
      QL(7,2)=QU7
      QL(7,3)=0.
      QL(7,4)=-QU5
      QL(7,5)=0.
      QL(7,6)=0.
      QL(7,7)=QU(2)
      QL(7,8)=0.

      QL(8,1)=0.
      QL(8,2)=GAMMA*QU(8)
      QL(8,3)=0.
      QL(8,4)=0.
      QL(8,5)=0.
      QL(8,6)=0.
      QL(8,7)=0.
      QL(8,8)=QU(2)

      RETURN
      END

      SUBROUTINE MATR_AY_ig(gamma,QU,BS,QL)
      Implicit NONE

      ENZO_REAL QU(8),QL(8,8),BS(3)  	
      ENZO_REAL GAMMA,QU5,QU6,QU7  	

      QU5=QU(5)+BS(1)  
      QU6=QU(6)+BS(2)  
      QU7=QU(7)+BS(3)  
	  
      QL(1,1)=QU(3)
      QL(1,2)=0.
      QL(1,3)=QU(1)
      QL(1,4)=0.
      QL(1,5)=0.
      QL(1,6)=0.
      QL(1,7)=0.
      QL(1,8)=0.
	
      QL(2,1)=0.
      QL(2,2)=QU(3)
      QL(2,3)=0.
      QL(2,4)=0.
      QL(2,5)=-QU6/QU(1)
      QL(2,6)=0.
      QL(2,7)=0.
      QL(2,8)=0.

      QL(3,1)=0.
      QL(3,2)=0.
      QL(3,3)=QU(3)
      QL(3,4)=0.
      QL(3,5)=QU5/QU(1)
      QL(3,6)=0.
      QL(3,7)=QU7/QU(1)
      QL(3,8)=1./QU(1)

      QL(4,1)=0.
      QL(4,2)=0.
      QL(4,3)=0.
      QL(4,4)=QU(3)
      QL(4,5)=0.
      QL(4,6)=0.
      QL(4,7)=-QU6/QU(1)
      QL(4,8)=0.

      QL(5,1)=0.
      QL(5,2)=-QU6
      QL(5,3)=QU5
      QL(5,4)=0.
      QL(5,5)=QU(3)
      QL(5,6)=0.
      QL(5,7)=0.
      QL(5,8)=0.

      QL(6,1)=0.
      QL(6,2)=0.
      QL(6,3)=0.
      QL(6,4)=0.
      QL(6,5)=0.
      QL(6,6)=QU(3)
      QL(6,7)=0.
      QL(6,8)=0.

      QL(7,1)=0.
      QL(7,2)=0.
      QL(7,3)=QU7
      QL(7,4)=-QU6
      QL(7,5)=0.
      QL(7,6)=0.
      QL(7,7)=QU(3)
      QL(7,8)=0.

      QL(8,1)=0.
      QL(8,2)=0.
      QL(8,3)=GAMMA*QU(8)
      QL(8,4)=0.
      QL(8,5)=0.
      QL(8,6)=0.
      QL(8,7)=0.
      QL(8,8)=QU(3)

      RETURN
      END

      SUBROUTINE MATR_AZ_ig(gamma,QU,BS,QL)
      Implicit NONE

      ENZO_REAL QU(8),QL(8,8),BS(3)  	
      ENZO_REAL GAMMA,QU5,QU6,QU7  	

      QU5=QU(5)+BS(1)  
      QU6=QU(6)+BS(2)  
      QU7=QU(7)+BS(3)  
	  
      QL(1,1)=QU(4)
      QL(1,2)=0.
      QL(1,3)=0.
      QL(1,4)=QU(1)
      QL(1,5)=0.
      QL(1,6)=0.
      QL(1,7)=0.
      QL(1,8)=0.
	
      QL(2,1)=0.
      QL(2,2)=QU(4)
      QL(2,3)=0.
      QL(2,4)=0.
      QL(2,5)=-QU7/QU(1)
      QL(2,6)=0.
      QL(2,7)=0.
      QL(2,8)=0.

      QL(3,1)=0.
      QL(3,2)=0.
      QL(3,3)=QU(4)
      QL(3,4)=0.
      QL(3,5)=0.
      QL(3,6)=-QU7/QU(1)
      QL(3,7)=0.
      QL(3,8)=0.

      QL(4,1)=0.
      QL(4,2)=0.
      QL(4,3)=0.
      QL(4,4)=QU(4)
      QL(4,5)=QU5/QU(1)
      QL(4,6)=QU6/QU(1)
      QL(4,7)=0.
      QL(4,8)=1./QU(1)


      QL(5,1)=0.
      QL(5,2)=-QU7
      QL(5,3)=0.
      QL(5,4)=QU5
      QL(5,5)=QU(4)
      QL(5,6)=0.
      QL(5,7)=0.
      QL(5,8)=0.

      QL(6,1)=0.
      QL(6,2)=0.
      QL(6,3)=-QU7
      QL(6,4)=QU6
      QL(6,5)=0.
      QL(6,6)=QU(4)
      QL(6,7)=0.
      QL(6,8)=0.

      QL(7,1)=0.
      QL(7,2)=0.
      QL(7,3)=0.
      QL(7,4)=0.
      QL(7,5)=0.
      QL(7,6)=0.
      QL(7,7)=QU(4)
      QL(7,8)=0.

      QL(8,1)=0.
      QL(8,2)=0.
      QL(8,3)=0.
      QL(8,4)=GAMMA*QU(8)
      QL(8,5)=0.
      QL(8,6)=0.
      QL(8,7)=0.
      QL(8,8)=QU(4)

      RETURN
      END

      SUBROUTINE DUDW_ig(gamma,QU,BS,QL)
      Implicit NONE

      ENZO_REAL QU(8),QL(8,8),BS(3)  	
      ENZO_REAL GAMMA,QU5,QU6,QU7  	

      QU5=QU(5)+BS(1)  
      QU6=QU(6)+BS(2)  
      QU7=QU(7)+BS(3)  
	  
      QL(1,1)=1.
      QL(1,2)=0.
      QL(1,3)=0.
      QL(1,4)=0.
      QL(1,5)=0.
      QL(1,6)=0.
      QL(1,7)=0.
      QL(1,8)=0.
	
      QL(2,1)=QU(2)
      QL(2,2)=QU(1)
      QL(2,3)=0.
      QL(2,4)=0.
      QL(2,5)=0.
      QL(2,6)=0.
      QL(2,7)=0.
      QL(2,8)=0.

      QL(3,1)=QU(3)
      QL(3,2)=0.
      QL(3,3)=QU(1)
      QL(3,4)=0.
      QL(3,5)=0.
      QL(3,6)=0.
      QL(3,7)=0.
      QL(3,8)=0.

      QL(4,1)=QU(4)
      QL(4,2)=0.
      QL(4,3)=0.
      QL(4,4)=QU(1)
      QL(4,5)=0.
      QL(4,6)=0.
      QL(4,7)=0.
      QL(4,8)=0.

      QL(5,1)=0.
      QL(5,2)=0.
      QL(5,3)=0.
      QL(5,4)=0.
      QL(5,5)=1.
      QL(5,6)=0.
      QL(5,7)=0.
      QL(5,8)=0.

      QL(6,1)=0.
      QL(6,2)=0.
      QL(6,3)=0.
      QL(6,4)=0.
      QL(6,5)=0.
      QL(6,6)=1.
      QL(6,7)=0.
      QL(6,8)=0.

      QL(7,1)=0.
      QL(7,2)=0.
      QL(7,3)=0.
      QL(7,4)=0.
      QL(7,5)=0.
      QL(7,6)=0.
      QL(7,7)=1.
      QL(7,8)=0.

      QL(8,1)=(QU(2)**2+QU(3)**2+QU(4)**2)/2.
      QL(8,2)=QU(1)*QU(2)
      QL(8,3)=QU(1)*QU(3)
      QL(8,4)=QU(1)*QU(4)
      QL(8,5)=QU5
      QL(8,6)=QU6
      QL(8,7)=QU7
      QL(8,8)=1./(GAMMA-1.)

      RETURN
      END

      SUBROUTINE AMPLTD_ig(QL,QD,SL)
      Implicit NONE

      Integer L,M
      ENZO_REAL QL(8,8),QD(8),SL(8)
      
      DO M=1,8
         SL(M)=0.   
         DO L=1,8
            SL(M)=SL(M)+QL(M,L)*QD(L)
         ENDDO
      ENDDO

      RETURN
      END 

      Function Xminmod_ig(x,y)
      
      ENZO_REAL x,y,Xminmod
	
      Xminmod = 0.5*(sign(1.,x) + sign(1.,y))*min(abs(x),abs(y))

      Return
      End

      Function Xmedian_ig(x,y,z)
      
      ENZO_REAL x,y,z,Xminmod,Xmedian

      External Xminmod

      Xmedian = x + Xminmod(y-x,z-x)

      Return
      End

      SUBROUTINE W5RECM_ig(VP1,VP2,VP3,VP4,VP5,VL,VR)
      Implicit NONE
      
      ENZO_REAL IS1,IS2,IS3,DELTA,WS1,WS2,WS3,WJ1,WJ2,WJ3
      ENZO_REAL AJ1,AJ2,AJ3,AJS,WM1,WM2,WM3,VL,VR
      ENZO_REAL VP1,VP2,VP3,VP4,VP5
	        
      DELTA = 1.0e-10
	  
      WS1 = 0.1
      WS2 = 0.6
      WS3 = 0.3 


      WJ1 = (11.*VP3 - 7.*VP2 + 2.*VP1)/6.
      WJ2 = (2. *VP4 + 5.*VP3 -      VP2)/6.
      WJ3 = (     -VP5 + 5.*VP4 + 2.*VP3)/6.

      IS1 = 13./12.*(VP3 - 2.*VP2 + VP1)**2 
     &	  + (3.*VP3 - 4.*VP2 + VP1)**2/4.
      IS2 = 13./12.*(VP4 - 2.*VP3 + VP2)**2 
     &	  + (VP4 - VP2)**2/4.
      IS3 = 13./12.*(VP5 - 2.*VP4 + VP3)**2 
     &	  + (VP5 - 4.*VP4 + 3.*VP3)**2/4.

      AJ1 = WS1/(IS1+DELTA)**2   
      AJ2 = WS2/(IS2+DELTA)**2   
      AJ3 = WS3/(IS3+DELTA)**2   

      AJS = AJ1 + AJ2 +AJ3 

      WM1 = AJ1/AJS
      WM2 = AJ2/AJS
      WM3 = AJ3/AJS

      AJ1 = WM1*(WS1 + WS1*WS1 - 3.*WS1*WM1 + WM1*WM1)/(WS1*WS1 
     &    + WM1*(1. - 2.*WS1)) 
      AJ2 = WM2*(WS2 + WS2*WS2 - 3.*WS2*WM2 + WM2*WM2)/(WS2*WS2 
     &    + WM2*(1. - 2.*WS2)) 
      AJ3 = WM3*(WS3 + WS3*WS3 - 3.*WS3*WM3 + WM3*WM3)/(WS3*WS3 
     &    + WM3*(1. - 2.*WS3)) 

      AJS = AJ1 + AJ2 +AJ3 

      WM1 = AJ1/AJS
      WM2 = AJ2/AJS
      WM3 = AJ3/AJS

      VR = WM1*WJ1 + WM2*WJ2 + WM3*WJ3 

      WJ1 = (11.*VP3 - 7.*VP4 + 2.*VP5)/6.
      WJ2 = (2. *VP2 + 5.*VP3 -      VP4)/6.
      WJ3 = (     -VP1 + 5.*VP2 + 2.*VP3)/6.

      IS1 = 13./12.*(VP3 - 2.*VP4 + VP5)**2 
     &	  + (3.*VP3 - 4.*VP4 + VP5)**2/4.
	  IS2 = 13./12.*(VP2 - 2.*VP3 + VP4)**2 
     &	  + (VP2 - VP4)**2/4.
	  IS3 = 13./12.*(VP1 - 2.*VP2 + VP3)**2 
     &	  + (VP1 - 4.*VP2 + 3.*VP3)**2/4.

      AJ1 = WS1/(IS1+DELTA)**2   
      AJ2 = WS2/(IS2+DELTA)**2   
      AJ3 = WS3/(IS3+DELTA)**2   

      AJS = AJ1 + AJ2 +AJ3 

      WM1 = AJ1/AJS
      WM2 = AJ2/AJS
      WM3 = AJ3/AJS

      AJ1 = WM1*(WS1 + WS1*WS1 - 3.*WS1*WM1 + WM1*WM1)/(WS1*WS1 
     &    + WM1*(1. - 2.*WS1)) 
      AJ2 = WM2*(WS2 + WS2*WS2 - 3.*WS2*WM2 + WM2*WM2)/(WS2*WS2 
     &    + WM2*(1. - 2.*WS2)) 
      AJ3 = WM3*(WS3 + WS3*WS3 - 3.*WS3*WM3 + WM3*WM3)/(WS3*WS3 
     &    + WM3*(1. - 2.*WS3)) 

      AJS = AJ1 + AJ2 +AJ3 

      WM1 = AJ1/AJS
      WM2 = AJ2/AJS
      WM3 = AJ3/AJS

      VL = WM1*WJ1 + WM2*WJ2 + WM3*WJ3 

      RETURN
      END
	
      Function VLeer_ig(a,b)

      ENZO_REAL a,b,VLeer

      if(a*a+b*b.ne.0)then
         VLeer=max(a*b,0.)*(a+b)/(a*a+b*b)

      else
         VLeer=0.
      endif 

      Return
      End	
	  
      SUBROUTINE PRIM_ig(gamma,QU,QP)
      Implicit NONE

      ENZO_REAL QU(8),QP(8)
      ENZO_REAL gamma,VKV,BKV

      QP(1)=QU(1)
      QP(2)=QU(2)/QU(1)
      QP(3)=QU(3)/QU(1)
      QP(4)=QU(4)/QU(1)
      QP(5)=QU(5)
      QP(6)=QU(6)
      QP(7)=QU(7)

      VKV=QU(2)**2+QU(3)**2+QU(4)**2
      BKV=QU(5)**2+QU(6)**2+QU(7)**2

      QP(8)=(GAMMA-1.)*(QU(8)-VKV/QU(1)/2.-BKV/2.)

      RETURN
      END

      SUBROUTINE CONS_ig(gamma,QP,QU)
      Implicit NONE

      ENZO_REAL QU(8),QP(8)
      ENZO_REAL gamma,VKV,BKV

      QU(1)=QP(1)
      QU(2)=QP(2)*QP(1)
      QU(3)=QP(3)*QP(1)
      QU(4)=QP(4)*QP(1)
      QU(5)=QP(5)
      QU(6)=QP(6)
      QU(7)=QP(7)

      VKV=QP(2)**2+QP(3)**2+QP(4)**2
      BKV=QP(5)**2+QP(6)**2+QP(7)**2
    
      QU(8)=QP(8)/(GAMMA-1.)+QP(1)*VKV/2.+BKV/2.

      RETURN
      END
