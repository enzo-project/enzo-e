c     $Id$
c     See LICENSE_CELLO file for license and copyright information

      SUBROUTINE VECTLRx(QU,QL,QR,U)
      Implicit NONE

      REAL QU(7),QL(7,7),QR(7,7),U(7)          
      REAL DLT,X22,SK1,BVX,BVY,BVZ,VAX,VAY,VAZ
      REAL VAH,CKV,A1,A2,CSR,CSV,CWH,DSCV,VF2
      REAL VS2,CF,CS,CA,BXX,BYZ,SBB,BMG,GPP,GPB
      REAL ASS,AS,AF,BY,BZ,BSGN,CD2

      a1=1.
      a2=1.
      
      DLT=1.E-12
      X22=1.e0/sqrt(2.e0)
      
      SK1=sqrt(QU(1))
      BVX=QU(5)/SK1
      BVY=QU(6)/SK1
      BVZ=QU(7)/SK1
      VAX=BVX**2
      VAY=BVY**2
      VAZ=BVZ**2
      VAH=VAX+VAY+VAZ
      CKV=a2
      CSR=sqrt(CKV)
      CSV=CKV+VAH
      CWH=CKV-VAH
      DSCV=CSV**2-4.e0*CKV*VAX
      IF(DSCV.LT.0.) DSCV=0.
      DSCV=sqrt(DSCV)
      VF2=(CSV+DSCV)/2.e0
      VS2=(CSV-DSCV)/2.e0
      IF(VS2.LT.0.) VS2=0.
      CF=sqrt(VF2)
      CS=sqrt(VS2)
      CA=abs(BVX)

      BXX=QU(5)**2
      BYZ=QU(6)**2+QU(7)**2
      SBB=sqrt(BYZ)
      BMG=BXX+BYZ
      GPP=a2*QU(1)
      GPB=abs(GPP-BXX)

      IF(BYZ.GT.DLT*BMG) THEN

         ASS=sqrt(CWH**2+4.e0*CKV*(VAY+VAZ))
         AS=sqrt(abs((1.e0-CWH/ASS)/2.e0))
         AF=sqrt(abs((1.e0+CWH/ASS)/2.e0))

         BY=QU(6)/SBB
         BZ=QU(7)/SBB


      ELSE

         IF(GPB.GT.DLT*GPP) THEN
            AS=sqrt(abs((1.e0-SIGN(1.e0,CWH))/2.e0))
            AF=sqrt(abs((1.e0+SIGN(1.e0,CWH))/2.e0))
            CS=sqrt(abs((CKV+VAX-SIGN(1.e0,CWH)*CWH)/2.e0))
            CF=sqrt(abs((CKV+VAX+SIGN(1.e0,CWH)*CWH)/2.e0))
         ELSE
            AS=X22
            AF=X22
            CS=CSR
            CF=CSR  
         ENDIF

         BY=X22
         BZ=X22

      ENDIF

      BSGN=SIGN(1.e0,QU(5))
      CD2=2.e0*CKV
      
      QL(1,1)=AF/2.e0
      QL(1,2)=-QU(1)*AF*CF/CD2
      QL(1,3)=QU(1)*AS*CS*BY*BSGN/CD2
      QL(1,4)=QU(1)*AS*CS*BZ*BSGN/CD2
      QL(1,5)=0.
      QL(1,6)=SK1*CSR*AS*BY/CD2
      QL(1,7)=SK1*CSR*AS*BZ/CD2
      
      QL(2,1)=0.
      QL(2,2)=0.
      QL(2,3)=-QU(1)*BZ*BSGN/2.e0
      QL(2,4)=QU(1)*BY*BSGN/2.e0
      QL(2,5)=0.
      QL(2,6)=-SK1*BZ/2.e0
      QL(2,7)=SK1*BY/2.e0

      QL(3,1)=AS/2.e0
      QL(3,2)=-QU(1)*AS*CS/CD2
      QL(3,3)=-QU(1)*AF*CF*BY*BSGN/CD2
      QL(3,4)=-QU(1)*AF*CF*BZ*BSGN/CD2
      QL(3,5)=0.
      QL(3,6)=-SK1*CSR*AF*BY/CD2
      QL(3,7)=-SK1*CSR*AF*BZ/CD2

      QL(4,1)=0.
      QL(4,2)=0.
      QL(4,3)=0.
      QL(4,4)=0.
      QL(4,5)=1.e0
      QL(4,6)=0.
      QL(4,7)=0.

      QL(5,1)=AS/2.e0
      QL(5,2)=QU(1)*AS*CS/CD2
      QL(5,3)=QU(1)*AF*CF*BY*BSGN/CD2
      QL(5,4)=QU(1)*AF*CF*BZ*BSGN/CD2
      QL(5,5)=0.
      QL(5,6)=-SK1*CSR*AF*BY/CD2
      QL(5,7)=-SK1*CSR*AF*BZ/CD2

      QL(6,1)=0.
      QL(6,2)=0.
      QL(6,3)=-QU(1)*BZ*BSGN/2.e0
      QL(6,4)=QU(1)*BY*BSGN/2.e0
      QL(6,5)=0.
      QL(6,6)=SK1*BZ/2.e0
      QL(6,7)=-SK1*BY/2.e0

      QL(7,1)=AF/2.e0
      QL(7,2)=QU(1)*AF*CF/CD2
      QL(7,3)=-QU(1)*AS*CS*BY*BSGN/CD2
      QL(7,4)=-QU(1)*AS*CS*BZ*BSGN/CD2
      QL(7,5)=0.
      QL(7,6)=SK1*CSR*AS*BY/CD2
      QL(7,7)=SK1*CSR*AS*BZ/CD2
      

      QR(1,1)=AF
      QR(1,2)=0.
      QR(1,3)=AS
      QR(1,4)=0.
      QR(1,5)=AS
      QR(1,6)=0.
      QR(1,7)=AF

      QR(2,1)=-AF*CF/QU(1)
      QR(2,2)=0.
      QR(2,3)=-AS*CS/QU(1)
      QR(2,4)=0.
      QR(2,5)=AS*CS/QU(1)
      QR(2,6)=0.
      QR(2,7)=AF*CF/QU(1)

      QR(3,1)=AS*CS*BY*BSGN/QU(1)
      QR(3,2)=-BZ*BSGN/QU(1)
      QR(3,3)=-AF*CF*BY*BSGN/QU(1)
      QR(3,4)=0.
      QR(3,5)=AF*CF*BY*BSGN/QU(1)
      QR(3,6)=-BZ*BSGN/QU(1)
      QR(3,7)=-AS*CS*BY*BSGN/QU(1)

      QR(4,1)=AS*CS*BZ*BSGN/QU(1)
      QR(4,2)=BY*BSGN/QU(1)
      QR(4,3)=-AF*CF*BZ*BSGN/QU(1)
      QR(4,4)=0.
      QR(4,5)=AF*CF*BZ*BSGN/QU(1)
      QR(4,6)=BY*BSGN/QU(1)
      QR(4,7)=-AS*CS*BZ*BSGN/QU(1)

      QR(5,1)=0.
      QR(5,2)=0.
      QR(5,3)=0.
      QR(5,4)=1.e0
      QR(5,5)=0.
      QR(5,6)=0.
      QR(5,7)=0.

      QR(6,1)=AS*BY*CSR/SK1
      QR(6,2)=-BZ/SK1
      QR(6,3)=-AF*BY*CSR/SK1
      QR(6,4)=0.
      QR(6,5)=-AF*BY*CSR/SK1
      QR(6,6)=BZ/SK1
      QR(6,7)=AS*BY*CSR/SK1

      QR(7,1)=AS*BZ*CSR/SK1
      QR(7,2)=BY/SK1
      QR(7,3)=-AF*BZ*CSR/SK1
      QR(7,4)=0.
      QR(7,5)=-AF*BZ*CSR/SK1
      QR(7,6)=-BY/SK1
      QR(7,7)=AS*BZ*CSR/SK1


      U(1)=QU(2)-CF 
      U(2)=QU(2)-CA
      U(3)=QU(2)-CS
      U(4)=QU(2)
      U(5)=QU(2)+CS
      U(6)=QU(2)+CA
      U(7)=QU(2)+CF 

      RETURN
      END

      SUBROUTINE VECTLRy(QU,QL,QR,U)
      Implicit NONE

      REAL QU(7),QL(7,7),QR(7,7),U(7)          
      REAL DLT,X22,SK1,BVX,BVY,BVZ,VAX,VAY,VAZ
      REAL VAH,CKV,A1,A2,CSR,CSV,CWH,DSCV,VF2
      REAL VS2,CF,CS,CA,BYY,BXZ,SBB,BMG,GPP,GPB
      REAL ASS,AS,AF,BX,BZ,BSGN,CD2
      
      a1=1.
      a2=1.

      DLT=1.E-12
      X22=1.e0/sqrt(2.e0)
      
      SK1=sqrt(QU(1))
      BVX=QU(5)/SK1
      BVY=QU(6)/SK1
      BVZ=QU(7)/SK1
      VAX=BVX**2
      VAY=BVY**2
      VAZ=BVZ**2
      VAH=VAX+VAY+VAZ
      CKV=A2
      CSR=sqrt(CKV)
      CSV=CKV+VAH
      CWH=CKV-VAH
      DSCV=CSV**2-4.e0*CKV*VAY
      IF(DSCV.LT.0.) DSCV=0.
      DSCV=sqrt(DSCV)
      VF2=(CSV+DSCV)/2.e0
      VS2=(CSV-DSCV)/2.e0
      IF(VS2.LT.0.) VS2=0.
      CF=sqrt(VF2)
      CS=sqrt(VS2)
      CA=abs(BVY)

      BYY=QU(6)**2
      BXZ=QU(5)**2+QU(7)**2
      SBB=sqrt(BXZ)
      BMG=BYY+BXZ
      GPP=A2*QU(1)
      GPB=abs(GPP-BYY)

      IF(BXZ.GT.DLT*BMG) THEN

         ASS=sqrt(CWH**2+4.e0*CKV*(VAX+VAZ))
         AS=sqrt(abs((1.e0-CWH/ASS)/2.e0))
         AF=sqrt(abs((1.e0+CWH/ASS)/2.e0))

         BX=QU(5)/SBB
         BZ=QU(7)/SBB

      ELSE

         IF(GPB.GT.DLT*GPP) THEN
            AS=sqrt(abs((1.e0-SIGN(1.e0,CWH))/2.e0))
            AF=sqrt(abs((1.e0+SIGN(1.e0,CWH))/2.e0))
            CS=sqrt(abs((CKV+VAY-SIGN(1.e0,CWH)*CWH)/2.e0))
            CF=sqrt(abs((CKV+VAY+SIGN(1.e0,CWH)*CWH)/2.e0))
         ELSE
            AS=X22
            AF=X22
            CS=CSR
            CF=CSR  
         ENDIF

         BX=X22
         BZ=X22

      ENDIF

      BSGN=SIGN(1.e0,QU(6))
      CD2=2.e0*CKV
      
      QL(1,1)=AF/2.e0
      QL(1,2)=QU(1)*AS*CS*BX*BSGN/CD2
      QL(1,3)=-QU(1)*AF*CF/CD2
      QL(1,4)=QU(1)*AS*CS*BZ*BSGN/CD2
      QL(1,5)=SK1*CSR*AS*BX/CD2
      QL(1,6)=0.
      QL(1,7)=SK1*CSR*AS*BZ/CD2
      
      QL(2,1)=0.
      QL(2,2)=-QU(1)*BZ*BSGN/2.e0
      QL(2,3)=0.
      QL(2,4)=QU(1)*BX*BSGN/2.e0
      QL(2,5)=-SK1*BZ/2.e0
      QL(2,6)=0.
      QL(2,7)=SK1*BX/2.e0

      QL(3,1)=AS/2.e0
      QL(3,2)=-QU(1)*AF*CF*BX*BSGN/CD2
      QL(3,3)=-QU(1)*AS*CS/CD2
      QL(3,4)=-QU(1)*AF*CF*BZ*BSGN/CD2
      QL(3,5)=-SK1*CSR*AF*BX/CD2
      QL(3,6)=0.
      QL(3,7)=-SK1*CSR*AF*BZ/CD2

      QL(4,1)=0.
      QL(4,2)=0.
      QL(4,3)=0.
      QL(4,4)=0.
      QL(4,5)=0.
      QL(4,6)=1.e0
      QL(4,7)=0.

      QL(5,1)=AS/2.e0
      QL(5,2)=QU(1)*AF*CF*BX*BSGN/CD2
      QL(5,3)=QU(1)*AS*CS/CD2
      QL(5,4)=QU(1)*AF*CF*BZ*BSGN/CD2
      QL(5,5)=-SK1*CSR*AF*BX/CD2
      QL(5,6)=0.
      QL(5,7)=-SK1*CSR*AF*BZ/CD2

      QL(6,1)=0.
      QL(6,2)=-QU(1)*BZ*BSGN/2.e0
      QL(6,3)=0.
      QL(6,4)=QU(1)*BX*BSGN/2.e0
      QL(6,5)=SK1*BZ/2.e0
      QL(6,6)=0.
      QL(6,7)=-SK1*BX/2.e0

      QL(7,1)=AF/2.e0
      QL(7,2)=-QU(1)*AS*CS*BX*BSGN/CD2
      QL(7,3)=QU(1)*AF*CF/CD2
      QL(7,4)=-QU(1)*AS*CS*BZ*BSGN/CD2
      QL(7,5)=SK1*CSR*AS*BX/CD2
      QL(7,6)=0.
      QL(7,7)=SK1*CSR*AS*BZ/CD2
      
      QR(1,1)=AF
      QR(1,2)=0.
      QR(1,3)=AS
      QR(1,4)=0.
      QR(1,5)=AS
      QR(1,6)=0.
      QR(1,7)=AF

      QR(2,1)=AS*CS*BX*BSGN/QU(1)
      QR(2,2)=-BZ*BSGN/QU(1)
      QR(2,3)=-AF*CF*BX*BSGN/QU(1)
      QR(2,4)=0.
      QR(2,5)=AF*CF*BX*BSGN/QU(1)
      QR(2,6)=-BZ*BSGN/QU(1)
      QR(2,7)=-AS*CS*BX*BSGN/QU(1)

      QR(3,1)=-AF*CF/QU(1)
      QR(3,2)=0.
      QR(3,3)=-AS*CS/QU(1)
      QR(3,4)=0.
      QR(3,5)=AS*CS/QU(1)
      QR(3,6)=0.
      QR(3,7)=AF*CF/QU(1)

      QR(4,1)=AS*CS*BZ*BSGN/QU(1)
      QR(4,2)=BX*BSGN/QU(1)
      QR(4,3)=-AF*CF*BZ*BSGN/QU(1)
      QR(4,4)=0.
      QR(4,5)=AF*CF*BZ*BSGN/QU(1)
      QR(4,6)=BX*BSGN/QU(1)
      QR(4,7)=-AS*CS*BZ*BSGN/QU(1)

      QR(5,1)=AS*BX*CSR/SK1
      QR(5,2)=-BZ/SK1
      QR(5,3)=-AF*BX*CSR/SK1
      QR(5,4)=0.
      QR(5,5)=-AF*BX*CSR/SK1
      QR(5,6)=BZ/SK1
      QR(5,7)=AS*BX*CSR/SK1

      QR(6,1)=0.
      QR(6,2)=0.
      QR(6,3)=0.
      QR(6,4)=1.e0
      QR(6,5)=0.
      QR(6,6)=0.
      QR(6,7)=0.

      QR(7,1)=AS*BZ*CSR/SK1
      QR(7,2)=BX/SK1
      QR(7,3)=-AF*BZ*CSR/SK1
      QR(7,4)=0.
      QR(7,5)=-AF*BZ*CSR/SK1
      QR(7,6)=-BX/SK1
      QR(7,7)=AS*BZ*CSR/SK1


      U(1)=QU(3)-CF 
      U(2)=QU(3)-CA
      U(3)=QU(3)-CS
      U(4)=QU(3)
      U(5)=QU(3)+CS
      U(6)=QU(3)+CA
      U(7)=QU(3)+CF 

      RETURN
      END

      SUBROUTINE VECTLRz(QU,QL,QR,U)
      Implicit NONE

      REAL QU(7),QL(7,7),QR(7,7),U(7)          
      REAL DLT,X22,SK1,BVX,BVY,BVZ,VAX,VAY,VAZ
      REAL VAH,CKV,A1,A2,CSR,CSV,CWH,DSCV,VF2
      REAL VS2,CF,CS,CA,BYY,BXZ,SBB,BMG,GPP,GPB
      REAL ASS,AS,AF,BX,BY,BSGN,CD2
      
      a1=1.
      a2=1.
      
      DLT=1.E-12
      X22=1.e0/sqrt(2.e0)
      
      SK1=sqrt(QU(1))
      BVX=QU(5)/SK1
      BVY=QU(6)/SK1
      BVZ=QU(7)/SK1
      VAX=BVX**2
      VAY=BVY**2
      VAZ=BVZ**2
      VAH=VAX+VAY+VAZ
      CKV=A2
      CSR=sqrt(CKV)
      CSV=CKV+VAH
      CWH=CKV-VAH
      DSCV=CSV**2-4.e0*CKV*VAZ
      IF(DSCV.LT.0.) DSCV=0.
      DSCV=sqrt(DSCV)
      VF2=(CSV+DSCV)/2.e0
      VS2=(CSV-DSCV)/2.e0
      IF(VS2.LT.0.) VS2=0.
      CF=sqrt(VF2)
      CS=sqrt(VS2)
      CA=abs(BVZ)

      BYY=QU(7)**2
      BXZ=QU(5)**2+QU(6)**2
      SBB=sqrt(BXZ)
      BMG=BYY+BXZ
      GPP=A2*QU(1)
      GPB=abs(GPP-BYY)

      IF(BXZ.GT.DLT*BMG) THEN

         ASS=sqrt(CWH**2+4.e0*CKV*(VAX+VAY))
         AS=sqrt(abs((1.e0-CWH/ASS)/2.e0))
         AF=sqrt(abs((1.e0+CWH/ASS)/2.e0))

         BX=QU(5)/SBB
         BY=QU(6)/SBB

      ELSE

         IF(GPB.GT.DLT*GPP) THEN
            AS=sqrt(abs((1.e0-SIGN(1.e0,CWH))/2.e0))
            AF=sqrt(abs((1.e0+SIGN(1.e0,CWH))/2.e0))
            CS=sqrt(abs((CKV+VAZ-SIGN(1.e0,CWH)*CWH)/2.e0))
            CF=sqrt(abs((CKV+VAZ+SIGN(1.e0,CWH)*CWH)/2.e0))
         ELSE
            AS=X22
            AF=X22
            CS=CSR
            CF=CSR  
         ENDIF

         BX=X22
         BY=X22

      ENDIF

      BSGN=SIGN(1.e0,QU(7))
      CD2=2.e0*CKV
      
      QL(1,1)=AF/2.e0
      QL(1,2)=QU(1)*AS*CS*BX*BSGN/CD2
      QL(1,3)=QU(1)*AS*CS*BY*BSGN/CD2
      QL(1,4)=-QU(1)*AF*CF/CD2
      QL(1,5)=SK1*CSR*AS*BX/CD2
      QL(1,6)=SK1*CSR*AS*BY/CD2
      QL(1,7)=0.
      
      QL(2,1)=0.
      QL(2,2)=-QU(1)*BY*BSGN/2.e0
      QL(2,3)=QU(1)*BX*BSGN/2.e0
      QL(2,4)=0.
      QL(2,5)=-SK1*BY/2.e0
      QL(2,6)=SK1*BX/2.e0
      QL(2,7)=0.

      QL(3,1)=AS/2.e0
      QL(3,2)=-QU(1)*AF*CF*BX*BSGN/CD2
      QL(3,3)=-QU(1)*AF*CF*BY*BSGN/CD2
      QL(3,4)=-QU(1)*AS*CS/CD2
      QL(3,5)=-SK1*CSR*AF*BX/CD2
      QL(3,6)=-SK1*CSR*AF*BY/CD2
      QL(3,7)=0.

      QL(4,1)=0.
      QL(4,2)=0.
      QL(4,3)=0.
      QL(4,4)=0.
      QL(4,5)=0.
      QL(4,6)=0.
      QL(4,7)=1.e0

      QL(5,1)=AS/2.e0
      QL(5,2)=QU(1)*AF*CF*BX*BSGN/CD2
      QL(5,3)=QU(1)*AF*CF*BY*BSGN/CD2
      QL(5,4)=QU(1)*AS*CS/CD2
      QL(5,5)=-SK1*CSR*AF*BX/CD2
      QL(5,6)=-SK1*CSR*AF*BY/CD2
      QL(5,7)=0.

      QL(6,1)=0.
      QL(6,2)=-QU(1)*BY*BSGN/2.e0
      QL(6,3)=QU(1)*BX*BSGN/2.e0
      QL(6,4)=0.
      QL(6,5)=SK1*BY/2.e0
      QL(6,6)=-SK1*BX/2.e0
      QL(6,7)=0.

      QL(7,1)=AF/2.e0
      QL(7,2)=-QU(1)*AS*CS*BX*BSGN/CD2
      QL(7,3)=-QU(1)*AS*CS*BY*BSGN/CD2
      QL(7,4)=QU(1)*AF*CF/CD2
      QL(7,5)=SK1*CSR*AS*BX/CD2
      QL(7,6)=SK1*CSR*AS*BY/CD2
      QL(7,7)=0.
      

      QR(1,1)=AF
      QR(1,2)=0.
      QR(1,3)=AS
      QR(1,4)=0.
      QR(1,5)=AS
      QR(1,6)=0.
      QR(1,7)=AF

      QR(2,1)=AS*CS*BX*BSGN/QU(1)
      QR(2,2)=-BY*BSGN/QU(1)
      QR(2,3)=-AF*CF*BX*BSGN/QU(1)
      QR(2,4)=0.
      QR(2,5)=AF*CF*BX*BSGN/QU(1)
      QR(2,6)=-BY*BSGN/QU(1)
      QR(2,7)=-AS*CS*BX*BSGN/QU(1)

      QR(3,1)=AS*CS*BY*BSGN/QU(1)
      QR(3,2)=BX*BSGN/QU(1)
      QR(3,3)=-AF*CF*BY*BSGN/QU(1)
      QR(3,4)=0.
      QR(3,5)=AF*CF*BY*BSGN/QU(1)
      QR(3,6)=BX*BSGN/QU(1)
      QR(3,7)=-AS*CS*BY*BSGN/QU(1)

      QR(4,1)=-AF*CF/QU(1)
      QR(4,2)=0.
      QR(4,3)=-AS*CS/QU(1)
      QR(4,4)=0.
      QR(4,5)=AS*CS/QU(1)
      QR(4,6)=0.
      QR(4,7)=AF*CF/QU(1)

      QR(5,1)=AS*BX*CSR/SK1
      QR(5,2)=-BY/SK1
      QR(5,3)=-AF*BX*CSR/SK1
      QR(5,4)=0.
      QR(5,5)=-AF*BX*CSR/SK1
      QR(5,6)=BY/SK1
      QR(5,7)=AS*BX*CSR/SK1

      QR(6,1)=AS*BY*CSR/SK1
      QR(6,2)=BX/SK1
      QR(6,3)=-AF*BY*CSR/SK1
      QR(6,4)=0.
      QR(6,5)=-AF*BY*CSR/SK1
      QR(6,6)=-BX/SK1
      QR(6,7)=AS*BY*CSR/SK1

      QR(7,1)=0.
      QR(7,2)=0.
      QR(7,3)=0.
      QR(7,4)=1.e0
      QR(7,5)=0.
      QR(7,6)=0.
      QR(7,7)=0.


      U(1)=QU(4)-CF 
      U(2)=QU(4)-CA
      U(3)=QU(4)-CS
      U(4)=QU(4)
      U(5)=QU(4)+CS
      U(6)=QU(4)+CA
      U(7)=QU(4)+CF 

      RETURN
      END


      SUBROUTINE VECTEGx(QU,U)
      Implicit NONE

      REAL QU(7),U(7)          
      REAL DLT,SK1,BVX,BVY,BVZ,VAX,VAY,VAZ
      REAL VAH,CKV,A1,A2,CSR,CSV,CWH,DSCV,VF2
      REAL VS2,CF,CS,VA,BXX,BYZ,SBB,BMG,GPP,GPB

      a1=1.
      a2=1.          

      DLT=1.E-12

      SK1=sqrt(QU(1))
      BVX=QU(5)/SK1
      BVY=QU(6)/SK1
      BVZ=QU(7)/SK1
      VAX=BVX**2
      VAY=BVY**2
      VAZ=BVZ**2
      VAH=VAX+VAY+VAZ
      CKV=a2
      CSR=sqrt(CKV)
      CSV=CKV+VAH
      CWH=CKV-VAH
      DSCV=CSV**2-4.e0*CKV*VAX
      IF(DSCV.LT.0.) DSCV=0.
      DSCV=sqrt(DSCV)
      VF2=(CSV+DSCV)/2.e0
      VS2=(CSV-DSCV)/2.e0
      IF(VS2.LT.0.) VS2=0.
      CF=sqrt(VF2)
      CS=sqrt(VS2)
      VA=abs(BVX)

      BXX=QU(5)**2
      BYZ=QU(6)**2+QU(7)**2
      SBB=sqrt(BYZ)
      BMG=BXX+BYZ
      GPP=a2*QU(1)
      GPB=abs(GPP-BXX)

      IF(BYZ.LE.DLT*BMG) THEN

         IF(GPB.GT.DLT*GPP) THEN
            CS=sqrt(abs((CKV+VAX-SIGN(1.e0,CWH)*CWH)/2.e0))
            CF=sqrt(abs((CKV+VAX+SIGN(1.e0,CWH)*CWH)/2.e0))
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

      RETURN
      END

      SUBROUTINE VECTEGy(QU,U)
      Implicit NONE

      REAL QU(7),U(7)          
      REAL DLT,SK1,BVX,BVY,BVZ,VAX,VAY,VAZ
      REAL VAH,CKV,A1,A2,CSR,CSV,CWH,DSCV,VF2
      REAL VS2,CF,CS,VA,BYY,BXZ,SBB,BMG,GPP,GPB
      
      a1=1.
      a2=1.
      
      DLT=1.E-12

      SK1=sqrt(QU(1))
      BVX=QU(5)/SK1
      BVY=QU(6)/SK1
      BVZ=QU(7)/SK1
      VAX=BVX**2
      VAY=BVY**2
      VAZ=BVZ**2
      VAH=VAX+VAY+VAZ
      CKV=a2
      CSR=sqrt(CKV)
      CSV=CKV+VAH
      CWH=CKV-VAH
      DSCV=CSV**2-4.e0*CKV*VAY
      IF(DSCV.LT.0.) DSCV=0.
      DSCV=sqrt(DSCV)
      VF2=(CSV+DSCV)/2.e0
      VS2=(CSV-DSCV)/2.e0
      IF(VS2.LT.0.) VS2=0.
      CF=sqrt(VF2)
      CS=sqrt(VS2)
      VA=abs(BVY)

      BYY=QU(6)**2
      BXZ=QU(5)**2+QU(7)**2
      SBB=sqrt(BXZ)
      BMG=BYY+BXZ
      GPP=a2*QU(1)
      GPB=abs(GPP-BYY)

      IF(BXZ.LE.DLT*BMG) THEN

         IF(GPB.GT.DLT*GPP) THEN
            CS=sqrt(abs((CKV+VAY-SIGN(1.e0,CWH)*CWH)/2.e0))
            CF=sqrt(abs((CKV+VAY+SIGN(1.e0,CWH)*CWH)/2.e0))
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

      RETURN
      END

      SUBROUTINE VECTEGz(QU,U)
      Implicit NONE

      REAL QU(7),U(7)          
      REAL DLT,SK1,BVX,BVY,BVZ,VAX,VAY,VAZ,VA
      REAL VAH,CKV,A1,A2,CSR,CSV,CWH,DSCV,VF2
      REAL VS2,CF,CS,BXX,BYZ,SBB,BMG,GPP,GPB

      a1=1.
      a2=1.  

      DLT=1.E-12

      SK1=sqrt(QU(1))
      BVX=QU(5)/SK1
      BVY=QU(6)/SK1
      BVZ=QU(7)/SK1
      VAX=BVX**2
      VAY=BVY**2
      VAZ=BVZ**2
      VAH=VAX+VAY+VAZ
      CKV=a2
      CSR=sqrt(CKV)
      CSV=CKV+VAH
      CWH=CKV-VAH
      DSCV=CSV**2-4.e0*CKV*VAZ
      IF(DSCV.LT.0.) DSCV=0.
      DSCV=sqrt(DSCV)
      VF2=(CSV+DSCV)/2.e0
      VS2=(CSV-DSCV)/2.e0
      IF(VS2.LT.0.) VS2=0.
      CF=sqrt(VF2)
      CS=sqrt(VS2)
      VA=abs(BVZ)

      BXX=QU(7)**2
      BYZ=QU(5)**2+QU(6)**2
      SBB=sqrt(BYZ)
      BMG=BXX+BYZ
      GPP=a2*QU(1)
      GPB=abs(GPP-BXX)

      IF(BYZ.LE.DLT*BMG) THEN

         IF(GPB.GT.DLT*GPP) THEN
            CS=sqrt(abs((CKV+VAZ-SIGN(1.e0,CWH)*CWH)/2.e0))
            CF=sqrt(abs((CKV+VAZ+SIGN(1.e0,CWH)*CWH)/2.e0))
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

      RETURN
      END

      SUBROUTINE FLUXx(QU,F)
      Implicit NONE
      
      REAL QU(7),F(7)
      REAL A1,A2,BKV,PXX 

      a1=1.
      a2=1.  

      BKV=QU(5)**2+QU(6)**2+QU(7)**2

      PXX=a2*qu(1)+BKV/2.e0

      F(1)=qu(1)*qu(2)
      F(2)=qu(1)*qu(2)**2+PXX-QU(5)**2
      F(3)=qu(1)*qu(2)*QU(3)-QU(5)*QU(6)
      F(4)=qu(1)*qu(2)*QU(4)-QU(5)*QU(7)
      F(5)=0.
      F(6)=qu(2)*QU(6)-QU(3)*QU(5)
      F(7)=qu(2)*QU(7)-QU(4)*QU(5)

      RETURN
      END

      SUBROUTINE FLUXy(QU,F)
      Implicit NONE
      
      REAL QU(7),F(7)
      REAL A1,A2,BKV,PXX 

      A1=1.
      A2=1.  

      BKV=QU(5)**2+QU(6)**2+QU(7)**2

      PXX=A2*qu(1)+BKV/2.e0

      F(1)=qu(1)*qu(3)
      F(2)=qu(1)*qu(2)*QU(3)-QU(6)*QU(5)
      F(3)=qu(1)*qu(3)**2+PXX-QU(6)**2
      F(4)=qu(1)*qu(3)*QU(4)-QU(6)*QU(7)
      F(5)=QU(3)*QU(5)-qu(2)*QU(6)
      F(6)=0.
      F(7)=QU(3)*QU(7)-qu(4)*QU(6)

      RETURN
      END

      SUBROUTINE FLUXz(QU,F)
      Implicit NONE
      
      REAL QU(7),F(7)
      REAL A1,A2,BKV,PXX 

      A1=1.
      A2=1.  

      BKV=QU(5)**2+QU(6)**2+QU(7)**2

      PXX=A2*qu(1)+BKV/2.e0

      F(1)=qu(1)*qu(4)
      F(2)=qu(1)*qu(2)*QU(4)-QU(5)*QU(7)
      F(3)=qu(1)*qu(3)*QU(4)-QU(6)*QU(7)
      F(4)=qu(1)*qu(4)**2+PXX-QU(7)**2
      F(5)=QU(4)*QU(5)-qu(2)*QU(7)
      F(6)=QU(4)*QU(6)-qu(3)*QU(7)
      F(7)=0.

      RETURN
      END

      SUBROUTINE FPMR(Y,QL,QR,QP,QC) 
      Implicit NONE
      
      Integer m
      REAL QL(7),QR(7),QP(7),QC(7)
      REAL DQ,Q6,Y

      DO M=1,7
         DQ=QR(M)-QL(M)
         Q6=6.0*(QP(M)-(QR(M)+QL(M))/2.e0)
         QC(M)=QL(M)+Y/2.e0*(DQ+(1.e0-2.e0/3.*Y)*Q6)
      ENDDO

      RETURN
      END
      
      SUBROUTINE FPML(Y,QL,QR,QP,QC) 
      Implicit NONE
      
      Integer m
      REAL QL(7),QR(7),QP(7),QC(7)
      REAL DQ,Q6,y

      DO M=1,7
         DQ=QR(M)-QL(M)
         Q6=6.0*(QP(M)-(QR(M)+QL(M))/2.e0)
         QC(M)=QR(M)-Y/2.e0*(DQ-(1.e0-2.e0/3.*Y)*Q6)
      ENDDO

      RETURN
      END
      
      SUBROUTINE FQML(Y,QL,QR,QP,QC) 
      Implicit NONE
      
      Integer m
      REAL QL(7),QR(7),QP(7),QC(7)
      REAL DQ,Q6,Y

      DO M=1,7
         DQ=QR(M)-QL(M)
         Q6=6.0*(QP(M)-(QR(M)+QL(M))/2.e0)
         QC(M)=QL(M)+Y*(DQ+(1.e0-Y)*Q6)
      ENDDO

      RETURN
      END

      subroutine averx(qul,qur,qsr)
      Implicit NONE
      
      REAL qul(7),qur(7),qsr(7)
      REAL RL,RR,RS

      rl=sqrt(qul(1))
      rr=sqrt(qur(1))
      rs=rl+rr
      
      qsr(1)=rl*rr
      qsr(2)=(rl*qul(2)+rr*qur(2))/rs
      qsr(3)=(rl*qul(3)+rr*qur(3))/rs
      qsr(4)=(rl*qul(4)+rr*qur(4))/rs
      qsr(5)=(qul(5)+qur(5))/2.e0
      qsr(6)=(qul(6)/rl+qur(6)/rr)/rs*qsr(1)
      qsr(7)=(qul(7)/rl+qur(7)/rr)/rs*qsr(1)

      return
      end

      subroutine avery(qul,qur,qsr)
      Implicit NONE
      
      REAL qul(7),qur(7),qsr(7)
      REAL RL,RR,RS
      
      rl=sqrt(qul(1))
      rr=sqrt(qur(1))
      rs=rl+rr
      
      qsr(1)=rl*rr
      qsr(2)=(rl*qul(2)+rr*qur(2))/rs
      qsr(3)=(rl*qul(3)+rr*qur(3))/rs
      qsr(4)=(rl*qul(4)+rr*qur(4))/rs
      qsr(5)=(qul(5)/rl+qur(5)/rr)/rs*qsr(1)
      qsr(6)=(qul(6)+qur(6))/2.e0
      qsr(7)=(qul(7)/rl+qur(7)/rr)/rs*qsr(1)

      return
      end

      subroutine averz(qul,qur,qsr)
      Implicit NONE
      
      REAL qul(7),qur(7),qsr(7)
      REAL RL,RR,RS
      
      rl=sqrt(qul(1))
      rr=sqrt(qur(1))
      rs=rl+rr
      
      qsr(1)=rl*rr
      qsr(2)=(rl*qul(2)+rr*qur(2))/rs
      qsr(3)=(rl*qul(3)+rr*qur(3))/rs
      qsr(4)=(rl*qul(4)+rr*qur(4))/rs
      qsr(5)=(qul(5)/rl+qur(5)/rr)/rs*qsr(1)
      qsr(6)=(qul(6)/rl+qur(6)/rr)/rs*qsr(1)
      qsr(7)=(qul(7)+qur(7))/2.e0

      return
      end


      Subroutine QDD6(Nx,Qvr,Qvl,Qrx,Qlx,Qp)
      Implicit NONE
      Integer m,nx,i
      REAL Qvr(Nx,7),Qvl(Nx,7),Qp(Nx,7) 
      REAL Qrx(Nx,7),Qlx(Nx,7) 
      REAL DQ,Q6

      Do I=1,Nx        
         
         DO M=1,7

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
      
      SUBROUTINE MATR_AX(QU,QL)
      Implicit NONE

      REAL QU(7),QL(7,7)          
      REAL A1,A2

      A1=1.
      A2=1.  

      QL(1,1)=QU(2)
      QL(1,2)=QU(1)
      QL(1,3)=0.
      QL(1,4)=0.
      QL(1,5)=0.
      QL(1,6)=0.
      QL(1,7)=0.
      
      QL(2,1)=a2/QU(1)
      QL(2,2)=QU(2)
      QL(2,3)=0.
      QL(2,4)=0.
      QL(2,5)=0.
      QL(2,6)=QU(6)/QU(1)
      QL(2,7)=QU(7)/QU(1)

      QL(3,1)=0.
      QL(3,2)=0.
      QL(3,3)=QU(2)
      QL(3,4)=0.
      QL(3,5)=0.
      QL(3,6)=-QU(5)/QU(1)
      QL(3,7)=0.

      QL(4,1)=0.
      QL(4,2)=0.
      QL(4,3)=0.
      QL(4,4)=QU(2)
      QL(4,5)=0.
      QL(4,6)=0.
      QL(4,7)=-QU(5)/QU(1)

      QL(5,1)=0.
      QL(5,2)=0.
      QL(5,3)=0.
      QL(5,4)=0.
      QL(5,5)=QU(2)
      QL(5,6)=0.
      QL(5,7)=0.

      QL(6,1)=0.
      QL(6,2)=QU(6)
      QL(6,3)=-QU(5)
      QL(6,4)=0.
      QL(6,5)=0.
      QL(6,6)=QU(2)
      QL(6,7)=0.

      QL(7,1)=0.
      QL(7,2)=QU(7)
      QL(7,3)=0.
      QL(7,4)=-QU(5)
      QL(7,5)=0.
      QL(7,6)=0.
      QL(7,7)=QU(2)

      RETURN
      END

      SUBROUTINE MATR_AY(QU,QL)
      Implicit NONE

      REAL QU(7),QL(7,7)          
      REAL A1,A2

      A1=1.
      A2=1.  
      
      QL(1,1)=QU(3)
      QL(1,2)=0.
      QL(1,3)=QU(1)
      QL(1,4)=0.
      QL(1,5)=0.
      QL(1,6)=0.
      QL(1,7)=0.
      
      QL(2,1)=0.
      QL(2,2)=QU(3)
      QL(2,3)=0.
      QL(2,4)=0.
      QL(2,5)=-QU(6)/QU(1)
      QL(2,6)=0.
      QL(2,7)=0.

      QL(3,1)=a2/QU(1)
      QL(3,2)=0.
      QL(3,3)=QU(3)
      QL(3,4)=0.
      QL(3,5)=QU(5)/QU(1)
      QL(3,6)=0.
      QL(3,7)=QU(7)/QU(1)

      QL(4,1)=0.
      QL(4,2)=0.
      QL(4,3)=0.
      QL(4,4)=QU(3)
      QL(4,5)=0.
      QL(4,6)=0.
      QL(4,7)=-QU(6)/QU(1)

      QL(5,1)=0.
      QL(5,2)=-QU(6)
      QL(5,3)=QU(5)
      QL(5,4)=0.
      QL(5,5)=QU(3)
      QL(5,6)=0.
      QL(5,7)=0.

      QL(6,1)=0.
      QL(6,2)=0.
      QL(6,3)=0.
      QL(6,4)=0.
      QL(6,5)=0.
      QL(6,6)=QU(3)
      QL(6,7)=0.

      QL(7,1)=0.
      QL(7,2)=0.
      QL(7,3)=QU(7)
      QL(7,4)=-QU(6)
      QL(7,5)=0.
      QL(7,6)=0.
      QL(7,7)=QU(3)

      RETURN
      END

      SUBROUTINE MATR_AZ(QU,QL)
      Implicit NONE

      REAL QU(7),QL(7,7)          
      REAL A1,A2

      A1=1.
      A2=1.  

      QL(1,1)=QU(4)
      QL(1,2)=0.
      QL(1,3)=0.
      QL(1,4)=QU(1)
      QL(1,5)=0.
      QL(1,6)=0.
      QL(1,7)=0.
      
      QL(2,1)=0.
      QL(2,2)=QU(4)
      QL(2,3)=0.
      QL(2,4)=0.
      QL(2,5)=-QU(7)/QU(1)
      QL(2,6)=0.
      QL(2,7)=0.

      QL(3,1)=0.
      QL(3,2)=0.
      QL(3,3)=QU(4)
      QL(3,4)=0.
      QL(3,5)=0.
      QL(3,6)=-QU(7)/QU(1)
      QL(3,7)=0.

      QL(4,1)=a2/QU(1)
      QL(4,2)=0.
      QL(4,3)=0.
      QL(4,4)=QU(4)
      QL(4,5)=QU(5)/QU(1)
      QL(4,6)=QU(6)/QU(1)
      QL(4,7)=0.

      QL(5,1)=0.
      QL(5,2)=-QU(7)
      QL(5,3)=0.
      QL(5,4)=QU(5)
      QL(5,5)=QU(4)
      QL(5,6)=0.
      QL(5,7)=0.

      QL(6,1)=0.
      QL(6,2)=0.
      QL(6,3)=-QU(7)
      QL(6,4)=QU(6)
      QL(6,5)=0.
      QL(6,6)=QU(4)
      QL(6,7)=0.

      QL(7,1)=0.
      QL(7,2)=0.
      QL(7,3)=0.
      QL(7,4)=0.
      QL(7,5)=0.
      QL(7,6)=0.
      QL(7,7)=QU(4)

      RETURN
      END

      SUBROUTINE DUDW(Qu,QL)
      Implicit NONE

      REAL QU(7),QL(7,7)          
      
      QL(1,1)=1.
      QL(1,2)=0.
      QL(1,3)=0.
      QL(1,4)=0.
      QL(1,5)=0.
      QL(1,6)=0.
      QL(1,7)=0.
      
      QL(2,1)=QU(2)
      QL(2,2)=QU(1)
      QL(2,3)=0.
      QL(2,4)=0.
      QL(2,5)=0.
      QL(2,6)=0.
      QL(2,7)=0.

      QL(3,1)=QU(3)
      QL(3,2)=0.
      QL(3,3)=QU(1)
      QL(3,4)=0.
      QL(3,5)=0.
      QL(3,6)=0.
      QL(3,7)=0.

      QL(4,1)=QU(4)
      QL(4,2)=0.
      QL(4,3)=0.
      QL(4,4)=QU(1)
      QL(4,5)=0.
      QL(4,6)=0.
      QL(4,7)=0.

      QL(5,1)=0.
      QL(5,2)=0.
      QL(5,3)=0.
      QL(5,4)=0.
      QL(5,5)=1.
      QL(5,6)=0.
      QL(5,7)=0.

      QL(6,1)=0.
      QL(6,2)=0.
      QL(6,3)=0.
      QL(6,4)=0.
      QL(6,5)=0.
      QL(6,6)=1.
      QL(6,7)=0.

      QL(7,1)=0.
      QL(7,2)=0.
      QL(7,3)=0.
      QL(7,4)=0.
      QL(7,5)=0.
      QL(7,6)=0.
      QL(7,7)=1.

      RETURN
      END

      SUBROUTINE AMPLTD(QL,QD,SL)
      Implicit NONE

      Integer m,l
      REAL QL(7,7),QD(7),SL(7)
      
      DO M=1,7
         SL(M)=0.   
         DO L=1,7
            SL(M)=SL(M)+QL(M,L)*QD(L)
         ENDDO
      ENDDO

      RETURN
      END 

      Function Xminmod(x,y)
      
      REAL x,y,Xminmod
      
      Xminmod = 0.5*(sign(1.e0,x) + sign(1.e0,y))*min(abs(x),abs(y))

      Return
      End

      Function Xmedian(x,y,z)
      
      REAL x,y,z,Xminmod,Xmedian

      External Xminmod

      Xmedian = x + Xminmod(y-x,z-x)

      Return
      End

      SUBROUTINE W5RECM(VP1,VP2,VP3,VP4,VP5,VL,VR)
      Implicit NONE
      
      REAL IS1,IS2,IS3,DELTA,WS1,WS2,WS3,WJ1,WJ2,WJ3
      REAL AJ1,AJ2,AJ3,AJS,WM1,WM2,WM3,VL,VR
      REAL VP1,VP2,VP3,VP4,VP5
      
      DELTA = 1.0d-20
      
      WS1 = 0.1E0
      WS2 = 0.6E0
      WS3 = 0.3E0 


      WJ1 = (11.E0*VP3 - 7.E0*VP2 + 2.E0*VP1)/6.E0
      WJ2 = (2.E0 *VP4 + 5.E0*VP3 -      VP2)/6.E0
      WJ3 = (     -VP5 + 5.E0*VP4 + 2.E0*VP3)/6.E0

      IS1 = 13.E0/12.E0*(VP3 - 2.E0*VP2 + VP1)**2 
     1     + (3.E0*VP3 - 4.E0*VP2 + VP1)**2/4.E0
      IS2 = 13.E0/12.E0*(VP4 - 2.E0*VP3 + VP2)**2 
     1     + (VP4 - VP2)**2/4.E0
      IS3 = 13.E0/12.E0*(VP5 - 2.E0*VP4 + VP3)**2 
     1     + (VP5 - 4.E0*VP4 + 3.E0*VP3)**2/4.E0

      AJ1 = WS1/(IS1+DELTA)**2   
      AJ2 = WS2/(IS2+DELTA)**2   
      AJ3 = WS3/(IS3+DELTA)**2   

      AJS = AJ1 + AJ2 +AJ3 

      WM1 = AJ1/AJS
      WM2 = AJ2/AJS
      WM3 = AJ3/AJS

      AJ1 = WM1*(WS1 + WS1*WS1 - 3.E0*WS1*WM1 + WM1*WM1)/(WS1*WS1 
     1     + WM1*(1.E0 - 2.E0*WS1)) 
      AJ2 = WM2*(WS2 + WS2*WS2 - 3.E0*WS2*WM2 + WM2*WM2)/(WS2*WS2 
     1     + WM2*(1.E0 - 2.E0*WS2)) 
      AJ3 = WM3*(WS3 + WS3*WS3 - 3.E0*WS3*WM3 + WM3*WM3)/(WS3*WS3 
     1     + WM3*(1.E0 - 2.E0*WS3)) 

      AJS = AJ1 + AJ2 +AJ3 

      WM1 = AJ1/AJS
      WM2 = AJ2/AJS
      WM3 = AJ3/AJS

      VR = WM1*WJ1 + WM2*WJ2 + WM3*WJ3 

      WJ1 = (11.E0*VP3 - 7.E0*VP4 + 2.E0*VP5)/6.E0
      WJ2 = (2.E0 *VP2 + 5.E0*VP3 -      VP4)/6.E0
      WJ3 = (     -VP1 + 5.E0*VP2 + 2.E0*VP3)/6.E0

      IS1 = 13.E0/12.E0*(VP3 - 2.E0*VP4 + VP5)**2 
     1     + (3.E0*VP3 - 4.E0*VP4 + VP5)**2/4.E0
      IS2 = 13.E0/12.E0*(VP2 - 2.E0*VP3 + VP4)**2 
     1     + (VP2 - VP4)**2/4.E0
      IS3 = 13.E0/12.E0*(VP1 - 2.E0*VP2 + VP3)**2 
     1     + (VP1 - 4.E0*VP2 + 3.E0*VP3)**2/4.E0

      AJ1 = WS1/(IS1+DELTA)**2   
      AJ2 = WS2/(IS2+DELTA)**2   
      AJ3 = WS3/(IS3+DELTA)**2   

      AJS = AJ1 + AJ2 +AJ3 

      WM1 = AJ1/AJS
      WM2 = AJ2/AJS
      WM3 = AJ3/AJS

      AJ1 = WM1*(WS1 + WS1*WS1 - 3.E0*WS1*WM1 + WM1*WM1)/(WS1*WS1 
     1     + WM1*(1.E0 - 2.E0*WS1)) 
      AJ2 = WM2*(WS2 + WS2*WS2 - 3.E0*WS2*WM2 + WM2*WM2)/(WS2*WS2 
     1     + WM2*(1.E0 - 2.E0*WS2)) 
      AJ3 = WM3*(WS3 + WS3*WS3 - 3.E0*WS3*WM3 + WM3*WM3)/(WS3*WS3 
     1     + WM3*(1.E0 - 2.E0*WS3)) 

      AJS = AJ1 + AJ2 +AJ3 

      WM1 = AJ1/AJS
      WM2 = AJ2/AJS
      WM3 = AJ3/AJS

      VL = WM1*WJ1 + WM2*WJ2 + WM3*WJ3 

      RETURN
      END
      
      Function VLeer(a,b)

      REAL a,b,VLeer

      if(a*a+b*b.ne.0)then
         VLeer=max(a*b,0.)*(a+b)/(a*a+b*b)

      else
         VLeer=0.
      endif 

      Return
      End        
      
      SUBROUTINE PRIM(QU,QP)

      REAL QU(7),QP(7)

      QP(1)=QU(1)
      QP(2)=QU(2)/QU(1)
      QP(3)=QU(3)/QU(1)
      QP(4)=QU(4)/QU(1)
      QP(5)=QU(5)
      QP(6)=QU(6)
      QP(7)=QU(7)

      RETURN
      END

      SUBROUTINE CONS(QP,QU)

      REAL QU(7),QP(7)

      QU(1)=QP(1)
      QU(2)=QP(2)*QP(1)
      QU(3)=QP(3)*QP(1)
      QU(4)=QP(4)*QP(1)
      QU(5)=QP(5)
      QU(6)=QP(6)
      QU(7)=QP(7)

      RETURN
      END  
