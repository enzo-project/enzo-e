c     $Id$
c     See LICENSE_ENZO file for license and copyright information

#include "fortran.h"

      Subroutine HLLDx(UL,UR,F,U) !SL and SR from min/max
      Implicit NONE
      
      ENZO_REAL U(7),F(7),UL(7),UR(7)

      Integer i
      ENZO_REAL Uhll(7),Fhll(7)
      ENZO_REAL FL(7),FR(7)
      ENZO_REAL SL,SR,SL_rot,SR_rot
      ENZO_REAL cf,ca           ! fast and Alfven velocities
      ENZO_REAL Bx,By,Bz,BB
      ENZO_REAL u_L,v_L,w_L,u_R,v_R,w_R
      ENZO_REAL ptL, ptR
      ENZO_REAL u_ast
      ENZO_REAL UL_ast(7),UR_ast(7)

!temporary
      ENZO_REAL f2,f3,X,g3,g4,g5,g6
      ENZO_REAL S1L, S1R, S2L, S2R
      ENZO_REAL a,a1,a2,ds


      ENZO_REAL one
      parameter(one = 1.0)

      a1=1.
      a2=1.  

      a=a1
      
      u_L=UL(2)/UL(1)
      v_L=UL(3)/UL(1)
      w_L=UL(4)/UL(1)

      u_R=UR(2)/UR(1)
      v_R=UR(3)/UR(1)
      w_R=UR(4)/UR (1)
      
      Bx=(Ul(5)+Ur(5))/2.e0

      By=UL(6)
      Bz=UL(7)
      
      BB=Bx*Bx+By*By+Bz*Bz
      ca=abs(Bx)/sqrt(UL(1))
      f2=a*a+BB/UL(1)
      ds=f2*f2-4.e0*a*a*ca*ca
      if(ds.lt.0.) ds=0.
      cf=sqrt(0.5e0*(f2+sqrt(ds)))
      ptL=a*a*UL(1)+BB/2.e0

      S1L=u_L-cf
      S1R=u_L+cf
      
      By=UR(6)
      Bz=UR(7)
      
      BB=Bx*Bx+By*By+Bz*Bz
      ca=abs(Bx)/sqrt(UR(1))
      f2=a*a+BB/UR(1)
      ds=f2*f2-4.e0*a*a*ca*ca
      if(ds.lt.0.) ds=0.
      cf=sqrt(0.5e0*(f2+sqrt(ds)))
      
      ptR=a*a*UR(1)+BB/2.e0

      S2L=u_R-cf
      S2R=u_R+cf

      SL=min(S1L,S2L)
      SR=max(S1R,S2R)

      if (SL.ge.0.) then

         F(1)=UL(2)
         F(2)=UL(2)*u_L+ptL-Ul(5)*Ul(5)
         F(3)=UL(3)*u_L-Ul(5)*Ul(6)
         F(4)=UL(4)*u_L-Ul(5)*Ul(7)
         F(5)=0.
         F(6)=UL(6)*u_L-Ul(5)*v_L
         F(7)=UL(7)*u_L-Ul(5)*w_L

         U(1)=UL(1)
         U(2)=UL(2)
         U(3)=UL(3)
         U(4)=UL(4)
         U(5)=UL(5)
         U(6)=UL(6)
         U(7)=UL(7)

         goto 8888
      endif

      if (SR.le.0.) then
         F(1)=UR(2)
         F(2)=UR(2)*u_R+ptR-Ur(5)*Ur(5)
         F(3)=UR(3)*u_R-Ur(5)*Ur(6)
         F(4)=UR(4)*u_R-Ur(5)*Ur(7)
         F(5)=0.
         F(6)=UR(6)*u_R-Ur(5)*v_R
         F(7)=UR(7)*u_R-Ur(5)*w_R
         U(1)=UR(1)
         U(2)=UR(2)
         U(3)=UR(3)
         U(4)=UR(4)
         U(5)=UR(5)
         U(6)=UR(6)
         U(7)=UR(7)
         goto 8888
      endif

      Fl(1)=UL(2)
      Fl(2)=UL(2)*u_L+ptL-Ul(5)*Ul(5)
      Fl(3)=UL(3)*u_L-Ul(5)*Ul(6)
      Fl(4)=UL(4)*u_L-Ul(5)*Ul(7)
      Fl(5)=0.
      Fl(6)=UL(6)*u_L-Ul(5)*v_L
      Fl(7)=UL(7)*u_L-Ul(5)*w_L

      Fr(1)=UR(2)
      Fr(2)=UR(2)*u_R+ptR-Ur(5)*Ur(5)
      Fr(3)=UR(3)*u_R-Ur(5)*Ur(6)
      Fr(4)=UR(4)*u_R-Ur(5)*Ur(7)
      Fr(5)=0.
      Fr(6)=UR(6)*u_R-Ur(5)*v_R
      Fr(7)=UR(7)*u_R-Ur(5)*w_R

      Do i=1,7
         Uhll(i)=(SR*UR(i)-SL*UL(i)-FR(i)+FL(i))/(SR-SL)
         Fhll(i)=(SR*FL(i)-SL*FR(i)+SL*SR*(UR(i)-UL(i)))/(SR-SL)
      enddo

      u_ast=Fhll(1)/Uhll(1)
      ca=abs(Bx)/sqrt(Uhll(1))
      
      SL_rot=u_ast-ca
      SR_rot=u_ast+ca

      UL_ast(1)=Uhll(1)
      UL_ast(2)=Uhll(2)
      
      f3=(SL-SL_rot)*(SL-SR_rot)

      if (f3.ne.0.) then
         UL_ast(3)=Uhll(1)*v_L-Bx*UL(6)*(u_ast-u_L)/f3
         UL_ast(4)=Uhll(1)*w_L-Bx*UL(7)*(u_ast-u_L)/f3
         Ul_ast(5)=bx
         UL_ast(6)=UL(6)*(UL(1)*((SL-u_L)**2)-Bx*Bx)/Uhll(1)/f3
         UL_ast(7)=UL(7)*(UL(1)*((SL-u_L)**2)-Bx*Bx)/Uhll(1)/f3
      else
         UL_ast(3)=UL(3)
         UL_ast(4)=UL(4)
         Ul_ast(5)=bx
         UL_ast(6)=UL(6)
         UL_ast(7)=UL(7)
      endif

      if ((SL.lt.0.).and.(0.le.SL_rot)) then
         Do i=1,7
            F(i)=FL(i)+SL*(UL_ast(i)-UL(i))
            U(i)=UL_ast(i)
         enddo
         goto 8888
      endif

      UR_ast(1)=Uhll(1)
      UR_ast(2)=Uhll(2)
      
      f3=(SR-SL_rot)*(SR-SR_rot)

      if (f3.ne.0.) then
         UR_ast(3)=Uhll(1)*v_R-Bx*UR(6)*(u_ast-u_R)/f3
         UR_ast(4)=Uhll(1)*w_R-Bx*UR(7)*(u_ast-u_R)/f3
         UR_ast(5)=bx
         UR_ast(6)=UR(6)*(UR(1)*((SR-u_R)**2)-Bx*Bx)/Uhll(1)/f3
         UR_ast(7)=UR(7)*(UR(1)*((SR-u_R)**2)-Bx*Bx)/Uhll(1)/f3
      else
         UR_ast(3)=UR(3)
         UR_ast(4)=UR(4)
         UR_ast(5)=bx
         UR_ast(6)=UR(6)
         UR_ast(7)=UR(7)
      endif

      if ((SR_rot.lt.0.).and.(0.le.SR)) then
         Do i=1,7
            F(i)=FR(i)+SR*(UR_ast(i)-UR(i))
            U(i)=UR_ast(i)
         enddo
         goto 8888
      endif

      X=sqrt(Uhll(1))*sign(one,Bx)

      g3=(UL_ast(3)+UR_ast(3)+X*(UR_ast(6)-UL_ast(6)))/2.e0
      g4=(UL_ast(4)+UR_ast(4)+X*(UR_ast(7)-UL_ast(7)))/2.e0
      g5=(UL_ast(6)+UR_ast(6)+(UR_ast(3)-UL_ast(3))/X)/2.e0
      g6=(UL_ast(7)+UR_ast(7)+(UR_ast(4)-UL_ast(4))/X)/2.e0


      F(1)=Fhll(1)
      F(2)=Fhll(2)
      F(3)=g3*u_ast-Bx*g5
      F(4)=g4*u_ast-Bx*g6
      F(5)=0.
      F(6)=g5*u_ast-Bx*g3/uhll(1)
      F(7)=g6*u_ast-Bx*g4/uhll(1)
      
      U(1)=Uhll(1)
      U(2)=Uhll(2)
      U(3)=g3
      U(4)=g4
      U(5)=bx        
      U(6)=g5
      U(7)=g6

c     if (.not.((SL_rot.le.0.).and.(0.le.SR_rot))) then
c     ERROR_MESSAGE

c     print*,'Error in HLLDS procedure! Press any key.' 
c     pause
c     endif

 8888 continue 

c     print*,'hlld finished'
c     pause

      End                       ! Subroutine HLLDx

      Subroutine HLLDy(UL,UR,F,U) !SL and SR from min/max
      
      ENZO_REAL U(7),F(7),UL(7),UR(7)

      Integer i
      ENZO_REAL Uhll(7),Fhll(7)
      ENZO_REAL FL(7),FR(7)
      ENZO_REAL SL,SR,SL_rot,SR_rot
      ENZO_REAL cf,ca           ! fast and Alfven velocities
      ENZO_REAL Bx,By,Bz,BB
      ENZO_REAL u_L,v_L,w_L,u_R,v_R,w_R
      ENZO_REAL ptL, ptR
      ENZO_REAL u_ast
      ENZO_REAL UL_ast(7),UR_ast(7)

!temporary
      ENZO_REAL f2,f3,X,g2,g4,g5,g7
      ENZO_REAL S1L, S1R, S2L, S2R
      ENZO_REAL a,a1,a2,ds


      ENZO_REAL one
      parameter(one = 1.0)

!     bug                common/aspd/a1,a2  
      a1=1.
      a2=1.

      a=a1
      
      u_L=UL(2)/UL(1)
      v_L=UL(3)/UL(1)
      w_L=UL(4)/UL(1)

      u_R=UR(2)/UR(1)
      v_R=UR(3)/UR(1)
      w_R=UR(4)/UR(1)
      
      By=(Ul(6)+Ur(6))/2.e0

      Bx=UL(5)
      Bz=UL(7)
      
      BB=Bx*Bx+By*By+Bz*Bz
      ca=abs(By)/sqrt(UL(1))
      f2=a*a+BB/UL(1)
      ds=f2*f2-4.e0*a*a*ca*ca
      if(ds.lt.0.) ds=0.
      cf=sqrt(0.5e0*(f2+sqrt(ds)))
      ptL=a*a*UL(1)+BB/2.e0

      S1L=v_L-cf
      S1R=v_L+cf
      
      Bx=UR(5)
      Bz=UR(7)
      
      BB=Bx*Bx+By*By+Bz*Bz
      ca=abs(By)/sqrt(UR(1))
      f2=a*a+BB/UR(1)
      ds=f2*f2-4.e0*a*a*ca*ca
      if(ds.lt.0.) ds=0.
      cf=sqrt(0.5e0*(f2+sqrt(ds)))
      ptR=a*a*UR(1)+BB/2.e0

      S2L=v_R-cf
      S2R=v_R+cf

      SL=min(S1L,S2L)
      SR=max(S1R,S2R)

      if (SL.ge.0.) then

         F(1)=UL(3)
         F(2)=UL(2)*v_L-Ul(5)*Ul(6)
         F(3)=UL(3)*v_L+ptl-Ul(6)*Ul(6)
         F(4)=UL(4)*v_L-Ul(6)*Ul(7)
         F(5)=UL(5)*v_L-Ul(6)*u_L
         F(6)=0.
         F(7)=UL(7)*v_L-Ul(6)*w_L

         U(1)=UL(1)
         U(2)=UL(2)
         U(3)=UL(3)
         U(4)=UL(4)
         U(5)=UL(5)
         U(6)=UL(6)
         U(7)=UL(7)

         goto 7777
      endif

      if (SR.le.0.) then

         F(1)=UR(3)
         F(2)=UR(2)*v_R-Ur(5)*UR(6)
         F(3)=UR(3)*v_R+ptR-UR(6)*UR(6)
         F(4)=UR(4)*v_R-UR(6)*UR(7)
         F(5)=UR(5)*v_R-UR(6)*u_R
         F(6)=0.
         F(7)=UR(7)*v_R-UR(6)*w_R

         U(1)=UR(1)
         U(2)=UR(2)
         U(3)=UR(3)
         U(4)=UR(4)
         U(5)=UR(5)
         U(6)=UR(6)
         U(7)=UR(7)

         goto 7777
      endif

      Fl(1)=UL(3)
      Fl(2)=UL(2)*v_L-Ul(5)*Ul(6)
      Fl(3)=UL(3)*v_L+ptl-Ul(6)*Ul(6)
      Fl(4)=UL(4)*v_L-Ul(6)*Ul(7)
      Fl(5)=UL(5)*v_L-Ul(6)*u_L
      Fl(6)=0.
      Fl(7)=UL(7)*v_L-Ul(6)*w_L

      Fr(1)=UR(3)
      Fr(2)=UR(2)*v_R-UR(5)*UR(6)
      Fr(3)=UR(3)*v_R+ptR-UR(6)*UR(6)
      Fr(4)=UR(4)*v_R-UR(6)*UR(7)
      Fr(5)=UR(5)*v_R-UR(6)*u_R
      Fr(6)=0.
      Fr(7)=UR(7)*v_R-UR(6)*w_R

      Do i=1,7
         Uhll(i)=(SR*UR(i)-SL*UL(i)-FR(i)+FL(i))/(SR-SL)
         Fhll(i)=(SR*FL(i)-SL*FR(i)+SL*SR*(UR(i)-UL(i)))/(SR-SL)
      enddo

      u_ast=Fhll(1)/Uhll(1)
      ca=abs(By)/sqrt(Uhll(1))
      
      SL_rot=u_ast-ca
      SR_rot=u_ast+ca

      UL_ast(1)=Uhll(1)
      UL_ast(3)=Uhll(3)
      
      f3=(SL-SL_rot)*(SL-SR_rot)

      if (f3.ne.0.) then
         UL_ast(2)=Uhll(1)*u_L-By*UL(5)*(u_ast-v_L)/f3
         UL_ast(4)=Uhll(1)*w_L-By*UL(7)*(u_ast-v_L)/f3
         UL_ast(5)=UL(5)*(UL(1)*((SL-v_L)**2)-By*By)/Uhll(1)/f3
         UL_ast(6)=by
         UL_ast(7)=UL(7)*(UL(1)*((SL-v_L)**2)-By*By)/Uhll(1)/f3
      else
         UL_ast(2)=UL(2)
         UL_ast(4)=UL(4)
         UL_ast(5)=UL(5)
         UL_ast(6)=by
         UL_ast(7)=UL(7)
      endif

      if ((SL.lt.0.).and.(0.le.SL_rot)) then
         Do i=1,7
            F(i)=FL(i)+SL*(UL_ast(i)-UL(i))
            U(i)=UL_ast(i)
         enddo
         goto 7777
      endif

      UR_ast(1)=Uhll(1)
      UR_ast(3)=Uhll(3)
      
      f3=(SR-SL_rot)*(SR-SR_rot)

      if (f3.ne.0.) then
         UR_ast(2)=Uhll(1)*u_R-By*UR(5)*(u_ast-v_R)/f3
         UR_ast(4)=Uhll(1)*w_R-By*UR(7)*(u_ast-v_R)/f3
         UR_ast(5)=UR(5)*(UR(1)*((SR-v_R)**2)-By*By)/Uhll(1)/f3
         UR_ast(6)=by
         UR_ast(7)=UR(7)*(UR(1)*((SR-v_R)**2)-By*By)/Uhll(1)/f3
      else
         UR_ast(2)=UR(2)
         UR_ast(4)=UR(4)
         UR_ast(5)=UR(5)
         UR_ast(6)=by
         UR_ast(7)=UR(7)
      endif

      if ((SR_rot.lt.0.).and.(0.le.SR)) then
         Do i=1,7
            F(i)=FR(i)+SR*(UR_ast(i)-UR(i))
            U(i)=UR_ast(i)
         enddo
         goto 7777
      endif

      X=sqrt(Uhll(1))*sign(one,By)
      
      g2=(UL_ast(2)+UR_ast(2)+X*(UR_ast(5)-UL_ast(5)))/2.e0
      g4=(UL_ast(4)+UR_ast(4)+X*(UR_ast(7)-UL_ast(7)))/2.e0
      g5=(UL_ast(5)+UR_ast(5)+(UR_ast(2)-UL_ast(2))/X)/2.e0
      g7=(UL_ast(7)+UR_ast(7)+(UR_ast(4)-UL_ast(4))/X)/2.e0

      F(1)=Fhll(1)
      F(2)=g2*u_ast-By*g5
      F(3)=Fhll(3)
      F(4)=g4*u_ast-By*g7
      F(5)=g5*u_ast-By*g2/uhll(1)
      F(6)=0.
      F(7)=g7*u_ast-By*g4/uhll(1)
      
      U(1)=Uhll(1)
      U(2)=g2
      U(3)=Uhll(3)
      U(4)=g4
      U(5)=g5        
      U(6)=by
      U(7)=g7

c     if (.not.((SL_rot.le.0.).and.(0.le.SR_rot))) then
c     ERROR_MESSAGE

c     print*,'Error in HLLDS procedure! Press any key.' 
c     pause
c     endif

7777        continue 

c      print*,'hlld finished'
c        pause

      End ! Subroutine HLLDy


        Subroutine HLLDz(UL,UR,F,U) !SL and SR from min/max
        
        ENZO_REAL U(7),F(7),UL(7),UR(7)

        Integer i
        ENZO_REAL Uhll(7),Fhll(7)
        ENZO_REAL FL(7),FR(7)
        ENZO_REAL SL,SR,SL_rot,SR_rot
        ENZO_REAL cf,ca  ! fast and Alfven velocities
        ENZO_REAL Bx,By,Bz,BB
        ENZO_REAL u_L,v_L,w_L,u_R,v_R,w_R
        ENZO_REAL ptL, ptR
        ENZO_REAL u_ast
        ENZO_REAL UL_ast(7),UR_ast(7)

        !temporary
        ENZO_REAL f2,f3,X,g2,g3,g5,g6
        ENZO_REAL S1L, S1R, S2L, S2R
        ENZO_REAL a,a1,a2,ds

        ENZO_REAL one
        parameter(one = 1.0)

!bug                common/aspd/a1,a2  
        a1=1.
        a2=1.

        a=a1
        
        u_L=UL(2)/UL(1)
        v_L=UL(3)/UL(1)
        w_L=UL(4)/UL(1)

        u_R=UR(2)/UR(1)
        v_R=UR(3)/UR(1)
        w_R=UR(4)/UR(1)
        
        Bz=(Ul(7)+Ur(7))/2.e0

        Bx=UL(5)
        By=UL(6)
                
        BB=Bx*Bx+By*By+Bz*Bz
        ca=abs(Bz)/sqrt(UL(1))
        f2=a*a+BB/UL(1)
        ds=f2*f2-4.e0*a*a*ca*ca
        if(ds.lt.0.) ds=0.
        cf=sqrt(0.5e0*(f2+sqrt(ds)))
        ptL=a*a*UL(1)+BB/2.e0

        S1L=w_L-cf
        S1R=w_L+cf
      
        Bx=UR(5)
        By=UR(6)
                
        BB=Bx*Bx+By*By+Bz*Bz
        ca=abs(Bz)/sqrt(UR(1))
        f2=a*a+BB/UR(1)
        ds=f2*f2-4.e0*a*a*ca*ca
        if(ds.lt.0.) ds=0.
        cf=sqrt(0.5e0*(f2+sqrt(ds)))
        
        ptR=a*a*UR(1)+BB/2.e0

        S2L=w_R-cf
        S2R=w_R+cf

        SL=min(S1L,S2L)
        SR=max(S1R,S2R)

        if (SL.ge.0.) then

        F(1)=UL(4)
        F(2)=UL(2)*w_L-Ul(5)*Ul(7)
        F(3)=UL(3)*w_L-Ul(6)*Ul(7)
        F(4)=UL(4)*w_L+ptl-Ul(7)*Ul(7)
        F(5)=UL(5)*w_L-Ul(7)*u_L
        F(6)=UL(6)*w_L-Ul(7)*v_L
        F(7)=0.

        U(1)=UL(1)
        U(2)=UL(2)
        U(3)=UL(3)
        U(4)=UL(4)
        U(5)=UL(5)
        U(6)=UL(6)
        U(7)=UL(7)

        goto 6666
        endif

        if (SR.le.0.) then

        F(1)=UR(4)
        F(2)=UR(2)*w_R-UR(5)*UR(7)
        F(3)=UR(3)*w_R-UR(6)*UR(7)
        F(4)=UR(4)*w_R+ptR-UR(7)*UR(7)
        F(5)=UR(5)*w_R-UR(7)*u_R
        F(6)=UR(6)*w_R-UR(7)*v_R
        F(7)=0.

        U(1)=UR(1)
        U(2)=UR(2)
        U(3)=UR(3)
        U(4)=UR(4)
        U(5)=UR(5)
        U(6)=UR(6)
        U(7)=UR(7)

        goto 6666
        endif

        Fl(1)=UL(4)
        Fl(2)=UL(2)*w_L-Ul(5)*Ul(7)
        Fl(3)=UL(3)*w_L-Ul(6)*Ul(7)
        Fl(4)=UL(4)*w_L+ptl-Ul(7)*Ul(7)
        Fl(5)=UL(5)*w_L-Ul(7)*u_L
        Fl(6)=UL(6)*w_L-Ul(7)*v_L
        Fl(7)=0.

        Fr(1)=UR(4)
        Fr(2)=UR(2)*w_R-UR(5)*UR(7)
        Fr(3)=UR(3)*w_R-UR(6)*UR(7)
        Fr(4)=UR(4)*w_R+ptR-UR(7)*UR(7)
        Fr(5)=UR(5)*w_R-UR(7)*u_R
        Fr(6)=UR(6)*w_R-UR(7)*v_R
        Fr(7)=0.

        Do i=1,7
        Uhll(i)=(SR*UR(i)-SL*UL(i)-FR(i)+FL(i))/(SR-SL)
        Fhll(i)=(SR*FL(i)-SL*FR(i)+SL*SR*(UR(i)-UL(i)))/(SR-SL)
        enddo

        u_ast=Fhll(1)/Uhll(1)
        ca=abs(Bz)/sqrt(Uhll(1))
        
        SL_rot=u_ast-ca
        SR_rot=u_ast+ca

        UL_ast(1)=Uhll(1)
        UL_ast(4)=Uhll(4)
        
        f3=(SL-SL_rot)*(SL-SR_rot)

        if (f3.ne.0.) then
        UL_ast(2)=Uhll(1)*u_L-Bz*UL(5)*(u_ast-w_L)/f3
        UL_ast(3)=Uhll(1)*v_L-Bz*UL(6)*(u_ast-w_L)/f3
        UL_ast(5)=UL(5)*(UL(1)*((SL-w_L)**2)-Bz*Bz)/Uhll(1)/f3
        UL_ast(6)=UL(6)*(UL(1)*((SL-w_L)**2)-Bz*Bz)/Uhll(1)/f3
        Ul_ast(7)=bz
        else
        UL_ast(2)=UL(2)
        UL_ast(3)=UL(3)
        UL_ast(5)=UL(5)
        UL_ast(6)=UL(6)
        Ul_ast(7)=bz
        endif

        if ((SL.lt.0.).and.(0.le.SL_rot)) then
        Do i=1,7
        F(i)=FL(i)+SL*(UL_ast(i)-UL(i))
        U(i)=UL_ast(i)
        enddo
        goto 6666
        endif

        UR_ast(1)=Uhll(1)
        UR_ast(4)=Uhll(4)
        
        f3=(SR-SL_rot)*(SR-SR_rot)

        if (f3.ne.0.) then
        UR_ast(2)=Uhll(1)*u_R-Bz*UR(5)*(u_ast-w_R)/f3
        UR_ast(3)=Uhll(1)*v_R-Bz*UR(6)*(u_ast-w_R)/f3
        UR_ast(5)=UR(5)*(UR(1)*((SR-w_R)**2)-Bz*Bz)/Uhll(1)/f3
        UR_ast(6)=UR(6)*(UR(1)*((SR-w_R)**2)-Bz*Bz)/Uhll(1)/f3
        UR_ast(7)=bz
        else
        UR_ast(2)=UR(2)
        UR_ast(3)=UR(3)
        UR_ast(5)=UR(5)
        UR_ast(6)=UR(6)
        UR_ast(7)=bz
        endif

        if ((SR_rot.lt.0.).and.(0.le.SR)) then
        Do i=1,7
        F(i)=FR(i)+SR*(UR_ast(i)-UR(i))
        U(i)=UR_ast(i)
        enddo
        goto 6666
        endif

        X=sqrt(Uhll(1))*sign(one,Bz)

        g2=(UL_ast(2)+UR_ast(2)+X*(UR_ast(5)-UL_ast(5)))/2.e0
        g3=(UL_ast(3)+UR_ast(3)+X*(UR_ast(6)-UL_ast(6)))/2.e0
        g5=(UL_ast(5)+UR_ast(5)+(UR_ast(2)-UL_ast(2))/X)/2.e0
        g6=(UL_ast(6)+UR_ast(6)+(UR_ast(3)-UL_ast(3))/X)/2.e0


        F(1)=Fhll(1)
        F(2)=g2*u_ast-Bz*g5
        F(3)=g3*u_ast-Bz*g6
        F(4)=Fhll(4)
        F(5)=g5*u_ast-Bz*g2/uhll(1)
        F(6)=g6*u_ast-Bz*g3/uhll(1)
        F(7)=0.
        
        U(1)=Uhll(1)
        U(2)=g2
        U(3)=g3
        U(4)=Uhll(4)
      U(5)=g5        
        U(6)=g6
        U(7)=bz

c        if (.not.((SL_rot.le.0.).and.(0.le.SR_rot))) then
c              ERROR_MESSAGE

c        print*,'Error in HLLDS procedure! Press any key.' 
c        pause
c        endif

6666        continue 

c      print*,'hlld finished'
c        pause

      End ! Subroutine HLLDz
