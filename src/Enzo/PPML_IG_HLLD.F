c     See LICENSE_PPML file for license and copyright information

#include "fortran.h"

	Subroutine HLLDx_ig(gamma,pl,pr,bs,f,us)
	Implicit NONE

	Integer m
	ENZO_REAL ql(8),qr(8),f(8),qlz(8),qrz(8),qzz(8),fl(8),fr(8),us(8)
        ENZO_REAL ul,ur,cfl,cfr,ptl,ptr,sl,sr,slu,sru,smu,sm,ptz
	ENZO_REAL rlz,slz,rrz,srz,gamma
        ENZO_REAL pl(8),pr(8),bs(3)

        call cons_ig(gamma,pl,ql)
        call cons_ig(gamma,pr,qr)
	call pressx_ig(gamma,ql,bs,ul,ptl,cfl)
	call pressx_ig(gamma,qr,bs,ur,ptr,cfr)

	sl = min(ul,ur) - max(cfl,cfr)
	sr = max(ul,ur) + max(cfl,cfr)

	if(sl.gt.0.)then

	   call fluxx_ig(gamma,pl,bs,f)

	   do m=1,8
	      us(m) = ql(m)
	   enddo

	   return

	endif

	if(sr.lt.0.)then

	   call fluxx_ig(gamma,pr,bs,f)   

	   do m=1,8
	      us(m) = qr(m)
	   enddo

	   return

	endif

	slu = sl - ul
	sru = sr - ur
	smu = sru*qr(1)-slu*ql(1)
	sm  = (sru*qr(1)*ur - slu*ql(1)*ul - ptr + ptl)/smu

	ptz = (sru*qr(1)*ptl - slu*ql(1)*ptr + ql(1)*qr(1)*sru*slu
     1        *(ur - ul))/smu


	rlz = ql(1)*slu/(sl - sm)
	slz = sm - abs(ql(5)+bs(1))/sqrt(rlz) 

	if(sl.le.0.and.slz.ge.0.)then
  
	   call fluxx_ig(gamma,pl,bs,fl)
	   call quzx_ig(gamma,sm,sl,cfl,ptz,ql,bs,qlz)  
	
	   do m=1,8
	      f(m) = fl(m) + sl*(qlz(m) - ql(m))
	      us(m) = qlz(m)
	   enddo
	    
	   return

	endif


	rrz = qr(1)*sru/(sr - sm)
	srz = sm + abs(qr(5)+bs(1))/sqrt(rrz) !+bs(1)

	if(srz.le.0.and.sr.ge.0.)then

	   call fluxx_ig(gamma,pr,bs,fr) 
	   call quzx_ig(gamma,sm,sr,cfr,ptz,qr,bs,qrz) 
	    
	   do m=1,8
	      f(m) = fr(m) + sr*(qrz(m) - qr(m))
	      us(m) = qrz(m)
	   enddo
	
	   return

	endif

	call quzx_ig(gamma,sm,sl,cfl,ptz,ql,bs,qlz)
	call quzx_ig(gamma,sm,sr,cfr,ptz,qr,bs,qrz)

	if(slz.le.0.and.sm.ge.0.)then

	   call fluxx_ig(gamma,pl,bs,fl)            

	   call quzzx_ig(1,gamma,sm,qlz,qrz,bs,qzz) 

	   do m=1,8
	      f(m) = fl(m) + slz*qzz(m) - (slz - sl)*qlz(m) - sl*ql(m)
	      us(m) = qzz(m)
	   enddo

	   return

	endif

	if(sm.le.0.and.srz.ge.0.)then

	   call fluxx_ig(gamma,pr,bs,fr)

	   call quzzx_ig(2,gamma,sm,qlz,qrz,bs,qzz)

	   do m=1,8
	      f(m) = fr(m) + srz*qzz(m) - (srz - sr)*qrz(m) - sr*qr(m)
	      us(m) = qzz(m)
	   enddo

	   return
	   
	endif

	return
	end


	Subroutine quzx_ig(gamma,sm,sa,cfa,ptz,qa,bs,qz)
	Implicit NONE

	ENZO_REAL qa(8),qz(8),qb(8),bs(3)
        ENZO_REAL gamma,ua,v2,b2,pg,sam,sap,slu,slm,psb,vb,vbz,rlz
	ENZO_REAL rsm,smb,smr,vba,eps,sm,sa,cfa,ptz
	ENZO_REAL qb5,qb6,qb7,bp,bp0,pta,rst

        call prim_ig(gamma,qa,qb)

	ua = qb(2)
	v2 = qb(2)**2 + qb(3)**2 + qb(4)**2
	b2 = qb(5)**2 + qb(6)**2 + qb(7)**2
	pg = qb(8)

	qb5 = qa(5)+bs(1)
	qb6 = qa(6)+bs(2)
	qb7 = qa(7)+bs(3)

	bp = qb5**2 + qb6**2 + qb7**2
	bp0 = bs(1)**2 + bs(2)**2 + bs(3)**2

	pta = pg + bp/2.-bp0/2. 

	sam = ua - cfa
	sap = ua + cfa

	slu = sa - ua
	slm = sa - sm
	psb = qb(1)*slu*slm - qb5**2

	if(abs(psb).le.1.e-12) then

	   qz(1) = qb(1)
	   qz(2) = ua
	   qz(3) = qb(3)
	   qz(4) = qb(4)
	   qz(5) = qb5
	   qz(6) = 0.
	   qz(7) = 0.

	   vb = ua*qb5
	   vbz= sm*qb5

	   qz(2) = qz(1)*qz(2)
	   qz(3) = qz(1)*qz(3)
	   qz(4) = qz(1)*qz(4)
c	qz(8)=(slu*qa(8) - pta*ua + ptz*sm + qa(5)*(vb - vbz))/(sa - sm)
	   qz(8) = qa(8)

	else
 
	   rlz = qb(1)*slu/slm
	   rsm = qb(1)*slu**2 - qb5**2
	   smb = (sm - ua)/psb
	   smr = rsm/psb

	   rst = qb(1)*slu*smb

	   qz(1) = rlz
	   qz(2) = sm
	   qz(3) = qb(3)-qb5*qb6*smb
	   qz(4) = qb(4)-qb5*qb7*smb
	   qz(5) = qb(5)
	   qz(6) = qb(6)*smr + bs(2)*rst
	   qz(7) = qb(7)*smr + bs(3)*rst

	   vba = qb(2)*qb(5) + qb(3)*qb(6) + qb(4)*qb(7)
	   vbz = qz(2)*qz(5) + qz(3)*qz(6) + qz(4)*qz(7)

	   qz(2) = qz(1)*qz(2)
	   qz(3) = qz(1)*qz(3)
	   qz(4) = qz(1)*qz(4)

	   qz(8) = (slu*qa(8) - pta*ua + ptz*sm + qb5*(vba - vbz))/slm

	endif

	return
	end

	Subroutine quzzx_ig(nt,gamma,sm,ql,qr,bs,q)
	Implicit NONE

        integer nt
	ENZO_REAL ql(8),qr(8),q(8),qlm(8),qrm(8),bs(3)
        ENZO_REAL sl,sr,xs,xt,r1,r5,vm,r8,sa,sn,tm,vz,sm,gamma
	ENZO_REAL o1

	o1 = 1.0d0

        call prim_ig(gamma,ql,qlm) 
        call prim_ig(gamma,qr,qrm) 

	sl = sqrt(qlm(1))
	sr = sqrt(qrm(1))
	xs = sl+sr
	xt = sqrt(qlm(1)*qrm(1))

	if(nt.eq.1)then

	   r1 = qlm(1)
	   r5 = qlm(5)
	   vm = qlm(2)*qlm(5) + qlm(3)*qlm(6) + qlm(4)*qlm(7)
	   r8 = ql(8)
	   sa = sl
	   sn = sign(o1,r5+bs(1)) 
	   tm = -1.

	else

	   r1 = qrm(1)
	   r5 = qrm(5)
	   vm = qrm(2)*qrm(5) + qrm(3)*qrm(6) + qrm(4)*qrm(7)
	   r8 = qr(8)
	   sa = sr
	   sn = sign(o1,r5+bs(1)) 
	   tm = 1.

	endif

	q(1) = r1
	q(2) = sm
	q(3) = (sl*qlm(3) + sr*qrm(3) + (qr(6) - ql(6))*sn)/xs
	q(4) = (sl*qlm(4) + sr*qrm(4) + (qr(7) - ql(7))*sn)/xs
	q(5) = r5
	q(6) = (sl*qr(6) + sr*ql(6) + (qrm(3) - qlm(3))*sn*xt)/xs
	q(7) = (sl*qr(7) + sr*ql(7) + (qrm(4) - qlm(4))*sn*xt)/xs

	vz = q(2)*q(5) + q(3)*q(6) + q(4)*q(7)
	
	q(2) = q(1)*q(2)
	q(3) = q(1)*q(3)
	q(4) = q(1)*q(4)
	
	q(8) = r8 + tm*sa*(vm-vz)*sn

	return
	end

	Subroutine pressx_ig(gamma,ql,bs,ul,pt,cfl)
	Implicit NONE

        ENZO_REAL ql(8),bs(3)  
	ENZO_REAL gamma,ul,ck,cb,ds,cfl,qb5,qb6,qb7,bp,bp0,vk,bm,pg,pt

	qb5 = ql(5)+bs(1)
	qb6 = ql(6)+bs(2)
	qb7 = ql(7)+bs(3)

        bp = qb5**2 + qb6**2 + qb7**2
        bp0 = bs(1)**2 + bs(2)**2 + bs(3)**2
	 
	vk = ql(2)**2 + ql(3)**2 + ql(4)**2
	bm = ql(5)**2 + ql(6)**2 + ql(7)**2
	pg = (gamma-1.)*(ql(8) - vk/ql(1)/2. - bm/2.)

	pt = pg + bp/2.-bp0/2.

	ul = ql(2)/ql(1)

	ck = gamma*pg
	cb = ck + bp
	ds = cb*cb - 4.*ck*qb5**2
	if(ds.le.0.) ds=0.
	cfl= sqrt((cb + sqrt(ds))/2./ql(1))

	return
	end

	Subroutine HLLDy_ig(gamma,pl,pr,bs,f,us)
	Implicit NONE

	Integer m
	ENZO_REAL ql(8),qr(8),f(8),qlz(8),qrz(8),qzz(8),fl(8),fr(8),us(8)
        ENZO_REAL ul,ur,cfl,cfr,ptl,ptr,sl,sr,slu,sru,smu,sm,ptz,gamma
	ENZO_REAL rlz,slz,rrz,srz
        ENZO_REAL pl(8),pr(8),bs(3)
 
        call cons_ig(gamma,pl,ql)
        call cons_ig(gamma,pr,qr)
	
	call pressy_ig(gamma,ql,bs,ul,ptl,cfl)
	call pressy_ig(gamma,qr,bs,ur,ptr,cfr)

	sl = min(ul,ur) - max(cfl,cfr)
	sr = max(ul,ur) + max(cfl,cfr)

	if(sl.gt.0.)then

	   call fluxy_ig(gamma,pl,bs,f)

	   do m=1,8
	      us(m) = ql(m)
	   enddo

	   return

	endif

	if(sr.lt.0.)then

	   call fluxy_ig(gamma,pr,bs,f)

	   do m=1,8
	      us(m) = qr(m)
	   enddo

	   return

	endif

	slu = sl - ul
	sru = sr - ur
	smu = sru*qr(1)-slu*ql(1)
	sm  = (sru*qr(1)*ur - slu*ql(1)*ul - ptr + ptl)/smu

	ptz = (sru*qr(1)*ptl - slu*ql(1)*ptr + ql(1)*qr(1)*sru*slu
     &        *(ur - ul))/smu


	rlz = ql(1)*slu/(sl - sm)
	slz = sm - abs(ql(6)+bs(2))/sqrt(rlz)

	if(sl.le.0.and.slz.ge.0.)then
  
	   call fluxy_ig(gamma,pl,bs,fl)
	   call quzy_ig(gamma,sm,sl,cfl,ptz,ql,bs,qlz)
	
	   do m=1,8
	      f(m) = fl(m) + sl*(qlz(m) - ql(m))
	      us(m) = qlz(m)
	   enddo
	    
	   return

	endif


	rrz = qr(1)*sru/(sr - sm)
	srz = sm + abs(qr(6)+bs(2))/sqrt(rrz)

	if(srz.le.0.and.sr.ge.0.)then

	   call fluxy_ig(gamma,pr,bs,fr)
	   call quzy_ig(gamma,sm,sr,cfr,ptz,qr,bs,qrz)
	    
	   do m=1,8
	      f(m) = fr(m) + sr*(qrz(m) - qr(m))
	      us(m) = qrz(m)
	   enddo
	   
	   return

	endif

	call quzy_ig(gamma,sm,sl,cfl,ptz,ql,bs,qlz)
	call quzy_ig(gamma,sm,sr,cfr,ptz,qr,bs,qrz)

	if(slz.le.0.and.sm.ge.0.)then

	   call fluxy_ig(gamma,pl,bs,fl)

	   call quzzy_ig(1,gamma,sm,qlz,qrz,bs,qzz)

	   do m=1,8
	      f(m) = fl(m) + slz*qzz(m) - (slz - sl)*qlz(m) - sl*ql(m)
	      us(m) = qzz(m)
	   enddo

	   return

	endif

	if(sm.le.0.and.srz.ge.0.)then

	   call fluxy_ig(gamma,pr,bs,fr)

	   call quzzy_ig(2,gamma,sm,qlz,qrz,bs,qzz)

	   do m=1,8
	      f(m) = fr(m) + srz*qzz(m) - (srz - sr)*qrz(m) - sr*qr(m)
	      us(m) = qzz(m)
	   enddo

	   return

	endif

	return
	end


	Subroutine quzy_ig(gamma,sm,sa,cfa,ptz,qa,bs,qz)
	Implicit NONE

	ENZO_REAL qa(8),qz(8),qb(8),bs(3)
	ENZO_REAL gamma,ua,v2,b2,pg,pta,sam,sap,slu,slm,psb,vb,vbz,rlz
	ENZO_REAL rsm,smb,smr,vba,eps,sm,sa,cfa,ptz
	ENZO_REAL qb5,qb6,qb7,bp,bp0,rst

        call prim_ig(gamma,qa,qb)

	ua = qb(3)
	v2 = qb(2)**2 + qb(3)**2 + qb(4)**2
	b2 = qb(5)**2 + qb(6)**2 + qb(7)**2
	pg = qb(8)

	qb5 = qa(5)+bs(1)
	qb6 = qa(6)+bs(2)
	qb7 = qa(7)+bs(3)

        bp = qb5**2 + qb6**2 + qb7**2
        bp0 = bs(1)**2 + bs(2)**2 + bs(3)**2
	 
	pta = pg + bp/2.-bp0/2.

	sam = ua - cfa
	sap = ua + cfa

	slu = sa - ua
	slm = sa - sm
	psb = qb(1)*slu*slm - qb6**2

	if(abs(psb).le.1.e-12) then

	   qz(1) = qb(1)
	   qz(2) = qb(2)
	   qz(3) = ua
	   qz(4) = qb(4)
	   qz(5) = 0.
	   qz(6) = qb6
	   qz(7) = 0.

	   vb = ua*qb6
	   vbz= sm*qb6

	   qz(2) = qz(1)*qz(2)
	   qz(3) = qz(1)*qz(3)
	   qz(4) = qz(1)*qz(4)
c	qz(8)=(slu*qa(8) - pta*ua + ptz*sm + qa(6)*(vb - vbz))/(sa - sm)
	   qz(8) = qa(8)

	else
 
	   rlz = qb(1)*slu/slm
	   rsm = qb(1)*slu**2 - qb6**2
	   smb = (sm - ua)/psb
	   smr = rsm/psb

	   rst = qb(1)*slu*smb

	   qz(1) = rlz
	   qz(2) = qb(2)-qb6*qb5*smb
	   qz(3) = sm
	   qz(4) = qb(4)-qb6*qb7*smb
	   qz(5) = qb(5)*smr + bs(1)*rst
	   qz(6) = qb(6)
	   qz(7) = qb(7)*smr + bs(3)*rst

	   vba = qb(2)*qb(5) + qb(3)*qb(6) + qb(4)*qb(7)
	   vbz = qz(2)*qz(5) + qz(3)*qz(6) + qz(4)*qz(7)

	   qz(2) = qz(1)*qz(2)
	   qz(3) = qz(1)*qz(3)
	   qz(4) = qz(1)*qz(4)

	   qz(8) = (slu*qa(8) - pta*ua + ptz*sm + qb6*(vba - vbz))/slm

	endif

	return
	end

	Subroutine quzzy_ig(nt,gamma,sm,ql,qr,bs,q)
	Implicit NONE

        integer nt
	ENZO_REAL ql(8),qr(8),q(8),qlm(8),qrm(8),bs(3)
        ENZO_REAL sl,sr,xs,xt,r1,r6,vm,r8,sa,sn,tm,vz,sm,gamma
	ENZO_REAL o1

	o1 = 1.0d0

        call prim_ig(gamma,ql,qlm) 
        call prim_ig(gamma,qr,qrm) 

	sl = sqrt(qlm(1))
	sr = sqrt(qrm(1))
	xs = sl+sr
	xt = sqrt(qlm(1)*qrm(1))

	if(nt.eq.1)then

	r1 = qlm(1)
	r6 = qlm(6)
	vm = qlm(2)*qlm(5) + qlm(3)*qlm(6) + qlm(4)*qlm(7)
	r8 = ql(8)
	sa = sl
	sn = sign(o1,r6+bs(2)) 
	tm = -1.

	else

	r1 = qrm(1)
	r6 = qrm(6)
	vm = qrm(2)*qrm(5) + qrm(3)*qrm(6) + qrm(4)*qrm(7)
	r8 = qr(8)
	sa = sr
	sn = sign(o1,r6+bs(2)) 
	tm = 1.

	endif

	q(1) = r1
	q(2) = (sl*qlm(2) + sr*qrm(2) + (qr(5) - ql(5))*sn)/xs
	q(3) = sm
	q(4) = (sl*qlm(4) + sr*qrm(4) + (qr(7) - ql(7))*sn)/xs
	q(5) = (sl*qr(5) + sr*ql(5) + (qrm(2) - qlm(2))*sn*xt)/xs
	q(6) = r6
	q(7) = (sl*qr(7) + sr*ql(7) + (qrm(4) - qlm(4))*sn*xt)/xs

	vz = q(2)*q(5) + q(3)*q(6) + q(4)*q(7)
	
	q(2) = q(1)*q(2)
	q(3) = q(1)*q(3)
	q(4) = q(1)*q(4)
	
	q(8) = r8 + tm*sa*(vm-vz)*sn

	return
	end

	Subroutine pressy_ig(gamma,ql,bs,ul,pt,cfl)
	Implicit NONE

	ENZO_REAL ql(8),bs(3)
	ENZO_REAL gamma,ul,ck,cb,ds,cfl,qb5,qb6,qb7,bp,bp0,vk,bm,pg,pt
	
	
	qb5 = ql(5)+bs(1)
	qb6 = ql(6)+bs(2)
	qb7 = ql(7)+bs(3)

        bp = qb5**2 + qb6**2 + qb7**2
        bp0 = bs(1)**2 + bs(2)**2 + bs(3)**2
	 
	vk = ql(2)**2 + ql(3)**2 + ql(4)**2
	bm = ql(5)**2 + ql(6)**2 + ql(7)**2
	pg = (gamma-1.)*(ql(8) - vk/ql(1)/2. - bm/2.)

	pt = pg + bp/2.-bp0/2.
	 
	ul = ql(3)/ql(1)

	ck = gamma*pg
	cb = ck + bp
	ds = cb*cb - 4.*ck*qb6**2
	if(ds.le.0.) ds=0.
	cfl= sqrt((cb + sqrt(ds))/2./ql(1))
	
	return
	end


	Subroutine HLLDz_ig(gamma,pl,pr,bs,f,us)
	Implicit NONE

	Integer m
	ENZO_REAL ql(8),qr(8),f(8),qlz(8),qrz(8),qzz(8),fl(8),fr(8),us(8)
	ENZO_REAL ul,ur,cfl,cfr,ptl,ptr,sl,sr,slu,sru,smu,sm,ptz
	ENZO_REAL rlz,slz,rrz,srz,gamma
	ENZO_REAL pl(8),pr(8),bs(3)
 
        call cons_ig(gamma,pl,ql)
        call cons_ig(gamma,pr,qr)
	   
	call pressz_ig(gamma,ql,bs,ul,ptl,cfl)
	call pressz_ig(gamma,qr,bs,ur,ptr,cfr)

	sl = min(ul,ur) - max(cfl,cfr)
	sr = max(ul,ur) + max(cfl,cfr)

	if(sl.gt.0.)then

	   call fluxz_ig(gamma,pl,bs,f)

	   do m=1,8
	      us(m) = ql(m)
	   enddo

	   return

	endif

	if(sr.lt.0.)then

	   call fluxz_ig(gamma,pr,bs,f)

	   do m=1,8
	      us(m) = qr(m)
	   enddo

	   return

	endif

	slu = sl - ul
	sru = sr - ur
	smu = sru*qr(1)-slu*ql(1)
	sm  = (sru*qr(1)*ur - slu*ql(1)*ul - ptr + ptl)/smu

	ptz = (sru*qr(1)*ptl - slu*ql(1)*ptr + ql(1)*qr(1)*sru*slu
     &        *(ur - ul))/smu


	rlz = ql(1)*slu/(sl - sm)
	slz = sm - abs(ql(7)+bs(3))/sqrt(rlz)

	if(sl.le.0.and.slz.ge.0.)then
  
	   call fluxz_ig(gamma,pl,bs,fl)
	   call quz_ig(gamma,sm,sl,cfl,ptz,ql,bs,qlz)
	
	   do m=1,8
	      f(m) = fl(m) + sl*(qlz(m) - ql(m))
	      us(m) = qlz(m)
	   enddo
	    
	   return
	   
	endif


	rrz = qr(1)*sru/(sr - sm)
	srz = sm + abs(qr(7)+bs(3))/sqrt(rrz)

	if(srz.le.0.and.sr.ge.0.)then

	   call fluxz_ig(gamma,pr,bs,fr)
	   call quz_ig(gamma,sm,sr,cfr,ptz,qr,bs,qrz)
	    
	   do m=1,8
	      f(m) = fr(m) + sr*(qrz(m) - qr(m))
	      us(m) = qrz(m)
	   enddo
	
	   return

	endif

	call quz_ig(gamma,sm,sl,cfl,ptz,ql,bs,qlz)
	call quz_ig(gamma,sm,sr,cfr,ptz,qr,bs,qrz)

	if(slz.le.0.and.sm.ge.0.)then

	   call fluxz_ig(gamma,pl,bs,fl)

	   call quzz_ig(1,gamma,sm,qlz,qrz,bs,qzz)

	   do m=1,8
	      f(m) = fl(m) + slz*qzz(m) - (slz - sl)*qlz(m) - sl*ql(m)
	      us(m) = qzz(m)
	   enddo

	   return

	endif

	if(sm.le.0.and.srz.ge.0.)then

	   call fluxz_ig(gamma,pr,bs,fr)

	   call quzz_ig(2,gamma,sm,qlz,qrz,bs,qzz)

	   do m=1,8
	      f(m) = fr(m) + srz*qzz(m) - (srz - sr)*qrz(m) - sr*qr(m)
	      us(m) = qzz(m)
	   enddo

	   return

	endif

	return
	end


	Subroutine quz_ig(gamma,sm,sa,cfa,ptz,qa,bs,qz)
	Implicit NONE

	ENZO_REAL qa(8),qz(8),qb(8),bs(3)
        ENZO_REAL gamma,ua,v2,b2,pg,sam,sap,slu,slm,psb,vb,vbz,rlz
	ENZO_REAL rsm,smb,smr,vba,eps,sm,sa,cfa,ptz
	ENZO_REAL qb5,qb6,qb7,bp,bp0,pta,rst

        call prim_ig(gamma,qa,qb)

	ua = qb(4)
	v2 = qb(2)**2 + qb(3)**2 + qb(4)**2
	b2 = qb(5)**2 + qb(6)**2 + qb(7)**2
	pg = qb(8)

	qb5 = qa(5)+bs(1)
	qb6 = qa(6)+bs(2)
	qb7 = qa(7)+bs(3)

        bp = qb5**2 + qb6**2 + qb7**2
        bp0 = bs(1)**2 + bs(2)**2 + bs(3)**2
	 
	pta = pg + bp/2.-bp0/2.

	sam = ua - cfa
	sap = ua + cfa

	slu = sa - ua
	slm = sa - sm
	psb = qb(1)*slu*slm - qb7**2

	if(abs(psb).le.1.e-12) then

	   qz(1) = qb(1)
	   qz(2) = qb(2)
	   qz(3) = qb(3)
	   qz(4) = ua
	   qz(5) = 0.
	   qz(6) = 0.
	   qz(7) = qb7

	   vb = ua*qb7
	   vbz= sm*qb7

	   qz(2) = qz(1)*qz(2)
	   qz(3) = qz(1)*qz(3)
	   qz(4) = qz(1)*qz(4)
c	qz(8)=(slu*qa(8) - pta*ua + ptz*sm + qa(7)*(vb - vbz))/(sa - sm)
	   qz(8) = qa(8)

	else
 
	   rlz = qb(1)*slu/slm
	   rsm = qb(1)*slu**2 - qb7**2
	   smb = (sm - ua)/psb
	   smr = rsm/psb
	   
	   rst = qb(1)*slu*smb
	
	   qz(1) = rlz
	   qz(2) = qb(2) - qb5*qb7*smb
	   qz(3) = qb(3) - qb6*qb7*smb
	   qz(4) = sm
	   qz(5) = qb(5)*smr + bs(1)*rst
	   qz(6) = qb(6)*smr + bs(2)*rst
	   qz(7) = qb(7)

	   vba = qb(2)*qb(5) + qb(3)*qb(6) + qb(4)*qb(7)
	   vbz = qz(2)*qz(5) + qz(3)*qz(6) + qz(4)*qz(7)

	   qz(2) = qz(1)*qz(2)
	   qz(3) = qz(1)*qz(3)
	   qz(4) = qz(1)*qz(4)
	
	   qz(8) = (slu*qa(8) - pta*ua + ptz*sm + qb7*(vba - vbz))/slm

	endif

	return
	end

	Subroutine quzz_ig(nt,gamma,sm,ql,qr,bs,q)
	Implicit NONE

	integer nt
	ENZO_REAL ql(8),qr(8),q(8),qlm(8),qrm(8),bs(3)
	ENZO_REAL sl,sr,xs,xt,r1,r7,vm,r8,sa,sn,tm,vz,sm,gamma
	ENZO_REAL o1

	o1 = 1.0d0

        call prim_ig(gamma,ql,qlm) 
        call prim_ig(gamma,qr,qrm) 

	sl = sqrt(qlm(1))
	sr = sqrt(qrm(1))
	xs = sl+sr
	xt = sqrt(qlm(1)*qrm(1))

	if(nt.eq.1)then

	   r1 = qlm(1)
	   r7 = qlm(7)
	   vm = qlm(2)*qlm(5) + qlm(3)*qlm(6) + qlm(4)*qlm(7)
	   r8 = ql(8)
	   sa = sl
	   sn = sign(o1,r7+bs(3)) 
	   tm = -1.

	else

	   r1 = qrm(1)
	   r7 = qrm(7)
	   vm = qrm(2)*qrm(5) + qrm(3)*qrm(6) + qrm(4)*qrm(7)
	   r8 = qr(8)
	   sa = sr
	   sn = sign(o1,r7+bs(3)) 
	   tm = 1.

	endif

	q(1) = r1
	q(2) = (sl*qlm(2) + sr*qrm(2) + (qr(5) - ql(5))*sn)/xs
	q(3) = (sl*qlm(3) + sr*qrm(3) + (qr(6) - ql(6))*sn)/xs
	q(4) = sm
	q(5) = (sl*qr(5) + sr*ql(5) + (qrm(2) - qlm(2))*sn*xt)/xs
	q(6) = (sl*qr(6) + sr*ql(6) + (qrm(3) - qlm(3))*sn*xt)/xs
	q(7) = r7

	vz = q(2)*q(5) + q(3)*q(6) + q(4)*q(7)
	
	q(2) = q(1)*q(2)
	q(3) = q(1)*q(3)
	q(4) = q(1)*q(4)
	
	q(8) = r8 + tm*sa*(vm-vz)*sn

	return
	end

	Subroutine pressz_ig(gamma,ql,bs,ul,pt,cfl)
	Implicit NONE

	ENZO_REAL ql(8),bs(3)
	ENZO_REAL gamma,ul,ck,cb,ds,cfl,qb5,qb6,qb7,bp,bp0,vk,bm,pg,pt

	qb5 = ql(5)+bs(1)
	qb6 = ql(6)+bs(2)
	qb7 = ql(7)+bs(3)

        bp = qb5**2 + qb6**2 + qb7**2
        bp0 = bs(1)**2 + bs(2)**2 + bs(3)**2
	 
	vk = ql(2)**2 + ql(3)**2 + ql(4)**2
	bm = ql(5)**2 + ql(6)**2 + ql(7)**2
	pg = (gamma-1.)*(ql(8) - vk/ql(1)/2. - bm/2.)

	pt = pg + bp/2.-bp0/2.
	 
	ul = ql(4)/ql(1)

	ck = gamma*pg
	cb = ck + bp
	ds = cb*cb - 4.*ck*qb7**2
	if(ds.le.0.) ds=0.
	cfl= sqrt((cb + sqrt(ds))/2./ql(1))

	return
	end

