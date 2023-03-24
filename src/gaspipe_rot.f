!
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
!     
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!     
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!     
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!     
      subroutine gaspipe_rot(node1,node2,nodem,nelem,lakon,kon,
     &        ipkon,nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,
     &        shcon,nshcon,rhcon,nrhcon,ntmat_,co,vold,mi,ttime,time,
     &        iaxial,iplausi)
!     
!     rotating pipe with friction and variable cross section
!     
      implicit none
!     
      logical identity,crit
      character*8 lakon(*)
      character*81 set(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,ielprop(*),
     &     idirf(*),index,kflag,nodef(*),inv,ipkon(*),kon(*),icase,ith,
     &     nshcon(*),nrhcon(*),ntmat_,mi(*),iaxial,iplausi
!
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(*),kappa,r,A,d,xl,
     &     T1,T2,Tt1,Tt2,pt1,pt2,cp,physcon(*),km1,dvi,kp1,kdkm1,
     &     reynolds,pi,lambda,kdkp1,rho,km1d2,ttime,time,pt2zpt1,ks,
     &     form_fact,pt2zpt1_c,Qred1_crit,Qred,phi,M1,M2,Qred1,co(3,*),
     &     shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*),vold(0:mi(2),*),
     &     bb,cc,ee1,ee2,dfdM1,dfdM2,M1_c,Qred2,za,zb,zc,Z1,Z2,r1,r2,
     &     term,omega,om2,d1,d2,A1,A2,alpha,beta,coef,M2_c,Qred2_crit
!
      ith=0
!
      if(kflag.eq.0) then
         identity=.true.
!     
         if(nactdog(2,node1).ne.0)then
            identity=.false.
         elseif(nactdog(2,node2).ne.0)then
            identity=.false.
         elseif(nactdog(1,nodem).ne.0)then
            identity=.false.
         endif
!     
      elseif(kflag.eq.1)then
         if(v(1,nodem).ne.0.d0) then
            xflow=v(1,nodem)
            return
         endif
!
         pi=4.d0*datan(1.d0)
!     
         index=ielprop(nelem)
         kappa=(cp/(cp-r))
!         
         A1=prop(index+1)
         A2=prop(index+2)
         xl=prop(index+3)
         ks=prop(index+4)
         form_fact=prop(index+5)
         d1=prop(index+6)
         if(form_fact.eq.1.d0) then
            d1=2.d0*dsqrt(A1/pi)
         endif
         d2=prop(index+7)
         if(form_fact.eq.1.d0) then
            d2=2.d0*dsqrt(A2/pi)
         endif
         r1=prop(index+8)
         r2=prop(index+9)
         omega=prop(index+10)
!
         pt1=v(2,node1)
         pt2=v(2,node2)
!
         Tt1=v(0,node1)-physcon(1)
         Tt2=v(0,node2)-physcon(1)
!
         A=(A1+A2)/2.d0
         d=(d1+d2)/2.d0
         if(r2.ge.r1) then
            om2=omega**2
         else
            om2=-omega**2
         endif
!     
!        calculation of the dynamic viscosity
!     
         if(dabs(dvi).lt.1d-30) then
            write(*,*) '*ERROR in gaspipe_fanno: '
            write(*,*) '       no dynamic viscosity defined'
            write(*,*) '       dvi= ',dvi
            call exit(201)
         endif        
!
!        assumed value for the reynolds number
!         
         reynolds=3.e3
!
         call friction_coefficient(xl,d,ks,reynolds,form_fact,
     &        lambda)
!
!        estimate of the flow using the incompressible relationships 
!        for a gas pipe
!
!        mean density
!
         rho=(pt1/(r*Tt1)+pt2/(r*Tt2))/2.d0
         term=(rho*(pt1-pt2)+rho**2*om2*(r2**2-r1**2)/2.d0)*
     &           2.d0*d/(lambda*xl)
!         
         if(term.ge.0.d0) then
            inv=1.d0
            if(v(1,nodem).le.0.d0) then
               xflow=A*dsqrt(term)
            else
               xflow=v(1,nodem)*iaxial
            endif
         else
!
!           if the term underneath the square root is negative,
!           lambda must have a negative sign, which means that the flow
!           direction has to be reversed
!            
            inv=-1.d0
            lambda=-lambda
            if(v(1,nodem).ge.0.d0) then
               xflow=-A*dsqrt(-term)
            else
               xflow=v(1,nodem)*iaxial
            endif
         endif
!
         kp1=kappa+1.d0
         km1=kappa-1.d0
         km1d2=km1/2.d0
!
!        check whether the flow does not exceed the critical one
!
         alpha=-4.d0*(d2-d1)/(xl*d)-(r1+r2)*om2*kp1/
     &        (cp*(Tt1+Tt2)*km1)
         if(alpha.eq.0.d0) then
            write(*,*) '*ERROR in gaspipe_rot: looks like the'
            write(*,*) '       cross section is constant and'
            write(*,*) '       the rotational speed is zero;'
            write(*,*) '       please use the gaspipe_fanno'
            write(*,*) '       element instead'
            call exit(201)
         endif
         beta=lambda*kappa/d
!
!        for subsonic flow:         
!        coef>0.d0 means that the Mach number increases from 1 to 2
!           (i.e. sonic conditions can only occur at 2)
!        coef<0.d0 means that the Mach number decreases from 1 to 2
!           (i.e. sonic conditions can only occur at 1)
!
         call ts_calc(xflow,Tt1,pt1,kappa,r,A1,T1,ith)
         M1=dsqrt(((Tt1/T1)-1.d0)/km1d2)
         call ts_calc(xflow,Tt2,pt2,kappa,r,A2,T2,ith)
         M2=dsqrt(((Tt2/T2)-1.d0)/km1d2)
         coef=alpha+beta*((M1+M2)/2.d0)**2
!
!        icase tells where sonic conditions can occur, if any.
!
         if(coef.ge.0.d0) then
            icase=2
         else
            icase=1
         endif
!
         za=1.d0/alpha
         zb=(alpha+beta)*beta/(alpha*(alpha*km1d2-beta))
         zc=kp1*km1d2/(2.d0*beta-alpha*km1)         
!
         if(omega.eq.0.d0) then
            call pt2zpt1_rot(pt2,pt1,kappa,r,xl,pt2zpt1_c,crit,icase,
     &        M1_c,M2_c,za,zb,zc,alpha,beta,Qred1_crit,Qred2_crit,A1,A2)
         else
            crit=.false.
         endif
!
!        check for critical flow only in the absence of rotational
!        speed
!         
         if((icase.eq.1).and.(omega.eq.0.d0)) then
c         if(icase.eq.1) then
!
!           decreasing Mach number from 1 to 2
!
            Qred2=dabs(xflow)*dsqrt(Tt2)/(A2*pt2)
!
!           check whether flow is critical
!           assigning the physcical correct sign to xflow
!
            if(crit) then
!     
!              the flow is set to half the critical value or
!              one of the pressures is adapted (depending on
!              which variable is unknown)
!
                  xflow=0.5d0*inv*Qred2_crit*pt2*A2/dsqrt(Tt2)
            elseif(Qred.gt.Qred2_crit) then
!     
!              the flow is set to half the critical value
!     
               xflow=0.5d0*inv*Qred2_crit*pt2*A2/dsqrt(Tt2)
            endif
         else
!
!           increasing Mach number from 1 to 2
! 
            Qred=dabs(xflow)*dsqrt(Tt1)/(A1*pt1)
            if(crit) then
!     
!              the flow is set to half the critical value or
!              one of the pressures is adapted (depending on
!              which variable is unknown)
!     
                  xflow=0.5d0*inv*Qred1_crit*pt1*A1/dsqrt(Tt1)
             elseif(Qred.gt.Qred1_crit) then
!     
!              the flow is set to half the critical value
!     
               xflow=0.5d0*inv*Qred1_crit*pt1*A1/dsqrt(Tt1)
            endif
         endif
!
      elseif(kflag.eq.2)then
!
         numf=5
!
         pi=4.d0*datan(1.d0)
!
         index=ielprop(nelem)
         kappa=(cp/(cp-r))
!
         A1=prop(index+1)
         A2=prop(index+2)
         xl=prop(index+3)
         ks=prop(index+4)
         form_fact=prop(index+5)
         d1=prop(index+6)
         if(form_fact.eq.1.d0) then
            d1=2.d0*dsqrt(A1/pi)
         endif
         d2=prop(index+7)
         if(form_fact.eq.1.d0) then
            d2=2.d0*dsqrt(A2/pi)
         endif
         r1=prop(index+8)
         r2=prop(index+9)
         omega=prop(index+10)
!
         pt1=v(2,node1)
         pt2=v(2,node2)
!
         Tt1=v(0,node1)-physcon(1)
         Tt2=v(0,node2)-physcon(1)
!
         xflow=v(1,nodem)*iaxial
!
         write(*,*) 'gaspipe_rot pt..: ',pt1,pt2,xflow,Tt1,Tt2
!
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=2
         idirf(5)=0
!
         nodef(1)=node1
         nodef(2)=node1
         nodef(3)=nodem
         nodef(4)=node2
         nodef(5)=node2
!     
         pt2zpt1=pt2/pt1
!         
         A=(A1+A2)/2.d0
         d=(d1+d2)/2.d0
         if(r2.ge.r1) then
            om2=omega**2
         else
            om2=-omega**2
         endif
!     
!     calculation of the dynamic viscosity
!     
         if(dabs(dvi).lt.1d-30) then
            write(*,*) '*ERROR in gaspipe_fanno: '
            write(*,*) '       no dynamic viscosity defined'
            write(*,*) '       dvi= ',dvi
            call exit(201)
         endif        
!     
         reynolds=dabs(xflow)*d/(dvi*A)
!     
!        calculation of the friction coefficient
!     
         call friction_coefficient(xl,d,ks,reynolds,form_fact,
     &        lambda)
!
         if(xflow.lt.0.d0) then
            lambda=-lambda
            inv=-1
         else
            inv=1
         endif
!
         kp1=kappa+1.d0
         km1=kappa-1.d0
         km1d2=km1/2.d0
!         
         alpha=-4.d0*(d2-d1)/(xl*d)-(r1+r2)*om2*kp1/
     &        (cp*(Tt1+Tt2)*km1)
         if(alpha.eq.0.d0) then
            write(*,*) '*ERROR in gaspipe_rot: looks like the'
            write(*,*) '       cross section is constant and'
            write(*,*) '       the rotational speed is zero;'
            write(*,*) '       please use the gaspipe_fanno'
            write(*,*) '       element instead'
            call exit(201)
         endif
         beta=lambda*kappa/d
!
!        for subsonic flow:         
!        coef>0.d0 means that the Mach number increases from 1 to 2
!           (i.e. sonic conditions can only occur at 2)
!        coef<0.d0 means that the Mach number decreases from 1 to 2
!           (i.e. sonic conditions can only occur at 1)
!
         call ts_calc(xflow,Tt1,pt1,kappa,r,A1,T1,ith)
         M1=dsqrt(((Tt1/T1)-1.d0)/km1d2)
         call ts_calc(xflow,Tt2,pt2,kappa,r,A2,T2,ith)
         M2=dsqrt(((Tt2/T2)-1.d0)/km1d2)
         coef=alpha+beta*((M1+M2)/2.d0)**2
         write(*,*) 'gaspipe_rot M1,M2',M1,M2
!
!        icase tells where sonic conditions can occur, if any.
!
         if(coef.ge.0.d0) then
            icase=2
         else
            icase=1
         endif
!
         za=1.d0/alpha
         zb=(alpha+beta)*beta/(alpha*(alpha*km1d2-beta))
         zc=kp1*km1d2/(2.d0*beta-alpha*km1)
         write(*,*) 'gaspipe_rot alpha beta',alpha,beta
         write(*,*) 'za zb zc',za,zb,zc
!
         if(omega.eq.0.d0) then
            call pt2zpt1_rot(pt2,pt1,kappa,r,xl,pt2zpt1_c,crit,icase,
     &        M1_c,M2_c,za,zb,zc,alpha,beta,Qred1_crit,Qred2_crit,A1,A2)
         else
            crit=.false.
         endif
         write(*,*) 'gaspipe_rot crit ',crit
!
!        check for critical flow only in the absence of rotational
!        speed
!         
         if((icase.eq.1).and.(omega.eq.0.d0)) then
!
!           decreasing Mach number from 1 to 2
!
            Qred2=dabs(xflow)*dsqrt(Tt2)/(A2*pt2)
!
!           check whether flow is critical
!           assigning the physcical correct sign to xflow
!
            if(crit) then
               write(*,*) '*WARNING in gaspipe_rot'
               write(*,*) '         critical conditions detected'
               write(*,*) '         in element ',nelem
               xflow=inv*Qred2_crit*A2*pt2/dsqrt(Tt2)
!
!              check whether flow has changed; if so, update v
!              for consistency
!
               if(dabs((xflow-iaxial*v(1,nodem))/xflow).gt.1.d-5) then
                  iplausi=0
c                  if(nactdog(1,nodem).ne.0) v(1,nodem)=xflow/iaxial
               endif
               M2=dsqrt(((Tt2/T2)-1.d0)/km1d2)
               M2=min(M2,0.999d0)
               if((alpha+beta*M2*M2)/(alpha+beta).lt.0.d0) then
                  M2=M2_c
               endif
!
!                 recalculate M1
!
               call ts_calc(xflow,Tt1,pt1,kappa,r,A1,T1,ith)
               M1=dsqrt(((Tt1/T1)-1.d0)/km1d2)
            else
               if(Qred2.gt.Qred2_crit) then
                  xflow=inv*Qred2_crit*A2*pt2/dsqrt(Tt2)
!
!                 check whether flow has changed; if so, update v
!                 for consistency
!
                  if(dabs((xflow-iaxial*v(1,nodem))/xflow).gt.1.d-5)
     &               then
                     iplausi=0
                     if(nactdog(1,nodem).ne.0) v(1,nodem)=xflow/iaxial
                  endif
!            
                  M2=M2_c
!
!                 recalculate M1
!
                  call ts_calc(xflow,Tt1,pt1,kappa,r,A1,T1,ith)
                  M1=dsqrt(((Tt1/T1)-1.d0)/km1d2)
c               else
c                  call ts_calc(xflow,Tt2,pt2,kappa,r,A2,T2,ith)
c                  M2=dsqrt(((Tt2/T2)-1.d0)/km1d2)
               endif
!     
c               Tt1=Tt2
c               call ts_calc(xflow,Tt1,pt1,kappa,r,A1,T1,ith)
c               M1=dsqrt(((Tt1/T1)-1.d0)/km1d2)
            endif
         elseif((icase.eq.2).and.(omega.eq.0.d0)) then
!
!           increasing Mach number from 1 to 2
!
            Qred1=dabs(xflow)*dsqrt(Tt1)/(A1*pt1)
!
!           check whether flow is critical
!           assigning the physcical correct sign to xflow
!
            if(crit) then
               write(*,*) '*WARNING in gaspipe_rot'
               write(*,*) '         critical conditions detected'
               write(*,*) '         in element ',nelem
               xflow=inv*Qred1_crit*A1*pt1/dsqrt(Tt1)
!
!              check whether flow has changed; if so, update v
!              for consistency
!
               if(dabs((xflow-iaxial*v(1,nodem))/xflow).gt.1.d-5) then
                  iplausi=0
c                  if(nactdog(1,nodem).ne.0) v(1,nodem)=xflow/iaxial
               endif
!
               M1=dsqrt(((Tt1/T1)-1.d0)/km1d2)
               M1=min(M1,0.999d0)
               if((alpha+beta)/(alpha+beta*M1*M1).lt.0.d0) then
                  M1=M1_c
               endif
!
!                 recalculate M2
!
               call ts_calc(xflow,Tt2,pt2,kappa,r,A2,T2,ith)
               M2=dsqrt(((Tt2/T2)-1.d0)/km1d2)
            else
               if(Qred1.gt.Qred1_crit) then
                  xflow=inv*Qred1_crit*A1*pt1/dsqrt(Tt1)
!
!                 check whether flow has changed; if so, update v
!                 for consistency
!
                  if(dabs((xflow-iaxial*v(1,nodem))/xflow).gt.1.d-5)
     &               then
                     iplausi=0
                     if(nactdog(1,nodem).ne.0) v(1,nodem)=xflow/iaxial
                  endif
!            
                  M1=M1_c
!
!                 recalculate M2
!
                  call ts_calc(xflow,Tt2,pt2,kappa,r,A2,T2,ith)
                  M2=dsqrt(((Tt2/T2)-1.d0)/km1d2)
c               else
c                  call ts_calc(xflow,Tt1,pt1,kappa,r,A1,T1,ith)
c                  M1=dsqrt(((Tt1/T1)-1.d0)/km1d2)
               endif
!     
c               Tt2=Tt1
               call ts_calc(xflow,Tt2,pt2,kappa,r,A2,T2,ith)
               M2=dsqrt(((Tt2/T2)-1.d0)/km1d2)
            endif
c         elseif(icase.eq.1) then
c         else
         endif
!
         bb=km1d2
         cc=-kp1/(2.d0*km1)
!     
!     definition of the coefficients
!
         if(icase.eq.1) then
!
!           decreasing Mach number from 1 to 2
!
            Z2=M2**2
            ee2=M2*(1.d0+bb*Z2)/(1.d0+bb*Z2*(1.d0+2.d0*cc))
            dfdM2=2.d0*M2*(za/Z2+zb/(alpha+beta*Z2)
     &                          +zc/(1.d0+km1d2*Z2))
!     
            if(.not.crit) then
!     
!              residual
!
               Z1=M1**2
               ee1=M1*(1.d0+bb*Z1)/(1.d0+bb*Z1*(1.d0+2.d0*cc))
               dfdM1=-2.d0*M1*(za/Z1+zb/(alpha+beta*Z1)
     &                              +zc/(1.d0+km1d2*Z1))
!
               f=dlog((Z2/Z1)**za
     &              *((alpha+beta*Z2)/(alpha+beta*Z1))**(zb/beta)
     &              *((1.d0+km1d2*Z2)/(1.d0+km1d2*Z1))**(zc/km1d2))
     &              -xl
!     
!              pressure node1
!     
               df(1)=-dfdM1*ee1/pt1
!     
!              temperature node1
!     
               df(2)=dfdM1*ee1/(2.d0*Tt1)
!     
!              mass flow
!     
               df(3)=(dfdM1*ee1+dfdM2*ee2)/xflow
!     
!              pressure node2
!     
               df(4)=-dfdM2*ee2/pt2
!     
!              temperature node2
!     
               df(5)=dfdM2*ee2/(2.d0*Tt2)
!               
            else
!
               f=dlog(Z2**za
     &              *((alpha+beta*Z2)/(alpha+beta))**(zb/beta)
     &              *((1.d0+km1d2*Z2)/(1.d0+km1d2))**(zc/km1d2))
     &              -xl
!     
!              pressure node1
!     
               df(1)=0.d0
!     
!              temperature node1
!     
               df(2)=0.d0
!     
!              mass flow
!     
               df(3)=(dfdM2*ee2)/xflow
!     
!              pressure node2
!     
               df(4)=-dfdM2*ee2/pt2
!     
!              temperature node2
!     
               df(5)=dfdM2*ee2/(2.d0*Tt2)
!
            endif
         elseif(icase.eq.2) then
!
!           increasing Mach number from 1 to 2
!
            Z1=M1**2
            ee1=M1*(1.d0+bb*Z1)/(1.d0+bb*Z1*(1.d0+2.d0*cc))
            dfdM1=-2.d0*M1*(za/Z1+zb/(alpha+beta*Z1)
     &                           +zc/(1.d0+km1d2*Z1))
!     
            if(.not.crit) then
!     
!              residual
!
               Z2=M2**2
               ee2=M2*(1.d0+bb*Z2)/(1.d0+bb*Z2*(1.d0+2.d0*cc))
               dfdM2=2.d0*M2*(za/Z2+zb/(alpha+beta*Z2)
     &                             +zc/(1.d0+km1d2*Z2))
!
               f=dlog((Z2/Z1)**za
     &              *((alpha+beta*Z2)/(alpha+beta*Z1))**(zb/beta)
     &              *((1.d0+km1d2*Z2)/(1.d0+km1d2*Z1))**(zc/km1d2))
     &              -xl
!     
!              pressure node1
!     
               df(1)=-dfdM1*ee1/pt1
!     
!              temperature node1
!     
               df(2)=dfdM1*ee1/(2.d0*Tt1)
!     
!              mass flow
!     
               df(3)=(dfdM1*ee1+dfdM2*ee2)/xflow
!     
!              pressure node2
!     
               df(4)=-dfdM2*ee2/pt2
!     
!              temperature node2
!     
               df(5)=dfdM2*ee2/(2.d0*Tt2)
!               
            else
!
               f=dlog((1.d0/Z1)**za
     &              *((alpha+beta)/(alpha+beta*Z1))**(zb/beta)
     &              *((1.d0+km1d2)/(1.d0+km1d2*Z1))**(zc/km1d2))
     &              -xl
!     
!              pressure node1
!     
               df(1)=-dfdM1*ee1/pt1
!     
!              temperature node1
!     
               df(2)=dfdM1*ee1/(2.d0*Tt1)
!     
!              mass flow
!     
               df(3)=dfdM1*ee1/xflow
!     
!              pressure node2
!     
               df(4)=0.d0
!     
!              temperature node2
!     
               df(5)=0.d0
!
            endif
         endif
!
!     output
!
      elseif(kflag.eq.3) then
!
         pi=4.d0*datan(1.d0)
!     
         kappa=(cp/(cp-r))
         km1=kappa-1.d0
         km1d2=km1/2.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         kdkp1=kappa/kp1
!         
         index=ielprop(nelem)
         A1=prop(index+1)
         A2=prop(index+2)    
         xl=prop(index+3)
!
         lambda=0.5d0
!
         ks=prop(index+4)
         form_fact=prop(index+5)
         d1=prop(index+6)
         if(form_fact.eq.1.d0) then
            d1=2.d0*dsqrt(A1/pi)
         endif
         d2=prop(index+7)
         if(form_fact.eq.1.d0) then
            d2=2.d0*dsqrt(A2/pi)
         endif
         r1=prop(index+8)
         r2=prop(index+9)
         omega=prop(index+10)
!
         pt1=v(2,node1)
         pt2=v(2,node2)
!     
         if(xflow.ge.0.d0) then
            inv=1
            xflow=v(1,nodem)*iaxial
            Tt1=v(0,node1)-physcon(1)
            Tt2=v(0,node2)-physcon(1)
         else
            inv=-1
            pt1=v(2,node2)
            pt2=v(2,node1)
            xflow=v(1,nodem)*iaxial
            Tt1=v(0,node2)-physcon(1)
            Tt2=v(0,node1)-physcon(1)
         endif
         call ts_calc(xflow,Tt1,pt1,kappa,r,A1,T1,ith)
c         Tt2=Tt1
         call ts_calc(xflow,Tt2,pt2,kappa,r,A2,T2,ith)
!     
         pt2zpt1=pt2/pt1
!         
         A=(A1+A2)/2.d0
         d=(d1+d2)/2.d0
!     
!     calculation of the dynamic viscosity 
!     
         if(dabs(dvi).lt.1d-30) then
            write(*,*) '*ERROR in gaspipe_fanno: '
            write(*,*) '       no dynamic viscosity defined'
            write(*,*) '       dvi= ',dvi
            call exit(201)
         endif
! 
         reynolds=dabs(xflow)*d/(dvi*A)
!
        if(reynolds.lt.1.d0) then
            reynolds= 1.d0
         endif
!     
!     definition of the friction coefficient for pure air
!     
         phi=1.d0
         call friction_coefficient(xl,d,ks,reynolds,form_fact,
     &        lambda)
!     
!     definition of the coefficients 
!     
         M1=dsqrt(((Tt1/T1)-1)/km1d2)
         M2=dsqrt(((Tt2/T2)-1)/km1d2)
!     
         write(1,*) ''
         write(1,55) ' from node ',node1,
     &        ' to node ', node2,':   air massflow rate = ',xflow
!     
         if(inv.eq.1) then
            write(1,53)'       Inlet node ',node1,' :    Tt1 = ',Tt1,
     &           ' , Ts1 = ',T1,'  , Pt1 = ',pt1,
     &           ' , M1 = ',M1
            write(1,*)'             Element ',nelem,lakon(nelem)
            write(1,57)'             dvi = ',dvi,' , Re = '
     &           ,reynolds
            write(1,58)'             PHI = ',phi,' , LAMBDA = ',
     &           lambda,
     &           ', LAMBDA*l/d = ',lambda*xl/d,' , ZETA_PHI = ',
     &           phi*lambda*xl/d
            write(1,53)'      Outlet node ',node2,' :    Tt2 = ',
     &           Tt2,
     &           ' , Ts2 = ',T2,'  , Pt2 = ',pt2,
     &           ' , M2 = ',M2
!     
         else if(inv.eq.-1) then
            write(1,53)'       Inlet node ',node2,':    Tt1= ',Tt1,
     &           ' , Ts1= ',T1,' , Pt1= ',pt1,
     &           ' , M1= ',M1
            write(1,*)'             Element ',nelem,lakon(nelem)
            write(1,57)'             dvi = ',dvi,' , Re = '
     &           ,reynolds
            write(1,58)'             PHI = ',phi,' , LAMBDA = ',
     &           lambda,
     &           ', LAMBDA*l/d = ',lambda*xl/d,' , ZETA_PHI = ',
     &           phi*lambda*xl/d
            write(1,53)'      Outlet node ',node1,' :    Tt2 = ',
     &           Tt2,
     &           ' , Ts2 = ',T2,'  , Pt2 =',pt2,
     &           ' , M2 = ',M2
         endif
      endif
!     
 55   format(1X,a,i6,a,i6,a,e11.4,a,e11.4)
 53   format(1X,a,i6,a,e11.4,a,e11.4,a,e11.4,a,e11.4)
 57   format(1X,a,e11.4,a,e11.4)
 58   format(1X,a,e11.4,a,e11.4,a,e11.4,a,e11.4)
!     
      xflow=xflow/iaxial
      df(3)=df(3)*iaxial
!     
      return
      end
