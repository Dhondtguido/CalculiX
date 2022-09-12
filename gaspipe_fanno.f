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
      subroutine gaspipe_fanno(node1,node2,nodem,nelem,lakon,kon,
     &        ipkon,nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,
     &        shcon,nshcon,rhcon,nrhcon,ntmat_,co,vold,mi,ttime,time,
     &        iaxial,iplausi)
!     
!     pipe with friction losses (Fanno Formulas) GAPF 
!
!     author: Yannick Muller
!     
      implicit none
!     
      logical identity,crit,wrongdir
      character*8 lakon(*)
      character*81 set(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(*),idirf(*),index,kflag,
     &     inv,ipkon(*),kon(*),icase,k_oil
     &     ,nshcon(*),nrhcon(*),ntmat_,mi(*),nodea,nodeb,
     &     nodec,iaxial,iplausi
!
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(*),kappa,r,A,d,l,
     &     T1,T2,Tt1,Tt2,pt1,pt2,cp,physcon(*),p2p1,km1,dvi,
     &     kp1,kdkm1,km1d2k,reynolds,pi,lambda,lld,xflow_air,
     &     C2,tdkp1,ttime,time,pt2zpt1,ks,form_fact,xflow_oil,
     &     pt2zpt1_c,Qred1_crit,Qred,phi,M1,M2,Qred1,co(3,*),
     &     shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*),vold(0:mi(2),*),
     &     radius,bb,cc,ee1,ee2,dfdM1,dfdM2,M1_c,Z1,Z2
!
!
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
!
         pi=4.d0*datan(1.d0)
!     
         index=ielprop(nelem)
         kappa=(cp/(cp-r))
         A=prop(index+1)
         d=prop(index+2)
         l=prop(index+3)
         ks=prop(index+4)
         if(lakon(nelem)(2:6).eq.'GAPFA') then
            icase=0
         elseif(lakon(nelem)(2:6).eq.'GAPFI') then
            icase=1
         endif
         form_fact=prop(index+5)
         xflow_oil=prop(index+6)
         k_oil=nint(prop(index+7))
!
         if(lakon(nelem)(7:8).eq.'FR') then
!
!           flexible radius
!            
            nodea=nint(prop(index+1))
            nodeb=nint(prop(index+2))
            radius=dsqrt((co(1,nodeb)+vold(1,nodeb)-
     &           co(1,nodea)-vold(1,nodea))**2)
!
            A=pi*radius**2
            d=2*radius
!
         elseif(lakon(nelem)(7:8).eq.'RL') then
!
!           flexible radius and length
!
            nodea=nint(prop(index+1))
            nodeb=nint(prop(index+2))
            nodec=nint(prop(index+3))
            radius=dsqrt((co(1,nodeb)+vold(1,nodeb)-
     &           co(1,nodea)-vold(1,nodea))**2)
            d=2*radius
            A=pi*radius**2
            l=dsqrt((co(2,nodec)+vold(2,nodec)-
     &           co(2,nodeb)-vold(2,nodeb))**2)
         endif
!
         pt1=v(2,node1)
         pt2=v(2,node2)
!
         if(pt1.ge.pt2) then
            inv=1
            Tt1=v(0,node1)-physcon(1)
            Tt2=v(0,node2)-physcon(1)
         else
            inv=-1
            pt1=v(2,node2)
            pt2=v(2,node1)
            Tt1=v(0,node2)-physcon(1)
            Tt2=v(0,node1)-physcon(1)
         endif
!
         p2p1=pt2/pt1
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         tdkp1=2.d0/kp1
         C2=tdkp1**kdkm1
! 
!        estimate of the flow using the orifice relationships
!        the flow is needed for Reynolds, Reynolds is needed
!        for the friction coefficient
!
         if(v(1,nodem)*inv.le.0.d0) then
            if(p2p1.gt.C2) then
               xflow=inv*pt1*A*dsqrt(2.d0*kdkm1*p2p1**(2.d0/kappa)
     &              *(1.d0-p2p1**(1.d0/kdkm1))/r)/dsqrt(Tt1)
            else
               xflow=inv*pt1*A*dsqrt(kappa/r)*tdkp1**(kp1/(2.d0*km1))/
     &              dsqrt(Tt1)
            endif
         else            
            xflow=v(1,nodem)*iaxial
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
         reynolds=dabs(xflow)*d/(dvi*A)
!
         call friction_coefficient(l,d,ks,reynolds,form_fact,
     &        lambda)
!
!        estimate of the flow using the incompressible relationships 
!        for a gas pipe (for zero initial flow or flow from low
!        to high pressure)
!
         if(v(1,nodem)*inv.le.0.d0) then
            xflow=inv*A*dsqrt(d/(lambda*l)*2*pt1/(r*Tt1)*(pt1-pt2))
         endif
!
         call pt2zpt1_crit(pt2,pt1,Tt1,lambda,kappa,r,l,d,
     &         pt2zpt1_c,Qred1_crit,crit,icase,M1_c)
!
         Qred=dabs(xflow)*dsqrt(Tt1)/(A*pt1)
!
!        correcting the mass flow if not physical 
!        (pt2/pt1 too small or Qred too big)
!
         if(crit) then
!
!           the flow is set to half the critical value
!
            xflow=0.5d0*inv*Qred1_crit*pt1*A/dsqrt(Tt1)
         elseif(Qred.gt.Qred1_crit) then
!
!           the flow is set to half the critical value
!
            xflow=0.5d0*inv*Qred1_crit*pt1*A/dsqrt(Tt1)
         endif
!
!        isothermal case: correcting the temperatures
!
         if(icase.eq.1) then
            call ts_calc(xflow,Tt1,pt1,kappa,r,A,T1,icase)
            call ts_calc(xflow,Tt2,pt2,kappa,r,A,T2,icase)
            if(inv.eq.1) then 
               v(3,node1)=T1
               v(3,node2)=T1
               if(nactdog(0,node2).eq.1) then
                  v(0,node2)=T1*(Tt2/T2)
               endif
            else
               v(3,node2)=T1
               v(3,node1)=T1
               if(nactdog(0,node1).eq.1) then
                  v(0,node1)=T1*(Tt2/T2) 
               endif
            endif
         endif
!
      elseif(kflag.eq.2)then
!
         numf=5
!
         pi=4.d0*datan(1.d0)
!
         kappa=(cp/(cp-r))
         km1=kappa-1.d0
         kp1=kappa+1.d0
!
         index=ielprop(nelem)
         A=prop(index+1)
         d=prop(index+2)
!
         l=prop(index+3)
         ks=prop(index+4)
         if(lakon(nelem)(2:6).eq.'GAPFA') then
            icase=0
         elseif(lakon(nelem)(2:6).eq.'GAPFI') then
            icase=1
         endif
         form_fact=prop(index+5)
         xflow_oil=prop(index+6)
         k_oil=nint(prop(index+7))
!
         if(lakon(nelem)(7:8).eq.'FR') then
!
!           flexible radius
!            
            nodea=nint(prop(index+1))
            nodeb=nint(prop(index+2))
            radius=dsqrt((co(1,nodeb)+vold(1,nodeb)-
     &           co(1,nodea)-vold(1,nodea))**2)
            A=pi*radius**2
            d=2*radius
!
         elseif(lakon(nelem)(7:8).eq.'RL') then
!
!           flexible radius and length
!
            nodea=nint(prop(index+1))
            nodeb=nint(prop(index+2))
            nodec=nint(prop(index+3))
            radius=dsqrt((co(1,nodeb)+vold(1,nodeb)-
     &           co(1,nodea)-vold(1,nodea))**2)
            d=2*radius
            A=pi*radius**2
            l=dsqrt((co(2,nodec)+vold(2,nodec)-
     &           co(2,nodeb)-vold(2,nodeb))**2)
         endif
!
         pt1=v(2,node1)
         pt2=v(2,node2)
         xflow=v(1,nodem)*iaxial
!
!        inv is the sign of the flow
!        xflow is replaced by its absolute value
!        wrongdir means that the flow goes from low
!        pressure to high pressure
!
         if(xflow.lt.0.d0) then
            inv=-1
         else
            inv=1
         endif
c         xflow=dabs(xflow)
         if((pt1-pt2)*inv.lt.0.d0) then
            wrongdir=.true.
         else
            wrongdir=.false.
         endif
!
!        the element is reoriented such that the mass flow
!        is directed from node 1 to node 2;
!        the pressure in node 1 may be more or less than
!        the pressure in node 2
!
         if(pt1.gt.pt2) then
!
            Tt1=v(0,node1)-physcon(1)
            call ts_calc(xflow,Tt1,pt1,kappa,r,A,T1,icase)
!     
            nodef(1)=node1
            nodef(2)=node1
            nodef(3)=nodem
            nodef(4)=node2
            nodef(5)=node2
         else
            pt1=v(2,node2)
            pt2=v(2,node1)
!
            Tt1=v(0,node2)-physcon(1)
            call ts_calc(xflow,Tt1,pt1,kappa,r,A,T1,icase)
!
            nodef(1)=node2
            nodef(2)=node2
            nodef(3)=nodem
            nodef(4)=node1
            nodef(5)=node1
         endif
!
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=2
         idirf(5)=0
!     
         pt2zpt1=pt2/pt1
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
         if(xflow_oil.ne.0d0) then
!     
!           two-phase-flow
!
            xflow_air=xflow
            call two_phase_flow(Tt1,pt1,T1,Tt2,pt2,T2,xflow_air,
     &           xflow_oil,nelem,lakon,kon,ipkon,ielprop,prop,
     &           v,dvi,cp,r,k_oil,phi,lambda,nshcon,nrhcon,
     &           shcon,rhcon,ntmat_,mi)
            lambda=lambda*phi
         else
!     
!           for pure air
!     
c            phi=1.d0
            call friction_coefficient(l,d,ks,reynolds,form_fact,
     &           lambda)
         endif
!     
!        calculating the critical conditions
!     
         call pt2zpt1_crit(pt2,pt1,Tt1,lambda,kappa,r,l,d,
     &         pt2zpt1_c,Qred1_crit,crit,icase,M1_c)
!
         if(wrongdir) lambda=-lambda
!
         Qred1=dabs(xflow)*dsqrt(Tt1)/(A*pt1)
!
!        check whether flow is critical
!        assigning the physcical correct sign to xflow
!
         if(crit) then
            xflow=inv*Qred1_crit*A*pt1/dsqrt(Tt1)
!
!           check whether flow has changed; if so, update v
!           for consistency
!
            if(dabs((xflow-iaxial*v(1,nodem))/xflow).gt.1.d-5) then
               iplausi=0
c               if(nactdog(1,nodem).ne.0) v(1,nodem)=xflow/iaxial
            endif
!            
            M1=dsqrt(2/km1*((Tt1/T1)-1.d0))
            if(icase.eq.0) then
               M1=min(M1,0.999d0)
            else
               M1=min(M1,0.999d0/dsqrt(kappa))
            endif
         else
            if(Qred1.gt.Qred1_crit) then
               xflow=inv*Qred1_crit*A*pt1/dsqrt(Tt1)
!
!              check whether flow has changed; if so, update v
!              for consistency
!
               if(dabs((xflow-iaxial*v(1,nodem))/xflow).gt.1.d-5) then
                  iplausi=0
                  if(nactdog(1,nodem).ne.0) v(1,nodem)=xflow/iaxial
               endif
!            
               M1=M1_c
            else
               M1=dsqrt(2/km1*((Tt1/T1)-1.d0))
            endif
!
!           determining M2: Tt2 -> T2 -> Tt2/T2 -> M2
!           the actual value of Tt2 is not relevant
!
            if(icase.eq.0) then
               Tt2=Tt1
            elseif(inv.eq.1) then
               Tt2=v(0,node2)-physcon(1)
            else
               Tt2=v(0,node1)-physcon(1)
            endif
            call ts_calc(xflow,Tt2,pt2,kappa,r,A,T2,icase)
            M2=dsqrt(2/km1*((Tt2/T2)-1.d0))
         endif
!
         Z1=M1**2
!         
         bb=km1/2.d0
         cc=-kp1/(2.d0*km1)
         ee1=M1*(1.d0+bb*Z1)/(1.d0+bb*Z1*(1.d0+2.d0*cc))
!     
!     definition of the coefficients
!
         lld=lambda*l/d
!
!        adiabatic case
!     
         if(icase.eq.0) then
!
            dfdM1=2.d0*(Z1-1.d0)/(kappa*M1**3*(1.d0+bb*Z1))
!     
            if(.not.crit) then
!     
!              residual
!
               Z2=M2**2
!               
               ee2=M2*(1.d0+bb*Z2)/(1.d0+bb*Z2*(1.d0+2.d0*cc))
               dfdM2=2.d0*(1.d0-Z2)/(kappa*M2**3*(1.d0+bb*Z2))
!
               f=(1.d0/Z1-1.d0/Z2)/kappa+kp1/(2.d0*kappa)*
     &           dlog(((1.d0+bb*Z2)*Z1)/
     &                ((1.d0+bb*Z1)*Z2))-lld
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
               df(3)=(dfdM1*ee1+dfdM2*ee2)/(xflow)
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
               f=(1.d0/Z1-1.d0)/kappa+kp1/(2.d0*kappa)*
     &           dlog(((1.d0+bb)*Z1)/
     &                ((1.d0+bb*Z1)))-lld
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
               df(3)=dfdM1*ee1/(xflow)
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
         elseif(icase.eq.1) then
!     
!           isothermal icase
!     
            dfdM1=2.d0*(kappa*Z1-1.d0)/(kappa*M1**3)
!
            if(.not.crit) then
!
               Z2=M2**2
!
               ee2=M2*(1.d0+bb*Z2)/(1.d0+bb*Z2*(1.d0+2.d0*cc))
               dfdM2=2.d0*(1.d0-kappa*Z2)/(kappa*M2**3)
!
!              redidual
!     
               f=(1.d0/Z1-1.d0/Z2)/kappa+dlog(Z1/Z2)-lld
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
               df(3)=(dfdM1*ee1+dfdM2*ee2)/(xflow)
c               df(3)=(dfdM1*ee1+dfdM2*ee2)/(inv*xflow)
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
!              residual
!  
               f=(1.d0/Z1-kappa)/kappa+dlog(kappa*Z1)-lld
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
               df(3)=dfdM1*ee1/(xflow)
c               df(3)=dfdM1*ee1/(inv*xflow)
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
         kp1=kappa+1.d0
         km1d2k=km1/(2.d0*kappa)
!         
         index=ielprop(nelem)
         A=prop(index+1)
         d=prop(index+2)    
         l=prop(index+3)
!
         lambda=0.5d0
!
         ks=prop(index+4)
         if(lakon(nelem)(2:6).eq.'GAPFA') then
            icase=0
         elseif(lakon(nelem)(2:6).eq.'GAPFI') then
            icase=1
         endif
         form_fact=prop(index+5)
         xflow_oil=prop(index+6)
         k_oil=nint(prop(index+7))
!
         pt1=v(2,node1)
         pt2=v(2,node2)
!     
         if(pt1.gt.pt2) then
            inv=1
            xflow=v(1,nodem)*iaxial
            Tt1=v(0,node1)-physcon(1)
            call ts_calc(xflow,Tt1,pt1,kappa,r,A,T1,icase)
            if(icase.eq.0) then
               Tt2=Tt1
            else
               T2=T1
            endif
!     
         else
            inv=-1
            pt1=v(2,node2)
            pt2=v(2,node1)
            xflow=v(1,nodem)*iaxial
            Tt1=v(0,node2)-physcon(1)
            call ts_calc(xflow,Tt1,pt1,kappa,r,A,T1,icase)
            if(icase.eq.0) then
               Tt2=Tt1
            else
               T2=T1
            endif
!     
         endif
!     
         pt2zpt1=pt2/pt1
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
!     definition of the friction coefficient for 2 phase flows and pure air
!     
         if(xflow_oil.ne.0d0) then
            call two_phase_flow(Tt1,pt1,T1,Tt2,pt2,T2,xflow,
     &           xflow_oil,nelem,lakon,kon,ipkon,ielprop,prop,
     &           v,dvi,cp,r,k_oil,phi,lambda,nshcon,nrhcon,
     &           shcon,rhcon,ntmat_,mi)
            lambda=lambda*phi
!     
!     for pure air
!     
         else
            phi=1.d0
            call friction_coefficient(l,d,ks,reynolds,form_fact,
     &           lambda)
         endif
!     
         call pt2zpt1_crit(pt2,pt1,Tt1,lambda,kappa,r,l,d,
     &        pt2zpt1_c,Qred1_crit,crit,icase,M1_c)
!     
!     definition of the coefficients 
!     
         M1=dsqrt(2/km1*((Tt1/T1)-1))
         if(crit) then
            if(icase.eq.0) then
               M2=1.d0
               T2=2.d0*Tt2/kp1
            else
               M2=1.d0/dsqrt(kappa)
               Tt2=T2*(1.d0+km1d2k)
            endif
         else
            if(icase.eq.0) then
               call ts_calc(xflow,Tt2,pt2,kappa,r,A,T2,icase)
            else
               Tt2=v(0,node2)-physcon(1)
               call tt_calc(xflow,Tt2,pt2,kappa,r,A,T2,icase)
            endif
            M2=dsqrt(2/km1*((Tt2/T2)-1))
         endif
!     
         write(1,*) ''
         write(1,55) ' from node ',node1,
     &           ' to node ', node2,':   air massflow rate = ',xflow,
     &           ' , oil massflow rate = ',xflow_oil
!     
         if(inv.eq.1) then
            write(1,53)'       Inlet node ',node1,' :    Tt1 = ',Tt1,
     &              ' , Ts1 = ',T1,'  , Pt1 = ',pt1,
     &              ' , M1 = ',M1
            write(1,*)'             Element ',nelem,lakon(nelem)
            write(1,57)'             dvi = ',dvi,' , Re = '
     &              ,reynolds
            write(1,58)'             PHI = ',phi,' , LAMBDA = ',
     &              lambda,
     &              ', LAMBDA*l/d = ',lambda*l/d,' , ZETA_PHI = ',
     &              phi*lambda*l/d
            write(1,53)'      Outlet node ',node2,' :    Tt2 = ',
     &              Tt2,
     &              ' , Ts2 = ',T2,'  , Pt2 = ',pt2,
     &              ' , M2 = ',M2
!    
         else if(inv.eq.-1) then
            write(1,53)'       Inlet node ',node2,':    Tt1= ',Tt1,
     &              ' , Ts1= ',T1,' , Pt1= ',pt1,
     &              ' , M1= ',M1
            write(1,*)'             Element ',nelem,lakon(nelem)
            write(1,57)'             dvi = ',dvi,' , Re = '
     &              ,reynolds
            write(1,58)'             PHI = ',phi,' , LAMBDA = ',
     &              lambda,
     &              ', LAMBDA*l/d = ',lambda*l/d,' , ZETA_PHI = ',
     &              phi*lambda*l/d
            write(1,53)'      Outlet node ',node1,' :    Tt2 = ',
     &              Tt2,
     &              ' , Ts2 = ',T2,'  , Pt2 =',pt2,
     &              ' , M2 = ',M2
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
