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
      subroutine moehring (node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,dvi,numf,set,mi,ttime,time,
     &     iaxial,iplausi)
!     
!     moehring element
!     This subroutines computes the evolution of the core swirl ratio 
!     for a disc stator system with either centrifugal or centripetal
!     flow.
!     Theoretical explanations can be found in
!     "Untersuchung dfes radialen Druckverlaufes und des übertragenen
!     drehmomentes im Radseitenraum von Kreiselpumpen bei glatter,
!     ebene Radwand und bei Anvendung von Rückenschaufeln"
!     Uwe Klaus Möhring , Dissertation, 
!     An der Üniversität Carolo-Wilhelmina zu Braunschweig 1976
!
!     author: Yannick Muller
!     
      implicit none
!     
      logical identity
      character*8 lakon(*)
      character*81 set(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(*),idirf(*),index,kflag,iaxial,
     &     ipkon(*),kon(*),kgas,key,neval,ier,limit,lenw,last,
     &     iwork2(400),node_up,node_down,mi(*),iplausi
!     
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(*),r,dvi,pi,
     &     R_min, R_max,cr, R_shroud,rsrmax,gap,swirl_up,
     &     pup,pdown,tup,tdown,kappa,cp,ttime,time,
     &     Rup,Rdown,K0,Kup,Cq,Re_phi,phi,lambda1, lambda2,
     &     a,b,epsabs,
     &     epsrel,result,abserr,work(1200),zk0,T1,T2,P1,P2,Pr,
     &     qred_crit,omega,rurd,C_p,Cm,Mr,f_k,f_t,f_p,f_m,f_cm,
     &     pdiff,pdiff_min,xflow_0,Cq_0,phi_0
!     
      external dKdX,dKdP,dkdT,dKdm,f_k,f_t,f_p,f_m,f_cm
!
!
!     
!      numf=4
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
         kappa=(cp/(cp-R))
         index=ielprop(nelem)
         qred_crit=dsqrt(kappa/R)*
     &     (1.d0+0.5d0*(kappa-1.d0))**(-0.5d0*(kappa+1.d0)/(kappa-1.d0))
!     
!     Because there is no explicit expression relating massflow
!     to pressure loss for möhrings
!     initial mass flow is set to arbitrarily
!     with consideration to flow direction
!     
         node1=kon(ipkon(nelem)+1)
         node2=kon(ipkon(nelem)+3)
         p1=v(2,node1)
         p2=v(2,node2)
         T1=v(0,node1)
         T2=v(0,node2)
!     
!     fictious cross section
!         A=1d-5
         if(p1.gt.p2) then
            xflow=1/dsqrt(T1)*P1*qred_crit*0.5d0
         else
            xflow=-1/dsqrt(T1)*P1*qred_crit*0.5d0
         endif
!     
      elseif(kflag.eq.2)then
!     
         numf=4
         index=ielprop(nelem) 
         pi=4.d0*datan(1.d0)
!     
!     Cr
         Cr=0.315d0
!     
!     minimal disc radius
         R_min=prop(index+1)
!     
!     maximal disc radius
         R_max=prop(index+2)
!     
!     R_min/R_max
         Rurd=R_min/R_max
!     
!     disc/stator gap
         gap=prop(index+3)
!     
!     shroud radius
         R_shroud=prop(index+4)
!     
!     R_schroud/R_max
         rsrmax=R_shroud/R_max      
!     
!     defining flow parameters
!     
!     massflow
         xflow=v(1,nodem)*iaxial
!     
!     upstream node
         node_up=nint(prop(index+5))
!     
!     downstream node
         node_down=nint(prop(index+6))
!   
!     centripetal
         if(lakon(nelem)(2:5).eq.'MRGP') then
            if(xflow.lt.0d0) then
               xflow=-xflow
            endif
!     
            Rup=R_max
            Rdown=R_min
!     
            nodef(1)=node_up
            nodef(2)=node_up
            nodef(3)=nodem
            nodef(4)=node_down
!     
!     centrifugal
         elseif(lakon(nelem)(2:5).eq.'MRGF') then
            if(xflow.gt.0d0) then
               xflow=-xflow
            endif
!     
            Rup=R_min
            Rdown=R_max
! 
!            if(xflow.gt.0) then
            nodef(1)=node_up
            nodef(2)=node_up
            nodef(3)=nodem
            nodef(4)=node_down
         endif
!     
!     upstream pressure
         pup=v(2,node_up)
!
!     downstream pressure
         pdown=v(2,node_down)
!
!     Upstream temperature
         Tup=v(0,node_up)
!
!     downstream temperature
         Tdown=v(0,node_down)
!     
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=2
!     
!     rotation related parameters
!     
!     rotation
!
         omega=prop(index+7)
!     
!     swirl at R_upstream
         swirl_up=prop(index+8)
!     
!     core swirl ratio when xflow=0
         K0=1.d0/(1.d0+(rsrmax**3.6d0*
     &        (rsrmax+4.6d0*gap/R_max))**(4.d0/7.d0))
!     
!     core swirl ratio at R_inlet
         Kup=swirl_up/(omega*Rup)
!     
!     dynamic_viscosity
         if(dabs(dvi).lt.1d-30) then
            write(*,*) '*ERROR in moehring: '
            write(*,*) '       no dynamic viscosity defined'
            write(*,*) '       dvi= ',dvi
            call exit(201)
         endif 
!     
!     defining common coefficients
!     
         Cq=xflow*R*Tup/(Pup*omega*(R_max)**3)
!     
         Re_phi=(omega*R_max**2*Pup)/(dvi*R*Tup)
!     
         phi=Cq*(Re_phi)**0.2d0
!
         zk0=(1-K0)/K0
!     
!     lambda1
         lambda1=(R_max-R_min)/dabs(R_max-R_min)*Pi*Cr/4*
     &        (dvi*R/(omega*R_max**2)**0.2d0*(omega*R_max**3)/R)
!     
!     lambda2
         lambda2=2d0*R/(omega**2*R_max**2)
!     
!*************************************************************************
!     integration of K(X),dKdp(X),dKdT(X),dKdm(X)

         limit=201
         lenw=5*limit
!     
!         if(lakon(nelem)(2:5).eq.'MRGF') then
!xflow.lt.0d0) then
!     
!     lower integration boundary
            a=rurd
!     
!     upper integration boundary
            b=1d0
!     
!         elseif(lakon(nelem)(2:5).eq.'MRGP') then
!     
!     lower integration boundary
!            a=1d0
!     
!     upper integration boundary
!            b=rurd
!         endif         
!     
!     absolute error
         epsabs=1d-7
!     
!     relative error
         epsrel=1d-7
!     
!     choice for local integration rule
         key=1
!
!     determining minimal pressure difference for xflow<<1
!
         if(lakon(nelem)(2:5).eq.'MRGF')then
            xflow_0=-0.00000003d-3
         elseif(lakon(nelem)(2:5).eq.'MRGP') then
            xflow_0=0.00000003d-3
         endif
!
         Cq_0=xflow_0*R*Tup/(Pup*omega*(R_max)**3)
!     
         phi_0=Cq_0*(Re_phi)**0.2d0
!
         call dqag(f_k,a,b,epsabs,epsrel,key,result,abserr,neval,ier,
     &        limit,lenw,last,iwork2,work,phi_0,lambda1,zk0,Pup,Tup,
     &        rurd,xflow_0,kup)
!
!     pressure coefficient
         C_p=2*result
         pdiff_min=C_p*Pup/(4*R*Tup)*omega**2*R_max**2
         pdiff=dabs(pdown-pup)
!
!     K(x)
!     
         call dqag(f_k,a,b,epsabs,epsrel,key,result,abserr,neval,ier,
     &        limit,lenw,last,iwork2,work,phi,lambda1,zk0,Pup,Tup,rurd,
     &        xflow,kup)
!        
!     residual
!
         if(lakon(nelem)(2:5).eq.'MRGF') then
            f=lambda2*(Pdown-Pup)/(Pdown+Pup)*Tup-result
         elseif(lakon(nelem)(2:5).eq.'MRGP') then
            f=lambda2*(Pup-Pdown)/(Pup+Pdown)*Tup-result
         endif
!     
!     pressure coefficient
         C_p=2*result
!     
!     moment coefficient
         call dqag(f_cm,a,b,epsabs,epsrel,key,result,abserr,neval,ier,
     &        limit,lenw,last,iwork2,work,phi,lambda1,zk0,Pup,Tup,rurd,
     &        xflow,kup)
!         
         Cm=0.5d0*Pi*Cr*Re_phi**(-0.2d0)*result
         Mr=0.5d0*Cm*(Pup/1000d0/(R*Tup))*omega**2*R_max**5
         Pr=Mr*omega
!     
!     pressure
!     
         call dqag(f_p,a,b,epsabs,epsrel,key,result,abserr,neval,ier,
     &        limit,lenw,last,iwork2,work,phi,lambda1,zk0,Pup,Tup,rurd,
     &        xflow,kup)
!     
!     partial derivative (upstream pressure)
! 
         if(lakon(nelem)(2:5).eq.'MRGF') then
            df(1)=-2*lambda2**Pdown/(Pdown+Pup)**2*Tup-result
         elseif(lakon(nelem)(2:5).eq.'MRGP') then
            df(1)=2*lambda2**Pdown/(Pdown+Pup)**2*Tup-result
         endif
!     
!     temperature
!     
         call dqag(f_t,a,b,epsabs,epsrel,key,result,abserr,neval,ier,
     &        limit,lenw,last,iwork2,work,phi,lambda1,zk0,Pup,Tup,rurd,
     &        xflow,kup)
!     
!     partial derivative (upstream temperature)       
         if(lakon(nelem)(2:5).eq.'MRGF') then
            df(2)=lambda2**(Pdown-Pup)/(Pdown+Pup)-result
         elseif(lakon(nelem)(2:5).eq.'MRGP') then  
            df(2)=lambda2**(Pup-Pdown)/(Pdown+Pup)-result
         endif
!     
!     mass flow
!     
         call dqag(f_m,a,b,epsabs,epsrel,key,result,abserr,neval,ier,
     &        limit,lenw,last,iwork2,work,phi,lambda1,zk0,Pup,Tup,rurd,
     &        xflow,kup)
!     
!     partial derivative (mass flow)             
         df(3)=-result
!     
!     pressure 
!     partial derivative (downstream pressure)
         if(lakon(nelem)(2:5).eq.'MRGF') then
            df(4)=2*lambda2**Pup/(Pdown+Pup)**2*Tup
         elseif(lakon(nelem)(2:5).eq.'MRGP') then
            df(4)=-2*lambda2**Pup/(Pdown+Pup)**2*Tup
         endif
!   
      endif
!     
      xflow=xflow/iaxial
      df(3)=df(3)*iaxial
!
      return
      end
!
******************************************************************************
      function f_k(x,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
!
      implicit none
      integer neq,idid,ipar,iwork(100),lrw,liw,j
      real*8 f_k,x,rpar(8),rtol,atol,y(1),info(15),
     &     Rurd,zk0,lambda1,Kup,xflow,pup,tup,phi,t,rwork(160)
!     
      external dKdX
!
!     storing the parameters
      rpar(1)=phi
      rpar(3)=zk0
!
!    relative error
         rtol=1.d-7
!     absolute error
         atol=1.d-7
!
!     initial value
         if(xflow.lt.0d0) then
            t=rurd
         else 
            t=1.d0
         endif

         neq=1
!
!     initialisation info field
         do j=1,15
            info(j)=0
         enddo
!     initial condition f(0)
!     core swirl ratio at Rup repectively Rdown depending 
!     on the type of element centrifugal or centripetal
!     
         y(1)= Kup 
!
         lrw=160
         liw=60
!
!     solving the differential equation Möhring 3.35
!     dK/dX=f(K(X))
!
         if(dabs(xflow).gt.1d-6) then
            call ddeabm(dKdX,neq,t,y,x,info,rtol,atol,idid,
     &        rwork,lrw,iwork,liw,rpar,ipar)
         else
            y(1)=1/(Zk0+1)
         endif
!
       f_k=y(1)**2*x
!
       return
       end
******************************************************************************
      function f_p(x,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
!
      implicit none
      integer neq,idid,ipar,iwork(100),lrw,liw,j
      real*8 f_p,x,rpar(8),rtol,atol,y(1),info(15),Rurd,
     &     zk0,lambda1,Kup,xflow,pup,tup,phi,t,rwork(160)
!     
      external dKdp
!     storing the parameters
      rpar(1)=phi
      rpar(2)=lambda1
      rpar(3)=zk0
      rpar(4)=Pup
      rpar(5)=Tup
      rpar(6)=rurd
      rpar(7)=xflow
      rpar(8)=kup
!

!    relative error
      rtol=1.d-7
!     absolute error
      atol=1.d-7
!
!     initial value
      if(xflow.lt.0d0) then
         t=rurd
      else
         t=1.d0
      endif
      neq=1
!     
!     initialisation info field
      do j=1,15
         info(j)=0
      enddo
!     initial condition f(0)
!     core swirl ratio at Rup
!     
      y(1)= 0d0    
!      
      lrw=160
      liw=60
!          
      call ddeabm(dKdp,neq,t,y,x,info,rtol,atol,idid,
     &     rwork,lrw,iwork,liw,rpar,ipar)
!     
      f_p=2*y(1)*x
!     
      return
      end
*****************************************************************************
      function f_t(x,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
!     
      implicit none
      integer neq,idid,ipar,iwork(100),lrw,liw,j
      real*8 f_t,x,rpar(8),rtol,atol,y(1),info(15),Rurd,
     &     zk0,lambda1,Kup,xflow, pup,tup,phi,t,rwork(160)
!     
      external dKdt
!
!     storing the parameters
      rpar(1)=phi
      rpar(2)=lambda1
      rpar(3)=zk0
      rpar(4)=Pup
      rpar(5)=Tup
      rpar(6)=rurd
      rpar(7)=xflow
      rpar(8)=kup
!             
!     relative error
      rtol=1.d-7
!     absolute error
      atol=1.d-7
!
!     initial value
      if(xflow.lt.0d0) then
         t=rurd
      else 
         t=1.d0
      endif
!      if(xflow.lt.0d0) then
!         t=rurd
!      else 
!         t=1.d0
!      endif     
      neq=1
!     
!     initialisation info field
      do j=1,15
         info(j)=0
      enddo
!     initial condition f(0)
!     core swirl ratio at Rup
!     
      y(1)= 0d0    
!     
       lrw=160
       liw=60
!     
      call ddeabm(dKdt,neq,t,y,x,info,rtol,atol,idid,
     &     rwork,lrw,iwork,liw,rpar,ipar)
!     
      f_t=2d0*y(1)*x
!     
      return
      end
******************************************************************************
      function f_m(x,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
!     
      implicit none
      integer neq,idid,ipar,j,iwork(100),lrw,liw
      real*8 f_m,x,rpar(8),rtol,atol,y(1),info(15),
     &     Rurd,zk0,lambda1,Kup,xflow,pup,tup,phi,t,rwork(160)
!     
      external dKdm
!
!     storing the parameters
      rpar(1)=phi
      rpar(2)=lambda1
      rpar(3)=zk0
      rpar(4)=Pup
      rpar(5)=Tup
      rpar(6)=rurd
      rpar(7)=xflow
      rpar(8)=kup
!     
!     relative error
      rtol=1.d-7
!     absolute error
      atol=1.d-7
! 
!     initial value
      if(xflow.lt.0d0) then
         t=rurd
      else 
         t=1.d0
      endif    
      neq=1
!     
!     initialisation info field
      do j=1,15
         info(j)=0
      enddo
!     initial condition f(0)
!     core swirl ratio at Rup
!     
      y(1)=0     
!     
      lrw=160
      liw=60
!     
!     solving the differential equation Möhring 3.35
!     dK/dX=f(K(X))
!     
      call ddeabm(dKdm,neq,t,y,x,info,rtol,atol,idid,
     &     rwork,lrw,iwork,liw,rpar,ipar)
!     
      f_m=2*y(1)*x
!     
      return
      end
******************************************************************************
      function f_cm(x,phi,lambda1,zk0,Pup,Tup,rurd,xflow,kup)
!     
      implicit none
      integer neq,idid,ipar,j,iwork(100),lrw,liw
      real*8 f_cm,x,rpar(8),rtol,atol,y(1),info(15),
     &     Rurd,zk0,lambda1,Kup,xflow,pup,tup,phi,t,rwork(160)
!     
      external dKdX
!  
!     storing the parameters
      rpar(1)=phi
      rpar(2)=lambda1
      rpar(3)=zk0
      rpar(4)=Pup
      rpar(5)=Tup
      rpar(6)=rurd
      rpar(7)=xflow
      rpar(8)=kup
!   
!     relative error
      rtol=1.d-7
!     absolute error
      atol=1.d-7
!     initial value
!
      if(xflow.lt.0d0) then
         t=rurd
      else
         t=1.d0
      endif
!
      neq=1
!     
!     initialisation info field
      do j=1,15
         info(j)=0
      enddo
!     initial condition f(0)
!     core swirl ratio at Rup
!     
      y(1)=Kup     
!     
      lrw=160
      liw=60
!     
!     solving the differential equation Möhring 3.35
!     dK/dX=f(K(X))
!     
      call ddeabm(dKdX,neq,t,y,x,info,rtol,atol,idid,
     &     rwork,lrw,iwork,liw,rpar,ipar)
!     
      f_cm=dabs(1-y(1))/(1-y(1))*dabs(1-y(1))**(1.75d0)
     &     *x**(3.6d0)
!     
      return
      end
