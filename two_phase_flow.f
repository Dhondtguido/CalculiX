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
      subroutine two_phase_flow(Tt1,pt1,T1,Tt2,pt2,T2,xflow_air,
     &     xflow_oil,nelem,lakon,kon,ipkon,ielprop,prop,v,
     &     dvi_air,cp,r,k_oil,phi,lambda,nshcon,nrhcon,shcon,
     &     rhcon,ntmat_,mi,iaxial)
!     
!    two phase flow correlations 
!
!     author: Yannick Muller
!     
      implicit none
!     
      character*8 lakon(*)
!     
      integer nelem,ielprop(*),index,mi(*),iaxial,
     &     ipkon(*),kon(*),icase,kgas,k_oil,mtlog,ier,nshcon(*),
     &     nrhcon(*),ntmat_,n1,n2,n11
!
!     note: Tt2 and T2 are used in company routines
!
      real*8 prop(*),v(0:mi(2),*),kappa,R,a,d,dl,
     &     T1,Tt1,pt1,pt2,cp,dvi_air,dvi_oil,
     &     reynolds,lambda,ks,form_fact,f,
     &     l_neg,xflow_air,xflow_oil,A1,A2,
     &     rho_air,rho_oil,nue_air,nue_oil,zeta,reynolds_h,mpg,
     &     xp,xpm2,xpmini,isothermal,dvi_h,zeta_h,auxphi,
     &     rad,theta,phi,phizeta,x,Tt2,T2,
     &     rho_q,p1,shcon(0:3,ntmat_,*),
     &     rhcon(0:1,ntmat_,*),cp_oil,r_oil
!
      parameter ( xpmini=1.E10)
!
!     this subroutine enables to take in account the existence of
!     2 phase flows (air /oil) in some flow elements.
!
!     lambda: friction coefficient solely due to air
!     phi: correction due to the presence of oil
!     (lambda_corrected=lambda*phi)
!    
!     the 2 following tables are used in Lockhart Martinelli Method.
!     See table p.44
!
      real*8 TX(17),TF(17)
      data TX
     &     /0.01d0,0.02d0,0.04d0,0.07d0,
     &      0.10d0,0.20d0,0.40d0,0.70d0,
     &      1.00d0,2.00d0,4.00d0,7.00d0,
     &      10.0d0,20.0d0,40.0d0,70.0d0,
     &      100.d0/
!
      data TF
     &    /1.28d0,1.37d0,1.54d0,1.71d0,
     &     1.85d0,2.23d0,2.83d0,3.53d0,
     &     4.20d0,6.20d0,9.50d0,13.7d0,
     &     17.5d0,29.5d0,51.5d0,82.0d0,
     &     111.d0/
!
      data n1 /1/
      data n2 /2/
      data n11 /11/
!
!
!
      index=ielprop(nelem)
!
      if(lakon(nelem)(2:5).eq.'GAPF') then
            A=prop(index+1)
            d=prop(index+2)
            dl=prop(index+3)
            ks=prop(index+4)
            form_fact=prop(index+5)
      endif
!      
      if(xflow_oil.eq.0.d0) then
         write(*,*) '*WARNING:in two_phase_flow'
         write(*,*) 'massflow oil for element',nelem,'in null'
         write(*,*) 'Calculation proceeds without oil correction'
         phi=1.d0
      endif
!
 
      xflow_air=dabs(xflow_air)
      kappa=Cp/(Cp-R)
!
!     First case:
!     the element is a restrictor of type
!     THICK-WALLED ORIFICE IN LARGE WALL (L/DH > 0.015)
!     I.E. IDL'CHIK (SECTION IV PAGE 144)!
!     and
!     Second case:
!     the element is a restrictor of type
!     SMOOTH BENDS B.H.R.A HANDBOOK (Miller)
!
!     Two phase flow correlations are taken from:
!     H.Zimmermann, A.Kammerer, R.Fischer and D. Rebhan
!     "Two phase flow correlations in Air/Oil systems of
!     Aero Engines."
!     ASME 91-GT-54     
!
      if((lakon(nelem)(2:7).eq.'RELOID').or.
     &     (lakon(nelem)(2:7).eq.'REBEMI')) then
!
         icase=0
!
         A1=prop(index+1)
         A2=prop(index+2)
         call ts_calc(xflow_air,Tt1,pt1,kappa,r,A1,T1,icase)
!
         d=dsqrt(A1*4/(4.d0*datan(1.d0)))
!
!     calculating the dynamic viscosity, the kinematic viscosity and
!     the density of air
!            
            kgas=0
!
         p1=pt1*(T1/Tt1)**(kappa/kappa-1)
         rho_air=p1/(R*T1)
         nue_air=dvi_air/rho_air
!     
!     calculating the dynamic viscosity, the kinematic viscosity and
!     the density of oil
!
         call materialdata_tg(k_oil,ntmat_,T1,shcon,nshcon,cp_oil,r_oil,
     &        dvi_oil,rhcon,nrhcon,rho_oil)
!
         if(xflow_oil.eq.0.d0) then
!     
!     pure air
            call zeta_calc(nelem,prop,ielprop,lakon,reynolds,zeta,
     &           isothermal,kon,ipkon,R,Kappa,v,mi,iaxial)
            lambda=zeta
            return
         else
!     
!     air/oil mixture for orifice or bend
!     For Bend see section 4.2.1
!     For orifices see 4.2.3
!
            mpg=xflow_air+xflow_oil
            xp=xflow_air/mpg
            if(mpg.gt.xflow_air*xpmini) then
               xpm2=xpmini**2
            else
               xpm2=(mpg/xflow_air)**2
            endif
!     
            rho_q=rho_oil/rho_air
!     
!     homogene dynamic viscosity (mass flow rate averaged)
!
            dvi_h=dvi_oil*dvi_air/((dvi_oil-dvi_air)*xp+dvi_air)
!     
!     homogene reynolds number
!
            reynolds_h=mpg*d/(A1*dvi_h)
!
            call zeta_calc(nelem,prop,ielprop,lakon,reynolds_h,zeta_h,
     &           isothermal,kon,ipkon,R,Kappa,v,mi,iaxial)
!     
!     orifice in a wall
!
            if(lakon(nelem)(2:7).eq.'RELOID') then
               
               auxphi=(1.d0+xp*(rho_q**(1.d0/6.d0)-1.d0))
     &               *(1.d0+xp*(rho_q**(5.d0/6.d0)-1.d0))
!     
!     bend
!
            elseif(lakon(nelem)(2:7).eq.'REBEMI') then
!     
!     radius of the bend
!
               rad=prop(index+4)
!
!     angle of the bend
!
               theta=prop(index+5)
!     
               f=(1.d0+2.2d0*theta/90.d0/(zeta_h*(2.d0+rad/d)))
     &              *xp*(1.d0-xp)+xp**2
!     
               auxphi=1.d0+(rho_q-1.d0)*f
            endif
!     
            phi=1/rho_q*auxphi*xpm2
            phizeta=zeta_h/rho_q*auxphi*xpm2
            lambda=zeta_h
!            
         endif   
!         
!     Third case:
!     the element is a pipe
!     the zeta coefficient is corrected according to 
!     Lockhart Martinelli Method
!     Reference: R.W. Lockhart and R.C. Martinelli
!                University of California, BErkeley, California
!                "Proposed correlation of data for 
!                 isothermal two-phase two-component
!                 flow in pipes"
!                 Chemical Engineering Progress vol.45, N°1
!
      elseif((lakon(nelem)(2:5).eq.'GAPF')
     &        .or.((lakon(nelem)(2:7).ne.'REBEMI')
     &        .and.(lakon(nelem)(2:7)).ne.'RELOID'))then
!
         if(lakon(nelem)(2:6).eq.'GAPFA') then
            icase=0
         elseif(lakon(nelem)(2:6).eq.'GAPFI') then
            icase=1
         else
            icase=0
         endif
!
         if((lakon(nelem)(2:3).eq.'RE').and.
     &        (lakon(nelem)(4:5).ne.'BR')) then
            a=min(prop(index+1),prop(index+2))
         endif
!
         call ts_calc(xflow_air,Tt1,pt1,kappa,r,a,T1,icase)
!
!     calculating kinematic viscosity and density for air 
!
         p1=pt1*(T1/Tt1)**(kappa/kappa-1)
         rho_air=p1/(R*T1)
         nue_air=dvi_air/rho_air
!         
!     calculation of the dynamic viscosity for oil
!
         call materialdata_tg(k_oil,ntmat_,T1,shcon,nshcon,cp_oil,r_oil,
     &        dvi_oil,rhcon,nrhcon,rho_oil)
!
         nue_oil=dvi_oil/rho_oil
!
!     Definition of the two phase flow modulus as defined in table 1
!
         x=dabs(xflow_oil/xflow_air)*(rho_air/rho_oil)**(0.553d0)
     &        *(nue_oil/nue_air)**(0.111d0)
!
         mtlog=17
!     Interpolating x in the table
         call onedint(TX,TF,mtlog,x,phi,n1,n2,n11,IER)
!
         if((lakon(nelem)(2:4).eq.'GAP'))then
!
!     Computing the friction coefficient
!           
            reynolds=dabs(xflow_air)*d/(dvi_air*a)
!         
            if(reynolds.lt.100.d0) then
               reynolds= 100.d0
            endif
!     
            call friction_coefficient(l_neg,d,ks,reynolds,form_fact,
     &           lambda)
         else 
            lambda=0
         endif
      endif
!
      return
      end
      
      
 
