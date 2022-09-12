!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
!     
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
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
      subroutine tee(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,r,physcon,numf,set,mi,ider,ttime,time,
     &     iaxial,iplausi)
!
!     A tee split element(zeta calculation according to Idel'chik)
!     Written by Yavor Dobrev
!
      implicit none
!
      logical identity
!
      character*81 set(*)
      character*8 lakon(*)
!
      integer
!     Number of the current element
     &nelem,
!     Array with the degrees of freedom
     &nactdog(0:3,*),
!     Number of the first node of the element
     &node1,
!     Number of the second node of the element
     &node2,
!     Number of the middle node of the element
     &nodem,
!     Number of the derivatives
     &numf,
!     Array with all nodes
     &kon(*),
!     Array, which shows where a given element starts in kon
     &ipkon(*),
!     Array which shows where the properties of a given element
!     begin in prop
     &ielprop(*),
!     Contains the nodes where the dependant varibles for this 
!     element are
     &nodef(*),
!     Contains the type of the dependant variables
     &idirf(*),
!     The starting index of the given element#s properties in prop
     &index,
!     Shows where we are in the program
     &kflag,
!     1 for normal, -1 for inverted flow
     &inv,
!     Shows how many types of DOF there are
     &mi(2),
!     This is the number of the previous element this element is
!     connected to
     &nelem1,
!     The middle node of the previous element
     &nodem1,
!     0 if a residual is to be calculated, 1 for derivatives
     &ider,chan_num,icase,iaxial,iplausi
!
      real*8 
!     Array with the properties of all elements
     &prop(*),
!     Array with the values of the unknown variables
     &v(0:mi(2),*),
!     Mass flow in the current element
     &xflow,
!     Residual value
     &f,
!     Array with the derivatives
     &df(*),
!     Isentropic exponent
     &kappa,
!     Gas constant
     &R,
!     Specific heat
     &cp,
!     Total temperatures and total pressures
     &Tt1,Tt2,pt1,pt2,
!     Array with physical constants
     &physcon(*),
!     Kappa terms
     &km1,kp1,kdkm1,tdkp1,
!     Pressure ratios
     &pt2pt1,pt2pt1_crit,
!     Areas
     &A,A1,A2,
!     Factor zeta_fac for influencing zeta
!     zeta_eff = zeta_fac*zeta
     &zeta_fac,
!     Variable for the function result
     &calc_residual_tee,
!     Mass flow in the previuos element
     &xflow1,
!     Mass flow in the current element(same as xflow)
     &xflow2,Ts0,pspt0,pspt2,M1,M2,Ts2,ttime,time,zeta
!
!
!
      index=ielprop(nelem)
!
      if(kflag.eq.0) then
!
!        Checking for degrees of freedom in the element
!
         identity=.true.
         if(nactdog(2,node1).ne.0)then
            identity=.false.
         elseif(nactdog(2,node2).ne.0)then
            identity=.false.
         elseif(nactdog(1,nodem).ne.0)then
            identity=.false.
!
         endif
!
      elseif(kflag.eq.1)then
         if(v(1,nodem).ne.0.d0) then
            xflow=v(1,nodem)
            return
         endif
!
!        Calculating the initial mass flow in the element
!
!        Some kappa terms
         kappa=(cp/(cp-R))
         kp1=kappa+1d0
         km1=kappa-1d0
         kdkm1=kappa/km1
         tdkp1=2.d0/kp1
!
!        Get the element properties
!        Cross sections
         if(nelem.eq.nint(prop(index+2))) then
            A=prop(index+5)
         elseif(nelem.eq.nint(prop(index+3)))then
            A=prop(index+6)
         endif
         zeta_fac=prop(index+11)
!
!        Pressures
         pt1=v(2,node1)
         pt2=v(2,node2)
!
!        Check for inverted flow
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
!        Pressure ratios
         pt2pt1=pt2/pt1
         pt2pt1_crit=tdkp1**kdkm1
!
         if(pt2pt1.gt.pt2pt1_crit) then
!           Subcritical case
            xflow=inv*pt1*A*dsqrt(2.d0*kdkm1*pt2pt1**(2.d0/kappa)
     &           *(1.d0-pt2pt1**(1.d0/kdkm1))/r)/dsqrt(Tt1)
         else
!           Critical case
            xflow=inv*pt1*A*dsqrt(kappa/r)*tdkp1**(kp1/(2.d0*km1))/
     &           dsqrt(Tt1)
         endif
!     
      elseif(kflag.eq.2)then
!
!        Calculation of residual/derivatives
!
!        Number of dependant variables(corresponding derivatives)
         numf=6
!
!        Calculate kappa
         kappa=(cp/(cp-R))
!
!        setting icase (always adiabatic)
         icase=0;
!
!        Inlet conditions are the same for both branches
!
!        Determining the previous element as
!        the incoming mass flow is defined there
         nelem1=nint(prop(index+1))
         nodem1=kon(ipkon(nelem1)+2)
!
!        Inlet conditions
         pt1=v(2,node1)
         Tt1=v(0,node1)-physcon(1)
         xflow1=v(1,nodem1)*iaxial
         A1=prop(index+4)
!
!        Outlet conditions
         Tt2=v(0,node2)
         xflow2=v(1,nodem)*iaxial
         pt2=v(2,node2)
         if(nelem.eq.nint(prop(index+2))) then
            A2=prop(index+5)
            zeta_fac=prop(index+11)
         elseif(nelem.eq.nint(prop(index+3))) then
            A2=prop(index+6)
            zeta_fac=prop(index+12)
         endif
!         zeta_fac = prop(index+11)
!
!        Set the node numbers for the degrees of freedom
         nodef(1)=node1
         nodef(2)=node1
         nodef(3)=nodem1
         nodef(4)=nodem
         nodef(5)=node2
         nodef(6)=node2
!
!        Sets the types of the degrees of freedom
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=1
         idirf(5)=2
         idirf(6)=0
!
!        Calculate residual or derivatives
         if(ider.eq.0) then
!           Residual
            f=calc_residual_tee(pt1,Tt1,xflow1,xflow2,pt2,
     &    Tt2,A1,A2,zeta_fac,kappa,R,ider,kflag,zeta)
         else
!           Derivatives
            call calc_ider_tee(df,pt1,Tt1,xflow1,xflow2,pt2,
     &    Tt2,A1,A2,zeta_fac,kappa,R,ider,kflag,zeta)
         endif
      elseif(kflag.eq.3) then
!
!        Element output
!
!        Calculate kappa
         kappa=(cp/(cp-R))
!
!        setting icase (always adiabatic)
         icase=0;
!
!        Inlet conditions are the same for both branches
!
!        Determining the previous element as
!        the incoming mass flow is defined there
         nelem1=nint(prop(index+1))
         nodem1=kon(ipkon(nelem1)+2)
!
!        Inlet conditions
         pt1=v(2,node1)
         Tt1=v(0,node1)-physcon(1)
         xflow1=v(1,nodem1)*iaxial
         A1=prop(index+4)
!
!        Outlet conditions
         Tt2=v(0,node2)
         xflow2=v(1,nodem)*iaxial
         pt2=v(2,node2)
         if(nelem.eq.nint(prop(index+2))) then
            A2=prop(index+5)
            chan_num=1
            zeta_fac=prop(index+11)
         elseif(nelem.eq.nint(prop(index+3))) then
            A2=prop(index+6)
            chan_num=2
            zeta_fac=prop(index+12)
         endif
!         zeta_fac = prop(index+11)
!
!        Write the main information about the element
!

!        Flow velocity at inlet
         call ts_calc(xflow1,Tt1,pt1,kappa,r,A1,Ts0,icase)
         pspt0=(Ts0/Tt1)**(kappa/(kappa-1))
!        Calculate Mach numbers
         call machpi(M1,pspt0,kappa,R)
         call ts_calc(xflow2,Tt2,pt2,kappa,r,A2,Ts2,icase)
!        Pressure ratio
         pspt2=(Ts2/Tt2)**(kappa/(kappa-1))
         call machpi(M2,pspt2,kappa,R)
!
         write(1,*) ''
         write(1,55) ' from node ',node1,
     &        ' to node ', node2,':   air massflow rate=' ,xflow
!     
         write(1,56)'       Inlet node ',node1,':    Tt1= ',Tt1,
     &        ' , Ts1= ',Ts0,' , Pt1= ',pt1,
     &        ', M1= ',M1
         write(1,*)'             Element ',nelem,lakon(nelem)
     &        ,', Branch ',chan_num
!     
 55      format(1x,a,i6,a,i6,a,e11.4,a)
 56      format(1x,a,i6,a,e11.4,a,e11.4,a,e11.4,a,e11.4)
!     
!     Set ider to calculate the residual
         ider=0
!     
!     Calculate the element one last time with enabled output
         f=calc_residual_tee(pt1,Tt1,xflow1,xflow2,pt2,
     &        Tt2,A1,A2,zeta_fac,kappa,R,ider,kflag,zeta)
!     
         write(1,56)'      Outlet node ',node2,':   Tt2= ',Tt2,
     &        ' , Ts2= ',Ts2,' , Pt2= ',pt2,
     &        ', M2= ',M2
!
      endif
!     
      xflow=xflow/iaxial
      df(3)=df(3)*iaxial
      df(4)=df(4)*iaxial
!     
      return
      end
