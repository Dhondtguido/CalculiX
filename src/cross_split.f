!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine cross_split(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,r,physcon,numf,set,mi,ider,ttime,time,
     &     iaxial,iplausi)
!
!     A cross split element(zeta calculation according to Idel'chik)
!     Written by Yavor Dobrev
!     For an explanation of the parameters see tee.f
!
      implicit none
!
      logical identity
!
      character*8 lakon(*)
      character*81 set(*)
!
      integer
     &nelem,nactdog(0:3,*),node1,node2,nodem,numf,kon(*),ipkon(*),
     &ielprop(*),nodef(*),idirf(*),kflag,inv,mi(2),nodem1,nelem1,
     &ichan_num,ider,icase,iaxial,index,iplausi
!
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(*),kappa,R,Tt1,Tt2,pt1,
     &pt2,cp,physcon(*),km1,kp1,kdkm1,pt2pt1,pt1pt2,pt2pt1_crit,
     &tdkp1,A,A1,A2,calc_residual_cross_split,dh1,dh2,alpha,xflow1,
     &xflow2,zeta_fac,A_s,Ts0,pspt0,pspt2,M1,M2,Ts2,ttime,time
!
!
!
      index=ielprop(nelem)
!
      if(kflag.eq.0) then
         identity=.true.
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
         kappa=(cp/(cp-R))
         kp1=kappa+1d0
         km1=kappa-1d0
         kdkm1=kappa/km1
         tdkp1=2.d0/kp1
!
         if(nelem.eq.nint(prop(index+2))) then
            A=prop(index+5)
         elseif(nelem.eq.nint(prop(index+3)))then
            A=prop(index+6)
         elseif(nelem.eq.nint(prop(index+4)))then
            A=prop(index+6)
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
         pt1pt2=pt1/pt2
         pt2pt1=1/pt1pt2
!
         pt2pt1_crit=tdkp1**kdkm1
!
         if(pt2pt1.gt.pt2pt1_crit) then
            xflow=inv*pt1*A*dsqrt(2.d0*kdkm1*pt2pt1**(2.d0/kappa)
     &           *(1.d0-pt2pt1**(1.d0/kdkm1))/r)/dsqrt(Tt1)
         else
            xflow=inv*pt1*A*dsqrt(kappa/r)*tdkp1**(kp1/(2.d0*km1))/
     &           dsqrt(Tt1)
         endif
!     
      elseif(kflag.eq.2)then
!     
         numf=6
!
         kappa=(cp/(cp-R))
!
!     setting icase (always adiabatic)
!     
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
         A1=prop(index+5)
         dh1=prop(index+9)
!
!        Set the node numbers for the degrees of freedom
         nodef(1)=node1
         nodef(2)=node1
         nodef(3)=nodem1
         nodef(4)=nodem
         nodef(5)=node2
         nodef(6)=node2
!
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=1
         idirf(5)=2
         idirf(6)=0
!
         if(nelem.eq.nint(prop(index+2))) then
            ichan_num=1
!
            Tt2=v(0,node2)
            xflow2=v(1,nodem)*iaxial
            pt2=v(2,node2)
            A2=A1
            A_s=prop(index+6)
            dh2=dh1
            alpha=0
            zeta_fac=prop(index+11)
!
         elseif(nelem.eq.nint(prop(index+3))) then
            ichan_num=2
!
            Tt2=v(0,node2)
            xflow2=v(1,nodem)*iaxial
            pt2=v(2,node2)
            A2=prop(index+6)
            dh2=prop(index+10)
            alpha=prop(index+7)
            zeta_fac=prop(index+12)
!
         elseif(nelem.eq.nint(prop(index+4))) then
            ichan_num=3
!
            Tt2=v(0,node2)
            xflow2=v(1,nodem)*iaxial
            pt2=v(2,node2)
            A1=prop(index+5)
            A2=prop(index+6)
            dh2=prop(index+10)
            alpha=prop(index+8)
            zeta_fac=prop(index+12)
!
         endif
!
         if(ider.eq.0)then
!           Residual
            f=calc_residual_cross_split(pt1,Tt1,xflow1,xflow2,pt2,
     &    Tt2,ichan_num,A1,A2,A_s,dh1,dh2,alpha,zeta_fac,
     &    kappa,R,ider,kflag)
         else
!           Derivatives
            call calc_ider_cross_split(df,pt1,Tt1,xflow1,xflow2,pt2,
     &    Tt2,ichan_num,A1,A2,A_s,dh1,dh2,alpha,zeta_fac,
     &    kappa,R,ider,kflag)
         endif
!
      elseif(kflag.eq.3) then
!
         kappa=(cp/(cp-R))
!
!         setting icase (always adiabatic)
!
         icase=0;
!
!        Inlet conditions are the same for both branches
!
!        Determining the previous element as
!        the incoming mass flow is defined there
!
         nelem1=nint(prop(index+1))
         nodem1=kon(ipkon(nelem1)+2)
!
!        Inlet conditions
!
         pt1=v(2,node1)
         Tt1=v(0,node1)-physcon(1)
         xflow1=v(1,nodem1)*iaxial
         A1=prop(index+5)
         dh1=prop(index+9)      
!
         if(nelem.eq.nint(prop(index+2))) then
            ichan_num=1
!
            Tt2=v(0,node2)
            xflow2=v(1,nodem)*iaxial
            pt2=v(2,node2)
            A2=A1
            A_s=prop(index+6)
            dh2=dh1
            alpha=0
            zeta_fac=prop(index+11)
!
         elseif(nelem.eq.nint(prop(index+3))) then
            ichan_num=2
!
            Tt2=v(0,node2)
            xflow2=v(1,nodem)*iaxial
            pt2=v(2,node2)
            A2=prop(index+6)
            dh2=prop(index+10)
            alpha=prop(index+7)
            zeta_fac=prop(index+12)
!
         elseif(nelem.eq.nint(prop(index+4))) then
            ichan_num=3
!
            Tt2=v(0,node2)
            xflow2=v(1,nodem)*iaxial
            pt2=v(2,node2)
            A1=prop(index+5)
            A2=prop(index+6)
            dh2=prop(index+10)
            alpha=prop(index+8)
            zeta_fac=prop(index+12)
!
         endif
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
     &        ' to node ', node2,':   air massflow rate= ',xflow
!     
         write(1,56)'       Inlet node ',node1,':    Tt1= ',Tt1,
     &        ', Ts1= ',Ts0,', Pt1= ',pt1,
     &        ', M1= ',M1
         write(1,*)'             Element ',nelem,lakon(nelem)
     &        ,', Branch ',ichan_num
!     
 55      format(1x,a,i6,a,i6,a,e11.4,a)
 56      format(1x,a,i6,a,e11.4,a,e11.4,a,e11.4,a,e11.4)
!     
!     Set ider to calculate the residual
         ider=0
!     
!     Calculate the element one last time with enabled output
         f=calc_residual_cross_split(pt1,Tt1,xflow1,xflow2,pt2,
     &        Tt2,ichan_num,A1,A2,A_s,dh1,dh2,alpha,zeta_fac,
     &        kappa,R,ider,kflag)
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
      
