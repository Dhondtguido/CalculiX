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
!     author: Yannick Muller
!
      real*8 function calc_residual_cross_split(pt1,Tt1,xflow1,xflow2,
     &pt2,Tt2,ichan_num,A1,A2,A_s,dh1,dh2,alpha,zeta_fac,
     &kappa,R,ider,iflag)
!
      implicit none
!
      integer icase,ichan_num,ider,icrit1,icrit2,iflag,ier
!
      real*8 
!     In- and Output
     &f,R,
!
!     Kappa stuff
     &kappa,km1,kp1,
!
     &pt1,pt2,Tt1,Tt2,xflow1,xflow2,
!
     &pt2_lim,
!
     &zeta,
!
     &A1,A2,
!
     &Ts0,Ts1,Ts2,dh1,dh2,alpha,Q_crit,pspt_crit,Q0,Q1,Q2,pspt0,
     &pspt1,pspt2,w1,w2,w1w2,w2w1,pi,z2d390,z1p090,z60,z90,hq,M1,M2,
     &zeta_fac,xflow_s,Q_s,Ts_s,pspt_s,w_s,wsw1,A_s,AsA1,VsV1
!
      real*8 Table_zeta(2,10)
!
      pi=4.d0*datan(1.d0)
!
      icrit1=0
      icrit2=0
!
!     setting icase (always adiabatic)
!     
      icase=0;
!
      km1=kappa-1.d0
      kp1=kappa+1.d0
      Q_crit=dsqrt(kappa/R)*
     &   (1+0.5d0*(kappa-1))**(-0.5d0*(kappa+1)/(kappa-1))
      pspt_crit=(2.d0/(KAPPA+1.d0))**(KAPPA/(KAPPA-1.d0))
!
      Q0=xflow1*dsqrt(Tt1)/pt1/A1
      Q1=xflow2*dsqrt(Tt1)/pt1/A2
      if(Q1.ge.Q_crit) then
         Q1=Q_crit
         icrit1=1
         write(*,*)'*WARNING in Cross Split:'
         write(*,*)'Critical conditions at 1'
      endif
      Q2=xflow2*dsqrt(Tt1)/pt2/A2
      if(Q2.ge.Q_crit) then
         Q2=Q_crit
         icrit2=1
         write(*,*)'*WARNING in Cross Split:'
         write(*,*)'Critical conditions at 2'
      endif
!
!     Flow velocity at inlet
      call ts_calc(xflow1,Tt1,pt1,kappa,r,A1,Ts0,icase)
      pspt0=(Ts0/Tt1)**(kappa/(kappa-1))
      call wpi(w1, pspt0, Q0, 
     &      dsqrt(Tt1),kappa,R)
!
!     Flow velocity at outlet
      call ts_calc(xflow2,Tt1,pt1,kappa,r,A2,Ts1,icase)
      pspt1=(Ts1/Tt1)**(kappa/(kappa-1))
      call wpi(w2, pspt1, Q1, 
     &      dsqrt(Tt2),kappa,R) 
!
      w2w1=w2/w1
      w1w2=w1/w2
!
!     Main branch
      if(ichan_num.eq.1) then
!
!          Zeta as in Calculix
         zeta=0.4d0*(1-W2W1)**2
!
         zeta=zeta*(W1W2)**2
!
!     First branch
      elseif((ichan_num.eq.2).or.(ichan_num.eq.3)) then            
         hq=dh2/dh1
         if(alpha.le.60.or.hq.le.2.d0/3.d0) then
            zeta=0.95d0*((W2W1-2d0*dcos(alpha*pi/180))
     &                 *W2W1+1.d0)
            zeta=zeta*(W1W2)**2
         else
            z2d390=0.95d0*((W2W1-2d0*dcos(90.d0*pi/180))
     &                 *W2W1+1.d0)
            z1p090=0.95d0*(0.34d0+W2W1**2)
            z90=z2d390+(3*hq-2.d0)*(z1p090-z2d390)
            Z60=0.95d0*((W2W1-2d0*dcos(60.d0*pi/180))
     &                 *W2W1+1.d0)
            zeta=z60+(alpha/30.d0-2.d0)*(z90-z60)
            zeta=zeta*(W1W2)**2
         endif
!
      endif
!
!     zeta_fac for side branches are all =1
!     main branch can be set by the user in ACC Designer
      zeta=zeta*zeta_fac
!
      if(icrit2.ne.1) then
         if(icrit1.ne.1) then
            f=pt2-pt1*pspt1**zeta
         else
            f=xflow2*dsqrt(Tt1)/pt1/A2-Q_crit
         endif
      else
         f=xflow2*dsqrt(Tt1)/pt2/A2-Q_crit
      endif
!
      if(iflag.eq.3) then
!
         write(1,57)'             zeta= ',zeta
 57      format(1x,a,f9.4)
!
      else if (iflag.eq.4) then
!
!        Calculate Mach numbers
         call machpi(M1,pspt0,kappa,R)
         call ts_calc(xflow2,Tt2,pt2,kappa,r,A2,Ts2,icase)
!        Pressure ratio
         pspt2=(Ts2/Tt2)**(kappa/(kappa-1))
         call machpi(M2,pspt2,kappa,R)
      
         write(1,80)'Inlet: Tt1= ',Tt1,
     &              ', pt1= ',pt1,', M1= ',M1
     
         write(1,77)'mass flow = ',xflow2,', kappa = ',kappa,
     &              ', zeta= ',zeta

         write(1,80)'Outlet: Tt2= ',Tt2,
     &              ', pt2= ',pt2,', M2= ',M2
     
 80   format(3x,a,f10.6,a,f10.2,a,f10.6)
 77   format(3x,a,f10.6,a,f10.2,a,f10.6)

      endif

!
      calc_residual_cross_split=f
!
      return
      end
