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
!     author: Yannick Muller
!     
      real*8 function calc_residual_wye(pt1,Tt1,xflow1,xflow2,
     &     pt2,Tt2,ichan_num,A1,A2,A_s,dh1,dh2,alpha,zeta_fac,kappa,
     &     R,ider,iflag,zeta)
!     
      implicit none
!     
      integer icase,ichan_num,icrit1,icrit2,ider,iflag
!     
      real*8 
!     In- and Output
     &     f,
     &     R,
!     
!     Kappa stuff
     &     kappa,
!     
     &     pt1,
     &     pt2,
     &     Tt1,
     &     Tt2,
     &     xflow1,
     &     xflow2,
!     
     &     zeta,
!     
     &     A1,
     &     A2,
     &     A_s,
!     
     &     Ts1,
     &     Ts2,
     &     Ts_s,
     &     dh1,
     &     dh2,
     &     alpha,
     &     Q_crit,
     &     pspt_crit,
     &     Q0,
     &     Q1,
     &     Q2,
     &     pspt0,
     &     pspt1,
     &     w1,
     &     w2,
     &     w1w2,
     &     w2w1,
     &     VsV1,
     &     pi,
     &     z2d390,
     &     z1p090,
     &     z60,
     &     z90,
     &     hq,
     &     Ts0,
     &     zeta_fac,
     &     M1,
     &     M2,
     &     pspt2
!     
      real*8 Table_A(2,11)
!     
      icrit1 = 0
      icrit2 = 0
!     
!     setting icase (always adiabatic)
      icase=0;
!     
      pi=4.d0*datan(1.d0)
      Q_crit = dsqrt(kappa/R)*
     &     (1+0.5d0*(kappa-1))**(-0.5d0*(kappa+1)/(kappa-1))
      pspt_crit = (2/(kappa+1)) ** (kappa/(kappa-1))
!     
      Q0 = xflow1*dsqrt(Tt1)/pt1/A1
      Q1 = xflow2*dsqrt(Tt1)/pt1/A2
      if(Q1.ge.Q_crit) then
        Q1 = Q_crit
        icrit1 = 1
        write(*,*)'*WARNING in Wye:'
        write(*,*)'Critical conditions at 1'
      endif
      Q2 = xflow2*dsqrt(Tt1)/pt2/A2
      if(Q2.ge.Q_crit) then
        Q2 = Q_crit
        icrit2 = 1
        write(*,*)'*WARNING in Wye:'
        write(*,*)'Critical conditions at 2'
      endif
!     
!     Flow velocity at inlet
      Ts0=Tt1
      
      call ts_calc(xflow1,Tt1,pt1,kappa,r,A1,Ts0,icase)
      pspt0 = (Ts0/Tt1)**(kappa/(kappa-1))
      call wpi(w1, pspt0, Q0, 
     &     dsqrt(Tt1),kappa,R) 
!     
!     Flow velocity at outlet
      call ts_calc(xflow2,Tt1,pt1,kappa,r,A2,Ts1,icase)
      pspt1 = (Ts1/Tt1)**(kappa/(kappa-1))
      call wpi(w2, pspt1, Q1, 
     &     dsqrt(Tt1),kappa,R) 
!     
      if(w2.eq.0.d0)then
        w1w2=1d30
      else
        w1w2=w1/w2
      endif
      if(w1.eq.0.d0)then
        w2w1=1d30
      else
        w2w1=w2/w1
      endif
!     
!     Zeta calculation
!     Main branch
      if(ichan_num.eq.1) then
!     Zeta as in Calculix and old ACC tool
        zeta=0.4d0*(1-W2W1)**2
        zeta=zeta*(W1W2)**2
!     
!     Branch
      elseif(ichan_num.eq.2) then            
        hq=dh2/dh1
!     
!     Interpolation as in CalculiX
        if((alpha.le.60.and.hq.le.1.d0).or.
     &       (alpha.le.90.d0.and.hq.le.2.d0/3.d0)) then
          zeta=0.95d0*((W2W1-2d0*dcos(alpha*pi/180))
     &         *W2W1+1.d0)
          zeta=zeta*(W1W2)**2
        elseif(alpha.le.90.d0.and.hq.le.1.d0) then
          z2d390=0.95d0*((W2W1-2d0*dcos(90.d0*pi/180))
     &         *W2W1+1.d0)
!     
          z1p090=0.95d0*(1.d0+0.3d0*W2W1**2)     
!     
          z90=z2d390+(3*hq-2.d0)*(z1p090-z2d390)
!     
          Z60=0.95d0*((W2W1-2d0*dcos(60.d0*pi/180))
     &         *W2W1+1.d0)
!     
          zeta=z60+(alpha/30.d0-2.d0)*(z90-z60)
          zeta=zeta*(W1W2)**2
        elseif(alpha.le.60.d0.and.hq.gt.1.d0) then
!     
!     extrapolation for hq>1 (not part of Idelchik
!     but the previous definition results sometimes 
!     into zeta<0 values
!     
          write(*,*)'WARNING in calc_residual_wye:'
          write(*,*)'   Branch element is 
     &outside valid range defined by Idelchik'
          zeta=0.95d0*((W2W1-2d0*dcos(alpha*pi/180))
     &         *W2W1+1.d0)
          zeta=zeta*(W1W2)**2 
        elseif(alpha.le.90.d0.and.hq.gt.1.d0) then
!     
!     extrapolation for hq>1 (not part of Idelchik
!     but the previous definition results sometimes 
!     into zeta<0 values
!     
          write(*,*)'WARNING in calc_residual_wye:'
          write(*,*)'   Branch element is 
     &outside valid range defined by Idelchik'
          Z60=0.95d0*((W2W1-2d0*dcos(60.d0*pi/180))
     &         *W2W1+1.d0)
!     
          z1p090=0.95d0*(1.d0+0.3d0*W2W1**2)    
!     
          zeta=Z60+(z1p090-Z60)/(90.d0-60.d0)*(alpha-60)     
!     
          zeta=zeta*(W1W2)**2
!     interpolation between exit loss zeta=1
!     and last point of definition range (hq=1)
          zeta=zeta+(1-zeta)/(50.d0-1.d0)*(hq-1)          
        else
!     
!     values are absolutely out of def. range 
!     
          write(*,*)'ERROR in wye.f:' 
          write(*,*)'   Branch element is 
     &outside valid range'
          call exit(201)  
        endif
      endif
!     
      zeta = zeta*zeta_fac
!     
      if(icrit2.ne.1) then
        if(icrit1.ne.1) then
          f = pt2 - pt1*pspt1**zeta
        else
          f = xflow2*dsqrt(Tt1)/pt1/A2-Q_crit
        endif
      else
        f = xflow2*dsqrt(Tt1)/pt2/A2-Q_crit
      endif
!     
      calc_residual_wye=f

      if(iflag.eq.4) then
!     Calculate Mach numbers
        call machpi(M1,pspt0,kappa,R)
        call ts_calc(xflow2,Tt2,pt2,kappa,r,A2,Ts2,icase)
!     Pressure ratio
        pspt2 = (Ts2/Tt2)**(kappa/(kappa-1))
        call machpi(M2,pspt2,kappa,R)
!        
 80     format(3x,a,f10.6,a,f10.2,a,f10.6)
 77     format(3x,a,f10.6,a,f10.2,a,f10.6)
      endif
!     
      return
      end
