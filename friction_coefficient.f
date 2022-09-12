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
!     This subroutine computes the friction coefficient of
!     the pipe flow for laminar and turbulent flow including the
!     transition region
!
!     author: Yannick Muller
!     
      subroutine friction_coefficient(l,d,ks,reynolds,form_fact,lambda)
!
!     INPUT:      
!
!     l           not used in internet-version of CalculiX
!     d           hydraulic diameter=4*area/(wetted perimeter)
!     ks          grain size
!     reynolds    Reynolds number 
!     form_fact   form factor for non-round sections and laminar flow
!                 (reference underneath)
!
!     INPUT:      
!
!     lambda: friction coefficient
!      
      implicit none
!     
      real*8 l,d,ks,reynolds,form_fact,lambda,alfa2,
     &     rey_turb_min,rey_lam_max,lzd,dd,ds,friction,dfriction,
     &     lambda_kr,lambda_turb,ksd
!
      rey_turb_min=4000
      rey_lam_max=2000
      lzd=l/d
      ksd=ks/d
!     
!     transition laminar turbulent domain
!     
      if((reynolds.gt.rey_lam_max).and.(reynolds.lt.rey_turb_min))then
!     
         lambda_kr=64.d0/rey_lam_max
!     
!     Solving the implicit White-Colebrook equation
!     1/dsqrt(friction)=-2*log10(2.51/(Reynolds*dsqrt(friction)+0.27*Ks))
!     
!     Using Haaland explicit relationship for the initial friction value
!     S.E. Haaland 1983 (Source en.Wikipwedia.org)
!     
         friction=(-1.8d0*dlog10(6.9d0/4000.d0+(ksd/3.7d0)**1.11d0))
     &        **(-2d0)
!     
         do
            ds=dsqrt(friction)
            dd=2.51d0/(4000.d0*ds)+0.27d0*ksd
            dfriction=(1.d0/ds+2.d0*dlog10(dd))*2.d0*friction*ds/
     &           (1.d0+2.51d0/(4000.d0*dd))
            if(dfriction.le.friction*1.d-3) then
               friction=friction+dfriction
               exit
            endif
            friction=friction+dfriction
         enddo
         lambda_turb=friction
         
!     
!     logarithmic interpolation in the trans laminar turbulent domain
!     
         lambda=lambda_kr*(lambda_turb/lambda_kr)
     &       **(log(reynolds/rey_lam_max)/log(rey_turb_min/rey_lam_max))
!     
!     laminar flow
!     using Couette-Poiseuille formula
!     the form factor for non round section can be found in works such as
!     Bohl,W
!     "Technische Strömungslehre Stoffeigenschaften von Flüssigkeiten und
!     Gasen, hydrostatik,aerostatik,incompressible Strömungen,
!     Strömungsmesstechnik
!     Vogel Würzburg Verlag 1980
!     
      elseif(reynolds.lt.rey_lam_max) then
         lambda=64.d0/reynolds
         lambda=form_fact*lambda
!     
!     turbulent
!     
      else
!     
!     Solving the implicit White-Colebrook equation
!     1/dsqrt(friction)=-2*log10(2.51/(Reynolds*dsqrt(friction)+0.27*Ks))
!     
!     Using Haaland explicit relationship for the initial friction value
!     S.E. Haaland 1983 (Source en.Wikipwedia.org)
!     
         friction=(-1.8d0*dlog10(6.9d0/reynolds+(ksd/3.7d0)
     &        **1.11d0))**(-2d0)
!     
         do
            ds=dsqrt(friction)
            dd=2.51d0/(reynolds*ds)+0.27d0*ksd
            dfriction=(1.d0/ds+2.d0*dlog10(dd))*2.d0*friction*ds/
     &           (1.d0+2.51d0/(reynolds*dd))
            if(dfriction.le.friction*1.d-3) then
               friction=friction+dfriction
               exit
            endif
            friction=friction+dfriction
         enddo
         lambda=friction
      endif
!     
      call interpol_alfa2(lzd,reynolds,alfa2)
!     
      return
!     
      end
      
      
      
