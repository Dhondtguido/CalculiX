!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine cd_ms_ms(p1,p2,T1,rad,d,xl,kappa,r,reynolds,u,vid,cd)
!
!     This subroutine enables to calculate the discharge coefficient for an 
!     orifice (shap edged , rotating..) following the results obtained 
!     by Mcgreehan and Schotsch
!     The decription of the method can be found in :
!     "Flow characteristics of long orifices with rotation and 
!     corner radiusing"
!     ASME 87-GT-162
!
!     author: Yannick Muller
!
      implicit none 
!
      real*8 p1,p2,T1,rad,d,xl,kappa,r,reynolds,u,cd,qlim,q,
     &     c1,c2,c3,fakt,aux,rzd,lkorr,qkorr,rv,vid
!
      qlim=10.d0
!
!     taking in account the influence of the Reynolds number
!
      cd=0.5885d0+372.d0/reynolds
      cd=min(cd,1.d0)
!
!     taking in account the edge radius
!
      rzd=rad/d
      aux=exp(-(3.5d0*rzd+5.5d0)*rzd)
      fakt=aux+0.008d0*(1.d0-aux)
      cd=1.d0-fakt*(1.d0-cd)
      cd=min(max(cd,0.d0),1.d0)
!
!     taking in account the lenght of the orifice
!
      lkorr=xl-rad
      q=lkorr/d
      qkorr=min(q,qlim)
      fakt=(1.d0+1.3d0*exp(-1.606d0*qkorr**2.d0))*
     &     (0.435d0+0.021d0*qkorr)/(2.3d0*0.435d0)
      cd=1.d0-fakt*(1.d0-cd)
      cd=min(max(cd,0.d0),1.d0)
!
!     taking in account the tangential velocity
!
      if(u.ne.0d0) then
         vid=dsqrt(2.d0*kappa/(kappa-1.d0)*r*T1*
     &        (1.d0-(p2/p1)**((kappa-1.d0)/kappa)))
         rv=u/vid*(0.6d0/cd)**3
         c1=exp(-rv**1.2d0)
         c2=0.5d0*rv**0.6d0*dsqrt(0.6d0/cd)
         c3=exp(-0.5d0*rv**0.9d0)
         cd=cd*(c1+c2*c3)
         cd=min(max(cd,0.d0),1.d0)
!
      endif
!
!     
      return
      end
