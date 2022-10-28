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
      subroutine cd_pk_ms(rad,d,xl,reynolds,p2,p1,beta,kappa,cd,u,
     &     T1,R)
!
!     This subroutines enable to calculate the compressible discharge 
!     coefficient for thin and long orifices with corner radiusing;
!
!     author: Yannick Muller
!     
      implicit none
!
      real*8 rad,d,xl,lqd,rqd,reynolds,p2,p1,p2p1,beta,beta_cor,kappa,
     &     cd,cdc_cl1,cdc_cl3,rldb,R,u,T1,c1,c2,
     &     c3,ms_cdr,rv,vid
!
      p2p1=p2/p1
      rqd=rad/d
      lqd=xl/d
      rldb=max(lqd,0.d0)
!
!     the method of cd calculation for a sharp edged aperture is only valid 
!     for beta comprised between 0 and 0.7 
!      
      if (beta.gt.0.7d0) then
         beta_cor=0.7d0
      else
         beta_cor=beta
      endif
!
!     differences between class1 or class2 or class3
!
      if (lqd.eq.rqd) then
!     
!     class1
!
         call pk_cdc_cl1(lqd,reynolds,p2p1,beta_cor,kappa,cdc_cl1)
         cd=cdc_cl1
      else
!     
!     class2 or class3 (clas2 is a sub class of class3 )
!     
         call pk_cdc_cl3(lqd,rqd,reynolds,p2p1,beta_cor,kappa,cdc_cl3)
         cd=cdc_cl3
      endif
!      
!     if rotating orifice with Mac Greehan & Scotch
!     The decription of the method can be found in :
!     "Flow characteristics of long orifices with rotation and 
!     corner radiusing" ASME 87-GT-16
!
!     rotating case eq 17
      
      if (u.ne.0) then
         vid=dsqrt(2.d0*kappa/(kappa-1.d0)*r*T1*
     &        (1.d0-(p2/p1)**((kappa-1.d0)/kappa)))
         rv=u/vid*(cd/0.6d0)**(-3d0)
         c1=exp(-rv**1.2d0)
         c2=0.5d0*rv**0.6d0*dsqrt(0.6d0/cd)
         c3=exp(-0.5d0*rv**0.9d0)
         ms_cdr=cd*(c1+c2*c3)
         cd=ms_cdr
         cd=min(max(cd,0.d0),1.d0)
      endif
!
      return
      end
      
