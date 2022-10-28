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
!cd compresssible for class3 orifices where l/d>0 and r/d>0
! typ b) with 0.28<l/d<0.5 eq 27
!
!     author: Yannick Muller
!
      subroutine pk_cdc_cl3b(lqd,rqd,reynolds,p2p1,beta,kappa,cdc_cl3b)
!
      implicit none
!
      real*8 lqd,rqd,reynolds,p2p1,beta,kappa,cdc_cl3b,cdi_rl,cdi_rl_b,
     &     cdi_rl1,cdi_rl2,cdc_cl3a,cdc_cl3d
!
!     eq.6
      call pk_cdi_rl(lqd,rqd,reynolds,beta,cdi_rl)
      cdi_rl_b=cdi_rl
!     eq.6 lqd=0.28
      call pk_cdi_rl(0.28d0,rqd,reynolds,beta,cdi_rl)
      cdi_rl1=cdi_rl
!     eq.6 lqd=0.5
      call pk_cdi_rl(0.5d0,rqd,reynolds,beta,cdi_rl)
      cdi_rl2=cdi_rl
!
!     as class 1 (class3 a) for lqd=0.28
      call pk_cdc_cl3a(0.28d0,rqd,reynolds,p2p1,beta,kappa,cdc_cl3a)
!     as class3 (class 3 d ) for lqd=0.5      
      call pk_cdc_cl3d(0.5d0,rqd,reynolds,p2p1,beta,cdc_cl3d)
!     eq 27 a (linear interpolation)
      cdc_cl3b=cdc_cl3a+(cdc_cl3d-cdc_cl3a)*(cdi_rl_b-cdi_rl1)
     &     /(cdi_rl2-cdi_rl1)
!      
      return
      end
