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
! cd compressible for class 3 orifice where, l/d>0 and r/d>0
!typ a) with 0 <= l/d<=0.28 (eq. 25 modified)
!
      subroutine pk_cdc_cl3a(lqd,rqd,reynolds,p2p1,beta,kappa,cdc_cl3a)
!
      implicit none
!      
      real*8 lqd,rqd,reynolds,p2p1,beta,kappa,cdc_cl3a,cdi_noz,
     &   cdi_rl,cdi_se,y0,yg,cdqcv_noz,cdqcv_rl
!
!     cd incompressible nozlle eq 4a 4b
      call pk_cdi_noz(reynolds,cdi_noz)
!     cd incompresible eq.6
      call pk_cdi_rl(lqd,rqd,reynolds,beta,cdi_rl)
!     cd incompressible sharp edge eq.3
      call pk_cdi_se(reynolds,beta,cdi_se)
!     y0,yg ,eq. 15-17, eq.18
      call pk_y0_yg(p2p1,beta,kappa,y0,yg)
!
      cdqcv_noz=cdi_noz/(0.0718d0*cdi_noz+0.9282d0)
      cdqcv_rl=cdi_rl/(0.0718d0*cdi_rl+0.9282d0)
!
!     eq.26 modified for class 3a
!
      cdc_cl3a=cdi_rl*((cdqcv_noz-cdqcv_rl)/(cdqcv_noz-cdi_se/0.971d0)
     &    *(y0/yg-1.d0)+1.d0)
!
      return
!
      end
