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
! cd compressible for class 3 orifices where l/d>0 and r/d>0
! type d) with 0.5<=l/d<=2 (eq. 27)
!
!     author: Yannick Muller
!
      subroutine pk_cdc_cl3d(lqd,rqd,reynolds,p2p1,beta,cdc_cl3d)
!
      implicit none
!
      real*8 lqd,rqd,reynolds,p2p1,beta,cdc_cl3d,cdi_rl,cdc_cl3_choked,
     &     jpsqpt,zeta
!
      cdc_cl3_choked=1.d0-(0.008d0+0.992d0*exp(-5.5d0*rqd
     &     -3.5d0*rqd**2.d0))*(1.d0-0.838d0)
!
      call pk_cdi_rl(lqd,rqd,reynolds,beta,cdi_rl)
!      
!     help function for eq 26
      if (p2p1.ge.1d0) then
         jpsqpt=1.d0
      elseif(p2p1.ge.0.1d0) then
         zeta=(1.d0-p2p1)/0.6d0
         jpsqpt=exp(-4.6d0*zeta**7d0-2.2d0*zeta**1.5d0)
      else
         jpsqpt=0.d0
      endif     
!
      cdc_cl3d=cdc_cl3_choked-jpsqpt*(cdc_cl3_choked-cdi_rl)
!
      return
!      
      end
