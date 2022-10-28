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
!     cd compresibble for class 3 orifices where l/d>0 and r/d>0
!
!     author: Yannick Muller
!
      subroutine pk_cdc_cl3(lqd,rqd,reynolds,p2p1,beta,kappa,cdc_cl3)
!     
      implicit none
!
      real*8  lqd,rqd,reynolds,p2p1,beta,kappa,cdc_cl3a,cdc_cl3b,
     &     cdc_cl3d,cdc_cl3
!     
      cdc_cl3a=0.d0
      cdc_cl3b=0.d0
      cdc_cl3d=0.d0
!
      if(lqd.le.0.28d0) then
         call pk_cdc_cl3a(lqd,rqd,reynolds,p2p1,beta,kappa,cdc_cl3a)
         cdc_cl3=cdc_cl3a
      elseif(lqd.le.0.5d0) then
         call pk_cdc_cl3b(lqd,rqd,reynolds,p2p1,beta,kappa,cdc_cl3b)
         cdc_cl3=cdc_cl3b
      else 
         call pk_cdc_cl3d(lqd,rqd,reynolds,p2p1,beta,cdc_cl3d)
         cdc_cl3=cdc_cl3d
!     
      endif
!     
      return
! 
      end
