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
!cd incompressible for long orifices (eq.6)
!
!     author: Yannick Muller
!
       subroutine pk_cdi_rl(lqd,rqd,reynolds,beta,cdi_rl)
!
       implicit none
!
      real*8 lqd,rqd,reynolds,beta,cdi_rl,rqd_cor,lrqd,cdi_r,glrqd
!
      rqd_cor=rqd
!
      if (rqd_cor.gt.lqd) then
         rqd_cor=lqd
      endif
!
      lrqd=lqd-rqd_cor
!      
      call pk_cdi_r(rqd_cor,reynolds,beta,cdi_r)
!      
      glrqd=(1d0+1.298d0*exp(-1.593d0*lrqd**2.33d0))
     &     *(0.435d0+0.021d0*lrqd)/(2.298d0*0.435d0)
!
      cdi_rl=1.d0-glrqd*(1.d0-cdi_r)
!
      return
!
      end
