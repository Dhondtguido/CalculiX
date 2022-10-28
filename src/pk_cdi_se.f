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
! cd incompressible for sharp edged orifices( eq.3)
!
!     author: Yannick Muller
!
      subroutine pk_cdi_se(reynolds,beta,cdi_se)
!
      implicit none
!
      real*8 reynolds,beta,cdi_se,reynolds_cor
!
      if(reynolds.eq.0d0) then
         reynolds_cor=1.d0
      else
         reynolds_cor=reynolds
      endif
!      
      cdi_se=0.5959d0+0.0312d0*beta**2.1d0-0.184d0*beta**8.d0
     &     +0.09d0*0.4333d0*beta**4.d0
     &     /(1.d0-beta**4.d0)-0.0337d0*0.47d0*beta**3.d0+91.71d0
     &     *(beta**1.75d0)/(reynolds_cor**0.75d0)
!      
      return
!
      end
