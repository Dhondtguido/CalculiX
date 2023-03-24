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
! cd inncompressible fro thin orifices with corner radiusing (eq 5)
!
!     author: Yannick Muller
!
      subroutine pk_cdi_r (rqd,reynolds,beta,cdi_r)
!
      implicit none
!
      real*8 rqd,reynolds,beta,cdi_r,frqd,cdi_se,cdi_noz
!      
      call pk_cdi_noz(reynolds,cdi_noz)
      call pk_cdi_se(reynolds,beta,cdi_se)
      
      frqd=0.008d0+0.992d0*exp(-5.5d0*rqd-3.5d0*rqd**2.d0)
!
      cdi_r=cdi_noz-frqd*(cdi_noz-cdi_se)
!      
      return 
      end
