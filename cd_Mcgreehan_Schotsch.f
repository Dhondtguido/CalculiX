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
!
!     this subroutine enables to calculate the basis incompressible
!     discharge coefficient 
!
!     "Flow Characteristics of long orifices with rotation and corner radiusing"
!     W.F. Mcgreehan and M.J. Schotsch
!     ASME 87-GT-162
!
!     author: Yannick Muller
!      
      subroutine cd_Mcgreehan_Schotsch(rzdh,bdh,reynolds,cdu)
!     
      implicit none
!     
      real*8 cdu,bdh,reynolds,cd_re,rzdh,cd_r
!     
      cd_re=0.5885d0+372d0/reynolds
!
!     the radius correction 
!
      cd_r=1-(0.008d0+0.992d0*exp(-5.5d0*rzdh-3.5d0*rzdh**2))
     &     *(1-cd_re)
!   
      cdu=1.d0-(1.d0-cd_r)*(1d0+1.3d0*exp(-1.606d0*(bdh*bdh)))
     &     *(0.435d0+0.021d0*bdh)
!     
      return 
!     
      end
      
