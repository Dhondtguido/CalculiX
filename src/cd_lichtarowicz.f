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
!     This subroutines enables to calculate the reynolds number correction after:
!     "Discharge coeffcients for incompressible non-cavitating flowthrough long orifices"
!     A. Lichtarowicz, R.K duggins and E. Markland
!     Journal  Mechanical Engineering Science , vol 7, No. 2, 1965
!
!     author: Yannick Muller
!
      subroutine cd_lichtarowicz(cd,cdu,reynolds,amod,bdh) 
!
      implicit none
!
      real*8 cdu,reynolds,amod,bdh,eps,A1,cd_diff,cd0,cd
!
      cd0=cdu
      cd_diff=1.d0
!
      do
!        
         if(cd_diff.lt.1.d-3) exit
!
         cd=cd0
         A1=20/(reynolds*dsqrt(1.d0-Amod**2))*(1.d0+2.25d0*bdh)
         eps=(0.005d0*bdh)/(1.d0+7.5d0*(log10(0.00015d0*reynolds*
     &        dsqrt(1.d0-Amod**2)/cd))**2)
               
         cd=((-1/cdu+eps)+dsqrt((1/cdu-eps)**2.d0+4.d0*A1))/(2*A1)
!
         cd_diff=dabs(cd-cd0)
!
         cd0=cd
!
       enddo
!     
         return
         end
