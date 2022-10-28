!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine randomval(randval,nev)              
!
!     Creation of normal ditributed unit variance gaussian       
!     random variables using the Box-Muller-Transformation
!
      implicit none
!
      integer nev,i
!
      real*8 randval(*),fac,v1,v2,rsq
!      
      call random_seed()
!
      do i=1,nev
         do
            call random_number(v1)
            call random_number(v2)
            rsq=v1**2+v2**2
            if((rsq.ge.1.d0).or.(rsq.le.0.d0)) cycle
            fac=sqrt(-2.d0*dlog(rsq)/rsq)  
            randval(i)=v1*fac
            exit
         enddo
      enddo
!
      return        
      end




