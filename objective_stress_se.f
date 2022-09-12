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
      subroutine objective_stress_se(nk,iobject,mi,dstn,objectset,
     &  ialnneigh,naneigh,nbneigh,stn,dksper)
!
!     calculates the sum of the square of the von Mises stress of a node
!     set
!
      implicit none
!
      character*81 objectset(5,*)
!
      integer nk,idir,iobject,mi(*),j,k,ialnneigh(*),naneigh,nbneigh
!
      real*8 stn(6,*),p,rho,xstress,dksper,dstn(6,*),mises,dmises 
!
!     reading rho and the mean stress for the Kreisselmeier-Steinhauser
!     function
!
      read(objectset(2,iobject)(41:60),'(f20.0)') rho
      read(objectset(2,iobject)(61:80),'(f20.0)') xstress
!
      dksper=0.d0
      do j=naneigh,nbneigh        
         k=ialnneigh(j)
!
!        Calculate unperturbed mises stress
!
         p=-(stn(1,k)+stn(2,k)+stn(3,k))/3.d0
         mises=dsqrt(1.5d0*((stn(1,k)+p)**2+
     &                      (stn(2,k)+p)**2+
     &                      (stn(3,k)+p)**2+
     &                 2.d0*(stn(4,k)**2+stn(5,k)**2+
     &                       stn(6,k)**2)))
!
!        Calculate perturbed mises stress
!
         p=-(dstn(1,k)+dstn(2,k)+dstn(3,k))/3.d0  
         dmises=dsqrt(1.5d0*((dstn(1,k)+p)**2+
     &                       (dstn(2,k)+p)**2+
     &                       (dstn(3,k)+p)**2+
     &                  2.d0*(dstn(4,k)**2+dstn(5,k)**2+
     &                        dstn(6,k)**2)))
!
!        Calculate delta mises stress
!
         dmises=dmises-mises
!
         dksper=dksper+dexp(rho*dsqrt(1.5d0*
     &   ((dstn(1,k)+p)**2+
     &   (dstn(2,k)+p)**2+
     &   (dstn(3,k)+p)**2+
     &   2.d0*(dstn(4,k)**2+
     &         dstn(5,k)**2+
     &         dstn(6,k)**2)))/xstress)*dmises
      enddo
!
      dksper=dksper/xstress
!
      return
      end

