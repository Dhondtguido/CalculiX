!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine inversewavspd(xi,et,c,rho,a)
!
!     Calculates the propagation wave speed in a material, up to its 21 
!     constants. Subroutine for calcmatwavsps.f

!     Based on the procedure in:
!     C. Lane. The Development of a 2D Ultrasonic Array Inspection 
!     for Single Crystal Turbine Blades.
!     Switzerland: Springer International Publishing, 2014.
!
!      CARLO MONJARAZ TEC (CMT)
!
!       INPUT:
!       
!       xi,et: values within a square domain between -1 and 1.
!              correlate in a unique way with 0<=phi<=pi and
!              -pi<=theta<=pi
!
!       elas: c(3,3,3,3) - The elasticity vector
!             
!       rho: double - Density of the material
!        
!        
!       OUTPUT:
!       
!       a: inverse of the group wave speed d(omega)/d(k)
!
      implicit none
!     
      integer i,j,k,l,ier,matz,ndim
!
      real*8 c(3,3,3,3),rho,xn(3),cm(3,3,3),a,xi,et,
     &       cmm(3,3),dd,al(3),alz(3,3),fv1(3),fv2(3),
     &       theta,phi,pi,p3(3),v(3),
     &       speed
!     
      pi=4.d0*datan(1.d0)
!
      theta=xi*pi
      phi=(et+1.d0)*pi/2.d0
!     
      do l=1,3
         v(l)=0.d0
         do k=1,3
            cmm(k,l)=0.d0
            do j=1,3
               cm(j,l,k)=0.d0
            enddo
         enddo
      enddo            
!     
      xn(1)=dcos(theta)*dsin(phi)
      xn(2)=dsin(theta)*dsin(phi)
      xn(3)=dcos(phi)
!     
!     find wave speed for this direction
!     
      do l=1,3
         do k=1,3
            do i=1,3
               do j=1,3
                  cm(l,k,i)=cm(l,k,i)+c(l,k,j,i)*xn(j)
               enddo
            enddo        
         enddo
      enddo
!     
      do k=1,3
         do i=1,3
            do l=1,3
               cmm(k,i)=cmm(k,i)+cm(l,k,i)*xn(l)
            enddo
         enddo        
      enddo
!     
      ndim=3
      matz=1
      ier=0
!
!     initialize variables for routine rs
!
      do j=1,3
         al(j)=0.d0
         fv1(j)=0.d0
         fv2(j)=0.d0
         do i=1,3
            alz(j,i)=0.d0
         enddo
      enddo
!     
      call rs(ndim,ndim,cmm,al,matz,alz,fv1,fv2,ier)
!     
!     normalize the eigenvectors
!     
      dd=dsqrt(alz(1,3)**2+alz(2,3)**2+alz(3,3)**2)
      p3(1)=alz(1,3)/dd
      p3(2)=alz(2,3)/dd
      p3(3)=alz(3,3)/dd
!     
      do l=1,3
         do k=1,3
            cmm(k,l)=0.d0
            do j=1,3
               cm(j,l,k)=0.d0
            enddo
         enddo
      enddo
!     
      do l=1,3
         do j=1,3
            do i=1,3
               do k=1,3
                  cm(l,j,i)=cm(l,j,i)+c(l,k,j,i)*xn(k);
               enddo
            enddo        
         enddo
      enddo
!     
      do  l=1,3
         do  j=1,3        
            do  i=1,3
               cmm(l,j)=cmm(l,j)+cm(l,j,i)*p3(i);        
            enddo
         enddo        
      enddo
!     
      do j=1,3
         do l=1,3
            v(j)=v(j)+cmm(l,j)*p3(l)
         enddo
      enddo
!     
      dd=dsqrt(v(1)**2+v(2)**2+v(3)**2) 
      speed=dd/dsqrt(rho*al(3))
!     
!     inverse wave speed
!     
      a=1.d0/speed
!
c      if(speed.gt.dsqrt(500000.d0/7.8d-9)) then
c      write(*,*) 'inversewavspd'
c      write(*,'(6(1x,e11.4))') xn(1),xn(2),xn(3),
c     &     (dsqrt(al(3)/rho)-speed)/speed,xn
c     &     (1)*p3(1)+xn(2)*p3(2)+xn(3)*p3(3)
c      endif
!     
      return
      end
      
