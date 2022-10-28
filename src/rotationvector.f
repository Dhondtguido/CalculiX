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
      subroutine rotationvector(a,v)
!
!     calculates rotation vector v from rotation matrix a
!
      implicit none
!
      real*8 a(3,3),q(4),v(3),theta,length,pi
!
!
!
!     based on: J.M.P. van Waveren, From Quaternion to Matrix and Back
!     February 27th 2005, Id Software, Inc.
!     and
!     Jay A. Farrel, Computation of the Quaternion from a Rotation 
!     Matrix, November 30, 2015, University of California, Riverside
!
      pi=4.d0*datan(1.d0)
!
      if ((a(1,1)+a(2,2)+a(3,3)).gt.0.d0) then
            q(1)=dsqrt(a(1,1)+a(2,2)+a(3,3)+1.d0)/2.d0
            q(2)=(a(3,2)-a(2,3))/(4.d0*q(1))
            q(3)=(a(1,3)-a(3,1))/(4.d0*q(1))
            q(4)=(a(2,1)-a(1,2))/(4.d0*q(1))
      else if (((a(1,1).gt.a(2,2)).and.(a(1,1).gt.a(3,3)))) then
            q(2)=dsqrt(a(1,1)-a(2,2)-a(3,3)+1.d0)/2.d0
            q(1)=(a(3,2)-a(2,3))/(4.d0*q(2))
            q(3)=(a(1,2)+a(2,1))/(4.d0*q(2))
            q(4)=(a(1,3)+a(3,1))/(4.d0*q(2))
      else if (a(2,2).gt.a(3,3)) then
            q(3)=dsqrt(-a(1,1)+a(2,2)-a(3,3)+1.d0)/2.d0
            q(1)=(a(1,3)-a(3,1))/(4.d0*q(3))
            q(2)=(a(1,2)+a(2,1))/(4.d0*q(3))
            q(4)=(a(2,3)+a(3,2))/(4.d0*q(3))
      else
            q(4)=dsqrt(-a(1,1)-a(2,2)+a(3,3)+1.d0)/2.d0
            q(1)=(a(2,1)-a(1,2))/(4.d0*q(4))
            q(2)=(a(1,3)+a(3,1))/(4.d0*q(4))
            q(3)=(a(2,3)+a(3,2))/(4.d0*q(4))            
      endif
!
      length=dsqrt(q(2)*q(2)+q(3)*q(3)+q(4)*q(4))
      theta=2.d0*acos(q(1))
!
!     if pi<theta<2*pi: reverse direction of rotation vector
!     and map angle in the range (0,pi)
!
      if (theta.gt.pi) then
            theta=2*pi-theta
            q(2)=-q(2)
            q(3)=-q(3)
            q(4)=-q(4)
      endif
!
      if (length.ne.0) then
            v(1)=theta*q(2)/length
            v(2)=theta*q(3)/length
            v(3)=theta*q(4)/length
      else
            v(1)=0.d0
            v(2)=0.d0
            v(3)=0.d0
      endif
!
!     if theta=pi:
!     if x<0 -> change sign of rotation vector
!     elseif x=0 and y<0 -> change sign of rotation vector
!     elseif x=0 and y=0 and z<0 -> change sign of rotation vector
!
      if (theta.eq.pi) then
            if(v(1).lt.0.d0) then
!                 +++ vs ---
!                 ++- vs --+
!                 +-+ vs -+-
!                 +-- vs -++
!                 +0+ vs -0-
!                 ++0 vs --0
!                 +0- vs -0+
!                 +-0 vs -+0
!                 +00 vs -00
                  v(1)=-v(1)
                  v(2)=-v(2)
                  v(3)=-v(3)
            elseif(v(1).eq.0.d0) then
                  if(v(2).lt.0.d0) then
!                       0+- vs 0-+
!                       0+0 vs 0-0
!                       0++ vs 0--
                        v(2)=-v(2)
                        v(3)=-v(3)
                  elseif(v(2).eq.0.d0) then
                        if(v(3).lt.0.d0) then
!                             00+ vs 00-
                              v(3)=-v(3)
                        endif
                  endif
            endif
      endif
!
c      write(*,*)'ROTATION VECTOR'
c      write(*,*)v(1),v(2),v(3)
!      
      return
!
      end
