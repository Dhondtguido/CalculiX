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
      subroutine rotationvectorinv(c,v)
!
!     calculates rotation matrix from rotation vector
!
      implicit none
!
      real*8 c(3,3),v(3),theta,ds,dc
!
!
!
      theta=dsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
!
      if (theta.eq.0.d0) then
            c(1,1)=1.d0
            c(1,2)=0.d0
            c(1,3)=0.d0
            c(2,1)=0.d0
            c(2,2)=1.d0
            c(2,3)=0.d0
            c(3,1)=0.d0
            c(3,2)=0.d0
            c(3,3)=1.d0
      else      
!
            dc=dcos(theta)
            ds=dsin(theta)
!
!     C-matrix from Guido Dhondt, The Finite Element
!     Method for Three-Dimensional Thermomechanical
!     Applications p 158
!     
            c(1,1)=dc+(1.d0-dc)*v(1)*v(1)/(theta*theta)
            c(1,2)=   (1.d0-dc)*v(1)*v(2)/(theta*theta)-ds*v(3)/theta
            c(1,3)=   (1.d0-dc)*v(1)*v(3)/(theta*theta)+ds*v(2)/theta
            c(2,1)=   (1.d0-dc)*v(2)*v(1)/(theta*theta)+ds*v(3)/theta
            c(2,2)=dc+(1.d0-dc)*v(2)*v(2)/(theta*theta)
            c(2,3)=   (1.d0-dc)*v(2)*v(3)/(theta*theta)-ds*v(1)/theta
            c(3,1)=   (1.d0-dc)*v(3)*v(1)/(theta*theta)-ds*v(2)/theta
            c(3,2)=   (1.d0-dc)*v(3)*v(2)/(theta*theta)+ds*v(1)/theta
            c(3,3)=dc+(1.d0-dc)*v(3)*v(3)/(theta*theta)
      endif
!
      return
!
      end
