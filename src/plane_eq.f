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
!     Subroutine plane_eq.f
!
!     Creates the plane equation from three known points. 
!     Gives the z-coordinate of the fourth point as an output.
!
!     x1,y1,z1: The coordinates of the first point
!     x2,y2,z2: The coordinates of the second point
!     x3,y3,z3: The coordinates of the third point
!     x0,y0: The x and y-coordinates for the fourth point
!     output: The z-coordinate according to the x0 and y0
!
!     by: Jaro Hokkanen
!
      subroutine plane_eq(x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,output)
!
      implicit none
!
      real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,output,
     &  a,b,c,d
!
      d=x1*y2*z3+y1*z2*x3+z1*x2*y3-x1*z2*y3-y1*x2*z3-z1*y2*x3
      if(d.ne.0.d0) then
         a=1.d0/d*(y2*z3+y1*z2+z1*y3-z2*y3-y1*z3-z1*y2)
      endif  
      if(d.ne.0.d0) then
         b=1.d0/d*(x1*z3+z2*x3+z1*x2-x1*z2-x2*z3-z1*x3)
      endif  
      if(d.ne.0.d0) then
         c=1.d0/d*(x1*y2+y1*x3+x2*y3-x1*y3-y1*x2-y2*x3)
      endif  
      if(d.ne.0.d0) then
         output=1.d0/c*(1.d0-a*x0-b*y0)
      else
         output=0.d0
      endif
      return
      end
