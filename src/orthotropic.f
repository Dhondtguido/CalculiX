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
      subroutine orthotropic(orthol,anisox)
!
!     expands the 9 orthotropic elastic constants into a
!     3x3x3x3 matrix
!
      implicit none
!
      real*8 orthol(9),anisox(3,3,3,3)
!
!
!
      anisox(1,1,1,1)=orthol(1)
      anisox(1,1,1,2)=0.d0
      anisox(1,1,1,3)=0.d0
      anisox(1,1,2,1)=0.d0
      anisox(1,1,2,2)=orthol(2)
      anisox(1,1,2,3)=0.d0
      anisox(1,1,3,1)=0.d0
      anisox(1,1,3,2)=0.d0
      anisox(1,1,3,3)=orthol(4)
      anisox(1,2,1,1)=0.d0
      anisox(1,2,1,2)=orthol(7)
      anisox(1,2,1,3)=0.d0
      anisox(1,2,2,1)=orthol(7)
      anisox(1,2,2,2)=0.d0
      anisox(1,2,2,3)=0.d0
      anisox(1,2,3,1)=0.d0
      anisox(1,2,3,2)=0.d0
      anisox(1,2,3,3)=0.d0
      anisox(1,3,1,1)=0.d0
      anisox(1,3,1,2)=0.d0
      anisox(1,3,1,3)=orthol(8)
      anisox(1,3,2,1)=0.d0
      anisox(1,3,2,2)=0.d0
      anisox(1,3,2,3)=0.d0
      anisox(1,3,3,1)=orthol(8)
      anisox(1,3,3,2)=0.d0
      anisox(1,3,3,3)=0.d0
      anisox(2,1,1,1)=0.d0
      anisox(2,1,1,2)=orthol(7)
      anisox(2,1,1,3)=0.d0
      anisox(2,1,2,1)=orthol(7)
      anisox(2,1,2,2)=0.d0
      anisox(2,1,2,3)=0.d0
      anisox(2,1,3,1)=0.d0
      anisox(2,1,3,2)=0.d0
      anisox(2,1,3,3)=0.d0
      anisox(2,2,1,1)=orthol(2)
      anisox(2,2,1,2)=0.d0
      anisox(2,2,1,3)=0.d0
      anisox(2,2,2,1)=0.d0
      anisox(2,2,2,2)=orthol(3)
      anisox(2,2,2,3)=0.d0
      anisox(2,2,3,1)=0.d0
      anisox(2,2,3,2)=0.d0
      anisox(2,2,3,3)=orthol(5)
      anisox(2,3,1,1)=0.d0
      anisox(2,3,1,2)=0.d0
      anisox(2,3,1,3)=0.d0
      anisox(2,3,2,1)=0.d0
      anisox(2,3,2,2)=0.d0
      anisox(2,3,2,3)=orthol(9)
      anisox(2,3,3,1)=0.d0
      anisox(2,3,3,2)=orthol(9)
      anisox(2,3,3,3)=0.d0
      anisox(3,1,1,1)=0.d0
      anisox(3,1,1,2)=0.d0
      anisox(3,1,1,3)=orthol(8)
      anisox(3,1,2,1)=0.d0
      anisox(3,1,2,2)=0.d0
      anisox(3,1,2,3)=0.d0
      anisox(3,1,3,1)=orthol(8)
      anisox(3,1,3,2)=0.d0
      anisox(3,1,3,3)=0.d0
      anisox(3,2,1,1)=0.d0
      anisox(3,2,1,2)=0.d0
      anisox(3,2,1,3)=0.d0
      anisox(3,2,2,1)=0.d0
      anisox(3,2,2,2)=0.d0
      anisox(3,2,2,3)=orthol(9)
      anisox(3,2,3,1)=0.d0
      anisox(3,2,3,2)=orthol(9)
      anisox(3,2,3,3)=0.d0
      anisox(3,3,1,1)=orthol(4)
      anisox(3,3,1,2)=0.d0
      anisox(3,3,1,3)=0.d0
      anisox(3,3,2,1)=0.d0
      anisox(3,3,2,2)=orthol(5)
      anisox(3,3,2,3)=0.d0
      anisox(3,3,3,1)=0.d0
      anisox(3,3,3,2)=0.d0
      anisox(3,3,3,3)=orthol(6)
!
      return
      end

