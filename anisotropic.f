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
      subroutine anisotropic(anisol,anisox)
!
!     expands the 21 anisotropic elastic constants into a
!     3x3x3x3 matrix
!
      implicit none
!
      real*8 anisol(21),anisox(3,3,3,3)
!
      anisox(1,1,1,1)=anisol(1)
      anisox(1,1,1,2)=anisol(7)
      anisox(1,1,1,3)=anisol(11)
      anisox(1,1,2,1)=anisol(7)
      anisox(1,1,2,2)=anisol(2)
      anisox(1,1,2,3)=anisol(16)
      anisox(1,1,3,1)=anisol(11)
      anisox(1,1,3,2)=anisol(16)
      anisox(1,1,3,3)=anisol(4)
      anisox(1,2,1,1)=anisol(7)
      anisox(1,2,1,2)=anisol(10)
      anisox(1,2,1,3)=anisol(14)
      anisox(1,2,2,1)=anisol(10)
      anisox(1,2,2,2)=anisol(8)
      anisox(1,2,2,3)=anisol(19)
      anisox(1,2,3,1)=anisol(14)
      anisox(1,2,3,2)=anisol(19)
      anisox(1,2,3,3)=anisol(9)
      anisox(1,3,1,1)=anisol(11)
      anisox(1,3,1,2)=anisol(14)
      anisox(1,3,1,3)=anisol(15)
      anisox(1,3,2,1)=anisol(14)
      anisox(1,3,2,2)=anisol(12)
      anisox(1,3,2,3)=anisol(20)
      anisox(1,3,3,1)=anisol(15)
      anisox(1,3,3,2)=anisol(20)
      anisox(1,3,3,3)=anisol(13)
      anisox(2,1,1,1)=anisol(7)
      anisox(2,1,1,2)=anisol(10)
      anisox(2,1,1,3)=anisol(14)
      anisox(2,1,2,1)=anisol(10)
      anisox(2,1,2,2)=anisol(8)
      anisox(2,1,2,3)=anisol(19)
      anisox(2,1,3,1)=anisol(14)
      anisox(2,1,3,2)=anisol(19)
      anisox(2,1,3,3)=anisol(9)
      anisox(2,2,1,1)=anisol(2)
      anisox(2,2,1,2)=anisol(8)
      anisox(2,2,1,3)=anisol(12)
      anisox(2,2,2,1)=anisol(8)
      anisox(2,2,2,2)=anisol(3)
      anisox(2,2,2,3)=anisol(17)
      anisox(2,2,3,1)=anisol(12)
      anisox(2,2,3,2)=anisol(17)
      anisox(2,2,3,3)=anisol(5)
      anisox(2,3,1,1)=anisol(16)
      anisox(2,3,1,2)=anisol(19)
      anisox(2,3,1,3)=anisol(20)
      anisox(2,3,2,1)=anisol(19)
      anisox(2,3,2,2)=anisol(17)
      anisox(2,3,2,3)=anisol(21)
      anisox(2,3,3,1)=anisol(20)
      anisox(2,3,3,2)=anisol(21)
      anisox(2,3,3,3)=anisol(18)
      anisox(3,1,1,1)=anisol(11)
      anisox(3,1,1,2)=anisol(14)
      anisox(3,1,1,3)=anisol(15)
      anisox(3,1,2,1)=anisol(14)
      anisox(3,1,2,2)=anisol(12)
      anisox(3,1,2,3)=anisol(20)
      anisox(3,1,3,1)=anisol(15)
      anisox(3,1,3,2)=anisol(20)
      anisox(3,1,3,3)=anisol(13)
      anisox(3,2,1,1)=anisol(16)
      anisox(3,2,1,2)=anisol(19)
      anisox(3,2,1,3)=anisol(20)
      anisox(3,2,2,1)=anisol(19)
      anisox(3,2,2,2)=anisol(17)
      anisox(3,2,2,3)=anisol(21)
      anisox(3,2,3,1)=anisol(20)
      anisox(3,2,3,2)=anisol(21)
      anisox(3,2,3,3)=anisol(18)
      anisox(3,3,1,1)=anisol(4)
      anisox(3,3,1,2)=anisol(9)
      anisox(3,3,1,3)=anisol(13)
      anisox(3,3,2,1)=anisol(9)
      anisox(3,3,2,2)=anisol(5)
      anisox(3,3,2,3)=anisol(18)
      anisox(3,3,3,1)=anisol(13)
      anisox(3,3,3,2)=anisol(18)
      anisox(3,3,3,3)=anisol(6)
!
      return
      end

