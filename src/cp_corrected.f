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
!     calculation of corrected cp      
!
!     author: Yannick Muller
!     
      subroutine cp_corrected(cp,Tg1,Tg2,cp_cor)
!
      implicit none
!
      real*8 cp,Tg1,Tg2,cp_cor
!
      Tg1=Tg1
      Tg2=Tg2
      cp_cor=cp
!
      return
      end
