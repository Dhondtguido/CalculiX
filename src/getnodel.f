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
!
!     subroutine to find the right node for different element types 
!     based on face number jface and node number ii
!
      integer function getnodel(ii,jface,nope)
!
!     author: Saskia Sitzmann
!
      implicit none
!
      integer ii,jface,nope,
     &        ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5)
!
      include "gauss.f"
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
!
      data ifacew1 /1,3,2,0,
     &             4,5,6,0,
     &             1,2,5,4,
     &             2,3,6,5,
     &             3,1,4,6/
!
!     nodes per face for quadratic wedge elements
!
      data ifacew2 /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             3,1,4,6,9,13,12,15/
!     
      getnodel=0
!            
      if((nope.eq.20).or.(nope.eq.8)) then
               getnodel=ifaceq(ii,jface)
      elseif((nope.eq.10).or.(nope.eq.4)) then
         getnodel=ifacet(ii,jface)
      elseif(nope.eq.6) then
         getnodel=ifacew1(ii,jface)
      else
         getnodel=ifacew2(ii,jface)
      endif
!     
      end
!
