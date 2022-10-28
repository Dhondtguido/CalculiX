!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
!     subroutine to find the local node number "getlocno" within the element
!     based on the local face number "jface" within the element
!     and the local node number "ii" within the face
!
      integer function getlocno(ii,jface,nope)
!      
!     author: Saskia Sitzmann
!      
      implicit none
!     
      integer ii,jface,nope,
     &     ifaceq(9,6),ifacet(7,4),ifacew(8,5)
!
!      
      include "gauss.f"
!     
      ifaceq=reshape((/4,3,2,1,11,10,9,12,21,
     &            5,6,7,8,13,14,15,16,22,
     &            1,2,6,5,9,18,13,17,23,
     &            2,3,7,6,10,19,14,18,24,
     &            3,4,8,7,11,20,15,19,25,
     &            4,1,5,8,12,17,16,20,26/),(/9,6/))
      ifacet=reshape((/1,3,2,7,6,5,11,
     &             1,2,4,5,9,8,12,
     &             2,3,4,6,10,9,13,
     &             1,4,3,8,10,7,14/),(/7,4/))
      ifacew=reshape((/1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/),(/8,5/))
!     
      getlocno=0
!     
      if((nope.eq.20).or.(nope.eq.8).or.
     &      (nope.eq.26).or.(nope.eq.11)) then
         getlocno=ifaceq(ii,jface)
      elseif((nope.eq.10).or.(nope.eq.4)) then
         getlocno=ifacet(ii,jface)
      else
         getlocno=ifacew(ii,jface)
      endif
!     
      end
!     
      
