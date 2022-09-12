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
      subroutine ufaceload(co,ipkon,kon,lakon,nboun,nodeboun,
     &  nelemload,sideload,nload,ne,nk)
!
!
!     INPUT:
!
!     co(0..3,1..nk)     coordinates of the nodes
!     ipkon(*)           element topology pointer into field kon
!     kon(*)             topology vector of all elements
!     lakon(*)           vector with elements labels
!     nboun              number of SPC's
!     nodeboun(*)        SPC node
!     nelemload(1..2,*)  1: elements faces of which are loaded
!                        2: nodes for environmental temperatures
!     sideload(*)        load label
!     nload              number of facial distributed loads
!     ne                 highest element number
!     nk                 highest node number
!
!     user routine called at the start of each step; possible use:
!     calculation of the area of sets of elements for
!     further use to calculate film or radiation coefficients.
!     The areas can be shared using common blocks.
!
      implicit none
!
      character*8 lakon(*)
      character*20 sideload(*)
!
      integer nelemload(2,*),nload,kon(*),ipkon(*),nk,ne,nboun,
     &  nodeboun(*)
!
      real*8 co(3,*)
!
!     enter code here
!
      return
      end
      


