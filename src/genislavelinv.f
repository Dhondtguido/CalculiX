!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
!   function calculation islavelinv
!
!  [out] islavelinv       (i)==0 if there is no slave node in the element, >0 otherwise
!  [in] jqtloc	        pointer into field irowtloc
!
      subroutine genislavelinv(islavelinv,jqtloc,
     &     lakon,ipkon,kon,ne,nasym)
!     
!     Author: Saskia Sitzmann
!     
      implicit none
!     
      character*8 lakon(*)
!     
      integer kon(*),islavelinv(*),i,j,ne,jqtloc(*),
     &     ipkon(*),konl(26),nope,node,indexe,nasym
!     
      do i=1,ne
c     if(islavelinv(i).lt.1)cycle
!     
         indexe=ipkon(i)
         if(lakon(i)(1:5).eq.'C3D8I') then
            nope=11
         elseif(lakon(i)(4:5).eq.'20') then
c     Bernhardi end
            nope=20
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:5).eq.'15') then
            nope=15
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         elseif((lakon(i)(1:2).eq.'ES').and.(lakon(i)(7:7).ne.'F')) then
!     
!     spring and contact spring elements (NO dashpot elements
!     = ED... elements)
!     
            read(lakon(i)(8:8),'(i1)') nope
            nope=nope+1
!     
!     local contact spring number
!     if friction is involved, the contact spring element
!     matrices are determined in mafillsmas.f
!     
            if(lakon(i)(7:7).eq.'C') then
               if(nasym.eq.1) cycle
               konl(nope+1)=kon(indexe+nope+1)
            endif
         else
            cycle
         endif
         do j=1,nope
            konl(j)=kon(indexe+j) 
         enddo
!     
         do j=1,nope
            node=konl(j)
cccccc error?????? (Guido 22 Nov 2019)
c            if(jqtloc(node+1)-jqtloc(node).gt.1)then
            if(jqtloc(node+1)-jqtloc(node).gt.0)then
               islavelinv(i)=1
            endif
            
         enddo
      enddo
!     
      return
      end
