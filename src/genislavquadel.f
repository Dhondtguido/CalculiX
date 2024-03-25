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
!     calculates islavquadel(i) for all elements i
!
!     islavquadel(i)=1 if i is a quadratic element containing
!                         at least one slave node
!                     else the value is 0.
!
!     the regular shape functions for elements i for which islavquadel(i)=1
!     have to be replaced by the tilde shape functions.
!
      subroutine genislavquadel(islavquadel,jqt,
     &     lakon,ipkon,kon,ne,nasym)
!     
!     Author: Saskia Sitzmann
!     
      implicit none
!     
      character*8 lakon(*)
!     
      integer kon(*),islavquadel(*),i,j,ne,jqt(*),nslavquadel,
     &     ipkon(*),konl(26),nope,node,indexe,nasym
!
      nslavquadel=0
      loop: do i=1,ne
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
            if(jqt(node+1)-jqt(node).gt.1)then
              islavquadel(i)=1
              nslavquadel=nslavquadel+1
              cycle loop
            endif
         enddo
      enddo loop
!     
      return
      end
