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
      subroutine assigndomtonodes(ne,lakon,ipkon,kon,ielmat,inomat,
     &  elcon,ncmat_,ntmat_,mi,ne2)
!
!     assigns the domain a node belongs to, to this node
!     (for electromagnetic calculations, only for nodes not
!      belonging to shells)
!
      implicit none
!
      character*8 lakon(*)
!
      integer i,j,nope,ne,imat,mi(*),ielmat(mi(3),*),ipkon(*),inomat(*),
     &  ncmat_,ntmat_,node,kon(*),indexe,ne2
!
      real*8 elcon(0:ncmat_,ntmat_,*)
!
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if(lakon(i)(7:7).eq.'L') cycle
         indexe=ipkon(i)
!
         if(lakon(i)(4:5).eq.'20') then
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
         else
            cycle
         endif
!
         imat=ielmat(1,i)
!
!        ne2 is the number of elements in domain 2 = A,V-domain
!
         if(int(elcon(2,1,imat)).eq.2) ne2=ne2+1
!
         do j=1,nope
            node=kon(indexe+j)
            if(inomat(node).ne.0) then
               if(inomat(node).ne.int(elcon(2,1,imat))) then
                  write(*,*) '*ERROR in assigndomtonodes: a node'
                  write(*,*) '       cannot belong to more than'
                  write(*,*) '       one domain'
                  call exit(201)
               else
                  cycle
               endif
            endif
            inomat(node)=int(elcon(2,1,imat))
         enddo
      enddo
!
      return
      end
