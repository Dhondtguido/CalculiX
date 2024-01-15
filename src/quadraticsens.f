
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2023 Guido Dhondt
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
      subroutine quadraticsens(ipkon,lakon,kon,nobject,dgdxglob,
     &   xinterpol,nnodes,ne,nk,nodedesiinv,objectset,nobjectstart)
!     
!     interpolation of the sensitivitites of the midnodes to the 
!     corner nodes - only valid for quadratic elements
!     
      implicit none
!     
      character*81 objectset(5,*)
      character*8 lakon(*)
!     
      integer i,ii,j,l,ielem,nodecor,nk,ne,
     &     nope,indexe,ipkon(*),konl(26),ifaceqc(2,20),
     &     ifacetc(2,10),ifacewc(2,15),kon(*),nnodes(nk),nobject,
     &     nodedesiinv(nk),istart,nobjectstart
!     
      real*8 dgdxglob(2,nk,nobject),xinterpol(nk,nobject)
!     
!     
!     
!     cornernodes next to the midnode for quadratic hex element
!     
      data ifaceqc /0,0,
     &     0,0,
     &     0,0,
     &     0,0,
     &     0,0,
     &     0,0,
     &     0,0,
     &     0,0,
     &     1,2,
     &     2,3,
     &     3,4,
     &     1,4,
     &     5,6,
     &     6,7,
     &     7,8,
     &     5,8,
     &     1,5,
     &     2,6,
     &     3,7,
     &     4,8/
!     
!     cornernodes next to the midnode for quadratic tet elements
!     
      data ifacetc /0,0,
     &     0,0,
     &     0,0,
     &     0,0,
     &     1,2,
     &     2,3,
     &     1,3,
     &     1,4,
     &     2,4,
     &     3,4/
!     
!     cornernodes next to the midnode for quadratic wedge elements
!     
      data ifacewc /0,0,
     &     0,0,
     &     0,0,
     &     0,0,
     &     0,0,
     &     0,0,
     &     1,2,
     &     2,3,
     &     1,3,
     &     4,5,
     &     5,6,
     &     4,6,
     &     1,4,
     &     2,5,
     &     3,6/
!     
!     Loop over all elements
!     
      do ielem=1,ne         
!     
        if(ipkon(ielem).lt.0) cycle
!     
!     Check if element is quadratic
!     
        if(lakon(ielem)(4:5).eq.'10') then
          nope=10
          istart=5
        elseif(lakon(ielem)(4:5).eq.'20') then
          nope=20
          istart=9
        elseif (lakon(ielem)(4:5).eq.'15') then  
          nope=15
          istart=7      
        else 
          cycle
        endif
        
        indexe=ipkon(ielem)        
!     
        do l=1,nope
          konl(l)=kon(indexe+l)
        enddo
!     
!     Loop over all midnodes of the element
!     
        do i=istart,nope
!     
          if(nodedesiinv(konl(i)).le.0) cycle
!     
!     Loop over the 2 neighbors of the midnode
!     
          do j=1,2
            if(lakon(ielem)(4:5).eq.'10') then
              nodecor=konl(ifacetc(j,i))
            elseif (lakon(ielem)(4:5).eq.'20') then
              nodecor=konl(ifaceqc(j,i))
            elseif (lakon(ielem)(4:5).eq.'15') then
              nodecor=konl(ifacewc(j,i))
            endif
            if(nodedesiinv(nodecor).eq.1) then
              nnodes(nodecor)=nnodes(nodecor)+1        
              do ii=1,nobject
                xinterpol(nodecor,ii)=xinterpol(nodecor,ii)+
     &               dgdxglob(1,konl(i),ii)
              enddo
            endif
          enddo
          nodedesiinv(konl(i))=-1
        enddo
      enddo
!     
      do i=1,nk
        if(nnodes(i).gt.0) then
          do j=1+nobjectstart,nobject
            if(objectset(1,j)(4:13).eq.'MEMBERSIZE') cycle
            if(objectset(1,j)(1:9).eq.'FIXGROWTH') cycle
            if(objectset(1,j)(1:12).eq.'FIXSHRINKAGE') cycle
            if(objectset(1,j)(1:9).eq.'PACKAGING') cycle
            dgdxglob(1,i,j)=xinterpol(i,j)/nnodes(i)
          enddo
        endif
      enddo  
!     
!     Correction of nodedesiinv
!     
      do i=1,nk
        if(nodedesiinv(i).eq.-1) then
          nodedesiinv(i)=1
        endif
      enddo
!     
      return
      end
