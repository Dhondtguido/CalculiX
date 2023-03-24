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
      subroutine findextsurface(nodface,ipoface,ne,ipkon,lakon,
     &  kon,konfa,ipkonfa,nk,lakonfa,nsurfs,ifreemax)
!
      implicit none
!
!     1) catalogues the external surface of the structure: see
!     comments further down
!     2) creates topology fields for the external faces:
!        ipkonfa(1...nsurfs): pointer into konfa
!        konfa(..): konfa(ipkonfa(i)+1,....ipkonfa(i+1)) contains
!                   the nodes belonging to face i
!        lakonfa(1...nsurfs): label (S3, S4, S6 or S8)
!
      character*8 lakon(*),lakonfa(*)
!
      integer ipkon(*),kon(*),ne,nodface(5,*),ipoface(*),nk,nopem,m,
     &  ithree,ifour,ifaceq(8,6),konfa(*),ipkonfa(*),nelemm,jfacem,
     &  ifacet(6,4),ifacew2(8,5),ifree,ifreenew,index,indexold,
     &  i,j,k,nodes(4),indexe,konl(26),nope,nsurfs,
     &  ifacew1(4,5),ifreemax
!
!
!
!     nodes belonging to the element faces
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
      data ifacew1 /1,3,2,0,
     &             4,5,6,0,
     &             1,2,5,4,
     &             2,3,6,5,
     &             3,1,4,6/
      data ifacew2 /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             3,1,4,6,9,13,12,15/
!
      ithree=3
      ifour=4
!
!     determining the external element faces of the mesh; 
!     the faces are catalogued by the three lowest nodes numbers
!     in ascending order. ipoface(i) points to a face for which
!     node i is the lowest node and nodface(1,ipoface(i)) and
!     nodface(2,ipoface(i)) are the next lower ones. 
!     nodface(3,ipoface(i)) contains the element number,
!     nodface(4,ipoface(i)) the face number and nodface(5,ipoface(i))
!     is a pointer to the next surface for which node i is the
!     lowest node; if there are no more such surfaces the pointer
!     has the value zero
!     An external element face is one which belongs to one element
!     only
!
      ifree=1
      ifreemax=1
      do i=1,6*ne-1
         nodface(5,i)=i+1
      enddo
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if(lakon(i)(1:3).ne.'C3D') cycle
         indexe=ipkon(i)
         if((lakon(i)(4:4).eq.'2').or.(lakon(i)(4:4).eq.'8')) then
            do j=1,6
               do k=1,4
                  nodes(k)=kon(indexe+ifaceq(k,j))
               enddo
               call insertsorti(nodes,ifour)
c               call isortii(nodes,iaux,ifour,kflag)
               indexold=0
               index=ipoface(nodes(1))
               do
!     
!     adding a surface which has not been 
!     catalogued so far
!     
                  if(index.eq.0) then
                     ifreemax=max(ifreemax,ifree)
                     ifreenew=nodface(5,ifree)
                     nodface(1,ifree)=nodes(2)
                     nodface(2,ifree)=nodes(3)
                     nodface(3,ifree)=i
                     nodface(4,ifree)=j
                     nodface(5,ifree)=ipoface(nodes(1))
                     ipoface(nodes(1))=ifree
                     ifree=ifreenew
                     exit
                  endif
!     
!     removing a surface which has already
!     been catalogued
!     
                  if((nodface(1,index).eq.nodes(2)).and.
     &                 (nodface(2,index).eq.nodes(3))) then
                     if(indexold.eq.0) then
                        ipoface(nodes(1))=nodface(5,index)
                     else
                        nodface(5,indexold)=nodface(5,index)
                     endif
                     nodface(5,index)=ifree
                     ifree=index
                     exit
                  endif
                  indexold=index
                  index=nodface(5,index)
               enddo
            enddo
         elseif((lakon(i)(4:4).eq.'4').or.(lakon(i)(4:5).eq.'10')) then
            do j=1,4
               do k=1,3
                  nodes(k)=kon(indexe+ifacet(k,j))
               enddo
               call insertsorti(nodes,ithree)
c               call isortii(nodes,iaux,ithree,kflag)
               indexold=0
               index=ipoface(nodes(1))
               do
!     
!     adding a surface which has not been 
!     catalogued so far
!     
                  if(index.eq.0) then
                     ifreemax=max(ifreemax,ifree)
                     ifreenew=nodface(5,ifree)
                     nodface(1,ifree)=nodes(2)
                     nodface(2,ifree)=nodes(3)
                     nodface(3,ifree)=i
                     nodface(4,ifree)=j
                     nodface(5,ifree)=ipoface(nodes(1))
                     ipoface(nodes(1))=ifree
                     ifree=ifreenew
                     exit
                  endif
!     
!     removing a surface which has already
!     been catalogued
!     
                  if((nodface(1,index).eq.nodes(2)).and.
     &                 (nodface(2,index).eq.nodes(3))) then
                     if(indexold.eq.0) then
                        ipoface(nodes(1))=nodface(5,index)
                     else
                        nodface(5,indexold)=nodface(5,index)
                     endif
                     nodface(5,index)=ifree
                     ifree=index
                     exit
                  endif
                  indexold=index
                  index=nodface(5,index)
               enddo
            enddo
         else
            do j=1,5
               if(j.le.2) then
                  do k=1,3
                     nodes(k)=kon(indexe+ifacew2(k,j))
                  enddo
                  call insertsorti(nodes,ithree)
c                  call isortii(nodes,iaux,ithree,kflag)
               else
                  do k=1,4
                     nodes(k)=kon(indexe+ifacew2(k,j))
                  enddo
                  call insertsorti(nodes,ifour)
c                  call isortii(nodes,iaux,ifour,kflag)
               endif
               indexold=0
               index=ipoface(nodes(1))
               do
!     
!     adding a surface which has not been 
!     catalogued so far
!     
                  if(index.eq.0) then
                     ifreemax=max(ifreemax,ifree)
                     ifreenew=nodface(5,ifree)
                     nodface(1,ifree)=nodes(2)
                     nodface(2,ifree)=nodes(3)
                     nodface(3,ifree)=i
                     nodface(4,ifree)=j
                     nodface(5,ifree)=ipoface(nodes(1))
                     ipoface(nodes(1))=ifree
                     ifree=ifreenew
                     exit
                  endif
!     
!     removing a surface which has already
!     been catalogued
!     
                  if((nodface(1,index).eq.nodes(2)).and.
     &                 (nodface(2,index).eq.nodes(3))) then
                     if(indexold.eq.0) then
                        ipoface(nodes(1))=nodface(5,index)
                     else
                        nodface(5,indexold)=nodface(5,index)
                     endif
                     nodface(5,index)=ifree
                     ifree=index
                     exit
                  endif
                  indexold=index
                  index=nodface(5,index)
               enddo
            enddo
         endif
      enddo
!     
!     Create fields konfa and ipkonfa for calculation of normals of shell
!     elements
!     
      nsurfs=0
      ifree=0
      do j=1,nk
!     
         if(ipoface(j).eq.0) cycle
         indexe=ipoface(j)
!     
         do
            nsurfs=nsurfs+1
            nelemm=nodface(3,indexe)
            jfacem=nodface(4,indexe)
!     
!     treatment of hexahedral elements
!     
            if(lakon(nelemm)(4:4).eq.'8') then
               nopem=4
               nope=8
            elseif(lakon(nelemm)(4:5).eq.'20') then
               nopem=8
               nope=20
!     
!     treatment of tethrahedral elements
!     
            elseif(lakon(nelemm)(4:5).eq.'10') then
               nopem=6
               nope=10
            elseif(lakon(nelemm)(4:4).eq.'4') then
               nopem=3
               nope=4
!     
!     treatment of wedge faces
!     
            elseif(lakon(nelemm)(4:4).eq.'6') then
               nope=6
               if(jfacem.le.2) then
                  nopem=3
               else
                  nopem=4
               endif
            elseif(lakon(nelemm)(4:5).eq.'15') then
               nope=15
               if(jfacem.le.2) then
                  nopem=6
               else
                  nopem=8
               endif
            endif 
!     
!     actual position of the nodes 
!     
            do k=1,nope
               konl(k)=kon(ipkon(nelemm)+k)
            enddo
!     
!     quadratic quad shell element
!     
            if((nope.eq.20).or.
     &           ((nope.eq.15).and.(jfacem.gt.2))) then    
               lakonfa(nsurfs)='S8'
               ipkonfa(nsurfs)=ifree
               do m=1,nopem
                  if(nope.eq.20) then
                     konfa(ifree+m)=konl(ifaceq(m,jfacem))
                  else
                     konfa(ifree+m)=konl(ifacew2(m,jfacem))
                  endif
               enddo
               ifree=ifree+nopem
!     
!     linear quad shell element
!     
            elseif((nope.eq.8).or.
     &              ((nope.eq.6).and.(jfacem.gt.2))) then
               lakonfa(nsurfs)='S4'
               ipkonfa(nsurfs)=ifree
               do m=1,nopem
                  if(nope.eq.8) then
                     konfa(ifree+m)=konl(ifaceq(m,jfacem))
                  else
                     konfa(ifree+m)=konl(ifacew1(m,jfacem))
                  endif
               enddo
               ifree=ifree+nopem
!     
!     quadratic tri shell element
!     
            elseif((nope.eq.10).or.
     &              ((nope.eq.15).and.(jfacem.le.2))) then    
               lakonfa(nsurfs)='S6'
               ipkonfa(nsurfs)=ifree
               do m=1,nopem
                  if(nope.eq.10) then
                     konfa(ifree+m)=konl(ifacet(m,jfacem))
                  else
                     konfa(ifree+m)=konl(ifacew2(m,jfacem))
                  endif
               enddo
               ifree=ifree+nopem
!     
!     linear tri shell element
!     
            elseif((nope.eq.4).or.
     &              ((nope.eq.6).and.(jfacem.le.2))) then
               lakonfa(nsurfs)='S3'
               ipkonfa(nsurfs)=ifree
               do m=1,nopem
                  if(nope.eq.4) then
                     konfa(ifree+m)=konl(ifacet(m,jfacem))
                  else
                     konfa(ifree+m)=konl(ifacew1(m,jfacem))
                  endif
               enddo
               ifree=ifree+nopem    
            endif
            indexe=nodface(5,indexe)
            if(indexe.eq.0) exit
         enddo      
      enddo
!     
      return
      end
