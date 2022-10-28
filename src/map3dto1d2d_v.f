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
      subroutine map3dto1d2d_v(yn,ipkon,inum,kon,lakon,nfield,nk,
     &  ne,nactdof)
!
!     interpolates basic degree of freedom nodal values 
!     (displacements, temperatures) to 1d/2d nodal locations
!
      implicit none
!
      logical quadratic
!
      character*8 lakon(*),lakonl
!
      integer ipkon(*),inum(*),kon(*),ne,indexe,nfield,nk,i,j,k,l,
     &  node3(8,3),node6(3,6),node8(3,8),node2d,node3d,indexe2d,ne1d2d,
     &  node3m(8,3),iflag,nactdof(nfield,*),jmax
!
      real*8 yn(nfield,*),ratioe(3)
!
!
!
      include "gauss.f"
!
      data node3 /1,4,8,5,12,20,16,17,9,11,15,13,
     &            0,0,0,0,2,3,7,6,10,19,14,18/
      data node3m /1,5,8,4,17,16,20,12,
     &             0,0,0,0,0,0,0,0,
     &             3,7,6,2,19,14,18,10/
      data node6 /1,13,4,2,14,5,3,15,6,7,0,10,8,0,11,9,0,12/
      data node8 /1,17,5,2,18,6,3,19,7,4,20,8,9,0,13,10,0,14,
     &      11,0,15,12,0,16/
      data ratioe /0.16666666666667d0,0.66666666666666d0,
     &     0.16666666666667d0/
      data iflag /2/
!
!     removing any results in 1d/2d nodes
!
      ne1d2d=0
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         lakonl=lakon(i)
         if((lakonl(7:7).eq.' ').or.(lakonl(7:7).eq.'I').or.
     &      (lakonl(1:1).ne.'C')) cycle
         ne1d2d=1
         indexe=ipkon(i)
c!
c!        inactivating the 3d expansion nodes of 1d/2d elements
c!
c         do j=1,20
c            inum(kon(indexe+j))=0
c         enddo
!
         if((lakonl(4:5).eq.'15').or.(lakonl(4:4).eq.'6')) then
            if(lakonl(4:5).eq.'15') then
               indexe2d=indexe+15
               jmax=6
            else
               indexe2d=indexe+6
               jmax=3
            endif
            do j=1,jmax
               node2d=kon(indexe2d+j)
               inum(node2d)=0
               do k=1,nfield
                  if(nactdof(k,node2d).le.0) yn(k,node2d)=0.d0
               enddo
            enddo
         elseif(lakonl(7:7).eq.'B') then
            if(lakonl(4:5).eq.'8I') then
               indexe2d=indexe+11
               jmax=2
            elseif(lakonl(4:5).eq.'8R') then
               indexe2d=indexe+8
               jmax=2
            elseif(lakonl(4:5).eq.'20') then
               indexe2d=indexe+20
               jmax=3
            endif
            do j=1,jmax
               node2d=kon(indexe2d+j)
               inum(node2d)=0
               do k=1,nfield
                  if(nactdof(k,node2d).le.0) yn(k,node2d)=0.d0
               enddo
            enddo
         else
            if(lakonl(4:5).eq.'8I') then
               indexe2d=indexe+11
               jmax=4
            elseif((lakonl(4:5).eq.'8R').or.(lakonl(4:5).eq.'8 ')) then
               indexe2d=indexe+8
               jmax=4
            elseif(lakonl(4:5).eq.'20') then
               indexe2d=indexe+20
               jmax=8
            endif
            do j=1,jmax
               node2d=kon(indexe2d+j)
               inum(node2d)=0
               do k=1,nfield
                  if(nactdof(k,node2d).le.0) yn(k,node2d)=0.d0
               enddo
            enddo
         endif
!
!        inactivating the 3d expansion nodes of 1d/2d elements
!
         do j=1,indexe2d-indexe
            inum(kon(indexe+j))=0
         enddo
!
      enddo
!
!     if no 1d/2d elements return
!
      if(ne1d2d.eq.0) return
!
!     interpolation of 3d results on 1d/2d nodes
!
      do i=1,ne
!     
         if(ipkon(i).lt.0) cycle
         lakonl=lakon(i)
         if((lakonl(7:7).eq.' ').or.(lakonl(7:7).eq.'I').or.
     &        (lakonl(1:1).ne.'C')) cycle
         indexe=ipkon(i)
!
!        check whether linear or quadratic element
!
         if((lakonl(4:4).eq.'6').or.(lakonl(4:4).eq.'8')) then
            quadratic=.false.
         else
            quadratic=.true.
         endif
!     
         if((lakonl(4:5).eq.'15').or.(lakonl(4:4).eq.'6')) then
            if(lakonl(4:5).eq.'15') then
               indexe2d=indexe+15
               jmax=6
            else
               indexe2d=indexe+6
               jmax=3
            endif
            do j=1,jmax
               node2d=kon(indexe2d+j)
               inum(node2d)=inum(node2d)-1
!     
!     taking the mean across the thickness
!     
               if((j.le.3).and.(quadratic)) then
!     
!     end nodes: weights 1/6,2/3 and 1/6
!     
                  do l=1,3
                     node3d=kon(indexe+node6(l,j))
                     do k=1,nfield
                        if(nactdof(k,node2d).le.0) yn(k,node2d)=
     &                       yn(k,node2d)+yn(k,node3d)*ratioe(l)
                     enddo
                  enddo
               else
!     
!     middle nodes: weights 1/2,1/2
!     
                  do l=1,3,2
                     node3d=kon(indexe+node6(l,j))
                     do k=1,nfield
                        if(nactdof(k,node2d).le.0) yn(k,node2d)=
     &                       yn(k,node2d)+yn(k,node3d)/2.d0
                     enddo
                  enddo
               endif
            enddo
         elseif(lakonl(7:7).eq.'B') then
            if(lakonl(4:5).eq.'8I') then
               indexe2d=indexe+11
               jmax=2
            elseif(lakonl(4:5).eq.'8R') then
               indexe2d=indexe+8
               jmax=2
            elseif(lakonl(4:5).eq.'20') then
               indexe2d=indexe+20
               jmax=3
            endif
!     
!     mean values for beam elements
!     
            do j=1,jmax
               node2d=kon(indexe2d+j)
!     
!     mean value of vertex values
!     
               do l=1,4
                  inum(node2d)=inum(node2d)-1
                  if(quadratic) then
                     node3d=kon(indexe+node3(l,j))
                  else
                     node3d=kon(indexe+node3(l,2*j-1))
                  endif
                  do k=1,nfield
                     if(nactdof(k,node2d).le.0) yn(k,node2d)=
     &                    yn(k,node2d)+yn(k,node3d)
                  enddo
               enddo
            enddo
         else
            if(lakonl(4:5).eq.'8I') then
               indexe2d=indexe+11
               jmax=4
            elseif((lakonl(4:5).eq.'8R').or.(lakonl(4:5).eq.'8 ')) then
               indexe2d=indexe+8
               jmax=4
            elseif(lakonl(4:5).eq.'20') then
               indexe2d=indexe+20
               jmax=8
            endif
            do j=1,jmax
               node2d=kon(indexe2d+j)
               inum(node2d)=inum(node2d)-1
!     
!     taking the mean across the thickness
!     
               if((j.le.4).and.(quadratic)) then
!     
!     end nodes: weights 1/6,2/3 and 1/6
!     
                  do l=1,3
                     node3d=kon(indexe+node8(l,j))
                     do k=1,nfield
                        if(nactdof(k,node2d).le.0) yn(k,node2d)=
     &                       yn(k,node2d)+yn(k,node3d)*ratioe(l)
                     enddo
                  enddo
               else
!     
!     middle nodes: weights 1/2,1/2
!     
                  do l=1,3,2
                     node3d=kon(indexe+node8(l,j))
                     do k=1,nfield
                        if(nactdof(k,node2d).le.0) yn(k,node2d)=
     &                       yn(k,node2d)+yn(k,node3d)/2.d0
                     enddo
                  enddo
               endif
            enddo
         endif
!     
      enddo
!     
!     taking the mean of nodal contributions coming from different
!     elements having the node in common
!     
      do i=1,nk
         if(inum(i).lt.0) then
            inum(i)=-inum(i)
            do j=1,nfield
               if(nactdof(j,i).le.0) yn(j,i)=yn(j,i)/inum(i)
            enddo
         endif
      enddo
!     
      return
      end
      
