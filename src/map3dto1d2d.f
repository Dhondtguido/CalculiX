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
      subroutine map3dto1d2d(yn,ipkon,inum,kon,lakon,nfield,nk,
     &  ne,cflag,co,vold,force,mi,ielprop,prop)
!
!     interpolates 3d field nodal values to 1d/2d nodal locations
!
!     the number of internal state variables is limited to 999
!     (cfr. array field)
!
      implicit none
!
      logical force,quadratic
!
      character*1 cflag
      character*8 lakon(*),lakonl
!
      integer ipkon(*),inum(*),kon(*),ne,indexe,nfield,nk,i,j,k,l,m,
     &  node3(8,3),node6(3,6),node8(3,8),node2d,node3d,indexe2d,ne1d2d,
     &  node3m(8,3),node(8),m1,m2,nodea,nodeb,nodec,iflag,mi(*),jmax,
     &  jinc,nope,mint3d,null,ielprop(*)
!
      real*8 yn(nfield,*),cg(3),p(3),pcg(3),t(3),xl(3,8),shp(7,8),
     &  xsj(3),e1(3),e2(3),e3(3),s(6),dd,xi,et,ze,co(3,*),xs(3,7),
     &  vold(0:mi(2),*),ratioe(3),weight,prop(*),xil,etl
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
      null=0
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
                  yn(k,node2d)=0.d0
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
                  yn(k,node2d)=0.d0
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
                  yn(k,node2d)=0.d0
               enddo
            enddo
         endif
!
!        inactivating the 3d expansion nodes of 1d/2d elements
!        in case forces are mapped this field is used to ensure
!        that the forces in the 3d-nodes are mapped only once onto the
!        2d-nodes
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
     &      (lakonl(1:1).ne.'C')) cycle
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
               if(.not.force) then
!
!                 taking the mean across the thickness
!
                  if((j.le.3).and.(quadratic)) then
!
!                    end nodes: weights 1/6,2/3 and 1/6
!
                     do l=1,3
                        node3d=kon(indexe+node6(l,j))
                        do k=1,nfield
                           yn(k,node2d)=yn(k,node2d)+
     &                                  yn(k,node3d)*ratioe(l)
                        enddo
                     enddo
                  else
!
!                    middle nodes: weights 1/2,1/2
!
                     do l=1,3,2
                        node3d=kon(indexe+node6(l,j))
                        do k=1,nfield
                           yn(k,node2d)=yn(k,node2d)+yn(k,node3d)/2.d0
                        enddo
                     enddo
                  endif
               else
!
!                 forces must be summed
!
!                 the contribution of each 3d-node should only be taken
!                 once. inum(node3d) is used as marker.
!
                  inum(node2d)=-1
!
                  if((j.le.3).and.(quadratic)) then
!
!                    end nodes of quadratic 2d-elements
!
                     do l=1,3
                        node3d=kon(indexe+node6(l,j))
                        if(inum(node3d).ne.0) cycle
                        inum(node3d)=1
                        do k=1,nfield
                           yn(k,node2d)=yn(k,node2d)+yn(k,node3d)
                        enddo
                     enddo
                  else
!
!                    middle nodes of quadratic 2d-elements
!                    or end nodes of linear 2d-elements
!
                     do l=1,3,2
                        node3d=kon(indexe+node6(l,j))
                        if(inum(node3d).ne.0) cycle
                        inum(node3d)=1
                        do k=1,nfield
                           yn(k,node2d)=yn(k,node2d)+yn(k,node3d)
                        enddo
                     enddo
                  endif
               endif
            enddo
         elseif(lakonl(7:7).eq.'B') then
            if(lakonl(4:5).eq.'8I') then
               indexe2d=indexe+11
               jmax=2
               jinc=1
               nope=4
            elseif(lakonl(4:5).eq.'8R') then
               indexe2d=indexe+8
               jmax=2
               jinc=1
               nope=4
            elseif(lakonl(4:5).eq.'20') then
               indexe2d=indexe+20
               jmax=3
               jinc=2
               nope=8
            endif
            if(cflag.ne.'M') then
!
!              mean values for beam elements
!
               do j=1,jmax
                  node2d=kon(indexe2d+j)
                  if(.not.force) then
!
!                    mean value of vertex values
!
                     do l=1,4
                        inum(node2d)=inum(node2d)-1
                        if(quadratic) then
                           node3d=kon(indexe+node3(l,j))
                        else
                           node3d=kon(indexe+node3(l,2*j-1))
                        endif
                        do k=1,nfield
                           yn(k,node2d)=yn(k,node2d)+yn(k,node3d)
                        enddo
                     enddo
                  else
!
!                    forces must be summed across the section
!
!                    the contribution of each 3d-node should only be taken
!                    once. inum(node3d) is used as marker.
!
                     inum(node2d)=-1
!
                     if((j.ne.2).and.(quadratic)) then
!
!                       end nodes of quadratic beam elements
!
                        do l=1,8
                           node3d=kon(indexe+node3(l,j))
                           if(inum(node3d).ne.0) cycle
                           inum(node3d)=1
                           do k=1,nfield
                              yn(k,node2d)=yn(k,node2d)+yn(k,node3d)
                           enddo
                        enddo
                     else
!
!                       middle nodes of quadratic beam elements or
!                       end nodes of linear beam elements
!
                        do l=1,4
                           if(quadratic) then
                              node3d=kon(indexe+node3(l,j))
                           else
                              node3d=kon(indexe+node3(l,2*j-1))
                           endif
                           if(inum(node3d).ne.0) cycle
                           inum(node3d)=1
                           do k=1,nfield
                              yn(k,node2d)=yn(k,node2d)+yn(k,node3d)
                           enddo
                        enddo
                     endif
                  endif
               enddo
            else
!
!              section forces for beam elements
!
               do j=1,jmax,jinc
                  node2d=kon(indexe2d+j)
                  inum(node2d)=inum(node2d)-1
!
!                 coordinates of the nodes belonging to the section
!
                  do l=1,nope
                     if(quadratic) then
                        node(l)=kon(indexe+node3m(l,j))
                     else
                        node(l)=kon(indexe+node3m(l,2*j-1))
                     endif
                     do m=1,3
                        xl(m,l)=co(m,node(l))+vold(m,node(l))
                     enddo
                  enddo
!
!                 center of gravity and unit vectors 1 and 2
!
                  do m=1,3
                     cg(m)=(xl(m,1)+xl(m,2)+xl(m,3)+xl(m,4))/4.d0
                     if(j.eq.1) then
                        e1(m)=(xl(m,1)+xl(m,4)-xl(m,2)-xl(m,3))
                     else
                        e1(m)=(xl(m,2)+xl(m,3)-xl(m,1)-xl(m,4))
                     endif
                     e2(m)=(xl(m,3)+xl(m,4)-xl(m,1)-xl(m,2))
                  enddo
!
!                 normalizing e1
!
                  dd=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
                  do m=1,3
                     e1(m)=e1(m)/dd
                  enddo
!
!                 making sure that e2 is orthogonal to e1
!
                  dd=e1(1)*e2(1)+e1(2)*e2(2)+e1(3)*e2(3)
                  do m=1,3
                     e2(m)=e2(m)-dd*e1(m)
                  enddo
!
!                 normalizing e2
!                  
                  dd=dsqrt(e2(1)*e2(1)+e2(2)*e2(2)+e2(3)*e2(3))
                  do m=1,3
                     e2(m)=e2(m)/dd
                  enddo
!
!                 e3 = e1 x e2 for j=3, e3 = e2 x e1 for j=1
!
                  if(j.eq.1) then
                     e3(1)=e2(2)*e1(3)-e1(2)*e2(3)
                     e3(2)=e2(3)*e1(1)-e1(3)*e2(1)
                     e3(3)=e2(1)*e1(2)-e1(1)*e2(2)
                  else
                     e3(1)=e1(2)*e2(3)-e2(2)*e1(3)
                     e3(2)=e1(3)*e2(1)-e2(3)*e1(1)
                     e3(3)=e1(1)*e2(2)-e2(1)*e1(2)
                  endif
!
!                 loop over the integration points (2x2 for
!                 rectangular or circular cross section, else
!                 section dependent)
!                  
                  if((lakonl(8:8).eq.'R').or.
     &               (lakonl(8:8).eq.'C')) then
                     mint3d=4
                  else
                     call beamintscheme(lakonl,mint3d,ielprop(i),
     &                    prop,null,xi,et,ze,weight)
!
!                    mint3d are the 3d integration points. It is
!                    assumed that along the beam axes 2 integration
!                    points are used, i.e. the number of 2d
!                    integration points is half the number of 3d
!                    integration points (ze is discarded).
!
                     mint3d=mint3d/2
                  endif
!
                  do l=1,mint3d
                     if(mint3d.eq.4) then
                        xi=gauss2d2(1,l)
                        et=gauss2d2(2,l)
                        weight=1.d0
                        if(quadratic) then
                           call shape8q(xi,et,xl,xsj,xs,shp,iflag)
                        else
                           call shape4q(xi,et,xl,xsj,xs,shp,iflag)
                        endif
                     else
                        call beamintscheme(lakonl,mint3d,ielprop(i),
     &                       prop,l,xi,et,ze,weight)
                        if(j.eq.1) then
                           xil=ze
                           etl=et
                           call shape8q(xil,etl,xl,xsj,xs,shp,iflag)
c                           write(*,*) 'l xi et',l,xil,etl
                        else
                           xil=ze
                           etl=-et
                           call shape8q(xil,etl,xl,xsj,xs,shp,iflag)
c                           write(*,*) 'l xi et',l,xil,etl
                        endif
                     endif
!
!                    local stress tensor
!
                     do m1=1,6
                        s(m1)=0.d0
                        do m2=1,nope
                           s(m1)=s(m1)+shp(4,m2)*yn(m1,node(m2))
c                           if(m1.eq.3)
c     &                     write(*,*) 'l,m1..',l,m1,m2,yn(m1,node(m2)),
c     &                         shp(4,m2)
                        enddo
c                        if(m1.eq.3) write(*,*) 'l,s(3)',l,s(3)
                     enddo
!
!                    local coordinates
!
                     do m1=1,3
                        p(m1)=0.d0
                        do m2=1,nope
                           p(m1)=p(m1)+shp(4,m2)*xl(m1,m2)
                        enddo
                        pcg(m1)=p(m1)-cg(m1)
                     enddo
!
!                    local stress vector on section
!
                     t(1)=(s(1)*xsj(1)+s(4)*xsj(2)+s(5)*xsj(3))*weight
                     t(2)=(s(4)*xsj(1)+s(2)*xsj(2)+s(6)*xsj(3))*weight
                     t(3)=(s(5)*xsj(1)+s(6)*xsj(2)+s(3)*xsj(3))*weight
c                     write(*,*) 'map3dto1d2d'
c                     write(*,*) 'element i,j,l',i,j,l
c                     write(*,*) 't',t(1),t(2),t(3)
c                     write(*,*) 'xsj',xsj(1),xsj(2),xsj(3)
c                     write(*,*) 's1-3',s(1),s(2),s(3)
c                     write(*,*) 's4-6',s(4),s(5),s(6)
!
!                    section forces
!
                     yn(1,node2d)=yn(1,node2d)+
     &                   (e1(1)*t(1)+e1(2)*t(2)+e1(3)*t(3))
                     yn(2,node2d)=yn(2,node2d)+
     &                   (e2(1)*t(1)+e2(2)*t(2)+e2(3)*t(3))
                     yn(3,node2d)=yn(3,node2d)+
     &                   (e3(1)*t(1)+e3(2)*t(2)+e3(3)*t(3))
!
!                    section moments
!
!                    about beam axis
!
                     yn(4,node2d)=yn(4,node2d)+
     &                     (e3(1)*pcg(2)*t(3)+e3(2)*pcg(3)*t(1)+
     &                     e3(3)*pcg(1)*t(2)-e3(3)*pcg(2)*t(1)-
     &                     e3(1)*pcg(3)*t(2)-e3(2)*pcg(1)*t(3))
!
!                    about 2-direction
!
                     yn(5,node2d)=yn(5,node2d)+
     &                     (e2(1)*pcg(2)*t(3)+e2(2)*pcg(3)*t(1)+
     &                     e2(3)*pcg(1)*t(2)-e2(3)*pcg(2)*t(1)-
     &                     e2(1)*pcg(3)*t(2)-e2(2)*pcg(1)*t(3))
!
!                    about 1-direction
!
                     yn(6,node2d)=yn(6,node2d)+
     &                     (e1(1)*pcg(2)*t(3)+e1(2)*pcg(3)*t(1)+
     &                     e1(3)*pcg(1)*t(2)-e1(3)*pcg(2)*t(1)-
     &                     e1(1)*pcg(3)*t(2)-e1(2)*pcg(1)*t(3))
!
!                    components 5 and 6 are switched in the frd
!                    format, so the final order is beam axis,
!                    1-direction and 2-direction, or s12, s23 and s31
!
                  enddo
               enddo
!
            endif
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
               if(.not.force) then
!
!                 taking the mean across the thickness
!
                  if((j.le.4).and.(quadratic)) then
!
!                    end nodes: weights 1/6,2/3 and 1/6
!
                     do l=1,3
                        node3d=kon(indexe+node8(l,j))
                        do k=1,nfield
                           yn(k,node2d)=yn(k,node2d)+
     &                                  yn(k,node3d)*ratioe(l)
                        enddo
                     enddo
                  else
!
!                    middle nodes: weights 1/2,1/2
!
                     do l=1,3,2
                        node3d=kon(indexe+node8(l,j))
                        do k=1,nfield
                           yn(k,node2d)=yn(k,node2d)+yn(k,node3d)/2.d0
                        enddo
                     enddo
                  endif
               else
!
!                 forces must be summed
!
!                 the contribution of each 3d-node should only be taken
!                 once. inum(node3d) is used as marker.
!
                  inum(node2d)=-1
!
                  if((j.le.4).and.(quadratic)) then
!
!                    end nodes of quadratic 2d-elements
!
                     do l=1,3
                        node3d=kon(indexe+node8(l,j))
                        if(inum(node3d).ne.0) cycle
                        inum(node3d)=1
                        do k=1,nfield
                           yn(k,node2d)=yn(k,node2d)+yn(k,node3d)
                        enddo
                     enddo
                  else
!
!                    middle nodes of quadratic 2d-elements
!                    or end nodes of linear 2d-elements
!
                     do l=1,3,2
                        node3d=kon(indexe+node8(l,j))
                        if(inum(node3d).ne.0) cycle
                        inum(node3d)=1
                        do k=1,nfield
                           yn(k,node2d)=yn(k,node2d)+yn(k,node3d)
                        enddo
                     enddo
                  endif
               endif
            enddo
         endif
!
      enddo
!
!     taking the mean of nodal contributions coming from different
!     elements having the node in common
!     restoring inum for the 3d-nodes to zero
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         lakonl=lakon(i)
         if((lakonl(7:7).eq.' ').or.(lakonl(7:7).eq.'I').or.
     &      (lakonl(1:1).ne.'C')) cycle
         indexe=ipkon(i)
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
               if(inum(node2d).lt.0) then
                  inum(node2d)=-inum(node2d)
                  do k=1,nfield
                     yn(k,node2d)=yn(k,node2d)/inum(node2d)
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
            do j=1,jmax
               node2d=kon(indexe2d+j)
               if(inum(node2d).lt.0) then
                  inum(node2d)=-inum(node2d)
                  do k=1,nfield
                     yn(k,node2d)=yn(k,node2d)/inum(node2d)
                  enddo
               endif
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
               if(inum(node2d).lt.0) then
                  inum(node2d)=-inum(node2d)
                  do k=1,nfield
                     yn(k,node2d)=yn(k,node2d)/inum(node2d)
                  enddo
               endif
            enddo
         endif
!     
!        inactivating the 3d expansion nodes of 1d/2d elements
!        in case forces are mapped this field is used to ensure
!        that the forces in the 3d-nodes are mapped only once onto the
!        2d-nodes
!
         do j=1,indexe2d-indexe
            inum(kon(indexe+j))=0
         enddo
!
      enddo
!
!     beam section forces in the middle nodes
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         lakonl=lakon(i)
         if((lakonl(7:7).eq.' ').or.(lakonl(7:7).eq.'I').or.
     &      (lakonl(1:1).ne.'C')) cycle
         indexe=ipkon(i)
!
!        not relevant for linear elements
!
         if((lakonl(4:4).eq.'6').or.(lakonl(4:4).eq.'8')) cycle
!
         if(lakonl(7:7).eq.'B') then
            indexe2d=indexe+20
            if(cflag.eq.'M') then
!
!              section forces in the middle node are the mean
!              of those in the end nodes
!
               nodea=kon(indexe2d+1)
               nodeb=kon(indexe2d+2)
               nodec=kon(indexe2d+3)
               inum(nodeb)=1
               do j=1,6
                  yn(j,nodeb)=yn(j,nodeb)+(yn(j,nodea)+yn(j,nodec))/2.d0
               enddo
!
            endif
         endif
!
      enddo
!
      return
      end
