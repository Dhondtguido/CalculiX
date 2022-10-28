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
      subroutine zienzhu(co,nk,kon,ipkon,lakon,ne,stn,
     & ipneigh,neigh,sti,mi)
!
!     modified zienkiewicz zhu pointwise error estimator
!
!     author: Sascha Merz
!
      implicit none
!
      character*8 lakon(*)
!
      integer maxmid
!
!     maximal number of midnodes surrounding a patch base node
!
      parameter(maxmid=400)
!
      integer kon(*),nk,ne,i,j,ipkon(*),indexe,nodebase,
     & k,ipneigh(*),neigh(2,*),ifree,node,ielem,ielem1,index,index1,m,
     & jj,mi(*),ii,ncount,nope,itypflag,inum(nk),nenei20(3,8),maxcommon,
     & icommon,idxnode,iecount,index2,lneigh8a(7,8),ipoints,icont,
     & nmids(maxmid),nelem(3),nelemv,iavflag,members(ne),ielem2,
     & linpatch,iterms,inodes(4),iaddelem,ielidx,nenei10(3,4),iscount,
     & idxnode1,maxcommon1
!
      real*8 co(3,*),stn(6,*),sti(6,mi(1),*),angle,scpav(6,nk),
     & angmax
!
      data lneigh8a /2,3,4,5,6,7,8,1,3,4,5,6,7,8,1,2,4,5,6,7,8,
     &               1,2,3,5,6,7,8,1,2,3,4,6,7,8,1,2,3,4,5,7,8,
     &               1,2,3,4,5,6,8,1,2,3,4,5,6,7/
!
      data nenei10 /5,7,8,5,6,9,6,7,10,8,9,10/
!
      data nenei20 /9,12,17,9,10,18,10,11,19,11,12,20,
     &              13,16,17,13,14,18,14,15,19,15,16,20/
!
      write(*,*) 'Estimating the stress errors'
      write(*,*)
!
!     initialization
!
      ifree=0
      do i=1,nk
         ipneigh(i)=0
         inum(i)=0
         do j=1,6
            scpav(j,i)=0.d0
         enddo
      enddo
!
!     build neighbour list
!
      do i=1,ne
         indexe=ipkon(i)
         if(lakon(i)(1:4).eq.'C3D4') then
            nope=4
         elseif(lakon(i)(1:4).eq.'C3D6') then
            nope=6
         elseif(lakon(i)(1:4).eq.'C3D8') then
            nope=8
         elseif(lakon(i)(1:5).eq.'C3D10') then
            nope=4
         elseif(lakon(i)(1:5).eq.'C3D15') then
            nope=6
         elseif(lakon(i)(1:5).eq.'C3D20') then
            nope=8
         else 
            cycle
         endif
         do j=1,nope
            node=kon(indexe+j)
            ifree=ifree+1
            neigh(2,ifree)=ipneigh(node)
            ipneigh(node)=ifree
            neigh(1,ifree)=i
         enddo
      enddo
!
      patches:do nodebase=1,nk
         if(ipneigh(nodebase).eq.0) cycle patches
         index=ipneigh(nodebase)
!
!        initialization
!
         do j=1,maxmid
            nmids(j)=0
         enddo
         do j=1,3
            nelem(j)=0
         enddo
         idxnode=nodebase
         linpatch=0;ipoints=0;iavflag=0;itypflag=0
!     
!        analyze neighbour structure
!     
         do
            ielem=neigh(1,index)
!     
            if(lakon(ielem)(1:5).eq.'C3D20') then
               nelem(1)=nelem(1)+1
            elseif(lakon(ielem)(1:5).eq.'C3D10') then
               nelem(2)=nelem(2)+1
            elseif(lakon(ielem)(1:4).eq.'C3D8') then
               nelem(3)=nelem(3)+1
            endif
!     
            if(neigh(2,index).eq.0)exit
            index=neigh(2,index)
         enddo
!     
!        itypflag 1: using hex20 estimator
!        itypflag 2: using tet10 estimator
!        itypflag 3: using hex8  estimator
!
         if(nelem(1).gt.0) then
            itypflag=1
         elseif(nelem(1).eq.0) then
            if(nelem(2).gt.0) then
               itypflag=2
            elseif(nelem(2).eq.0) then
               if(nelem(3).gt.0) then
                  itypflag=3
               else
                  itypflag=0
               endif
            endif
         endif
!
!        case 1: element type not supported
!            
         if(itypflag.eq.0) then
            write(*,*) '*WARINING in estimator: Elements of node',
     &           nodebase,' cannot be used for error estimation.'
            do j=1,6
               scpav(j,nodebase)=stn(j,nodebase)
            enddo
!
!................ case 2: using hex estimator...........................
!                                                                      !
!                                                                      !         
!
         elseif(itypflag.eq.1.or.itypflag.eq.3) then
!     
!           get spaceangle
!
            call angsum(lakon,kon,ipkon,neigh,ipneigh,
     &           co,nodebase,itypflag,angle)   
!
!           determine surface structure, if surface node
!
            if(angle.lt.12.56535d0) then
               call chksurf(lakon,kon,ipkon,neigh,ipneigh,co,
     &              itypflag,nodebase,icont,iscount,angmax)
            endif
!
            if(itypflag.eq.1) then
               nelemv=nelem(1)
            elseif(itypflag.eq.3) then
               nelemv=nelem(3)
            endif
!
!           if the elements neighbouring the node are not appropriate
!           to form a valid patch, then look for an alternative patch
!           base node.
!
            if(
     &         angle.lt.8.64d0
     &         .and.
     &         iscount.eq.nelemv
     &         .and.
     &         angmax.lt.3.d1       
     &         .or.
     &         nelemv.lt.5
     &        ) then
c               write(*,*) 'Node',nodebase,' is not valid as',
c     &              ' a patch base node. Seeking another',
c     &              ' patch base node.'
!     
!              node is on surface. seeking a node within volume having 
!              the most elements in common
!     
               index=ipneigh(nodebase)
!
!              index points to the location in neigh, where the element
!              number is given
!
               maxcommon=0
               elements: do
                  if(index.eq.0) exit elements
                  ielem=neigh(1,index)
!
                  if(.not.((lakon(ielem)(1:5).eq.'C3D20'
     &                 .and.itypflag.eq.1)
     &                 .or.(lakon(ielem)(1:4).eq.'C3D8'
     &                 .and.itypflag.eq.3))) then
                     index=neigh(2,index)
                     cycle elements
                  endif
!
!                 ielem: element number
!
                  indexe=ipkon(ielem)
!     
!                 find element # corresp. 2 global node #
!                  
                  do m=1,8
                     if(kon(indexe+m).eq.nodebase) exit
                  enddo
!
!                 for every corner node of a neighbouring
!                 element count the common elements
!
                  nodes: do j=1,7
                     k=lneigh8a(j,m)
!
!                    get global node # for neighbouring node
!                    and count the elements in common
!
                     node=kon(ipkon(ielem)+k)
!
!                    skip, if this node is on surface as well
!
                     call  angsum(lakon,kon,ipkon,neigh,ipneigh,
     &                    co,node,itypflag,angle)
                     if(angle.lt.12.56535d0) cycle nodes
!
!                    go through neighbouring elements to find common
!                    elements
!
                     icommon=0
                     index1=ipneigh(nodebase)
                     outer: do
                        if(index1.eq.0) exit outer
                        ielem1=neigh(1,index1)
                        if(.not.((lakon(ielem1)(1:5).eq.'C3D20'
     &                       .and.itypflag.eq.1)
     &                       .or.(lakon(ielem1)(1:4).eq.'C3D8'
     &                       .and.itypflag.eq.3))) then
                           index1=neigh(2,index1)
                           cycle outer
                        endif
                        index2=ipneigh(node)
                        inner: do
                           if(index2.eq.0) exit inner
                           ielem2=neigh(1,index2)
                           if(.not.((lakon(ielem2)(1:5).eq.'C3D20'
     &                          .and.itypflag.eq.1)
     &                          .or.(lakon(ielem2)(1:4).eq.'C3D8'
     &                          .and.itypflag.eq.3))) then
                              index2=neigh(2,index2)
                              cycle inner
                           endif
!
!                          check, if element is common
!
                           if(ielem1.eq.ielem2) then
                              icommon=icommon+1
                           endif
                           index2=neigh(2,index2)
                        enddo inner
                        index1=neigh(2,index1)
                     enddo outer
!
!                    check, if the last two nodes have the most
!                    common elements
!
                     if(icommon.gt.maxcommon) then
                        maxcommon=icommon
                        idxnode=node
                     endif
                  enddo nodes
!
                  index=neigh(2,index)
               enddo elements
! 
!              however, if there is a node on the surface, having more
!              elements in common and the original base node is not on
!              a free surface, then use this node as basenode.
!
               if(icont.eq.0) then
                  index=ipneigh(nodebase)
!
                  maxcommon1=0
                  elements1: do
                     if(index.eq.0) exit elements1
                     ielem=neigh(1,index)
!
                     if(.not.((lakon(ielem)(1:5).eq.'C3D20'
     &                    .and.itypflag.eq.1)
     &                    .or.(lakon(ielem)(1:4).eq.'C3D8'
     &                    .and.itypflag.eq.3))) then
                        index=neigh(2,index)
                        cycle elements1
                     endif
!
                     indexe=ipkon(ielem)
!     
                     do m=1,8
                        if(kon(indexe+m).eq.nodebase) exit
                     enddo
!     
                     nodes1: do j=1,7
                        k=lneigh8a(j,m)
!
                        node=kon(ipkon(ielem)+k)
!
                        icommon=0
                        index1=ipneigh(nodebase)
                        outer1: do
                           if(index1.eq.0) exit outer1
                           ielem1=neigh(1,index1)
                           if(.not.((lakon(ielem1)(1:5).eq.'C3D20'
     &                          .and.itypflag.eq.1)
     &                          .or.(lakon(ielem1)(1:4).eq.'C3D8'
     &                          .and.itypflag.eq.3))) then
                              index1=neigh(2,index1)
                              cycle outer1
                           endif
                           index2=ipneigh(node)
                           inner1: do
                              if(index2.eq.0) exit inner1
                              ielem2=neigh(1,index2)
                              if(.not.((lakon(ielem2)(1:5).eq.'C3D20'
     &                             .and.itypflag.eq.1)
     &                             .or.(lakon(ielem2)(1:4).eq.'C3D8'
     &                             .and.itypflag.eq.3))) then
                                 index2=neigh(2,index2)
                                 cycle inner1
                              endif
!
                              if(ielem1.eq.ielem2) then
                                 icommon=icommon+1
                              endif
                              index2=neigh(2,index2)
                           enddo inner1
                           index1=neigh(2,index1)
                        enddo outer1
!
                        if(icommon.gt.maxcommon1) then
                           maxcommon1=icommon
                           idxnode1=node
                        endif
                     enddo nodes1
!     
                     index=neigh(2,index)
                  enddo elements1
!
                  if(maxcommon1.gt.maxcommon) then
                     idxnode=idxnode1
                  endif
               endif
               index=ipneigh(idxnode)
               nelemv=0
               do
                  if(index.eq.0) exit
                  ielem=neigh(1,index)
                  if(.not.((lakon(ielem)(1:5).eq.'C3D20'
     &                 .and.itypflag.eq.1)
     &                 .or.(lakon(ielem)(1:4).eq.'C3D8'
     &                 .and.itypflag.eq.3))) then
                     index=neigh(2,index)
                     cycle
                  endif
                  nelemv=nelemv+1
                  index=neigh(2,index)
               enddo
!
!              verify
!
               call angsum(lakon,kon,ipkon,neigh,ipneigh,
     &              co,idxnode,itypflag,angle)   
               if(angle.lt.12.56535d0) then
                  call chksurf(lakon,kon,ipkon,neigh,ipneigh,co,
     &                 itypflag,idxnode,icont,iscount,angmax)
               endif
               if(
     &              angle.lt.8.64d0
     &              .and.
     &              iscount.eq.nelemv
     &              .and.
     &              angmax.lt.3.d1       
     &              .or.
     &              nelemv.lt.5
     &            ) then
                  write(*,*) '*WARNING in estimator: Patch not',
     &                 ' appropriate for patch recovery,'
                  write(*,*) '                       using',
     &                 ' average of sampling point values', nodebase
                  iavflag=1
               endif
!     
c               write(*,*) 'Alternative node is',idxnode
            endif
!
            index=ipneigh(idxnode)
!
!           determine patch elements (members), the number of
!           patch elements (linpatch) and the number of 
!           sampling points (ipoints) in the patch
!
            do
               if(index.eq.0) exit
               ielem=neigh(1,index)
                  if(.not.((lakon(ielem)(1:5).eq.'C3D20'
     &                 .and.itypflag.eq.1)
     &                 .or.(lakon(ielem)(1:4).eq.'C3D8'
     &                 .and.itypflag.eq.3))) then
                  index=neigh(2,index)
                  cycle
               endif
               if(lakon(ielem)(1:4).eq.'C3D8') then
                  if(lakon(ielem)(5:5).eq.'R') then
                     ipoints=ipoints+1
                  else
                     ipoints=ipoints+8
                  endif
               elseif(lakon(ielem)(1:5).eq.'C3D20') then
                  if(lakon(ielem)(6:6).eq.'R') then
                     ipoints=ipoints+8
                  else
                     ipoints=ipoints+27
                  endif
               endif
               linpatch=linpatch+1
               members(linpatch)=ielem
               index=neigh(2,index)
            enddo
!
            if(itypflag.eq.1) iterms=20
            if(itypflag.eq.3) iterms=4
!
!           evaluate patch for patch base node (nodebase)
!
            call patch(iterms,nodebase,sti,scpav,mi(1),kon,ipkon,
     &           ipoints,members,linpatch,co,lakon,iavflag)
            inum(nodebase)=1
!     
!           midnodes
!
            if(itypflag.eq.1) then
               k=1
               hexelements: do ielidx=1,linpatch
                  ielem=members(ielidx)
!     
!                 find out index of base node in kon
!
                  do m=1,8
                     if(nodebase.eq.kon(ipkon(ielem)+m)) exit
                     if(m.eq.8) then
                        cycle hexelements
                     endif
                  enddo
!
                  hexnodes: do j=1,3
!
!                    if node already in patch, skip
!
                     if(k.gt.1) then
                        do ii=1,k-1
                           if(nmids(ii)
     &                          .eq.kon(ipkon(ielem)+nenei20(j,m))) then
                              cycle hexnodes
                           endif
                        enddo
                     endif
                     nmids(k)=kon(ipkon(ielem)+nenei20(j,m))
                     inum(nmids(k))=inum(nmids(k))+1
                     call patch(iterms,nmids(k),sti,scpav,mi(1),kon,
     &                    ipkon,ipoints,members,linpatch,co,lakon,
     &                    iavflag)
                     k=k+1
                     if(k.gt.maxmid) then
                        write(*,*) '*ERROR in estimator: array size',
     &                       ' for midnodes exceeded'
                        exit hexelements
                     endif
                  enddo hexnodes
               enddo hexelements
            endif
!
!................ case 3: using tet estimator..........................
!                                                                      !
!                                                                      !         
         elseif(itypflag.eq.2) then
!
!
!     for the tetrahedral elements it could happen, that additional
!     elements have to be added to the patch. 
!     the information already obtained is the number of elements
!     neighbouring the node.
!     now a check should be undertaken, if the shape of the initial
!     patch is useful. 
!
!
!        case 1: the patch has a conical shape. this occurs, if
!        if all neighbouring elements of the 
!        patch base node have another vertice node in common 
!        (integration points are then not reasonably distributed).
!        the patch looks then like a cone or pyramid:
!
!                       x
!                      /|\
!                     / | \
!                    /  |  \
!                   /   |   \
!                  /    |    \
!        .........x_____o_____x.................
!
!        ...... surface
!    
!        x patch member node
!        o patch base node
!
!        a alternative patch base node has to be found, if the following
!        conditions occur:
!
!        1. the node is not a volume node and there is a conical
!           patch with a even number of patch members
!
!        2. the node has only one neighbouring element. then it is
!           most likely a corner node and a good patch will be created
!           if another edgenode of this element is used as patch base node
!         
!        3. the node is located at a free surface and the number of
!           neighbouring elements is odd and the elements have conical
!           shape
!
            call angsum(lakon,kon,ipkon,neigh,ipneigh,co,nodebase,
     &           itypflag,angle)
!
!           if the node is on the surface, then check, if the surfaces
!           adjacent to the node have normal vectors with an angle of
!           maximal 10 degree (free surface is assumed, otherwise edge
!           or corner)
!
            if(angle.lt.12.56535d0) then
               call chksurf(lakon,kon,ipkon,neigh,ipneigh,co,
     &              itypflag,nodebase,icont,iscount,angmax)
            endif
!
            nelemv=nelem(2)
!
            icommon=0
!
            if(
     &         angle.lt.12.56535d0
     &         .and.(
     &               mod(nelemv,2).eq.0
     &               .or.
     &               nelemv.eq.1
     &               )
     &          .or.(
     &               icont.eq.1
     &               .and.
     &               mod(nelemv,2).ne.0
     &              )
     &        ) then
!     
!              searching for another common vertex node
!
               index=ipneigh(nodebase)
               idxnode=0
               node=0
               conical: do
                  if(index.eq.0) exit
                  ielem=neigh(1,index)
                  if(.not.(lakon(ielem)(1:5).eq.'C3D10'
     &                 .and.itypflag.eq.2)) then
                     index=neigh(2,index)
                     cycle
                  endif
                  do j=1,4
                     node=kon(ipkon(ielem)+j)
                     if(node.eq.nodebase) cycle
!     
!                    is this node a member of all patch elements?
!
                     index1=ipneigh(nodebase)
                     ncount=0
                     do
                        if(index1.eq.0) exit
                        ielem1=neigh(1,index1)
                        if(.not.(lakon(ielem1)(1:5).eq.'C3D10'
     &                       .and.itypflag.eq.2)) then
                           index1=neigh(2,index1)
                           cycle
                        endif
                        do jj=1,4
                           if(idxnode.eq.kon(ipkon(ielem1)+jj)
     &                          .and.idxnode.ne.0) cycle
                           if(node.eq.kon(ipkon(ielem1)+jj))
     &                          ncount=ncount+1
                        enddo
                        index1=neigh(2,index1)
                     enddo
                     if(ncount.eq.nelemv) then
                        icommon=icommon+1
                        idxnode=node
                     endif
                  enddo
                  index=neigh(2,index)
               enddo conical
!     
               if(icommon.gt.0) then
c                  write(*,*) 'conical patch found, node',nodebase,
c     &                 ' using node',idxnode,' as base node instead'
!
!                 count elements
!
                  index=ipneigh(idxnode)
                  nelemv=0
                  do
                     if(index.eq.0) exit
                     ielem=neigh(1,index)
                     if(.not.(lakon(ielem)(1:5).eq.'C3D10'
     &                    .and.itypflag.eq.2)) then
                        index=neigh(2,index)
                        cycle
                     endif
                     nelemv=nelemv+1
                     index=neigh(2,index)
                  enddo
!               
               else
                  idxnode=nodebase
               endif   
            endif
!
!           add neighbouring elements to patch firstly
!
            linpatch=0
            index=ipneigh(idxnode)
            do
               if(index.eq.0) exit
               ielem=neigh(1,index)
               if(.not.(lakon(ielem)(1:5).eq.'C3D10'
     &              .and.itypflag.eq.2)) then
                  index=neigh(2,index)
                  cycle
               endif
               linpatch=linpatch+1
               members(linpatch)=ielem
               index=neigh(2,index)
            enddo
!     
            if(nelemv.gt.0.and.nelemv.lt.4) then
c               write(*,*) 'node',nodebase,' has',nelemv,
c     &              ' neighbouring elements'
!     
!              patch has to be extended
!     
!              add more elements to patch, where the elements
!              have to share a face with the initial patch elements
!
               index=ipneigh(idxnode)
               iecount=0
               do
                  if(index.eq.0) exit 
                  ielem=neigh(1,index)
                  if(.not.(lakon(ielem)(1:5).eq.'C3D10'
     &                 .and.itypflag.eq.2)) then
                     index=neigh(2,index)
                     cycle
                  endif
!
!                 check every element in the neighbour list
!                 if there is any element in the element list
!                 having three nodes in common (is a common
!                 surface). 
!               
!                 save element nodes
!     
                  do j=1,4
                     inodes(j)=kon(ipkon(ielem)+j)
                  enddo
!
!                 loop over each element in the model and
!                 check against patch element
!
                  tetloop: do k=1,ne
                     if(.not.(lakon(ielem)(1:5).eq.'C3D10'
     &                 .and.itypflag.eq.2)) then
                        cycle
                     endif
!
!                    skip counting, if element is already in patch
!
                     index1=ipneigh(idxnode)
                     do
                        if(index1.eq.0) exit
                        ielem1=neigh(1,index1)
                        if(.not.(lakon(ielem1)(1:5).eq.'C3D10'
     &                       .and.itypflag.eq.2)) then
                           index1=neigh(2,index1)
                           cycle
                        endif
                        if(ielem1.eq.k) cycle tetloop
                        index1=neigh(2,index1)
                     enddo
!     
                     ncount=0
                     do j=1,4
                        do jj=1,4
                           if(inodes(jj).eq.kon(ipkon(k)+j)) then
                              ncount=ncount+1
                           endif
                           if(ncount.gt.2) then
                              linpatch=linpatch+1
                              iecount=iecount+1
                              members(linpatch)=k
!
!                             keep number of 2nd found element in mind
!
                              if(nelemv+iecount.eq.2)
     &                             iaddelem=k
                              cycle tetloop
                           endif
                        enddo
                     enddo
                  enddo tetloop
                  index=neigh(2,index)
               enddo
!
!              if there are still just two elements, treat
!              those elements as initial patch and extend
!
               if(nelemv+iecount.eq.2) then
!
!                 search again all elements
!
                  do j=1,4
                     inodes(j)=kon(ipkon(iaddelem)+j)
                  enddo
                  loop3: do k=1,ne
                     if(.not.(lakon(k)(1:5).eq.'C3D10'
     &                 .and.itypflag.eq.2)) then
                        cycle loop3
                     endif
                     index=ipneigh(idxnode)
                     do
                        if(index.eq.0) exit
                        ielem=neigh(1,index)
                        if(.not.(lakon(ielem)(1:5).eq.'C3D10'
     &                       .and.itypflag.eq.2)) then
                           exit
                        endif
                        index=neigh(2,index)
                     enddo
                     if(ielem.eq.k.or.iaddelem.eq.k) cycle loop3
                     ncount=0
                     do j=1,4
                        do jj=1,4
                           if(inodes(jj).eq.kon(ipkon(k)+j))
     &                          ncount=ncount+1
                           if(ncount.gt.2) then
                              linpatch=linpatch+1
                              iecount=iecount+1
                              members(linpatch)=k
                              cycle loop3
                           endif
                        enddo
                     enddo
                  enddo loop3
               endif
               nelemv=nelemv+iecount
c               write(*,*) 'now node',nodebase,' has',nelemv,
c     &              ' elements'
!
!              verify
!               
               if(nelemv.lt.5) then
                  write(*,*) '*WARNING in estimator: Patch not',
     &                 ' appropriate for patch recovery,'
                  write(*,*) '                       using',
     &                 ' average of sampling point values:', nodebase
                  iavflag=1
               endif
            endif
!
!           nodal stresses for patch
!
!           determine degree of patch polynomial
!     
            ipoints=linpatch*4
            iterms=10
            if(ipoints.ge.17.and.ipoints.lt.21) then
               iterms=10
            elseif(ipoints.ge.21.and.ipoints.lt.63) then
               iterms=11
            elseif(ipoints.ge.63) then
               iterms=17
            endif
!
c            write(*,*) 'using',iterms,' terms for patch',nodebase
!     
!           patch for patch base node
!
            call patch(iterms,nodebase,sti,scpav,mi(1),kon,ipkon,
     &           ipoints,members,linpatch,co,lakon,iavflag)
            inum(nodebase)=1
!
!           midnodes
! 
            k=1
            tetelements: do ielidx=1,linpatch
               ielem=members(ielidx)
!
!                 find out index of base node in kon
!
               do m=1,4
                  if(nodebase.eq.kon(ipkon(ielem)+m)) exit
                  if(m.eq.4) then
                     cycle tetelements
                  endif
               enddo
!
               tetnodes: do j=1,3
!
!                    if node already in patch, skip
!
                  if(k.gt.1) then
                     do ii=1,k-1
                        if(nmids(ii)
     &                       .eq.kon(ipkon(ielem)+nenei10(j,m))) then
                           cycle tetnodes
                        endif
                     enddo
                  endif
                  nmids(k)=kon(ipkon(ielem)+nenei10(j,m))
                  inum(nmids(k))=inum(nmids(k))+1
                  call patch(iterms,nmids(k),sti,scpav,mi(1),kon,
     &                 ipkon,ipoints,members,linpatch,co,lakon,
     &                 iavflag)
                  k=k+1
                  if(k.gt.maxmid) then
                     write(*,*) '*ERROR in estimator: array size',
     &                    ' for midnodes exceeded'
                     exit tetelements
                  endif
               enddo tetnodes
            enddo tetelements
         endif
      enddo patches
!.......................................................................
!
      do i=1,nk
         if(inum(i).gt.0) then
            do j=1,6
               stn(j,i)=scpav(j,i)/inum(i)
            enddo
         else
            do j=1,6
               stn(j,i)=0.d0
            enddo
         endif
      enddo
!
      return
      end
