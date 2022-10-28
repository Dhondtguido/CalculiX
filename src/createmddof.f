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
      subroutine createmddof(imddof,nmddof,istartset,iendset,
     &            ialset,nactdof,ithermal,mi,imdnode,nmdnode,ikmpc,
     &            ilmpc,ipompc,nodempc,nmpc,
     &            imdmpc,nmdmpc,imdboun,nmdboun,ikboun,nboun,
     &            nset,ntie,tieset,set,lakon,kon,ipkon,labmpc,
     &            ilboun,filab,prlab,prset,nprint,ne,cyclicsymmetry)
!
!     creating a set imddof containing the degrees of freedom
!     selected by the user for modal dynamic calculations. The
!     solution will be calculated for these dof's only in order
!     to speed up the calculation.
!
      implicit none
!
      logical nodeslavsurf
!
      character*6 prlab(*)
      character*8 lakon(*)
      character*20 labmpc(*)
      character*81 tieset(3,*),rightset,set(*),slavset,noset,prset(*)
      character*87 filab(*)
!
      integer imddof(*),nmddof,nrset,istartset(*),iendset(*),mi(*),
     &  ialset(*),nactdof(0:mi(2),*),node,ithermal(*),j,k,l,
     &  ikmpc(*),ilmpc(*),ipompc(*),nodempc(3,*),nmpc,id,
     &  imdnode(*),nmdnode,imdmpc(*),nmdmpc,nprint,ipos,
     &  imdboun(*),nmdboun,ikboun(*),nboun,indexe1,indexe,islav,
     &  jface,nset,ntie,nnodelem,nope,nodef(8),nelem,nface,imast,
     &  ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),kon(*),
     &  ipkon(*),i,ilboun(*),nlabel,ne,cyclicsymmetry
!
!     nodes per face for hex elements
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
!
!     nodes per face for tet elements
!
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
!
!     nodes per face for linear wedge elements
!
      data ifacew1 /1,3,2,0,
     &             4,5,6,0,
     &             1,2,5,4,
     &             2,3,6,5,
     &             3,1,4,6/
!
!     nodes per face for quadratic wedge elements
!
      data ifacew2 /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             3,1,4,6,9,13,12,15/
!
      data nlabel /55/
!
!     if 1d/2d elements are part of the mesh, no node selection
!     is performed (because of the renumbering due to the
!     expansion node selection is excessively difficult)
!
      do i=1,ne
         if((lakon(i)(7:7).eq.'E').or.
     &      (lakon(i)(7:7).eq.'S').or.
     &      ((lakon(i)(7:7).eq.'A').and.(lakon(i)(1:1).eq.'C')).or.
     &      (lakon(i)(7:7).eq.'L').or.
     &      (lakon(i)(7:7).eq.'B')) then
            nmdnode=0
            nmddof=0
            nmdboun=0
            nmdmpc=0
            return
         endif
      enddo
!
!     storing the nodes, dofs, spcs and mpcs for which *NODE FILE 
!     or *EL FILE was selected
!
      do i=1,nlabel
!
!        CDIS,CSTR und CELS are not taken into account:
!        contact area is treated separately (no set can
!        be specified for CDIS, CSTR und CELS)
!
         if((i.eq.26).or.(i.eq.27)) cycle
!
         if(filab(i)(1:1).ne.' ') then
            read(filab(i)(7:87),'(a81)') noset
c            nrset=0
c            do k=1,nset
c               if(set(k).eq.noset) then
c                  nrset=k
c                  exit
c               endif
c            enddo
            call cident81(set,noset,nset,id)
            nrset=0
            if(id.gt.0) then
              if(noset.eq.set(id)) then
                nrset=id
              endif
            endif
!
!           if output for all nodes is selected, use
!           of imdnode is deactivated
!
            if(nrset.eq.0) then
               if(cyclicsymmetry.eq.1) then
                  write(*,*) '*ERROR in createmddof: in a cylic'
                  write(*,*) '       symmetric modal dynamic or'
                  write(*,*) '       steady static dynamics calculation'
                  write(*,*) '       a node set MUST be defined on each'
                  write(*,*) '       *NODE FILE, *NODE OUTPUT, *EL FILE'
                  write(*,*) '       or *ELEMENT OUTPUT card.'
                  write(*,*) '       Justification: in a steady state'
                  write(*,*) '       dynamics calculation with cyclic'
                  write(*,*) '       symmetry the segment is expanded'
                  write(*,*) '       into 360 °. Storing results for'
                  write(*,*) '       this expansion may lead to huge'
                  write(*,*) '       frd-files. Specifying a set can'
                  write(*,*) '       reduce this output.'
                  call exit(201)
               endif
               nmdnode=0
               nmddof=0
               nmdboun=0
               nmdmpc=0
               return
            endif
!
!           adding the nodes belonging to nrset
!
            do j=istartset(nrset),iendset(nrset)
               if(ialset(j).gt.0) then
                  node=ialset(j)
                  call addimd(imdnode,nmdnode,node)
                  if(ithermal(1).ne.2) then
                     do k=1,3
                        call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                   nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                   nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &                   ikboun,nboun,ilboun)
                     enddo
                  else
                     k=0
                     call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,ikboun,
     &                nboun,ilboun)
                  endif
               else
                  node=ialset(j-2)
                  do
                     node=node-ialset(j)
                     if(node.ge.ialset(j-1)) exit
                     call addimd(imdnode,nmdnode,node)
                     if(ithermal(1).ne.2) then
                        do k=1,3
                           call addimdnodedof(node,k,ikmpc,ilmpc,
     &                      ipompc,nodempc,nmpc,imdnode,nmdnode,imddof,
     &                      nmddof,nactdof,mi,imdmpc,nmdmpc,imdboun,
     &                      nmdboun,ikboun,nboun,ilboun)
                        enddo
                     else
                        k=0
                        call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                   nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                   nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &                   ikboun,nboun,ilboun)
                     endif
                  enddo
               endif
            enddo
!
         endif
      enddo
!
!     storing the nodes, dofs, spcs and mpcs for which *NODE PRINT 
!     was selected
!
      do i=1,nprint
         if((prlab(i)(1:4).eq.'U   ').or.
     &        (prlab(i)(1:4).eq.'NT  ').or.
     &        (prlab(i)(1:4).eq.'RF  ').or.
     &        (prlab(i)(1:4).eq.'RFL ').or.
     &        (prlab(i)(1:4).eq.'PS  ').or.
     &        (prlab(i)(1:4).eq.'PN  ').or.
     &        (prlab(i)(1:4).eq.'MF  ').or.
     &        (prlab(i)(1:4).eq.'V   ')) then
            noset=prset(i)
c            nrset=0
c            do k=1,nset
c               if(set(k).eq.noset) then
c                  nrset=k
c                  exit
c               endif
c            enddo
            call cident81(set,noset,nset,id)
            nrset=0
            if(id.gt.0) then
              if(noset.eq.set(id)) then
                nrset=id
              endif
            endif
!
!           adding the nodes belonging to nrset
!
            do j=istartset(nrset),iendset(nrset)
               if(ialset(j).gt.0) then
                  node=ialset(j)
                  call addimd(imdnode,nmdnode,node)
                  if(ithermal(1).ne.2) then
                     do k=1,3
                        call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                   nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                   nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &                   ikboun,nboun,ilboun)
                     enddo
                  else
                     k=0
                     call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,ikboun,
     &                nboun,ilboun)
                  endif
               else
                  node=ialset(j-2)
                  do
                     node=node-ialset(j)
                     if(node.ge.ialset(j-1)) exit
                     call addimd(imdnode,nmdnode,node)
                     if(ithermal(1).ne.2) then
                        do k=1,3
                           call addimdnodedof(node,k,ikmpc,ilmpc,
     &                      ipompc,nodempc,nmpc,imdnode,nmdnode,imddof,
     &                      nmddof,nactdof,mi,imdmpc,nmdmpc,imdboun,
     &                      nmdboun,ikboun,nboun,ilboun)
                        enddo
                     else
                        k=0
                        call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                   nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                   nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &                   ikboun,nboun,ilboun)
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
!
!     check whether all contact slave and master nodes (and corresponding
!     dofs, spcs and mpcs were selected
!
      do i=1,ntie
!     
!     check for contact conditions
!     'C' are active contact conditions
!     '-' are temporarily deactivated contact conditions
!     
         if((tieset(1,i)(81:81).eq.'C').or.
     &      (tieset(1,i)(81:81).eq.'-')) then
            rightset=tieset(3,i)
!     
!     determining the master surface
!     
c            do j=1,nset
c               if(set(j).eq.rightset) exit
c            enddo
            call cident81(set,rightset,nset,id)
            j=nset+1
            if(id.gt.0) then
              if(rightset.eq.set(id)) then
                j=id
              endif
            endif
            if(j.gt.nset) then
               write(*,*) '*ERROR in createmddof: master surface',
     &              rightset
               write(*,*) '       does not exist'
               call exit(201)
            endif
            imast=j
!     
            do j=istartset(imast),iendset(imast)
!     
               nelem=int(ialset(j)/10.d0)
               jface=ialset(j)-10*nelem
!     
               indexe=ipkon(nelem)
!     
               if(lakon(nelem)(4:4).eq.'2') then
                  nnodelem=8
                  nface=6
               elseif(lakon(nelem)(4:4).eq.'8') then
                  nnodelem=4
                  nface=6
               elseif(lakon(nelem)(4:5).eq.'10') then
                  nnodelem=6
                  nface=4
               elseif(lakon(nelem)(4:4).eq.'4') then
                  nnodelem=3
                  nface=4
               elseif(lakon(nelem)(4:5).eq.'15') then
                  if(jface.le.2) then
                     nnodelem=6
                  else
                     nnodelem=8
                  endif
                  nface=5
                  nope=15
               elseif(lakon(nelem)(4:4).eq.'6') then
                  if(jface.le.2) then
                     nnodelem=3
                  else
                     nnodelem=4
                  endif
                  nface=5
                  nope=6
               else
                  cycle
               endif
!     
!     determining the master nodes 
!     
               if(nface.eq.4) then
                  do k=1,nnodelem
                     nodef(k)=kon(indexe+ifacet(k,jface))
                  enddo
               elseif(nface.eq.5) then
                  if(nope.eq.6) then
                     do k=1,nnodelem
                        nodef(k)=kon(indexe+ifacew1(k,jface))
                     enddo
                  elseif(nope.eq.15) then
                     do k=1,nnodelem
                        nodef(k)=kon(indexe+ifacew2(k,jface))
                     enddo
                  endif
               elseif(nface.eq.6) then
                  do k=1,nnodelem
                     nodef(k)=kon(indexe+ifaceq(k,jface))
                  enddo
               endif
!
               do l=1,nnodelem
                  node=nodef(l)
                  call addimd(imdnode,nmdnode,node)
                  if(ithermal(1).ne.2) then
                     do k=1,3
                        call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                       nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                       nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &                       ikboun,nboun,ilboun)
                     enddo
                  endif
               enddo
            enddo
!
!           determining the slave nodes 
!
            slavset=tieset(2,i)
!
!           check whether facial slave surface; 
!
            ipos=index(slavset,' ')-1
! 
!           default for node-to-surface contact is
!           a nodal slave surface
!
            if(slavset(ipos:ipos).eq.'S') then
!
!              'S' means node-to-face contact, it does not mean
!              that the slave surface is node-based;
!              start with assuming a nodal surface (default),
!              if not present, check for a facial surface
!
               nodeslavsurf=.true.
            else
               nodeslavsurf=.false.
            endif
!
!           determining the slave surface 
!
c            do j=1,nset
c               if(set(j).eq.slavset) exit
c            enddo
            call cident81(set,slavset,nset,id)
            j=nset+1
            if(id.gt.0) then
              if(slavset.eq.set(id)) then
                j=id
              endif
            endif
            if(j.gt.nset) then
               do j=1,nset
                  if((set(j)(1:ipos-1).eq.slavset(1:ipos-1)).and.
     &                 (set(j)(ipos:ipos).eq.'T')) then
                     nodeslavsurf=.false.
                     exit
                  endif
               enddo
            endif
!
            islav=j
!
            if(nodeslavsurf) then
!
!              nodal slave surface
!
               do j=istartset(islav),iendset(islav)
                  if(ialset(j).gt.0) then
                     node=ialset(j)
                     call addimd(imdnode,nmdnode,node)
                     if(ithermal(1).ne.2) then
                        do k=1,3
                           call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                       nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                       nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &                       ikboun,nboun,ilboun)
                        enddo
                     endif
                  else
                     k=ialset(j-2)
                     do
                        k=k-ialset(j)
                        if(k.ge.ialset(j-1)) exit
                        node=k
                        call addimd(imdnode,nmdnode,node)
                        if(ithermal(1).ne.2) then
                           do k=1,3
                              call addimdnodedof(node,k,ikmpc,ilmpc,
     &                          ipompc,nodempc,nmpc,imdnode,nmdnode,
     &                          imddof,nmddof,nactdof,mi,imdmpc,nmdmpc,
     &                          imdboun,nmdboun,ikboun,nboun,ilboun)
                           enddo
                        endif
                     enddo
                  endif
               enddo
            else
!     
!             facial slave surface
!
               do j=istartset(islav),iendset(islav)
!     
                  nelem=int(ialset(j)/10.d0)
                  jface=ialset(j)-10*nelem
!     
                  indexe=ipkon(nelem)
!     
                  if(lakon(nelem)(4:4).eq.'2') then
                     nnodelem=8
                     nface=6
                  elseif(lakon(nelem)(4:4).eq.'8') then
                     nnodelem=4
                     nface=6
                  elseif(lakon(nelem)(4:5).eq.'10') then
                     nnodelem=6
                     nface=4
                  elseif(lakon(nelem)(4:4).eq.'4') then
                     nnodelem=3
                     nface=4
                  elseif(lakon(nelem)(4:5).eq.'15') then
                     if(jface.le.2) then
                        nnodelem=6
                     else
                        nnodelem=8
                     endif
                     nface=5
                     nope=15
                  elseif(lakon(nelem)(4:4).eq.'6') then
                     if(jface.le.2) then
                        nnodelem=3
                     else
                        nnodelem=4
                     endif
                     nface=5
                     nope=6
                  else
                     cycle
                  endif
!     
!     determining the slave nodes 
!     
                  if(nface.eq.4) then
                     do k=1,nnodelem
                        nodef(k)=kon(indexe+ifacet(k,jface))
                     enddo
                  elseif(nface.eq.5) then
                     if(nope.eq.6) then
                        do k=1,nnodelem
                           nodef(k)=kon(indexe+ifacew1(k,jface))
                        enddo
                     elseif(nope.eq.15) then
                        do k=1,nnodelem
                           nodef(k)=kon(indexe+ifacew2(k,jface))
                        enddo
                     endif
                  elseif(nface.eq.6) then
                     do k=1,nnodelem
                        nodef(k)=kon(indexe+ifaceq(k,jface))
                     enddo
                  endif
!     
                  do l=1,nnodelem
                     node=nodef(l)
                     call addimd(imdnode,nmdnode,node)
                     if(ithermal(1).ne.2) then
                        do k=1,3
                           call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                       nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                       nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &                       ikboun,nboun,ilboun)
                        enddo
                     endif
                  enddo
               enddo
            endif
!     
         endif
      enddo
!
!     adding nodes, dofs, spcs and mpcs belonging to nonlinear MPC's 
!     (why only dependent nodes?)
!      
      do i=1,nmpc
         if((labmpc(i)(1:20).ne.'                    ').and.
     &          (labmpc(i)(1:6).ne.'CYCLIC').and.
     &          (labmpc(i)(1:9).ne.'SUBCYCLIC')) then
            indexe1=ipompc(i)
            if(indexe1.eq.0) cycle
            node=nodempc(1,indexe1)
            call addimd(imdnode,nmdnode,node)
            if(ithermal(1).ne.2) then
               do k=1,3
                  call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                 nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                 nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &                 ikboun,nboun,ilboun)
               enddo
            else
               k=0
               call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &              nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &              nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,ikboun,
     &              nboun,ilboun)
            endif
         endif
      enddo
!
      return
      end




