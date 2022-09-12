!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine tiefaccont(lakon,ipkon,kon,ntie,tieset,nset,set,
     &     istartset,iendset,ialset,itiefac,islavsurf,islavnode,
     &     imastnode,nslavnode,nmastnode,nslavs,nmasts,ifacecount,
     &     iponoels,inoels,ifreenoels,mortar,ipoface,nodface,nk,xnoels)
!     
!     Catalogueing the slave faces (itieface, islavsurf)
!     the slave nodes (islavnode, nslavnode)
!     the slave faces to which the slave nodes
!     belong
!     the master nodes (imastnode, nmastnode; only
!     for surface-to-surface contact)
!     
!     Authors: Li,Yang; Rakotonanahary, Samoela; 
!     
      implicit none
!     
      logical nodeslavsurf
!     
      character*8 lakon(*)
      character*81 tieset(3,*),slavset,mastset,set(*)
!     
      logical exist
!     
      integer ntie,i,j,k,l,nset,istartset(*),iendset(*),ialset(*),
     &     ifaces,nelems,jfaces,ifacem,nelemm,nslavs,nmasts,jface,
     &     jfacem,indexe,nopes,nopem,ipkon(*),kon(*),id,nodef(9),
     &     ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),node,
     &     itiefac(2,*),islavsurf(2,*),islavnode(*),imastnode(*),
     &     nslavnode(ntie+1),nmastnode(ntie+1),ifacecount,islav,imast,
     &     ipos,index1,iponoels(*),inoels(2,*),ifreenoels,ifreenoelold,
     &     mortar,numbern,numberf,iface,kflag,nk,ipoface(*),
     &     nodface(5,*),nface,nelem,nope
!     
      real*8 xnoels(*)
!     
!     nslavnode: num of slave nodes
!     islavnode: all slave nodes, tie by tie, ordered within one tie constraint
!     nmastnode: num of master nodes
!     imastnode: all master nodes, tie by tie, ordered within one tie constraint
!     islavsurf: all slave faces
!     itiefac: pointer into field islavsurf
!     
!     nodes per face for hex elements
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
!     
!     nodes per face for tet elements
!     
      data ifacet /1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/
!     
!     nodes per face for linear wedge elements
!     
      data ifacew1 /1,3,2,0,
     &     4,5,6,0,
     &     1,2,5,4,
     &     2,3,6,5,
     &     3,1,4,6/
!     
!     nodes per face for quadratic wedge elements
!     
      data ifacew2 /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     3,1,4,6,9,13,12,15/
!     
      ifacecount=0
      nslavs=0
      nmasts=0
      ifreenoels=0
      nodeslavsurf=.false.
!     
!     counters for new fields islavsurf and itiefac
!     
      do i=1,ntie
!     
!     check for contact conditions
!     
        if((tieset(1,i)(81:81).eq.'C').or.
     &       (tieset(1,i)(81:81).eq.'-')) then
          slavset=tieset(2,i)
!     
!     check whether facial slave surface; 
!     
          ipos=index(slavset,' ')-1
!     
!     default for node-to-surface contact is
!     a nodal slave surface
!     
          if(slavset(ipos:ipos).eq.'S') then
            nodeslavsurf=.true.
          endif
!     
!     determining the slave surface 
!     
c          do j=1,nset
c            if(set(j).eq.slavset) exit
c          enddo
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
     &             (set(j)(ipos:ipos).eq.'T')) then
                nodeslavsurf=.false.
                exit
              endif
            enddo
          endif
!     
          islav=j
!     
          if((mortar.le.0).and.(nodeslavsurf)) then
c          if((mortar.eq.0).and.(nodeslavsurf)) then
!     
!     nodal slave surface and node-to-surface contact
!     
!     storing the slave nodes in islavnode (sorted)
!     
            nslavnode(i)=nslavs
            numbern=0
            do j=istartset(islav),iendset(islav)
              if(ialset(j).gt.0) then
                k=ialset(j)
                call nident(islavnode(nslavs+1),k,numbern,id)
                if(id.gt.0) then
                  if(islavnode(nslavs+id).eq.k) cycle
                endif
                numbern=numbern+1
                do l=numbern,id+2,-1
                  islavnode(nslavs+l)=islavnode(nslavs+l-1)
                enddo
                islavnode(nslavs+id+1)=k
              else
                k=ialset(j-2)
                do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  call nident(islavnode(nslavs+1),k,numbern,id)
                  if(id.gt.0) then
                    if(islavnode(nslavs+id).eq.k) cycle
                  endif
                  numbern=numbern+1
                  do l=numbern,id+2,-1
                    islavnode(nslavs+l)=islavnode(nslavs+l-1)
                  enddo
                  islavnode(nslavs+id+1)=k
                enddo
              endif
              nslavnode(i+1)=nslavnode(i)+numbern
            enddo
!     
!     check all external solid faces whether they contain
!     slave nodes
!     
!     islavsurf(1,*) contains the faces (ordered)
!     islavsurf(2,*) contains the position of the
!     original order
!     
            itiefac(1,i)=ifacecount+1
            numberf=0
!     
            do j=1,nk
              index1=ipoface(j)
              do
                if(index1.eq.0) exit
                iface=nodface(4,index1)
                do k=0,3
                  if(k.eq.0) then
                    node=j
                  else
                    node=nodface(k,index1)
                  endif
!     
!     check whether node belongs to slave surface
!     
                  call nident(islavnode(nslavs+1),node,numbern,
     &                 id)
                  if(id.gt.0) then
                    if(islavnode(nslavs+id).eq.node) then
                      call nident2(islavsurf(1,ifacecount+1),
     &                     iface,numberf,id)
!     
                      ipos=0
                      if(id.gt.0) then
                        if(islavsurf(1,ifacecount+id).eq.iface) 
     &                       then
!     original position (not ordered)
                          ipos=islavsurf(2,ifacecount+id)
                        endif
                      endif
!     
!     check whether new face
!     
                      if(ipos.eq.0) then
                        numberf=numberf+1
                        do l=ifacecount+numberf,ifacecount+id+2
     &                       ,-1
                          islavsurf(1,l)=islavsurf(1,l-1)
                          islavsurf(2,l)=islavsurf(2,l-1)
                        enddo
                        islavsurf(1,ifacecount+id+1)=iface
                        islavsurf(2,ifacecount+id+1)=numberf
!     original position (not ordered)
                        ipos=numberf
                      endif
!     
!     update info to which faces a slave node
!     belongs
!     
                      ifreenoelold=iponoels(node)
                      ifreenoels=ifreenoels+1
                      iponoels(node)=ifreenoels
                      inoels(1,ifreenoels)=ifacecount+ipos
!     
!     determining the area coefficient (= the
!     force coefficient corresponding to a 
!     uniform pressure on the face)
!     
                      nelem=int(iface/10.d0)
                      jface=iface-10*nelem
!     
                      indexe=ipkon(nelem)
                      if(lakon(nelem)(4:4).eq.'2') then
                        nopes=8
                        nface=6
                      elseif(lakon(nelem)(4:4).eq.'8') then
                        nopes=4
                        nface=6
                      elseif(lakon(nelem)(4:5).eq.'10') then
                        nopes=6
                        nface=4
                      elseif(lakon(nelem)(4:4).eq.'4') then
                        nopes=3
                        nface=4
                      elseif(lakon(nelem)(4:5).eq.'15') then
                        if(jface.le.2) then
                          nopes=6
                        else
                          nopes=8
                        endif
                        nface=5
                        nope=15
                      elseif(lakon(nelem)(4:4).eq.'6') then
                        if(jface.le.2) then
                          nopes=3
                        else
                          nopes=4
                        endif
                        nface=5
                        nope=6
                      else
                        cycle
                      endif
!     
!     determining the nodes of the face
!     
                      if(nface.eq.4) then
                        do l=1,nopes
                          nodef(l)=kon(indexe+ifacet(l,jface))
                          if(nodef(l).eq.node) then
                            if(nopes.eq.3) then
                              xnoels(ifreenoels)=1.d0/3.d0
!     
!     for 6-node faces the weights 0 and 1/3 are changed into
!     1/999 and 332/999 in order to avoid problems with the
!     zero-area-check in springforc_n2f*f and springstiff_n2f*f
!     
                            elseif(l.le.3) then
                              xnoels(ifreenoels)=1.d0/999.d0
                            else
                              xnoels(ifreenoels)
     &                             =332.d0/999.d0
                            endif
                          endif
                        enddo
                      elseif(nface.eq.5) then
                        if(nope.eq.6) then
                          do l=1,nopes
                            nodef(l)=
     &                           kon(indexe+ifacew1(l,jface))
                            if(nodef(l).eq.node) then
                              if(nopes.eq.3) then
                                xnoels(ifreenoels)=1.d0/3.d0
                              else
                                xnoels(ifreenoels)=1.d0/4.d0
                              endif
                            endif
                          enddo
                        elseif(nope.eq.15) then
                          do l=1,nopes
                            nodef(l)=
     &                           kon(indexe+ifacew2(l,jface))
                            if(nodef(l).eq.node) then
                              if(nopes.eq.6) then
                                if(l.le.3) then
                                  xnoels(ifreenoels)=
     &                                 1.d0/999.d0
                                else
                                  xnoels(ifreenoels)=
     &                                 332.d0/999.d0
                                endif
                              else
!     
!     for a 8-node face a distribution of 1/100 at the vertex
!     nodes and 24/100 at the midnodes (instead of -1/12 and 1/3)
!     
                                if(l.le.4) then
                                  xnoels(ifreenoels)=
     &                                 1.d0/100.d0
                                else
                                  xnoels(ifreenoels)=
     &                                 24.d0/100.d0
                                endif
                              endif
                            endif
                          enddo
                        endif
                      elseif(nface.eq.6) then
                        do l=1,nopes
                          nodef(l)=kon(indexe+ifaceq(l,jface))
                          if(nodef(l).eq.node) then
                            if(nopes.eq.4) then
                              xnoels(ifreenoels)=1.d0/4.d0
!     
!     for a 8-node face a distribution of 1/100 at the vertex
!     nodes and 24/100 at the midnodes (instead of -1/12 and 1/3)
!     
                            elseif(l.le.4) then
                              xnoels(ifreenoels)=1.d0/100.d0
                            else
                              xnoels(ifreenoels)=
     &                             24.d0/100.d0
                            endif
                            
                          endif
                        enddo
                      endif
!     
                      inoels(2,ifreenoels)=ifreenoelold
                    endif
                  endif
                enddo
                index1=nodface(5,index1)
              enddo
            enddo
!     
!     restoring the right order of islavsurf(1,*);
!     thereafter islavsurf(2,*) is obsolete for node-to-surface
!     contact
!     
            kflag=2
            call isorti(islavsurf(1,ifacecount+1),numberf,kflag)
!     
!     update ifacecount and itiefac
!     
            itiefac(2,i)=ifacecount+numberf
!     
            nslavs=nslavnode(i+1)
            ifacecount=itiefac(2,i)
!     
          else
!     
!     element face slave surface (node-to-surface or
!     surface-to-surface contact)
!     
            nslavnode(i)=nslavs
!     
            itiefac(1,i)=ifacecount+1
            do j=istartset(islav),iendset(islav)
              if(ialset(j).gt.0) then
!     
!     store the slave face in islavsurf               
!     
                ifacecount=ifacecount+1
                islavsurf(1,ifacecount)=ialset(j)
!     
!     store the nodes belonging to the slave face
!     in islavnode
!     
                ifaces=ialset(j)
                nelems=int(ifaces/10)
                jfaces=ifaces - nelems*10
                indexe=ipkon(nelems)
!     
                if(lakon(nelems)(4:5).eq.'20') then
                  nopes=8
                elseif(lakon(nelems)(4:4).eq.'8') then
                  nopes=4
                elseif(lakon(nelems)(4:5).eq.'10') then
                  nopes=6
                elseif(lakon(nelems)(4:4).eq.'4') then
                  nopes=3
                endif
!     
                if(lakon(nelems)(4:4).eq.'6') then
                  if(jfaces.le.2) then
                    nopes=3
                  else
                    nopes=4
                  endif
                endif
                if(lakon(nelems)(4:5).eq.'15') then
                  if(jfaces.le.2) then
                    nopes=6
                  else
                    nopes=8
                  endif
                endif   
!     
                do l=1,nopes
                  if((lakon(nelems)(4:4).eq.'2').or.
     &                 (lakon(nelems)(4:4).eq.'8')) then
                    node=kon(indexe+ifaceq(l,jfaces))
                  elseif((lakon(nelems)(4:4).eq.'4').or.
     &                   (lakon(nelems)(4:5).eq.'10')) then
                    node=kon(indexe+ifacet(l,jfaces))
                  elseif(lakon(nelems)(4:4).eq.'6') then
                    node=kon(indexe+ifacew1(l,jfaces))
                  elseif(lakon(nelems)(4:5).eq.'15') then
                    node=kon(indexe+ifacew2(l,jfaces))
                  endif
                  call nident(islavnode(nslavnode(i)+1),node,
     &                 nslavs-nslavnode(i),id)
                  exist=.FALSE.
                  if(id.gt.0) then
                    if(islavnode(nslavnode(i)+id).eq.node) then
                      exist=.TRUE.
                    endif
                  endif
                  if(.not.exist) then
                    nslavs=nslavs+1
                    do k=nslavs,nslavnode(i)+id+2,-1
                      islavnode(k)=islavnode(k-1)
                    enddo
                    islavnode(nslavnode(i)+id+1)=node
                  endif
!     
!     filling fields iponoels and inoels
!     
                  if(mortar.le.0) then
c                  if(mortar.eq.0) then
                    ifreenoels=ifreenoels+1
                    inoels(1,ifreenoels)=ifacecount
                    inoels(2,ifreenoels)=iponoels(node)
                    iponoels(node)=ifreenoels
!     
!     filling xnoels with the coefficient corresponding
!     to a constant pressure (sum over all nodes must
!     be 1)
!     
                    if(nopes.eq.3) then
                      xnoels(ifreenoels)=1.d0/3.d0
                    elseif(nopes.eq.4) then
                      xnoels(ifreenoels)=1.d0/4.d0
                    elseif(nopes.eq.6) then
                      if(l.le.3) then
                        xnoels(ifreenoels)=1.d0/999.d0
                      else
                        xnoels(ifreenoels)=332.d0/999.d0
                      endif
                    elseif(nopes.eq.8) then
!     
!     for a 8-node face a distribution of 1/100 at the vertex
!     nodes and 24/100 at the midnodes (instead of -1/12 and 1/3)
!     
                      if(l.le.4) then
                        xnoels(ifreenoels)=1.d0/100.d0
                      else
                        xnoels(ifreenoels)=24.d0/100.d0
                      endif
                    endif
c     inoels(2,ifreenoels)=ifreenoelold
                  endif
                enddo
!     
              endif
            enddo
            nslavnode(i+1)=nslavs
            itiefac(2,i)=ifacecount
          endif
!
!     for massless contact: marking the elements adjacent to
!     the slave faces with the slave face number
!
c   TODO CMT: label SLAVE faces for contact IF mass redistribution
c          if(mortar.eq.-1) then
c            do j=itiefac(1,i),itiefac(2,i)
c              iface=islavsurf(1,j)
c              nelem=int(iface/10.d0)
c              jface=iface-10*nelem
c              write(lakon(nelem)(7:7),'(i1)') jface
c            enddo
c          endif
!     
!     determining the master surface
!     
          mastset=tieset(3,i)
c          do j=1,nset
c            if(set(j).eq.mastset) exit
c          enddo
          call cident81(set,mastset,nset,id)
          j=nset+1
          if(id.gt.0) then
            if(mastset.eq.set(id)) then
              j=id
            endif
          endif
          if(j.gt.nset) then
            write(*,*) '*ERROR in tiefaccont: master surface'
            write(*,*) '       does not exist'
            call exit(201)
          endif
          imast=j
          nmastnode(i)=nmasts
!     
          do j=istartset(imast),iendset(imast)
!     
!     create imastnode, and nmastnode
!     
            ifacem=ialset(j)
            nelemm=int(ifacem/10)
            jfacem=ifacem - nelemm*10
            indexe=ipkon(nelemm)
!     
            if(lakon(nelemm)(4:5).eq.'20') then
              nopem=8
            elseif(lakon(nelemm)(4:4).eq.'8') then
              nopem=4
            elseif(lakon(nelemm)(4:5).eq.'10') then
              nopem=6
            elseif(lakon(nelemm)(4:4).eq.'4') then
              nopem=3
            endif
!     
            if(lakon(nelemm)(4:4).eq.'6') then
              if(jfacem.le.2) then
                nopem=3
              else
                nopem=4
              endif
            endif
            if(lakon(nelemm)(4:5).eq.'15') then
              if(jfacem.le.2) then
                nopem=6
              else
                nopem=8
              endif
            endif   
!     
            do l=1,nopem
              if((lakon(nelemm)(4:4).eq.'2').or.
     &             (lakon(nelemm)(4:4).eq.'8')) then
                node=kon(indexe+ifaceq(l,jfacem))
              elseif((lakon(nelemm)(4:4).eq.'4').or.
     &               (lakon(nelemm)(4:5).eq.'10')) then
                node=kon(indexe+ifacet(l,jfacem))
              elseif(lakon(nelemm)(4:4).eq.'6') then
                node=kon(indexe+ifacew1(l,jfacem))
              elseif(lakon(nelemm)(4:5).eq.'15') then
                node=kon(indexe+ifacew2(l,jfacem))
              endif
              call nident(imastnode(nmastnode(i)+1),node,
     &             nmasts-nmastnode(i),id)
              exist=.FALSE.
              if(id.gt.0) then
                if(imastnode(nmastnode(i)+id).eq.node) then
                  exist=.TRUE.
                endif
              endif
              if(exist) cycle
              nmasts=nmasts+1
              do k=nmasts,nmastnode(i)+id+2,-1
                imastnode(k)=imastnode(k-1)
              enddo
              imastnode(nmastnode(i)+id+1)=node
            enddo
!     
          enddo
          nmastnode(i+1)=nmasts
!     
        else
!     
!     no contact tie
!     
          nslavnode(i+1)=nslavnode(i)
          nmastnode(i+1)=nmastnode(i)
        endif 
      enddo      
!     
      return
      end
