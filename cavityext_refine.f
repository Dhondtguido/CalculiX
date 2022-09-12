!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998 Guido Dhondt
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
      subroutine cavityext_refine(kontet,ifatet,ifreetet,bc,ifac,itetfa,
     &     ifreefa,planfa,ipofa,cotet,ibase,node,iexternfa,ipoed,
     &     iedg,ifreeed,ipoeled,ieled,ifreele,nktet,netet_,
     &     ibasenewnodes,conewnodes,ipoeln,ieln,ifreeln,
     &     nnewnodes,iexternedg,iedtet,cotetorig,iedge,iexternnode,
     &     cg,height,iparentel)
!     
      implicit none
!     
!     ig(*): face queue
!     ige(*): cavity element to which face ig(*) belongs
!     ng: number of face queue elements
!     
!     iecav(i): elements belonging to the cavity
!     necav
!     
!     ikcav(i): nodes belonging to the cavity (including the border)
!     nkcav
!     
!     ifcav(i): border face of the cavity
!     nfcav
!     
      integer nktet,netet_,nnewnodes,ibasis,incav(4,netet_),isol,
     &     kontet(4,*),ifatet(4,*),ifreetet,ifac(4,*),itetfa(2,*),
     &     ifreefa,ipofa(*),nodes(4),iexternfa(*),ipoed(*),ifreeed,
     &     ikboun(nktet),iedg(3,*),irevert,ipoeled(*),ieled(2,*),
     &     ibase,node,i,ig(4*netet_),ige(netet_),iecav(netet_),ifs,
     &     ikcav(nktet),ipoeln(*),ieln(2,*),ifreeln,ibasenewnodes(*),
     &     ifcav(4*netet_),kflag,ng,necav,nkcav,nfcav,iface,iel,k,
     &     ichange,id,idum(1),ielnew,j,ifanew,ielement,l,inode,ifreele,
     &     inewel(netet_),
     &     iexternedg(*),iedtet(6,*),iedge,iexternnode(*),node1,node2,
     &     n1,n2,n3,indexe,indexeold,index,indexold,jface,imasteredg,
     &     iparentel(*)
!     
      real*8 bc(4,*),planfa(4,*),cotet(3,*),cotetorig(3),cg(3,*),
     &     dd,x,y,z,conewnodes(3,*),volume,height,hh,hei(4),surf,
     &     cotetcpy(3)
!     
!     
!     
!     start: cavity consists of base element only
!     
      ng=4
      do i=1,4
        ig(i)=abs(ifatet(i,ibase))
        ige(i)=ibase
      enddo
!     
      kflag=2
      call isortii(ig,ige,ng,kflag)
!     
      necav=1
      iecav(necav)=ibase
!     
      nkcav=4
      do i=1,4
        ikcav(i)=kontet(i,ibase)
      enddo
      kflag=1
      call insertsorti(ikcav,nkcav)
!     
      nfcav=0
!     
!     add elements through adjacency
!     
      loop3: do
        iface=ig(1)
        iel=ige(1)
        if(iface.eq.0) exit
        do i=1,2
          ielnew=itetfa(i,iface)
!     
          if(ielnew.eq.0) then
!     
!     outside face: new cavity face
!     remove iface from face queue
!     
            nfcav=nfcav+1
            ifcav(nfcav)=iface
            do k=1,ng-1
              ig(k)=ig(k+1)
              ige(k)=ige(k+1)
            enddo
            ig(ng)=0
            ige(ng)=0
            ng=ng-1
            cycle loop3
          endif
!     
          if(ielnew.ne.iel) exit
        enddo
        iel=ielnew
        dd=dsqrt((cotet(1,node)-bc(1,iel))**2+
     &       (cotet(2,node)-bc(2,iel))**2+
     &       (cotet(3,node)-bc(3,iel))**2)
!     
        if(dd.ge.bc(4,iel)) then
!     
!     new cavity face detected
!     
          nfcav=nfcav+1
c     write(*,*) 'cavity add2',ibase,nfcav
          ifcav(nfcav)=iface
!     
!     remove face iface from face queue
!     
          do k=1,ng-1
            ig(k)=ig(k+1)
            ige(k)=ige(k+1)
          enddo
          ig(ng)=0
          ige(ng)=0
          ng=ng-1
        else
!     
!     new cavity element detected
!     
          necav=necav+1
          iecav(necav)=iel
!     
!     updating the cavity nodes
!     
          do j=1,4
            inode=kontet(j,iel)
            call nident(ikcav,inode,nkcav,id)
            if(id.ne.0) then
              if(ikcav(id).eq.inode) cycle
            endif
!     
!     new cavity node detected
!     
            do k=nkcav,id+1,-1
              ikcav(k+1)=ikcav(k)
            enddo
            ikcav(id+1)=inode
            nkcav=nkcav+1
          enddo
!     
!     remove face iface from face queue
!     
          do k=1,ng-1
            ig(k)=ig(k+1)
            ige(k)=ige(k+1)
          enddo
          ig(ng)=0
          ige(ng)=0
          ng=ng-1
!     
!     updating the cavity faces
!     
          do j=1,4
            ifanew=abs(ifatet(j,iel))
            if(ifanew.eq.iface) cycle
            call nident(ig,ifanew,ng,id)
            if(id.gt.0) then
              if(ig(id).eq.ifanew) then
!     
!     remove face ifanew from face queue
!     
                do k=id,ng-1
                  ig(k)=ig(k+1)
                  ige(k)=ige(k+1)
                enddo
                ig(ng)=0
                ige(ng)=0
                ng=ng-1
                cycle
              endif
            endif
!     
!     add face ifanew to face queue
!     
            do k=ng,id+1,-1
              ig(k+1)=ig(k)
              ige(k+1)=ige(k)
            enddo
            ig(id+1)=ifanew
            ige(id+1)=iel
            ng=ng+1
          enddo
        endif
!     
      enddo loop3
!     
!     check of the cavity
!     
      kflag=1
      call isortii(iecav,idum,necav,kflag)
      call isortii(ifcav,idum,nfcav,kflag)
!     
!     check for nodes within the cavity
!     
      do
        if(necav.eq.1) exit
        do j=1,nkcav
          ikboun(j)=1
        enddo
        do j=1,nfcav
          iface=ifcav(j)
          do k=1,3
            inode=ifac(k,iface)
            call nident(ikcav,inode,nkcav,id)
            ikboun(id)=0
          enddo
        enddo
        do j=1,nkcav
          if(ikboun(j).eq.1) exit
        enddo
        if(j.gt.nkcav) exit
!     
!     node inside cavity
!     
        write(*,*) '*WARNING: NODE WITHIN CAVITY'
        inode=ikcav(j)
        loop1: do j=1,necav
          ielement=iecav(j)
          if(ielement.eq.ibase) cycle
          do k=1,4
            if(kontet(k,ielement).eq.inode) exit loop1
          enddo
        enddo loop1
!     
!     remove ielement
!     
        call nident(iecav,ielement,necav,id)
        do k=id,necav-1
          iecav(k)=iecav(k+1)
        enddo
        iecav(necav)=0
        necav=necav-1
!     
!     remove/add faces of element
!     
        do l=1,4
          iface=abs(ifatet(l,ielement))
          call nident(ifcav,iface,nfcav,id)
          if(id.ne.0) then
            if(ifcav(id).eq.iface) then
!     
!     remove face from ifcav
!     
              do k=id,nfcav-1
                ifcav(k)=ifcav(k+1)
              enddo
              ifcav(nfcav)=0
              nfcav=nfcav-1
              cycle
            endif
          endif
!     
!     add face to ifcav
!     
          do k=nfcav,id+1,-1
            ifcav(k+1)=ifcav(k)
          enddo
          ifcav(id+1)=iface
          nfcav=nfcav+1
        enddo
      enddo
!     
!     check for too smal elements (too small height)
!     
!     restore the surface coordinates of the node
!     
      do i=1,3
        cotetcpy(i)=cotet(i,node)
        cotet(i,node)=cotetorig(i)
      enddo
!     
      do
        if(necav.eq.1) exit
        ichange=0
        do j=1,nfcav
          iface=ifcav(j)
!     
          do k=1,3
            nodes(k)=ifac(k,iface)
          enddo
          nodes(4)=node
!     
!     calculate the volume of the tetrahedron obtained by
!     connecting node with the face "iface"
!     
          call calcvol(nodes(1),nodes(2),nodes(3),nodes(4),cotet,
     &         volume)
          volume=dabs(volume)
!     
!     calculate the area of the faces of the tetrahedron
!     
          call calcsurf(nodes(1),nodes(2),nodes(3),cotet,hei(1))
          call calcsurf(nodes(2),nodes(3),nodes(4),cotet,hei(2))
          call calcsurf(nodes(3),nodes(4),nodes(1),cotet,hei(3))
          call calcsurf(nodes(4),nodes(1),nodes(2),cotet,hei(4))
!     
!     calculate the height
!     
          do k=1,4
            hei(k)=volume/hei(k)
          enddo
!     
!     look for the element inside the cavity bordering the
!     face
!     
          do k=1,2
            ielement=itetfa(k,iface)
            call nident(iecav,ielement,necav,id)
            if(id.ne.0) then
              if(iecav(id).eq.ielement) exit
            endif
          enddo
!     
!     the element is removed if
!     - the volume is not extremely small AND
!     - at least one height is very small AND
!     - it connects to an external face
!     
          if((dabs(volume).gt.1.d-12).and.
     &         ((hei(1).lt.height/100.d0).or.
     &         (hei(2).lt.height/100.d0).or.
     &         (hei(3).lt.height/100.d0).or.
     &         (hei(4).lt.height/100.d0)).and.
     &         (iexternfa(iface).gt.0).and.
     &         (ielement.ne.ibase)) then
c            write(*,*) '*INFO in cavityext_refine: element is removed'
c            write(*,*) '      from cavity'
          else
            cycle
          endif
!     
          ichange=1
!     
!     remove ielement
!     
          call nident(iecav,ielement,necav,id)
          do k=id,necav-1
            iecav(k)=iecav(k+1)
          enddo
          iecav(necav)=0
          necav=necav-1
!     
!     remove/add faces of element
!     
          do l=1,4
            iface=abs(ifatet(l,ielement))
            call nident(ifcav,iface,nfcav,id)
            if(id.ne.0) then
              if(ifcav(id).eq.iface) then
!     
!     remove face from ifcav
!     
                do k=id,nfcav-1
                  ifcav(k)=ifcav(k+1)
                enddo
                ifcav(nfcav)=0
                nfcav=nfcav-1
                cycle
              endif
            endif
!     
!     add face to ifcav
!     
            do k=nfcav,id+1,-1
              ifcav(k+1)=ifcav(k)
            enddo
            ifcav(id+1)=iface
            nfcav=nfcav+1
          enddo
!     
!     remove nodes from the cavity
!     
          loop4: do l=1,4
            inode=kontet(l,ielement)
            do k=1,necav
              iel=iecav(k)
              if(iel.eq.ielement) cycle
              do i=1,4
                if(inode.eq.kontet(i,iel)) cycle loop4
              enddo
            enddo
            call nident(ikcav,inode,nkcav,id)
            do k=id,nkcav-1
              ikcav(k)=ikcav(k+1)
            enddo
            ikcav(nkcav)=0
            nkcav=nkcav-1
          enddo loop4
          if(ichange.ne.0) exit
        enddo
        if(ichange.eq.0) exit
      enddo
!     
!     revert to the subsurface coordinates of the node at stake
!     
      do i=1,3
        cotet(i,node)=cotetcpy(i)
      enddo
!     
!     check for nonconvexity of the cavity
!     
      do
        if(necav.eq.1) exit
        ichange=0
        do j=1,nfcav
          iface=ifcav(j)
!     
!     look for the element inside the cavity bordering the
!     face
!     
          do k=1,2
            ielement=itetfa(k,iface)
            call nident(iecav,ielement,necav,id)
            if(id.ne.0) then
              if(iecav(id).eq.ielement) exit
            endif
          enddo
!     
!     the base element should not be removed
!     
          if(ielement.eq.ibase) cycle
!     
          do k=1,4
            if(abs(ifatet(k,ielement)).eq.iface) exit
          enddo
          if(ifatet(k,ielement)*(planfa(1,iface)*cotet(1,node)+
     &         planfa(2,iface)*cotet(2,node)+
     &         planfa(3,iface)*cotet(3,node)+
     &         planfa(4,iface)).ge.(1.d-12*iface)) cycle
          ichange=1
!     
!     remove ielement
!     
          call nident(iecav,ielement,necav,id)
          do k=id,necav-1
            iecav(k)=iecav(k+1)
          enddo
          iecav(necav)=0
          necav=necav-1
!     
!     remove/add faces of element
!     
          do l=1,4
            iface=abs(ifatet(l,ielement))
            call nident(ifcav,iface,nfcav,id)
            if(id.ne.0) then
              if(ifcav(id).eq.iface) then
!     
!     remove face from ifcav
!     
                do k=id,nfcav-1
                  ifcav(k)=ifcav(k+1)
                enddo
                ifcav(nfcav)=0
                nfcav=nfcav-1
                cycle
              endif
            endif
!     
!     add face to ifcav
!     
            do k=nfcav,id+1,-1
              ifcav(k+1)=ifcav(k)
            enddo
            ifcav(id+1)=iface
            nfcav=nfcav+1
          enddo
!     
!     remove nodes from the cavity
!     
          loop2: do l=1,4
            inode=kontet(l,ielement)
            do k=1,necav
              iel=iecav(k)
              if(iel.eq.ielement) cycle
              do i=1,4
                if(inode.eq.kontet(i,iel)) cycle loop2
              enddo
            enddo
            call nident(ikcav,inode,nkcav,id)
            do k=id,nkcav-1
              ikcav(k)=ikcav(k+1)
            enddo
            ikcav(nkcav)=0
            nkcav=nkcav-1
          enddo loop2
          if(ichange.ne.0) exit
        enddo
        if(ichange.eq.0) exit
      enddo
!     
!     restoring the unperturbed coordinates
!     
      do i=1,3
        cotet(i,node)=cotetorig(i)
      enddo
!     
!     determining the nodes belonging to the edge
!     (only if the edge is part of an originally external
!     edge, i.e. if iexternedg(iedge)=1)
!     
!     labeling the new node as external
!     
      if(iexternedg(iedge).gt.0) then
        node1=iedg(1,iedge)
        node2=iedg(2,iedge)
      else
        node1=0
        node2=0
      endif
      iexternnode(node)=1
      imasteredg=iexternedg(iedge)
c      write(*,*) 'start cavityext_refine',node,iedge,node1,node2,
c     &     iexternedg(iedge)
c      write(*,*) 'start cavityext_refine',node,iedg(1,12),iedg(2,12),
c     &     iexternedg(12)
!     
!     remove the elements and faces physically
!     
      do i=1,necav
!     
!     storing the nodes of the removed elements:
!     are needed in case the newly generated elements
!     have a really bad shape and the original mesh
!     has to be restored
!     
        incav(1,i)=kontet(1,iecav(i))
        incav(2,i)=kontet(2,iecav(i))
        incav(3,i)=kontet(3,iecav(i))
        incav(4,i)=kontet(4,iecav(i))
!     
!     remove element iecav(i)
!     
        call removetet_refine(kontet,ifatet,ifreetet,ifac,itetfa,
     &       ifreefa,ipofa,iecav(i),ipoeln,ieln,ifreeln,
     &       ipoeled,ieled,ifreele,iedtet,ipoed,iedg,
     &       ifreeed,iexternfa,iexternedg)
      enddo
!     
!     fill the cavity
!     
      irevert=0
      do i=1,nfcav
        iface=ifcav(i)
!     
!     right orientation of the tetrahedron; can be replaced
!     by a check on the final mesh, if this is faster
!     
        if(planfa(1,iface)*cotet(1,node)+
     &       planfa(2,iface)*cotet(2,node)+
     &       planfa(3,iface)*cotet(3,node)+
     &       planfa(4,iface).ge.0) then
          do k=1,3
            nodes(k)=ifac(k,iface)
          enddo
        else
          do k=1,3
            nodes(k)=ifac(4-k,iface)
          enddo
        endif
!     
        nodes(4)=node
!     
!     check the volume of the new element; if the volume is
!     zero (within the round-off accuracy) the element is not 
!     created (such an element is assumed to be on the
!     surface) and ifcav(i) gets a negative sign
!     
!     calculate the volume of the new element and the surface
!     of its faces
!     
        call calcvol(nodes(1),nodes(2),nodes(3),nodes(4),cotet,volume)
!     
        if(dabs(volume).lt.1.d-12) then
          if(iexternfa(iface).gt.0) then
            ifcav(i)=-ifcav(i)
            cycle
          endif
        endif
!     
        inewel(i)=ifreetet
        call generatetet_refine2(kontet,ifatet,ifreetet,bc,ifac,itetfa,
     &       ifreefa,planfa,ipofa,nodes,cotet,ipoeln,ieln,
     &       ifreeln,ipoed,ifreeed,iedg,ipoeled,ieled,ifreele,iedtet,
     &       cg)
        iparentel(inewel(i))=iparentel(ibase)
!     
        do j=1,4
          iface=abs(ifatet(j,inewel(i)))
          call calcsurf(ifac(1,iface),ifac(2,iface),ifac(3,iface),
     &         cotet,surf)
          hei(j)=volume/surf
        enddo
!     
!     check the volume and the height orthogonal to each face of
!     the new element
!     
        hh=height
!     
        if((volume.lt.hh**3/6000.d0).or.
     &       (hei(1).lt.hh/100.d0).or.
     &       (hei(2).lt.hh/100.d0).or.
     &       (hei(3).lt.hh/100.d0).or.
     &       (hei(4).lt.hh/100.d0)) then
          irevert=1
          nktet=nktet-1
c          write(*,*) '*INFO in cavityext_refine: bad element'
c          write(*,*) '      node is not inserted'
!     
!     remove the new elements generated so far
!     
          do j=1,i
            if(ifcav(j).lt.0) cycle
            call removetet_refine(kontet,ifatet,ifreetet,ifac,itetfa,
     &           ifreefa,ipofa,inewel(j),ipoeln,ieln,ifreeln,
     &           ipoeled,ieled,ifreele,iedtet,ipoed,iedg,ifreeed,
     &           iexternfa,iexternedg)
          enddo
!     
!     recreate the original elements (element numbers etc..
!     may not be the same as before)
!     
          do j=1,necav
            do k=1,4
              nodes(k)=incav(k,j)
            enddo
            inewel(j)=ifreetet
            call generatetet_refine2(kontet,ifatet,ifreetet,bc,ifac,
     &           itetfa,ifreefa,planfa,ipofa,nodes,cotet,ipoeln,ieln,
     &           ifreeln,ipoed,ifreeed,iedg,ipoeled,ieled,ifreele,
     &           iedtet,cg)
            iparentel(inewel(i))=iparentel(ibase)
          enddo
          exit
        endif
      enddo
!     
!     due to the remeshing of the cavity the base element numbers in
!     ibasenewnodes(*) have to be updated      
!     
!     sorting the original cavity elements
!     
      kflag=1
      call isortii(iecav,idum,necav,kflag)
!     
!     reassigns the basis element numbers for the nodes which still
!     need to be inserted
!     
      do j=1,nnewnodes
!     
!     check for nodes already inserted
!     
        if(ibasenewnodes(j).eq.0) cycle
!     
        ibasis=ibasenewnodes(j)
        call nident(iecav,ibasis,necav,id)
        if(id.gt.0) then
          if(iecav(id).eq.ibasis) then
!     
!     base element of new node j belonged to remeshed cavity
!     
            x=conewnodes(1,j)
            y=conewnodes(2,j)
            z=conewnodes(3,j)
!     
            if(irevert.eq.0) then
!     
!     truly new mesh: nfcav new elements         
!     
              do i=1,nfcav
                if(ifcav(i).lt.0) cycle
                iel=inewel(i)
                isol=1
!     
!     check the node lies on the "inside" of all
!     4 faces bordering the element
!     
                do k=1,4
                  ifs=ifatet(k,iel)
                  iface=abs(ifs)
                  if(((planfa(1,iface)*x+planfa(2,iface)*y+
     &                 planfa(3,iface)*z+planfa(4,iface))*ifs).lt.
     &                 -1.d-12*iface) then
                    isol=0
                    exit
                  endif
                enddo
                if(isol.eq.1) then
!     
!     new base element found
!     
                  ibasenewnodes(j)=iel
                  exit
                endif
              enddo
!     
              if(isol.eq.0) then
c                write(*,*) '*ERROR in cavity: node belongs to'
c                write(*,*) '       zero elements'
c     call exit(201)
                ibasenewnodes(j)=0
              endif
            else
!     
!     no new mesh, only numbers may have changed:
!     necav elemetns
!     
              do i=1,necav
                iel=inewel(i)
                isol=1
!     
!     check the node lies on the "inside" of all
!     4 faces bordering the element
!     
                do k=1,4
                  ifs=ifatet(k,iel)
                  iface=abs(ifs)
                  if(((planfa(1,iface)*x+planfa(2,iface)*y+
     &                 planfa(3,iface)*z+planfa(4,iface))*ifs).lt.
     &                 -1.d-12*iface) then
                    isol=0
                    exit
                  endif
                enddo
                if(isol.eq.1) then
!     
!     new base element found
!     
                  ibasenewnodes(j)=iel
                  exit
                endif
              enddo
!     
              if(isol.eq.0) then
c                write(*,*) '*ERROR in cavity: node belongs to'
c                write(*,*) '       zero elements'
c     call exit(201)
                ibasenewnodes(j)=0
              endif
            endif
          endif
        endif
      enddo
!     
      if(irevert.eq.1) return
!     
!     labeling the newly generated faces and edges as external or
!     internal
!     
      do i=1,nfcav
        if(ifcav(i).lt.0) cycle
        do j=1,4
          iface=abs(ifatet(j,inewel(i)))
          if(itetfa(2,iface).eq.0) then
c     iexternfa(iface)=1
!     
!     only if the face is new iexternfa has to be determined
!     
            if(iexternfa(iface).le.0) then
!     
!     n1<n2 are two nodes belonging to the face and < node
!     
              n1=ifac(1,iface)
              n2=ifac(2,iface)
c              write(*,*) 'external face ',n1,n2,node
!     
              do k=1,nfcav
                if(ifcav(k).gt.0) cycle
                jface=abs(ifcav(k))
                if(((n1.eq.ifac(1,jface)).or.
     &               (n1.eq.ifac(2,jface))).and.
     &               ((n2.eq.ifac(2,jface)).or.
     &               (n2.eq.ifac(3,jface)))) then
                  iexternfa(iface)=iexternfa(jface)
c     write(*,*) 'cavityext heritage',jface,iface
                  exit
                endif
              enddo
!     
              if(iexternfa(iface).eq.0) then
                write(*,*) '*ERROR in cavityext_refine: external'
                write(*,*) '       face has no parent'
                call exit(201)
              endif
*** start change  23.01.2020          
c            endif
*** end change  23.01.2020          
!     
              do k=1,3
                n1=k
                n2=k+1
                if(n2.gt.3) n2=1
                n1=ifac(n1,iface)
                n2=ifac(n2,iface)
                if(n1.gt.n2) then
                  n3=n1
                  n1=n2
                  n2=n3
                endif
c                write(*,*) 'cavityext_refine edge ',n1,n2
c                write(*,*) 'node1 node node2',node1,node,node2
                indexe=ipoed(n1)
                do
                  if(indexe.eq.0) then
                    write(*,*) '*ERROR in cavityext_refine'
                    write(*,*) '       iedge data base is corrupt'
                    call exit(201)
                  endif
                  if(iedg(2,indexe).eq.n2) then
                    if(n1.eq.node) then
                      if((n2.eq.node1).or.(n2.eq.node2)) then
c     iexternedg(indexe)=iexternedg(iedge)
                        iexternedg(indexe)=imasteredg
c                        write(*,*) 'cavityext_refine',n1,n2,
c     &                       iexternedg(indexe)
                      else
                        if(iexternedg(indexe).eq.0)
     &                       iexternedg(indexe)=-1
c                        write(*,*) 'cavityext_refine',n1,n2,
c     &                       iexternedg(indexe)
                      endif
                    elseif(n2.eq.node) then
                      if((n1.eq.node1).or.(n1.eq.node2)) then
c     iexternedg(indexe)=iexternedg(iedge)
                        iexternedg(indexe)=imasteredg
c                        write(*,*) 'cavityext_refine',n1,n2,
c     &                       iexternedg(indexe)
                      else
                        if(iexternedg(indexe).eq.0)
     &                       iexternedg(indexe)=-1
c                        write(*,*) 'cavityext_refine',n1,n2,
c     &                       iexternedg(indexe)
                      endif
                    else
                      if(iexternedg(indexe).eq.0)
     &                     iexternedg(indexe)=-1
c                      write(*,*) 'cavityext_refine',n1,n2,
c     &                     iexternedg(indexe)
                    endif
                    exit
                  endif
                  indexe=iedg(3,indexe)
                enddo
              enddo
***   start change  23.01.2020          
            endif
*** end change 23.01.2020
          endif
        enddo
      enddo
!     
!     removing the orphan faces and edges
!     
      do i=1,nfcav
        if(ifcav(i).gt.0) cycle
        iface=abs(ifcav(i))
        do k=1,3
          n1=k
          n2=k+1
          if(n2.gt.3) n2=1
          n1=ifac(n1,iface)
          n2=ifac(n2,iface)
          if(n1.gt.n2) then
            n3=n1
            n1=n2
            n2=n3
          endif
          indexe=ipoed(n1)
!     
          if(iedg(2,indexe).eq.n2) then
!     
!     if the edge does not belong to any element it is
!     removed
!     
            if(ipoeled(indexe).eq.0) then
              iexternedg(indexe)=0
              ipoed(n1)=iedg(3,indexe)
              iedg(3,indexe)=ifreeed
              ifreeed=indexe
            endif
          else
            do
              indexeold=indexe
              indexe=iedg(3,indexe)
              if(indexe.eq.0) exit
              if(iedg(2,indexe).eq.n2) then
                if(ipoeled(indexe).eq.0) then
!     
!     if the edge does not belong to any element
!     it is removed
!     
                  iexternedg(indexe)=0
                  iedg(3,indexeold)=iedg(3,indexe)
                  iedg(3,indexe)=ifreeed
                  ifreeed=indexe
                endif
                exit
              endif
            enddo
          endif
        enddo
!     
!     removing the face
!     
        inode=ifac(1,iface)
        index=ipofa(inode)
        if(index.eq.iface) then
          iexternfa(iface)=0
          ipofa(inode)=ifac(4,index)
          ifac(4,index)=ifreefa
          ifreefa=index
        else
          do
            indexold=index
            index=ifac(4,index)
            if(index.eq.0) then
              write(*,*) '*ERROR in cavityext_refine: face to be'
              write(*,*) '       deleted is not catalogued'
              write(*,*) '       in field ifac'
              write(*,*) (ifac(j,iface),j=1,4)
              call exit(201)
            endif
            if(index.eq.iface) then
              iexternfa(iface)=0
              ifac(4,indexold)=ifac(4,index)
              ifac(4,index)=ifreefa
              ifreefa=index
              exit
            endif
          enddo
        endif
      enddo
!     
      return
      end
      



      
