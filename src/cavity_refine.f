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
      subroutine cavity_refine(kontet,ifatet,ifreetet,bc,ifac,itetfa,
     &     ifreefa,planfa,ipofa,cotet,ibase,node,iexternfa,ipoed,
     &     iedg,ifreeed,ipoeled,ieled,ifreele,nktet,netet_,
     &     ibasenewnodes,conewnodes,ipoeln,ieln,ifreeln,
     &     nnewnodes,iexternedg,iedtet,cg,height,iparentel)
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
     &     id,idum(1),ielnew,j,ifanew,ielement,l,inode,ifreele,
     &     iexternedg(*),iedtet(6,*),ichange,inewel(netet_),
     &     iparentel(*)
!     
      real*8 bc(4,*),planfa(4,*),cotet(3,*),cg(3,*),height,hh,
     &     dd,x,y,z,conewnodes(3,*),volume,quality
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
c     write(*,*) 'cavity add1',ibase,nfcav
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
     &     (cotet(2,node)-bc(2,iel))**2+
     &     (cotet(3,node)-bc(3,iel))**2)
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
        write(*,*) '*WARNING in cavity_refine: node within cavity'
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
c     write(*,*) '*WARNING: CAVITY NOT CONVEX'
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
c     write(*,*) 'cavity subtract2',ibase,nfcav
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
c     write(*,*) 'cavity add4',ibase,nfcav
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
c     write(*,*) i,nfcav
c     write(*,*) ifcav(i)
c     write(*,*) 'cavity1 ',i,ifcav(i)
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
c     write(*,*) 'cavity2 ',i,ifcav(i)
        nodes(4)=node
        inewel(i)=ifreetet
        call generatetet_refine2(kontet,ifatet,ifreetet,bc,ifac,itetfa,
     &       ifreefa,planfa,ipofa,nodes,cotet,ipoeln,ieln,
     &       ifreeln,ipoed,ifreeed,iedg,ipoeled,ieled,ifreele,iedtet,
     &       cg)
        iparentel(inewel(i))=iparentel(ibase)
!     
!     calculate the element quality (aspect ratio)
!     NOTE: this is not the same as meshquality.f; in meshquality.f
!     the quality is calculated for the existing element/mesh but
!     here the element does not yet exist; therefore we cannot
!     use meshquality.f
!     
        call meshqualitycavity(nodes(1),nodes(2),nodes(3),nodes(4),
     &       cotet,quality,volume)
!     
        hh=height
!     
!     check the volume and the quality of the new element
!     
        if((volume.lt.hh**3/6000.d0).or.
     &       (quality.gt.10.d0)) then
          irevert=1
          nktet=nktet-1
c          write(*,*) '*INFO in cavity_refine: bad element'
c          write(*,*) '      node is not inserted'
!     
!     remove the new elements generated so far
!     
          do j=1,i
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
c     write(*,*) 'cavity 6',i,ifcav(i)
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
c                write(*,*)
c     &               '*ERROR in cavity_refine: node belongs to'
c                write(*,*) '       zero elements'
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
c     write(*,*)
c     &                   '*ERROR in cavity_refine: node belongs to'
c     write(*,*) '       zero elements'
                ibasenewnodes(j)=0
              endif
            endif
          endif
        endif
      enddo
c     write(*,*) 'END OF CAVITY.F'
!     
      return
      end
      



      
