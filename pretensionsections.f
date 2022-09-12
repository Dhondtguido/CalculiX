!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine pretensionsections(inpc,textpart,ipompc,nodempc,
     &     coefmpc,nmpc,nmpc_,mpcfree,nk,ikmpc,ilmpc,
     &     labmpc,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,lakon,
     &     kon,ipkon,set,nset,istartset,iendset,ialset,co,ics,dcs,t0,
     &     ithermal,ne,ier)
!     
!     reading the input deck: *PRE-TENSION SECTION
!     
      implicit none
!     
      logical twod,quad,normal
!     
      character*1 inpc(*)
      character*8 lakon(*)
      character*20 labmpc(*)
      character*81 surface,set(*)
      character*132 textpart(16)
!     
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,istep,istat,
     &     n,i,j,key,nk,node,ifacequad(3,4),ifacetria(3,3),npt,
     &     mpcfreeold,ikmpc(*),ilmpc(*),id,idof,iline,ipol,inl,
     &     ipoinp(2,*),inp(3,*),ipoinpc(0:*),irefnode,lathyp(3,6),inum,
     &     jn,jt,jd,iside,nelem,jface,nopes,nface,nodef(8),nodel(8),
     &     ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),indexpret,
     &     k,ipos,nkold,nope,m,kon(*),ipkon(*),indexe,iset,nset,idir,
     &     istartset(*),iendset(*),ialset(*),index1,ics(2,*),mpcpret,
     &     mint,iflag,ithermal(*),ielem,three,in(3),node1,node2,isign,
     &     ndep,nind,kflag,ne,nkref,noderef,ier
!     
      real*8 coefmpc(*),xn(3),xt(3),xd(3),dd,co(3,*),dcs(*),area,
     &     areanodal(8),xl2(3,8),xi,et,weight,shp2(7,8),t0(*),
     &     xs2(3,2),xsj2(3),xsj,yn(3),r,xnl(3)
!     
      include "gauss.f"
!     
!     latin hypercube positions in a 3 x 3 matrix
!     
      data lathyp /1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1/
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
!     nodes per face for quad elements
!     
      data ifacequad /1,2,5,
     &     2,3,6,
     &     3,4,7,
     &     4,1,8/
!     
!     nodes per face for tria elements
!     
      data ifacetria /1,2,4,
     &     2,3,5,
     &     3,1,6/
!     
!     flag for shape functions
!     
      data iflag /2/
!     
      if(istep.gt.0) then
        write(*,*) 
     &       '*ERROR reading *PRE-TENSION SECTION: *EQUATION should'
        write(*,*) '       be placed before all step definitions'
        ier=1
        return
      endif
!     
      ielem=0
      do i=1,81
        surface(i:i)=' '
      enddo
!     
      do i=2,n
        if(textpart(i)(1:8).eq.'SURFACE=') then
          if(ielem.ne.0) then
            write(*,*) '*ERROR reading PRE-TENSION SECTION:'
            write(*,*) '       ELEMENT and SURFACE are'
            write(*,*) '       mutually exclusive'
            call inputerror(inpc,ipoinpc,iline,
     &           "*PRE-TENSION SECTION%",ier)
            return
          endif
          surface=textpart(i)(9:88)
          ipos=index(surface,' ')
          surface(ipos:ipos)='T'
        elseif(textpart(i)(1:5).eq.'NODE=') then
          read(textpart(i)(6:15),'(i10)',iostat=istat) irefnode
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*PRE-TENSION SECTION%",ier)
            return
          endif
          if((irefnode.gt.nk).or.(irefnode.le.0)) then
            write(*,*) '*ERROR reading *PRE-TENSION SECTION:'
            write(*,*) '       node ',irefnode,' is not defined'
            ier=1
            return
          endif
        elseif(textpart(i)(1:8).eq.'ELEMENT=') then
          if(surface(1:1).ne.' ') then
            write(*,*) '*ERROR reading PRE-TENSION SECTION:'
            write(*,*) '       ELEMENT and SURFACE are'
            write(*,*) '       mutually exclusive'
            call inputerror(inpc,ipoinpc,iline,
     &           "*PRE-TENSION SECTION%",ier)
            return
          endif
          read(textpart(i)(9:18),'(i10)',iostat=istat) ielem
          if(istat.gt.0) then
            write(*,*) '*ERROR reading PRE-TENSION SECTION:'
            write(*,*) '       cannot read element number'
            call inputerror(inpc,ipoinpc,iline,
     &           "*PRE-TENSION SECTION%",ier)
            return
          endif
          if((ielem.gt.ne).or.(ielem.le.0)) then
            write(*,*) '*ERROR reading PRE-TENSION SECTION:'
            write(*,*) '       element',ielem,'is not defined'
            call inputerror(inpc,ipoinpc,iline,
     &           "*PRE-TENSION SECTION%",ier)
            return
          endif
        else
          write(*,*) 
     &         '*WARNING reading *PRE-TENSION SECTION: parameter not recog
     &nized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*PRE-TENSION SECTION%")
        endif
      enddo
!     
!     checking whether the surface exists and is an element face
!     surface
!     
      if(surface(1:1).ne.' ') then
c     iset=0
c     do i=1,nset
c     if(set(i).eq.surface) then
c     iset=i
c     exit
c     endif
c     enddo
        call cident81(set,surface,nset,id)
        iset=0
        if(id.gt.0) then
          if(surface.eq.set(id)) then
            iset=id
          endif
        endif
        if(iset.eq.0) then
          write(*,*) 
     &        '*ERROR reading *PRE-TENSION SECTION: nonexistent surface'
          write(*,*) '       or surface consists of nodes'
          call inputerror(inpc,ipoinpc,iline,
     &         "*PRE-TENSION SECTION%",ier)
          return
        endif
      elseif(ielem.gt.0) then
        if(lakon(ielem)(1:3).ne.'B31') then
          write(*,*) '*ERROR reading PRE-TENSION SECTION:'
          write(*,*) '       element',ielem,' is not a linear'
          write(*,*) '       beam element'
          ier=1
          return
        endif
      else
        write(*,*) '*ERROR reading PRE-TENSION SECTION:'
        write(*,*) '       either the parameter SURFACE or the'
        write(*,*) '       parameter ELEMENT must be used'
        call inputerror(inpc,ipoinpc,iline,
     &       "*PRE-TENSION SECTION%",ier)
        return
      endif
!     
!     reading the normal vector and normalizing
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      normal=.true.
      do i=1,3
        read(textpart(i)(1:20),'(f20.0)',iostat=istat) xn(i)
        if(istat.gt.0) then
          normal=.false.
          exit
        endif
      enddo
      if(normal) then
        dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
        do i=1,3
          xn(i)=xn(i)/dd
        enddo
      endif
!     
!     beam pretension
!     
      if(ielem.gt.0) then
        indexe=ipkon(ielem)
!     
!     removing the beam element (is not needed any more; if not removed
!     the beam should have a small cross section compared to the
!     pre-tensioned structure in order to get the same results as in 
!     ABAQUS
!     
        ipkon(ielem)=-1
!     
        node1=kon(indexe+1)
        node2=kon(indexe+2)
!     
!     if the normal is not defined by the user:
!     normal = connection from node 1 to node 2
!     
        if(.not.normal) then
          do j=1,3
            xn(j)=co(j,node2)-co(j,node1)
          enddo
          dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
          do i=1,3
            xn(i)=xn(i)/dd
          enddo
          write(*,*) '*INFO: pre-tension normal for element',ielem
          write(*,*) '       ',xn(1),xn(2),xn(3)
          write(*,*) 
        endif
!     
        do i=1,3
          yn(i)=dabs(xn(i))
          in(i)=i
        enddo
!     
!     sorting yn in descending order
!     
        three=3
        kflag=-2
        call dsort(yn,in,three,kflag)
!     
!     equation: u_2*n-u_1*n+d=0
!     
!     looking for the largest coefficient
!     
        do i=1,3
          if(yn(i).lt.1.d-10) cycle
!     
!     default dependent and independent nodes
!     
          ndep=node2
          nind=node1
          isign=1
!     
          idof=8*(node2-1)+in(i)
          call nident(ikmpc,idof,nmpc,id)
          if(id.ne.0) then
            if(ikmpc(id).eq.idof) then
              idof=8*(node1-1)+in(i)
              call nident(ikmpc,idof,nmpc,id)
              if(id.ne.0) then
                if(ikmpc(id).eq.idof) then
                  cycle
                endif
              endif
              ndep=node1
              nind=node2
              isign=-1
            endif
          endif
!     
          nmpc=nmpc+1
          if(nmpc.gt.nmpc_) then
            write(*,*) 
     &           '*ERROR reading *PRE-TENSION SECTION: increase nmpc_'
            ier=1
            return
          endif
          labmpc(nmpc)='PRETENSION          '
          ipompc(nmpc)=mpcfree
!     
!     updating ikmpc and ilmpc
!     
          do j=nmpc,id+2,-1
            ikmpc(j)=ikmpc(j-1)
            ilmpc(j)=ilmpc(j-1)
          enddo
          ikmpc(id+1)=idof
          ilmpc(id+1)=nmpc
!     
c     jn=in(i)
          idir=in(i)
!     
          if(dabs(xn(idir)).gt.1.d-10) then
            nodempc(1,mpcfree)=ndep
            nodempc(2,mpcfree)=idir
            coefmpc(mpcfree)=isign*xn(idir)
            mpcfree=nodempc(3,mpcfree)
!     
            nodempc(1,mpcfree)=nind
            nodempc(2,mpcfree)=idir
            coefmpc(mpcfree)=-isign*xn(idir)
            mpcfree=nodempc(3,mpcfree)
          endif
!     
          idir=idir+1
          if(idir.eq.4) idir=1
          if(dabs(xn(idir)).gt.1.d-10) then
            nodempc(1,mpcfree)=ndep
            nodempc(2,mpcfree)=idir
            coefmpc(mpcfree)=isign*xn(idir)
            mpcfree=nodempc(3,mpcfree)
!     
            nodempc(1,mpcfree)=nind
            nodempc(2,mpcfree)=idir
            coefmpc(mpcfree)=-isign*xn(idir)
            mpcfree=nodempc(3,mpcfree)
          endif
!     
          idir=idir+1
          if(idir.eq.4) idir=1
          if(dabs(xn(idir)).gt.1.d-10) then
            nodempc(1,mpcfree)=ndep
            nodempc(2,mpcfree)=idir
            coefmpc(mpcfree)=isign*xn(idir)
            mpcfree=nodempc(3,mpcfree)
!     
            nodempc(1,mpcfree)=nind
            nodempc(2,mpcfree)=idir
            coefmpc(mpcfree)=-isign*xn(idir)
            mpcfree=nodempc(3,mpcfree)
          endif
!     
!     inhomogeneous term
!     
          nodempc(1,mpcfree)=irefnode
          nodempc(2,mpcfree)=1
          coefmpc(mpcfree)=1.d0
          mpcfreeold=mpcfree
          mpcfree=nodempc(3,mpcfree)
          nodempc(3,mpcfreeold)=0
!     
          if(normal) call getnewline(inpc,textpart,istat,n,key,iline,
     &         ipol,inl,ipoinp,inp,ipoinpc)
!     
          return
        enddo
        write(*,*) '*ERROR reading *PRE-TENSION SECTION'
        write(*,*) '       all DOFS of the beam elements'
        write(*,*) '       have been used previously'
        ier=1
        return
      endif
!     
!     pre-tension section is a surface
!     
!     check whether normal was defined by the user  
!     if not, calculate it as the mean normal on the surface
!     
      if(.not.normal) then
        do j=1,3
          xn(j)=0.d0
        enddo
        do k=1,iendset(iset)-istartset(iset)+1
          iside=ialset(istartset(iset)+k-1)
          nelem=int(iside/10.d0)
          indexe=ipkon(nelem)
          jface=iside-10*nelem
          call calcnormal(nelem,jface,lakon,co,xnl,indexe,kon)
          do j=1,3
            xn(j)=xn(j)+xnl(j)
          enddo
        enddo
        dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
        do i=1,3
          xn(i)=xn(i)/dd
        enddo
        write(*,*) '*INFO: pre-tension normal for surface ',
     &       surface(1:ipos-1)
        write(*,*) '       ',xn(1),xn(2),xn(3)
        write(*,*) 
      endif
!     
!     finding a unit vector xt perpendicular to the normal vector
!     using a unit vector in x or in y 
!     
      if(dabs(xn(1)).lt.0.95d0) then
        xt(1)=1.d0-xn(1)*xn(1)
        xt(2)=-xn(1)*xn(2)
        xt(3)=-xn(1)*xn(3)
      else
        xt(1)=-xn(2)*xn(1)
        xt(2)=1.d0-xn(2)*xn(2)
        xt(3)=-xn(2)*xn(3)
      endif
      dd=dsqrt(xt(1)*xt(1)+xt(2)*xt(2)+xt(3)*xt(3))
      do i=1,3
        xt(i)=xt(i)/dd
      enddo
!     
!     xd=xn x xt
!     
      xd(1)=xn(2)*xt(3)-xn(3)*xt(2)
      xd(2)=xn(3)*xt(1)-xn(1)*xt(3)
      xd(3)=xn(1)*xt(2)-xn(2)*xt(1)
!     
!     generating a Latin hypercube
!     checking which DOF's of xn, xt and xd are nonzero
!     
      do inum=1,6
        if((dabs(xn(lathyp(1,inum))).gt.1.d-3).and.
     &       (dabs(xt(lathyp(2,inum))).gt.1.d-3).and.
     &       (dabs(xd(lathyp(3,inum))).gt.1.d-3)) exit
      enddo
      jn=lathyp(1,inum)
      jt=lathyp(2,inum)
      jd=lathyp(3,inum)
!     
!     generating the MPCs
!     
      indexpret=0
      nkold=nk
      m=iendset(iset)-istartset(iset)+1
!     
!     number of distinct pre-strain nodes for the present keyword
!     
      npt=0
      area=0.d0
!     
!     loop over all element faces belonging to the surface
!     
!     twod=.true. means that the face belongs to a plane stress,
!     plane strain or axisymmetric element
!     
!     quad=.true. means that the face belongs to a 4-sided
!     plane stress, plane strain, axisymmetric or shell element
!     
      do k=1,m
        twod=.false.
        quad=.false.
        iside=ialset(istartset(iset)+k-1)
        nelem=int(iside/10.d0)
        indexe=ipkon(nelem)
        jface=iside-10*nelem
!     
!     nodes: #nodes in the face
!     nface: # of faces belonging to the element
!     the nodes are stored in nodef(*)
!     
        if(lakon(nelem)(4:4).eq.'2') then
          nopes=8
          nface=6
        elseif(lakon(nelem)(3:4).eq.'D8') then
          nopes=4
          nface=6
        elseif(lakon(nelem)(4:5).eq.'10') then
          nopes=6
          nface=4
        elseif(lakon(nelem)(3:4).eq.'D4') then
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
        elseif(lakon(nelem)(3:4).eq.'D6') then
          if(jface.le.2) then
            nopes=3
          else
            nopes=4
          endif
          nface=5
          nope=6
        elseif((lakon(nelem)(2:2).eq.'8').or.
     &         (lakon(nelem)(4:4).eq.'8')) then
!     
!     8-node 2-D elements
!     
          nopes=3
          nface=4
          quad=.true.
          jface=jface-2
          if(lakon(nelem)(4:4).eq.'8') then
            twod=.true.
c     jface=jface-2
          endif
        elseif((lakon(nelem)(2:2).eq.'6').or.
     &         (lakon(nelem)(4:4).eq.'6')) then
!     
!     6-node 2-D elements
!     
          nopes=3
          nface=3
          jface=jface-2
          if(lakon(nelem)(4:4).eq.'6') then
            twod=.true.
c     jface=jface-2
          endif
        elseif((lakon(nelem)(2:2).eq.'4').or.
     &         (lakon(nelem)(4:4).eq.'4')) then
!     
!     4-node 2-D elements
!     
          nopes=2
          nface=4
          quad=.true.
          jface=jface-2
          if(lakon(nelem)(4:4).eq.'4') then
            twod=.true.
c     jface=jface-2
          endif
        elseif((lakon(nelem)(2:2).eq.'3').or.
     &         (lakon(nelem)(4:4).eq.'3')) then
!     
!     3-node 2-D elements
!     
          nopes=2
          nface=3
          jface=jface-2
          if(lakon(nelem)(4:4).eq.'3') then
            twod=.true.
c     jface=jface-2
          endif
        else
          cycle
        endif
!     
!     determining the nodes of the face
!     
        if(nface.eq.3) then
          do i=1,nopes
            nodef(i)=kon(indexe+ifacetria(i,jface))
            nodel(i)=ifacetria(i,jface)
          enddo
        elseif(nface.eq.4) then
          if(quad) then
            do i=1,nopes
              nodef(i)=kon(indexe+ifacequad(i,jface))
              nodel(i)=ifacequad(i,jface)
            enddo
          else
            do i=1,nopes
              nodef(i)=kon(indexe+ifacet(i,jface))
              nodel(i)=ifacet(i,jface)
            enddo
          endif
        elseif(nface.eq.5) then
          if(nope.eq.6) then
            do i=1,nopes
              nodef(i)=kon(indexe+ifacew1(i,jface))
              nodel(i)=ifacew1(i,jface)
            enddo
          elseif(nope.eq.15) then
            do i=1,nopes
              nodef(i)=kon(indexe+ifacew2(i,jface))
              nodel(i)=ifacew2(i,jface)
            enddo
          endif
        elseif(nface.eq.6) then
          do i=1,nopes
            nodef(i)=kon(indexe+ifaceq(i,jface))
            nodel(i)=ifaceq(i,jface)
          enddo
        endif
!     
!     loop over the nodes belonging to the face   
!     ics(1,*): pretension node
!     ics(2,*): corresponding partner node
!     dcs(*): area corresponding to pretension node   
!     
        do i=1,nopes
          node=nodef(i)
          call nident2(ics,node,npt,id)
          if(id.gt.0) then
            if(ics(1,id).eq.node) then
!     
!     node was already treated: replacing the node
!     by the partner node
!     
              kon(indexe+nodel(i))=ics(2,id)
              cycle
            endif
          endif
!     
!     generating a partner node
!     
          nk=nk+1
!     
!     coordinates for the new node
!     
          do j=1,3
            co(j,nk)=co(j,node)
          enddo
!     
!     updating the topology
!     
          kon(indexe+nodel(i))=nk
!     
!     updating ics
!     
          npt=npt+1
          do j=npt,id+2,-1
            ics(1,j)=ics(1,j-1)
            ics(2,j)=ics(2,j-1)
c     dcs(j)=dcs(j-1)
          enddo
          ics(1,id+1)=node
          ics(2,id+1)=nk
c     dcs(id+1)=0.d0
!     
!     first MPC perpendicular to the normal direction
!     
          idof=8*(nk-1)+jt
          call nident(ikmpc,idof,nmpc,id)
!     
          nmpc=nmpc+1
          if(nmpc.gt.nmpc_) then
            write(*,*) 
     &           '*ERROR reading *PRE-TENSION SECTION: increase nmpc_'
            ier=1
            return
          endif
          ipompc(nmpc)=mpcfree
          labmpc(nmpc)='                    '
!     
!     updating ikmpc and ilmpc
!     
          do j=nmpc,id+2,-1
            ikmpc(j)=ikmpc(j-1)
            ilmpc(j)=ilmpc(j-1)
          enddo
          ikmpc(id+1)=idof
          ilmpc(id+1)=nmpc
!     
          idir=jt
          if(dabs(xt(idir)).gt.1.d-10) then
            nodempc(1,mpcfree)=nk
            nodempc(2,mpcfree)=idir
            coefmpc(mpcfree)=-xt(idir)
c     write(*,*) 'pretensionsections ',nk,idir,-xt(idir)
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
          endif
!     
          idir=idir+1
          if(idir.eq.4) idir=1
          if(dabs(xt(idir)).gt.1.d-10) then
            nodempc(1,mpcfree)=nk
            nodempc(2,mpcfree)=idir
            coefmpc(mpcfree)=-xt(idir)
c     write(*,*) 'pretensionsections ',nk,idir,-xt(idir)
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
          endif
!     
          idir=idir+1
          if(idir.eq.4) idir=1
          if(dabs(xt(idir)).gt.1.d-10) then
            nodempc(1,mpcfree)=nk
            nodempc(2,mpcfree)=idir
            coefmpc(mpcfree)=-xt(idir)
c     write(*,*) 'pretensionsections ',nk,idir,-xt(idir)
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
          endif
!     
          idir=jt
          if(dabs(xt(idir)).gt.1.d-10) then
            nodempc(1,mpcfree)=node
            nodempc(2,mpcfree)=jt
            coefmpc(mpcfree)=xt(idir)
c     write(*,*) 'pretensionsections ',node,jt,xt(idir)
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
          endif
!     
          idir=idir+1
          if(idir.eq.4) idir=1
          if(dabs(xt(idir)).gt.1.d-10) then
            nodempc(1,mpcfree)=node
            nodempc(2,mpcfree)=idir
            coefmpc(mpcfree)=xt(idir)
c     write(*,*) 'pretensionsections ',node,idir,xt(idir)
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
          endif
!     
          idir=idir+1
          if(idir.eq.4) idir=1
          if(dabs(xt(idir)).gt.1.d-10) then
            nodempc(1,mpcfree)=node
            nodempc(2,mpcfree)=idir
            coefmpc(mpcfree)=xt(idir)
c     write(*,*) 'pretensionsections ',node,idir,xt(idir)
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
          endif
          nodempc(3,mpcfreeold)=0
!     
!     second MPC perpendicular to the normal direction
!     
          if(.not.twod) then
            idof=8*(nk-1)+jd
            call nident(ikmpc,idof,nmpc,id)
!     
            nmpc=nmpc+1
            if(nmpc.gt.nmpc_) then
              write(*,*) 
     &             '*ERROR reading *PRE-TENSION SECTION: increase nmpc_'
              ier=1
              return
            endif
            labmpc(nmpc)='                    '
            ipompc(nmpc)=mpcfree
!     
!     updating ikmpc and ilmpc
!     
            do j=nmpc,id+2,-1
              ikmpc(j)=ikmpc(j-1)
              ilmpc(j)=ilmpc(j-1)
            enddo
            ikmpc(id+1)=idof
            ilmpc(id+1)=nmpc
!     
            idir=jd
            if(dabs(xd(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=nk
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=-xd(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xd(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=nk
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=-xd(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xd(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=nk
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=-xd(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=jd
            if(dabs(xd(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=node
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=xd(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xd(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=node
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=xd(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xd(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=node
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=xd(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
            nodempc(3,mpcfreeold)=0
          endif
!     
!     MPC in normal direction
!     
!     check whether initialized
!     
          if(indexpret.eq.0) then
            idof=8*(nk-1)+jn
            call nident(ikmpc,idof,nmpc,id)
!     
            nmpc=nmpc+1
            if(nmpc.gt.nmpc_) then
              write(*,*) 
     &             '*ERROR reading *PRE-TENSION SECTION: increase nmpc_'
              ier=1
              return
            endif
            labmpc(nmpc)='PRETENSION          '
            ipompc(nmpc)=mpcfree
            mpcpret=nmpc
!     
!     taking the first node as reference node for the condition:
!     "the displacement in pre-tension direction in all nodes
!     should be the same as for the reference node" (pre-tension
!     surfaces stay parallel to each other)
!     
            nkref=nk
            noderef=node
!     
!     updating ikmpc and ilmpc
!     
            do j=nmpc,id+2,-1
              ikmpc(j)=ikmpc(j-1)
              ilmpc(j)=ilmpc(j-1)
            enddo
            ikmpc(id+1)=idof
            ilmpc(id+1)=nmpc
c     else
c     nodempc(3,indexpret)=mpcfree
          endif
!     
!     MPC's specifying that the pre-tension surfaces should stay
!     parallel
!     
          if(indexpret.ne.0) then
            idof=8*(nk-1)+jn
            call nident(ikmpc,idof,nmpc,id)
!     
            nmpc=nmpc+1
            if(nmpc.gt.nmpc_) then
              write(*,*) 
     &             '*ERROR reading *PRE-TENSION SECTION: increase nmpc_'
              ier=1
              return
            endif
            ipompc(nmpc)=mpcfree
            labmpc(nmpc)='                    '
!     
!     updating ikmpc and ilmpc
!     
            do j=nmpc,id+2,-1
              ikmpc(j)=ikmpc(j-1)
              ilmpc(j)=ilmpc(j-1)
            enddo
            ikmpc(id+1)=idof
            ilmpc(id+1)=nmpc
!     
            idir=jn
            if(dabs(xn(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=nk
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=-xn(idir)
c     write(*,*) 'pretensionsections ',nk,idir,-xn(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xn(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=nk
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=-xn(idir)
c     write(*,*) 'pretensionsections ',nk,idir,-xn(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xn(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=nk
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=-xn(idir)
c     write(*,*) 'pretensionsections ',nk,idir,-xn(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=jn
            if(dabs(xn(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=node
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=xn(idir)
c     write(*,*) 'pretensionsections ',node,idir,xn(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xn(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=node
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=xn(idir)
c     write(*,*) 'pretensionsections ',node,idir,xn(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xn(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=node
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=xn(idir)
c     write(*,*) 'pretensionsections ',node,idir,xn(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
!     corresponding terms in the reference node
!     
            idir=jn
            if(dabs(xn(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=nkref
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=xn(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xn(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=nkref
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=xn(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xn(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=nkref
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=xn(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=jn
            if(dabs(xn(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=noderef
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=-xn(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xn(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=noderef
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=-xn(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xn(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=noderef
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=-xn(idir)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
            nodempc(3,mpcfreeold)=0
!     
            nodempc(3,indexpret)=mpcfree
c     endif
          else
!     
!     governing pre-tension equation
!     
            idir=jn
            if(dabs(xn(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=nk
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=-xn(idir)
c     write(*,*) 'pretensionsectionsf1 ',nk,idir,-xn(idir)
              indexpret=mpcfree
              mpcfree=nodempc(3,mpcfree)
!     
              nodempc(1,mpcfree)=node
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=xn(idir)
c     write(*,*) 'pretensionsectionsf2 ',node,idir,xn(idir)
              indexpret=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xn(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=nk
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=-xn(idir)
c     write(*,*) 'pretensionsectionsf3 ',nk,idir,-xn(idir)
              indexpret=mpcfree
              mpcfree=nodempc(3,mpcfree)
!     
              nodempc(1,mpcfree)=node
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=xn(idir)
c     write(*,*) 'pretensionsectionsf4 ',node,idir,xn(idir)
              indexpret=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
!     
            idir=idir+1
            if(idir.eq.4) idir=1
            if(dabs(xn(idir)).gt.1.d-10) then
              nodempc(1,mpcfree)=nk
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=-xn(idir)
c     write(*,*) 'pretensionsectionsf5 ',nk,idir,-xn(idir)
              indexpret=mpcfree
              mpcfree=nodempc(3,mpcfree)
!     
              nodempc(1,mpcfree)=node
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=xn(idir)
c     write(*,*) 'pretensionsectionsf6 ',node,idir,xn(idir)
              indexpret=mpcfree
              mpcfree=nodempc(3,mpcfree)
            endif
          endif
!     
!     thermal MPC
!     
          if(ithermal(2).gt.0) then
            idof=8*(nk-1)
            call nident(ikmpc,idof,nmpc,id)
!     
            nmpc=nmpc+1
            if(nmpc.gt.nmpc_) then
              write(*,*) 
     &             '*ERROR reading *PRE-TENSION SECTION: increase nmpc_'
              ier=1
              return
            endif
!     
!     the label for purely mechanical calculations is
!     necessary in order to copy t1 at the newly created
!     nodes in file calinput.f
!     
            if(ithermal(2).gt.1) then
              labmpc(nmpc)='                    '
            else
              labmpc(nmpc)='THERMALPRET         '
            endif
!     
            ipompc(nmpc)=mpcfree
!     
!     updating ikmpc and ilmpc
!     
            do j=nmpc,id+2,-1
              ikmpc(j)=ikmpc(j-1)
              ilmpc(j)=ilmpc(j-1)
            enddo
            ikmpc(id+1)=idof
            ilmpc(id+1)=nmpc
!     
            nodempc(1,mpcfree)=nk
            nodempc(2,mpcfree)=0
            coefmpc(mpcfree)=1.d0
            mpcfree=nodempc(3,mpcfree)
!     
            nodempc(1,mpcfree)=node
            nodempc(2,mpcfree)=0
            coefmpc(mpcfree)=-1.d0
!     
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
            nodempc(3,mpcfreeold)=0
          endif
!     
!     initial value for new node
!     
          if(ithermal(1).gt.0) then
            t0(nk)=t0(node)
          endif
!     
        enddo
!     
!     calculating the area of the face and its contributions
!     to the facial nodes
!     
!     number of integration points
!     
        if(lakon(nelem)(3:5).eq.'D8R') then
          mint=1
        elseif(lakon(nelem)(3:4).eq.'D8') then
          mint=4
        elseif(lakon(nelem)(4:6).eq.'20R') then
          mint=4
        elseif(lakon(nelem)(4:4).eq.'2') then
          mint=9
        elseif(lakon(nelem)(4:5).eq.'10') then
          mint=3
        elseif(lakon(nelem)(4:4).eq.'4') then
          mint=1
        elseif(lakon(nelem)(3:4).eq.'D6') then
          mint=1
        elseif(lakon(nelem)(4:5).eq.'15') then
          if(jface.le.2) then
            mint=3
          else
            mint=4
          endif
!     
!     faces of 2-D elements    
!     
        elseif((lakon(nelem)(3:3).eq.'R').or.
     &         (lakon(nelem)(5:5).eq.'R')) then
!     
!     reduced integration
!     
          if(nopes.eq.2) then
            mint=1
          else
            mint=2
          endif
        else
!     
!     full integration
!     
          if(nopes.eq.2) then
            mint=2
          else
            mint=3
          endif
        endif
!     
        do i=1,nopes
c     areanodal(i)=0.d0
          do j=1,3
            xl2(j,i)=co(j,nodef(i))
          enddo
        enddo
!     
        do m=1,mint
          if((lakon(nelem)(3:5).eq.'D8R').or.
     &         ((lakon(nelem)(3:4).eq.'D6').and.(nopes.eq.4))) then
            xi=gauss2d1(1,m)
            et=gauss2d1(2,m)
            weight=weight2d1(m)
          elseif((lakon(nelem)(3:4).eq.'D8').or.
     &           (lakon(nelem)(4:6).eq.'20R').or.
     &           ((lakon(nelem)(4:5).eq.'15').and.
     &           (nopes.eq.8))) then
            xi=gauss2d2(1,m)
            et=gauss2d2(2,m)
            weight=weight2d2(m)
          elseif(lakon(nelem)(4:4).eq.'2') then
            xi=gauss2d3(1,m)
            et=gauss2d3(2,m)
            weight=weight2d3(m)
          elseif((lakon(nelem)(4:5).eq.'10').or.
     &           ((lakon(nelem)(4:5).eq.'15').and.
     &           (nopes.eq.6))) then
            xi=gauss2d5(1,m)
            et=gauss2d5(2,m)
            weight=weight2d5(m)
          elseif((lakon(nelem)(4:4).eq.'4').or.
     &           ((lakon(nelem)(3:4).eq.'D6').and.
     &           (nopes.eq.3))) then
            xi=gauss2d4(1,m)
            et=gauss2d4(2,m)
            weight=weight2d4(m)
!     
!     faces of 2-D elements    
!     
          elseif((lakon(nelem)(3:3).eq.'R').or.
     &           (lakon(nelem)(5:5).eq.'R')) then
!     
!     reduced integration
!     
            if(nopes.eq.2) then
              xi=gauss1d1(1,m)
              weight=weight1d1(m)
            else
              xi=gauss1d2(1,m)
              weight=weight1d2(m)
            endif
          else
!     
!     full integration
!     
            if(nopes.eq.2) then
              xi=gauss1d2(1,m)
              weight=weight1d2(m)
            else
              xi=gauss1d3(1,m)
              weight=weight1d3(m)
            endif
          endif
!     
c     write(*,*) 'pretension12 ',nopes
          if(nopes.eq.8) then
            call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
          elseif(nopes.eq.4) then
            call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
c     write(*,*) 'pretension11 ',xsj2(1)
          elseif(nopes.eq.6) then
            call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
          elseif((nopes.eq.3).and.(.not.twod).and.
     &           (lakon(nelem)(1:1).ne.'S')) then
            call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
          elseif(nopes.eq.3) then
!     
!     3-node line
!     
            call shape3l(xi,xl2,xsj2,xs2,shp2,iflag)
          else
!     
!     2-node line
!     
            call shape2l(xi,xl2,xsj2,xs2,shp2,iflag)
          endif
!     
!     calculating the total area and nodal area
!     
          if(.not.twod) then
            xsj=weight*dsqrt(xsj2(1)**2+xsj2(2)**2+xsj2(3)**2)
          elseif(lakon(nelem)(1:3).eq.'CAX') then
!     
!     radial distance is taken into account for the area
!     
            r=0.d0
            do i=1,nopes
              r=r+shp2(4,i)*xl2(1,i)
            enddo
c     xsj=weight*xsj2(1)*r
            xsj=weight*dsqrt(xsj2(1)**2+xsj2(2)**2+xsj2(3)**2)*r
          else
c     xsj=weight*xsj2(1)
            xsj=weight*dsqrt(xsj2(1)**2+xsj2(2)**2+xsj2(3)**2)
c     write(*,*) 'pretension10',weight
c     write(*,*) 'pretension10',xsj
          endif
          area=area+xsj
c     write(*,*) 'pretension10 ',area
c     do i=1,nopes
c     areanodal(i)=areanodal(i)+xsj*shp2(4,i)
c     enddo
!     
        enddo
!     
!     inserting the nodal area into field dcs
!     
        do i=1,nopes
          node=nodef(i)
          call nident2(ics,node,npt,id)
c     dcs(id)=dcs(id)+areanodal(i)
        enddo
!     
      enddo
!     
      nodempc(3,indexpret)=mpcfree
      nodempc(1,mpcfree)=irefnode
      nodempc(2,mpcfree)=1
      coefmpc(mpcfree)=area
c     write(*,*) 'pretensionsectionsf77 ',area
      mpcfreeold=mpcfree
      mpcfree=nodempc(3,mpcfree)
      nodempc(3,mpcfreeold)=0
!     
!     changing the coefficients of the pretension MPC
!     
      index1=ipompc(mpcpret)
      do
        node=nodempc(1,nodempc(3,index1))
        call nident2(ics,node,npt,id)
        do j=1,2
c     coefmpc(index1)=coefmpc(index1)*dcs(id)
          coefmpc(index1)=coefmpc(index1)*area
          index1=nodempc(3,index1)
        enddo
        if(nodempc(1,index1).eq.irefnode) exit
      enddo
!     
      if(normal) call getnewline(inpc,textpart,istat,n,key,iline,
     &     ipol,inl,ipoinp,inp,ipoinpc)
!     
      return
      end

