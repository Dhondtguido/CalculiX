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
      subroutine gentiedmpc(tieset,ntie,itietri,ipkon,kon,
     &     lakon,set,istartset,iendset,ialset,cg,straight,
     &     koncont,co,xo,yo,zo,x,y,z,nx,ny,nz,nset,
     &     ifaceslave,istartfield,iendfield,ifield,
     &     ipompc,nodempc,coefmpc,nmpc,nmpctied,mpcfree,ikmpc,ilmpc,
     &     labmpc,ithermal,tietol,nef,ncont,imastop,ikboun,nboun,kind,
     &     jobnamef)
!     
!     generates MPC's for the slave tied contact nodes
!     
      implicit none
!     
      character*1 kind
      character*8 lakon(*)
      character*20 labmpc(*)
      character*81 tieset(3,*),slavset,set(*)
      character*132 jobnamef
      character*256 fn
!     
      integer ntie,nset,istartset(*),iendset(*),ialset(*),iwrite,
     &     itietri(2,ntie),ipkon(*),kon(*),koncont(4,*),node,
     &     neigh(1),iflag,kneigh,i,j,k,l,isol,itri,ll,kflag,n,nx(*),
     &     ny(*),nz(*),nstart,ifaceq(8,6),ifacet(6,4),nboun,
     &     ifacew1(4,5),ifacew2(8,5),nelem,jface,indexe,imastop(3,*),
     &     nnodelem,nface,nope,nodef(8),idof,kstart,kend,jstart,id,
     &     jend,ifield(*),istartfield(*),iendfield(*),ifaceslave(*),
     &     ipompc(*),nodempc(3,*),nmpc,nmpctied,mpcfree,ikmpc(*),
     &     ilmpc(*),ithermal(*),nef,ncont,mpcfreeold,m,id1,ikboun(*),
     &     itriold,itrinew,ntriangle,ntriangle_,itriangle(100),
     &     ilen
!     
      real*8 cg(3,*),straight(16,*),co(3,*),p(3),
     &     dist,xo(*),yo(*),zo(*),x(*),y(*),z(*),pl(3,9),
     &     ratio(9),xi,et,coefmpc(*),tietol(4,*),tolloc
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
!     opening a file to store the nodes which are not connected
!     
      iwrite=0
      ilen=index(jobnamef,' ')-1
      fn=jobnamef(1:ilen)//'_WarnNodeMissTiedContact.nam'
      open(40,file=fn,status='unknown')
      write(40,*) '*NSET,NSET=WarnNodeMissTiedContact'
!     
      nmpctied=nmpc
!     
!     calculating a typical element size
!     
      tolloc=0.d0
      do i=1,ncont
        tolloc=tolloc+dabs(straight(1,i)*cg(1,i)+
     &       straight(2,i)*cg(2,i)+
     &       straight(3,i)*cg(3,i)+
     &       straight(4,i))
      enddo
      tolloc=0.025d0*tolloc/ncont
!     
!     maximum number of neighboring master triangles for a slave node
!     
      kflag=2
!     
      do i=1,ntie
        if(tieset(1,i)(81:81).ne.kind) cycle
!     
!     determining for which dofs MPC's have to be generated
!     
        if(kind.eq.'T') then
!     
!     thermomechanical ties or CFD
!     
          if(nef.gt.0) then
            if(ithermal(2).le.1) then
              kstart=1
              kend=4
            else
              kstart=0
              kend=4
            endif
          else
            if(ithermal(2).le.1) then
              kstart=1
              kend=3
            elseif(ithermal(2).eq.2) then
              kstart=0
              kend=0
            else
              kstart=0
              kend=3
            endif
          endif
        elseif(kind.eq.'E') then
!     
!     electromagnetic ties between the domains
!     
          if((tieset(1,i)(1:2).eq.'12').or.
     &         (tieset(1,i)(1:2).eq.'13').or.
     &         (tieset(1,i)(1:2).eq.'23')) then
            kstart=1
            kend=3
          elseif((tieset(1,i)(1:2).eq.'21').or.
     &           (tieset(1,i)(1:2).eq.'31')) then
            kstart=5
            kend=5
          endif
        endif
!     
        iflag=0
        kneigh=1
        slavset=tieset(2,i)
!     
!     default tolerance if none is specified 
!     
        if(tietol(1,i).lt.1.d-10) tietol(1,i)=tolloc
!     
!     determining the slave set
!     
        if(ifaceslave(i).eq.0) then
c     do j=1,nset
c     if(set(j).eq.slavset) then
c     exit
c     endif
c     enddo
          call cident81(set,slavset,nset,id)
          j=nset+1
          if(id.gt.0) then
            if(slavset.eq.set(id)) then
              j=id
            endif
          endif
          jstart=istartset(j)
          jend=iendset(j)
        else
          jstart=istartfield(i)
          jend=iendfield(i)
        endif
!     
        nstart=itietri(1,i)-1
        n=itietri(2,i)-nstart
        if(n.lt.kneigh) kneigh=n
        do j=1,n
          xo(j)=cg(1,nstart+j)
          x(j)=xo(j)
          nx(j)=j
          yo(j)=cg(2,nstart+j)
          y(j)=yo(j)
          ny(j)=j
          zo(j)=cg(3,nstart+j)
          z(j)=zo(j)
          nz(j)=j
        enddo
        call dsort(x,nx,n,kflag)
        call dsort(y,ny,n,kflag)
        call dsort(z,nz,n,kflag)
!     
        do j=jstart,jend
          if(((ifaceslave(i).eq.0).and.(ialset(j).gt.0)).or.
     &         (ifaceslave(i).eq.1)) then
!     
            if(ifaceslave(i).eq.0) then
              node=ialset(j)
            else
              node=ifield(j)
            endif
!     
            do k=1,3
              p(k)=co(k,node)
            enddo
!     
!     determining the kneigh neighboring master contact
!     triangle centers of gravity
!     
            call near3d(xo,yo,zo,x,y,z,nx,ny,nz,p(1),p(2),p(3),
     &           n,neigh,kneigh)
!     
            isol=0
!     
            itriold=0
            itri=neigh(1)+itietri(1,i)-1
            ntriangle=0
            ntriangle_=100
!     
            loop1: do
            do l=1,3
              ll=4*l-3
              dist=straight(ll,itri)*p(1)+
     &             straight(ll+1,itri)*p(2)+
     &             straight(ll+2,itri)*p(3)+
     &             straight(ll+3,itri)
c     if(dist.gt.1.d-6) then
c     if(dist.gt.tietol(1,i)) then
              if(dist.gt.tolloc) then
                itrinew=imastop(l,itri)
                if(itrinew.eq.0) then
c     write(*,*) '**border reached'
                  isol=-1
                  exit loop1
                elseif(itrinew.eq.itriold) then
c     write(*,*) '**solution in between triangles'
                  isol=itri
                  exit loop1
                else
                  call nident(itriangle,itrinew,ntriangle,id)
                  if(id.gt.0) then
                    if(itriangle(id).eq.itrinew) then
c     write(*,*) '**circular path; no solution'
                      isol=-2
                      exit loop1
                    endif
                  endif
                  ntriangle=ntriangle+1
                  if(ntriangle.gt.ntriangle_) then
c     write(*,*) '**too many iterations'
                    isol=-3
                    exit loop1
                  endif
                  do k=ntriangle,id+2,-1
                    itriangle(k)=itriangle(k-1)
                  enddo
                  itriangle(id+1)=itrinew
                  itriold=itri
                  itri=itrinew
                  cycle loop1
                endif
              elseif(l.eq.3) then
c     write(*,*) '**regular solution'
                isol=itri
                exit loop1
              endif
            enddo
          enddo loop1
!     
!     if an opposite triangle is found: check the distance
!     perpendicular to the triangle
!     
          if(isol.gt.0) then
            dist=dabs(straight(13,itri)*p(1)+
     &           straight(14,itri)*p(2)+
     &           straight(15,itri)*p(3)+
     &           straight(16,itri))
            if(dist.gt.tietol(1,i)) isol=0
          endif
!     
          if(isol.le.0) then
!     
!     no MPC is generated
!     
            write(*,*) '*WARNING in gentiedmpc: no tied MPC'
            write(*,*) '         generated for node ',node
            if(isol.eq.0) then
              write(*,*)
     &             '         opposite master face found but'
              write(*,*) '         too far away; distance: ',dist
              write(*,*) '         tolerance: ',tietol(1,i)
            else
              write(*,*) '         no opposite master face'
              write(*,*) '         found; in-face tolerance: ',
     &             tolloc
            endif
            write(40,*) node
            iwrite=1
          else
!     
            nelem=int(koncont(4,itri)/10.d0)
            jface=koncont(4,itri)-10*nelem
!     
            indexe=ipkon(nelem)
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
!     determining the nodes of the face
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
!     attaching the node with coordinates in p
!     to the face
!     
            do k=1,nnodelem
              do l=1,3
                pl(l,k)=co(l,nodef(k))
              enddo
            enddo
            call attach_2d(pl,p,nnodelem,ratio,dist,xi,et)
!     
!     adjusting the coordinates of the node (only
!     if the user did not specify ADJUST=NO on the
!     *TIE card)
!     
            if(tietol(2,i).gt.0.d0) then
              do k=1,3
                co(k,node)=p(k)
              enddo
            endif
!     
!     generating MPC's
!     
            do l=kstart,kend
              idof=8*(node-1)+l
              call nident(ikmpc,idof,nmpc,id)
              if(id.gt.0) then
                if(ikmpc(id).eq.idof) then
                  write(*,*) '*WARNING in gentiedmpc:'
                  write(*,*) '         DOF ',l,' of node ',
     &                 node,' is not active;'
                  write(*,*) '         no tied constraint ',
     &                 'is generated'
                  write(40,*) node
                  iwrite=1
                  cycle
                endif
              endif
!     
              call nident(ikboun,idof,nboun,id1)
              if(id1.gt.0) then
                if(ikboun(id1).eq.idof) then
                  write(*,*) '*WARNING in gentiedmpc:'
                  write(*,*) '         DOF ',l,' of node ',
     &                 node,' is not active;'
                  write(*,*) '         no tied constraint ',
     &                 'is generated'
                  write(40,*) node
                  iwrite=1
                  cycle
                endif
              endif
!     
              nmpc=nmpc+1
              labmpc(nmpc)='                    '
              ipompc(nmpc)=mpcfree
!     
!     updating ikmpc and ilmpc
!     
              do m=nmpc,id+2,-1
                ikmpc(m)=ikmpc(m-1)
                ilmpc(m)=ilmpc(m-1)
              enddo
              ikmpc(id+1)=idof
              ilmpc(id+1)=nmpc
!     
              nodempc(1,mpcfree)=node
              nodempc(2,mpcfree)=l
              coefmpc(mpcfree)=1.d0
              mpcfree=nodempc(3,mpcfree)
              if(mpcfree.eq.0) then
                write(*,*)
     &               '*ERROR in gentiedmpc: increase memmpc_'
                call exit(201)
              endif
              do k=1,nnodelem
                nodempc(1,mpcfree)=nodef(k)
                nodempc(2,mpcfree)=l
                coefmpc(mpcfree)=-ratio(k)
                mpcfreeold=mpcfree
                mpcfree=nodempc(3,mpcfree)
                if(mpcfree.eq.0) then
                  write(*,*)
     &                 '*ERROR in gentiedmpc: increase memmpc_'
                  call exit(201)
                endif
              enddo
              nodempc(3,mpcfreeold)=0
            enddo
!     
          endif
!     
        else
          node=ialset(j-2)
          do
            node=node-ialset(j)
            if(node.ge.ialset(j-1)) exit
!     
            do k=1,3
              p(k)=co(k,node)
            enddo
!     
!     determining the kneigh neighboring master contact
!     triangle centers of gravity
!     
            call near3d(xo,yo,zo,x,y,z,nx,ny,nz,p(1),p(2),p(3),
     &           n,neigh,kneigh)
!     
            isol=0
!     
            itriold=0
            itri=neigh(1)+itietri(1,i)-1
            ntriangle=0
            ntriangle_=100
!     
            loop2: do
            do l=1,3
              ll=4*l-3
              dist=straight(ll,itri)*p(1)+
     &             straight(ll+1,itri)*p(2)+
     &             straight(ll+2,itri)*p(3)+
     &             straight(ll+3,itri)
c     if(dist.gt.1.d-6) then
c     if(dist.gt.tietol(1,i)) then
              if(dist.gt.tolloc) then
                itrinew=imastop(l,itri)
                if(itrinew.eq.0) then
c     write(*,*) '**border reached'
                  isol=-1
                  exit loop2
                elseif(itrinew.eq.itriold) then
c     write(*,*) '**solution in between triangles'
                  isol=itri
                  exit loop2
                else
                  call nident(itriangle,itrinew,ntriangle,id)
                  if(id.gt.0) then
                    if(itriangle(id).eq.itrinew) then
c     write(*,*) '**circular path; no solution'
                      isol=-2
                      exit loop2
                    endif
                  endif
                  ntriangle=ntriangle+1
                  if(ntriangle.gt.ntriangle_) then
c     write(*,*) '**too many iterations'
                    isol=-3
                    exit loop2
                  endif
                  do k=ntriangle,id+2,-1
                    itriangle(k)=itriangle(k-1)
                  enddo
                  itriangle(id+1)=itrinew
                  itriold=itri
                  itri=itrinew
                  cycle loop2
                endif
              elseif(l.eq.3) then
c     write(*,*) '**regular solution'
                isol=itri
                exit loop2
              endif
            enddo
          enddo loop2
!     
!     if an opposite triangle is found: check the distance
!     perpendicular to the triangle
!     
          if(isol.gt.0) then
            dist=dabs(straight(13,itri)*p(1)+
     &           straight(14,itri)*p(2)+
     &           straight(15,itri)*p(3)+
     &           straight(16,itri))
            if(dist.gt.tietol(1,i)) isol=0
          endif
!     
!     check whether distance is larger than tietol(1,i):
!     no element is generated
!     
          if(isol.le.0) then
!     
!     no MPC is generated
!     
            write(*,*) '*WARNING in gentiedmpc: no tied MPC'
            write(*,*) '         generated for node ',node
            if(isol.eq.0) then
              write(*,*) '         master face too far away'
              write(*,*) '         distance: ',dist
              write(*,*) '         tolerance: ',tietol(1,i)
            else
              write(*,*) '         no corresponding master face'
              write(*,*) '         found; tolerance: ',
     &             tolloc
c     &                                tietol(1,i)
            endif
            write(40,*) node
            iwrite=1
          else
!     
            nelem=int(koncont(4,itri)/10.d0)
            jface=koncont(4,itri)-10*nelem
!     
            indexe=ipkon(nelem)
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
!     determining the nodes of the face
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
!     attaching the node with coordinates in p
!     to the face
!     
            do k=1,nnodelem
              do l=1,3
                pl(l,k)=co(l,nodef(k))
              enddo
            enddo
            call attach_2d(pl,p,nnodelem,ratio,dist,xi,et)
!     
!     adjusting the coordinates of the node (only
!     if the user did not specify ADJUST=NO on the
!     *TIE card)
!     
            if(tietol(2,i).gt.0.d0) then
              do k=1,3
                co(k,node)=p(k)
              enddo
            endif
!     
!     generating MPC's
!     
            do l=kstart,kend
              idof=8*(node-1)+l
              call nident(ikmpc,idof,nmpc,id)
              if(id.gt.0) then
                if(ikmpc(id).eq.idof) then
                  write(*,*) '*WARNING in gentiedmpc:'
                  write(*,*) '         DOF ',l,' of node ',
     &                 node,' is not active;'
                  write(*,*) '         no tied constraint ',
     &                 'is generated'
                  write(40,*) node
                  iwrite=1
                  cycle
                endif
              endif
!     
              call nident(ikboun,idof,nboun,id1)
              if(id1.gt.0) then
                if(ikboun(id1).eq.idof) then
                  write(*,*) '*WARNING in gentiedmpc:'
                  write(*,*) '         DOF ',l,' of node ',
     &                 node,' is not active;'
                  write(*,*) '         no tied constraint ',
     &                 'is generated'
                  write(40,*) node
                  iwrite=1
                  cycle
                endif
              endif
!     
              nmpc=nmpc+1
              labmpc(nmpc)='                    '
              ipompc(nmpc)=mpcfree
!     
!     updating ikmpc and ilmpc
!     
              do m=nmpc,id+2,-1
                ikmpc(m)=ikmpc(m-1)
                ilmpc(m)=ilmpc(m-1)
              enddo
              ikmpc(id+1)=idof
              ilmpc(id+1)=nmpc
!     
              nodempc(1,mpcfree)=node
              nodempc(2,mpcfree)=l
              coefmpc(mpcfree)=1.d0
              mpcfree=nodempc(3,mpcfree)
              if(mpcfree.eq.0) then
                write(*,*)
     &               '*ERROR in gentiedmpc: increase memmpc_'
                call exit(201)
              endif
              do k=1,nnodelem
                nodempc(1,mpcfree)=nodef(k)
                nodempc(2,mpcfree)=l
                coefmpc(mpcfree)=-ratio(k)
                mpcfreeold=mpcfree
                mpcfree=nodempc(3,mpcfree)
                if(mpcfree.eq.0) then
                  write(*,*)
     &                 '*ERROR in gentiedmpc: increase memmpc_'
                  call exit(201)
                endif
              enddo
              nodempc(3,mpcfreeold)=0
            enddo
          endif
!     
        enddo
      endif
      enddo
      enddo
!     
!     number of tied MPC's
!     
      nmpctied=nmpc-nmpctied
!     
      if(iwrite.eq.1) then
        write(*,*) '*INFO in gentiedmpc:'
        write(*,*) '      failed nodes are stored in file'
        write(*,*) '      ',fn(1:ilen+28)
        write(*,*) '      This file can be loaded into'
        write(*,*) '      an active cgx-session by typing'
        write(*,*) 
     &       '      read ',fn(1:ilen+28),' inp'
        write(*,*)
        close(40)
      else
        close(40,status='delete')
      endif
!     
      return
!     
      end

