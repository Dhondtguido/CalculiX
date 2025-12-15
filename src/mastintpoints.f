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
!
!     Determining the location of the integration points in the
!     master surfaces opposite of slave surface sideloadloc(2:2)
!     of element nelemloadloc. This location depends on the triangulation of 
!     the master surface. For the master surface loads the local
!     coordinates and the integration weight is stored in pmastload,
!     the slave surface is stored in imastload
!
      subroutine mastintpoints(ipkon,kon,lakon,straight,
     &     nintpoint,koncont,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,
     &     imastop,mi,ncont,ipe,ime,
     &     nelemload,sideload,nload,nload_,imastload,pmastload,
     &     nelemloadloc,sideloadloc)
!
!     Based on slavintpoints.f
!
      implicit none
!
      character*8 lakon(*)
      character*20 sideload(*),sideloadloc
!     
      integer nintpoint,imastop(3,*),ncont,nelemload(2,*),nload,
     &     ipkon(*),kon(*),koncont(4,*),node,neigh(10),iflag,kneigh,
     &     j,k,l,ii,itri,nx(*),ny(*),nz(*),nelemm,jfacem,nload_,
     &     indexe,nopesm,nope,ifaces,nelems,jfaces,mi(*),
     &     m,nopes,konl(20),id,maface(8),nmaface,imastload(2,*),
     &     mafacecorner(8,8),line,iactiveline(2,3*ncont),
     &     icoveredmelem(3*ncont),nactiveline,ipe(*),ime(4,*),k1,j1,
     &     ncoveredmelem,nodem(8),nodepg(8),nelemloadloc,
     &     il,nodel,getnodel,ifacem,idummy,nopem,npg,k2,j2
!     
      real*8 straight(16,*),co(3,*),vold(0:mi(2),*),pmastload(3,*),
     &     xo(*),yo(*),zo(*),x(*),y(*),z(*),
     &     xl2m(3,8),xl2s(3,8),xlpg(3,8),
     &     pmiddle(3),xl2sr(3,8),xl3sp(3,8),xl3s(3,8),
     &     dd,areaslav,al,xn(3),slavstraight(36),err,xquad(2,8),
     &     xtri(2,6),xi,et,xsj2(3),xs2(3,2),shp2(7,8),anglesm
!
      data iflag /2/
!
      data xquad /-1.d0,-1.d0,
     &     1.d0,-1.d0,
     &     1.d0,1.d0,
     &     -1.d0,1.d0,
     &     0.d0,-1.d0,
     &     1.d0,0.d0,
     &     0.d0,1.d0,
     &     -1.d0,0.d0/
!     
      data xtri /0.d0,0.d0,
     &     1.d0,0.d0,
     &     0.d0,1.d0,
     &     0.5d0,0.d0,
     &     0.5d0,0.5d0,
     &     0.d0,0.5d0/
!     
      kneigh=1
      err=0.1d0
      areaslav=0.d0
!     
!     Research of the contact integration points
!     
      nelems=nelemloadloc
      if(ipkon(nelems).lt.0) then
        write(*,*) '*WARNING in slavintpoints'
        write(*,*) '         element ',nelems,' on slave contact'
        write(*,*) '         surface does not exist'
        return
      endif
      read(sideloadloc(2:2),'(i1)') jfaces
      ifaces=10*nelems+jfaces
!     
!     get nope,nopes
!     
      call faceinfo(nelems,jfaces,lakon,nope,nopes,idummy)
!     
!     actual position of the nodes belonging to the
!     slave surface
!     
      do j=1,nope
         konl(j)=kon(ipkon(nelems)+j)
      enddo
!     
      do m=1,nopes
         do j=1,3
            nodel=getnodel(m,jfaces,nope)
            xl2s(j,m)=co(j,konl(nodel))+
     &           vold(j,konl(nodel))     
         enddo
      enddo  
!     
!     slightly reducing the size of the slave surface in
!     an aleatoric way
!     
      do j=1,3
         pmiddle(j)=0.d0
         do m=1,nopes
            pmiddle(j)=pmiddle(j)+xl2s(j,m)
         enddo
         pmiddle(j)=pmiddle(j)/nopes
      enddo
!
!     sort vertices for quadratic elements in ccw direction
!     
      if((nopes.eq.3).or.(nopes.eq.4))then
         do j=1,nopes
            do k=1,3
               xl3s(k,j)=xl2s(k,j)
            enddo
         enddo
      else
         do j=1,int(nopes/2.d0)
            do k=1,3
               xl3s(k,2*j-1)=xl2s(k,j)
               xl3s(k,2*j)=xl2s(k,(int(nopes/2.d0))+j)           
            enddo
         enddo
      endif
!     
!     calculate the mean normal vector xn on the slave face
!     +
!     determine the equations of the slave face
!     (mean)plane and of the planes perpendicular to its 
!     piecewise linear approximation 
!     
      if(nopes.eq.3) then
         call straighteq3d(xl3s,slavstraight)
         do k=1,3
            xn(k)=slavstraight(4*nopes+k)
         enddo               
      else
         do k=1,3
            xn(k)=0.d0
         enddo
!     
         do m=1,nopes
            if((nopes.eq.4).or.(nopes.eq.8))then
               xi=xquad(1,m)
               et=xquad(2,m)
            else
               xi=xtri(1,m)
               et=xtri(2,m)
            endif
            if(nopes.eq.8)then
               call shape8q(xi,et,xl2s,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.4)then
               call shape4q(xi,et,xl2s,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.6)then
               call shape6tri(xi,et,xl2s,xsj2,xs2,shp2,iflag)
            else
               call shape3tri(xi,et,xl2s,xsj2,xs2,shp2,iflag)
            endif   
            dd=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)
     &          +xsj2(3)*xsj2(3))
            xsj2(1)=xsj2(1)/dd
            xsj2(2)=xsj2(2)/dd
            xsj2(3)=xsj2(3)/dd
!     
            do k=1,3
               xn(k)=xn(k)
     &              +xsj2(k)
            enddo
         enddo 
!     
!     normalizing the mean normal on the slave surface
!     
         dd=dsqrt(xn(1)**2+xn(2)**2+xn(3)**2)
         do k=1,3
            xn(k)=xn(k)/dd
         enddo           
         call approxplane(xl3s,slavstraight,xn,nopes)
      endif
!     
!     Project slave nodes to meanplane, needed for Sutherland-Hodgman
!     
      do j=1,nopes
         al=-xn(1)*xl3s(1,j)-xn(2)*
     &        xl3s(2,j)-xn(3)*xl3s(3,j)-
     &        slavstraight(nopes*4+4)
         if(nopes.ne.3)then
            do k=1,3
               xl3sp(k,j)=xl3s(k,j)+al*xn(k)
            enddo
         else
            do k=1,3
               xl3sp(k,j)=xl3s(k,j)
            enddo
         endif
      enddo
!     
!     determine the master faces corresponding to the corner
!     nodes
!     
      nmaface=0
      do j=1,8
         maface(j)=0
         do k=1,8
            mafacecorner(j,k)=0
         enddo
      enddo
!     
      do j=1,3
         do m=1,nopes
            xl2sr(j,m)=xl2s(j,m)-2*err*(xl2s(j,m)-pmiddle(j))
         enddo
      enddo
      do j=1,nopes
         call neartriangle_load(xl2sr(1,j),xn,xo,yo,zo,x,y,z,nx,ny,nz,
     &        neigh,kneigh,ncont,straight,imastop,itri)
         nodel=getnodel(j,jfaces,nope)
         node=konl(nodel) 
         if(itri.eq.0) cycle
         ifacem=koncont(4,itri)
!     
         call nident(maface,ifacem,nmaface,id)
         if(id.gt.0) then
            if(maface(id).eq.ifacem) then
               mafacecorner(j,id)=1
               cycle
            endif
         endif
!     
!     triangle was not covered yet: add to stack
!     
!     angle criteria 
!
         anglesm=xn(1)*straight(13,itri)
     &        +xn(2)*straight(14,itri)
     &        +xn(3)*straight(15,itri)
!     
         if(anglesm.lt.-0.7d0)then
!     angle between surface normals between 135 and 225 degree
!     angle between surfaces between 0 and 45 degree
            nmaface=nmaface+1
            do k=nmaface,id+2,-1
               maface(k)=maface(k-1)
               do m=1,j-1
                  mafacecorner(m,k)=mafacecorner(m,k-1)
               enddo
            enddo
            maface(id+1)=ifacem
            mafacecorner(j,id+1)=1
            do m=1,j-1
               mafacecorner(m,id+1)=0
            enddo 
         endif              
      enddo
!     
      nactiveline=0
!     
!     treating the corner elements first
!     
      ncoveredmelem=0
      do j=1,nmaface
         ifacem=maface(j)
         nelemm=int(ifacem/10.d0)
         jfacem=ifacem-10*nelemm
!
!        add master element to covered stack
!
         call nident(icoveredmelem,nelemm,ncoveredmelem,id)
         if((id.ne.0).and.(icoveredmelem(id).eq.nelemm)) then
!     master element was already treated
            cycle
         else
!     add master element to covered elements
            ncoveredmelem=ncoveredmelem+1
            do ii=ncoveredmelem,id+2,-1
               icoveredmelem(ii)=icoveredmelem(ii-1)
            enddo
            icoveredmelem(id+1)=nelemm
         endif
         call faceinfo(nelemm,jfacem,lakon,nopem,
     &        nopesm,idummy)     
!     
!     determining the nodes of the face
!    
         do j1=1,nopem
            konl(j1)=kon(ipkon(nelemm)+j1)
         enddo
         do k1=1,nopesm
            nodel=getnodel(k1,jfacem,nopem)
            nodem(k1)=konl(nodel)
            do j1=1,3
               xl2m(j1,k1)=co(j1,konl(nodel))+
     &              vold(j1,konl(nodel))
            enddo
         enddo 
         dd=dsqrt(xn(1)**2+xn(2)**2+xn(3)**2)
!     
!     divide master element into konvex subelements
!
!     npg: number of polygon nodes     
!     nodepg: node numbers of polygon nodes
!     xlpg: coordinates of polygon nodes
!     
         if((nopesm.eq.3).or.(nopesm.eq.4))then
!     
!     no loop needed
!     
            npg=nopesm
            do k2=1,nopesm
               nodepg(k2)=nodem(k2)
               do j2=1,3
                  xlpg(j2,k2)=xl2m(j2,k2)
               enddo
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!     
         elseif(nopesm.eq.6)then
!
!     tri6 surface is divided into 4 tri3 polygons
!     
!     1. triangle
!
            npg=3
            nodepg(1)=nodem(1)
            nodepg(2)=nodem(4)
            nodepg(3)=nodem(6)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,1)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,4)
            enddo
            do j2=1,3
               xlpg(j2,3)=xl2m(j2,6)
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!     
!     2. triangle
!     
            npg=3
            nodepg(1)=nodem(4)
            nodepg(2)=nodem(2)
            nodepg(3)=nodem(5)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,4)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,2)
            enddo
            do j2=1,3
               xlpg(j2,3)=xl2m(j2,5)
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!     
!     3. triangle
!     
            npg=3
            nodepg(1)=nodem(5)
            nodepg(2)=nodem(3)
            nodepg(3)=nodem(6)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,5)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,3)
            enddo
            do j2=1,3
               xlpg(j2,3)=xl2m(j2,6)
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!     
!     4. triangle
!     
            npg=3
            nodepg(1)=nodem(4)
            nodepg(2)=nodem(5)
            nodepg(3)=nodem(6)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,4)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,5)
            enddo
            do j2=1,3
               xlpg(j2,3)=xl2m(j2,6)
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!
         elseif(nopesm.eq.8)then
!     
!     quad8 surface is divided into 4 tri3 + 1 quad4 polygons
!     
!     1. triangle
!
            npg=3
            nodepg(1)=nodem(1)
            nodepg(2)=nodem(5)
            nodepg(3)=nodem(8)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,1)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,5)
            enddo
            do j2=1,3
               xlpg(j2,3)=xl2m(j2,8)
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!     
!     2. triangle
!     
            npg=3
            nodepg(1)=nodem(5)
            nodepg(2)=nodem(2)
            nodepg(3)=nodem(6)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,5)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,2)
            enddo
            do j2=1,3
               xlpg(j2,3)=xl2m(j2,6)
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!     
!     3. triangle
!     
            npg=3
            nodepg(1)=nodem(6)
            nodepg(2)=nodem(3)
            nodepg(3)=nodem(7)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,6)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,3)
            enddo
            do j2=1,3
               xlpg(j2,3)=xl2m(j2,7)
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!     
!     4. triangle
!     
            npg=3
            nodepg(1)=nodem(7)
            nodepg(2)=nodem(4)
            nodepg(3)=nodem(8)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,7)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,4)
            enddo
            do j2=1,3
               xlpg(j2,3)=xl2m(j2,8)
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!     
!     quad
!     
            npg=4
            nodepg(1)=nodem(5)
            nodepg(2)=nodem(6)
            nodepg(3)=nodem(7)
            nodepg(4)=nodem(8)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,5)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,6)
            enddo
            do j2=1,3
              xlpg(j2,3)=xl2m(j2,7)
           enddo
           do j2=1,3
              xlpg(j2,4)=xl2m(j2,8)
           enddo
           call treatmasterface_load(
     &          nopes,slavstraight,xn,xl2m,xl3sp,
     &          ipe,ime,iactiveline,nactiveline,
     &          ifacem,nintpoint,imastload,pmastload,
     &          xlpg,npg,nodepg,areaslav,nopem)
        endif
      enddo
!     
!     corners of the Slave surface have already been treated
!     
      do j=1,nopes
         mafacecorner(j,1)=0
      enddo
!     
!     retrieving all triangles by neighborhood search
!     
      do
         line=iactiveline(1,1)
         if(nactiveline.eq.0) exit
         if(koncont(4,ime(2,line)).eq.iactiveline(2,1)) then
            itri=imastop(ime(3,line),ime(2,line))
         else
            itri=ime(2,line)
         endif
!
!     check whether still in contact tie
!
         if((itri.gt.ncont).or.(itri.lt.0))then
            itri=0
         endif                            
         
         if(itri.eq.0) then
            nactiveline=nactiveline-1
            do il=1,nactiveline
               do k=1,2
                  iactiveline(k,il)=iactiveline(k,il+1)
               enddo
            enddo
            cycle
         endif
         
!     
         ifacem=koncont(4,itri)
         nelemm=int(koncont(4,itri)/10.d0)
         jfacem=koncont(4,itri)-10*nelemm     
!     
!     add master element to covered stack
!     
         call nident(icoveredmelem,nelemm,ncoveredmelem,id)
         if((id.gt.0).and.(icoveredmelem(id).eq.nelemm)) then
!
!     master element was already treated
!
            nactiveline=nactiveline-1
            do il=1,nactiveline
               do k=1,2
                  iactiveline(k,il)=iactiveline(k,il+1)
               enddo
            enddo
            cycle
         else
!
!     add master element to covered elements
!
            ncoveredmelem=ncoveredmelem+1
            do ii=ncoveredmelem,id+2,-1
               icoveredmelem(ii)=icoveredmelem(ii-1)
            enddo
            icoveredmelem(id+1)=nelemm
         endif 
         indexe=ipkon(nelemm)
         call faceinfo(nelemm,jfacem,lakon,nopem,
     &        nopesm,idummy)
!     
!     determining the nodes of the face
!     
         do j1=1,nopem
            konl(j1)=kon(ipkon(nelemm)+j1)
         enddo
         do k1=1,nopesm
            nodel=getnodel(k1,jfacem,nopem)
            nodem(k1)=konl(nodel)
            do j1=1,3
               xl2m(j1,k1)=co(j1,konl(nodel))+
     &              vold(j1,konl(nodel))
            enddo
         enddo
!     
!     divide master element into konvex subelements
!     
         if((nopesm.eq.3).or.(nopesm.eq.4))then
!     
!     no loop needed
!     
            npg=nopesm
            do k2=1,nopesm
               nodepg(k2)=nodem(k2)
               do j2=1,3
                  xlpg(j2,k2)=xl2m(j2,k2)
               enddo
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!     
         elseif(nopesm.eq.6)then
!
!     tri6 surface is divided into 4 tri3 polygons
!     
!     1. triangle
!
            npg=3
            nodepg(1)=nodem(1)
            nodepg(2)=nodem(4)
            nodepg(3)=nodem(6)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,1)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,4)
            enddo
            do j2=1,3
               xlpg(j2,3)=xl2m(j2,6)
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!     
!     2. triangle
!     
            npg=3
            nodepg(1)=nodem(4)
            nodepg(2)=nodem(2)
            nodepg(3)=nodem(5)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,4)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,2)
            enddo
            do j2=1,3
               xlpg(j2,3)=xl2m(j2,5)
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!     
!     3. triangle
!     
            npg=3
            nodepg(1)=nodem(5)
            nodepg(2)=nodem(3)
            nodepg(3)=nodem(6)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,5)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,3)
            enddo
            do j2=1,3
               xlpg(j2,3)=xl2m(j2,6)
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!     
!     4. triangle
!     
            npg=3
            nodepg(1)=nodem(4)
            nodepg(2)=nodem(5)
            nodepg(3)=nodem(6)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,4)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,5)
            enddo
            do j2=1,3
               xlpg(j2,3)=xl2m(j2,6)
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!
         elseif(nopesm.eq.8)then
!     
!     quad8 surface is divided into 4 tri3 + 1 quad4 polygons
!     
!     1. triangle
!
            npg=3
            nodepg(1)=nodem(1)
            nodepg(2)=nodem(5)
            nodepg(3)=nodem(8)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,1)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,5)
            enddo
            do j2=1,3
               xlpg(j2,3)=xl2m(j2,8)
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!     
!     2. triangle
!     
            npg=3
            nodepg(1)=nodem(5)
            nodepg(2)=nodem(2)
            nodepg(3)=nodem(6)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,5)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,2)
            enddo
            do j2=1,3
               xlpg(j2,3)=xl2m(j2,6)
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!     
!     3. triangle
!     
            npg=3
            nodepg(1)=nodem(6)
            nodepg(2)=nodem(3)
            nodepg(3)=nodem(7)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,6)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,3)
            enddo
            do j2=1,3
               xlpg(j2,3)=xl2m(j2,7)
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!     
!     4. triangle
!     
            npg=3
            nodepg(1)=nodem(7)
            nodepg(2)=nodem(4)
            nodepg(3)=nodem(8)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,7)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,4)
            enddo
            do j2=1,3
               xlpg(j2,3)=xl2m(j2,8)
            enddo
            call treatmasterface_load(
     &           nopes,slavstraight,xn,xl2m,xl3sp,
     &           ipe,ime,iactiveline,nactiveline,
     &           ifacem,nintpoint,imastload,pmastload,
     &           xlpg,npg,nodepg,areaslav,nopem)
!     
!     quad
!     
            npg=4
            nodepg(1)=nodem(5)
            nodepg(2)=nodem(6)
            nodepg(3)=nodem(7)
            nodepg(4)=nodem(8)
            do j2=1,3
               xlpg(j2,1)=xl2m(j2,5)
            enddo
            do j2=1,3
               xlpg(j2,2)=xl2m(j2,6)
            enddo
            do j2=1,3
              xlpg(j2,3)=xl2m(j2,7)
           enddo
           do j2=1,3
              xlpg(j2,4)=xl2m(j2,8)
           enddo
           call treatmasterface_load(
     &          nopes,slavstraight,xn,xl2m,xl3sp,
     &          ipe,ime,iactiveline,nactiveline,
     &          ifacem,nintpoint,imastload,pmastload,
     &          xlpg,npg,nodepg,areaslav,nopem)
        endif
      enddo
!     
      return
      end
