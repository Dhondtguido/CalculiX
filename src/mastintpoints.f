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
     &     imastop,mi,ncont,ipe,ime,nelemload,sideload,nload,nload_,
     &     imastload,pmastload,nelemloadloc,sideloadloc,
     &     xload,iamload,nam)
!
!     Based on slavintpoints.f
!
      implicit none
!
      character*8 lakon(*)
      character*20 sideload(*),sideloadloc
!     
      integer nintpoint,imastop(3,*),ncont,nelemload(2,*),nload,
     &     ipkon(*),kon(*),koncont(4,*),neigh(10),iflag,kneigh,i,
     &     j,k,ii,itri,nx(*),ny(*),nz(*),nelemm,jfacem,nload_,
     &     indexe,nopesm,ifaces,nelems,jfaces,mi(*),
     &     m,nopes,konl(20),id,maface(8),nmaface,imastload(2,*),
     &     mafacecorner(8,8),line,iactiveline(2,3*ncont),
     &     icoveredmelem(3*ncont),nactiveline,ipe(*),ime(4,*),k1,j1,
     &     ncoveredmelem,nodem(8),nodepg(8),nelemloadloc,
     &     il,nodel,getnodel,ifacem,idummy,nopem,npg,k2,j2,
     &     iamload(2,*),nam
!     
      real*8 straight(16,*),co(3,*),vold(0:mi(2),*),pmastload(3,*),
     &     xo(*),yo(*),zo(*),x(*),y(*),z(*),xload(2,*),p12(3),
     &     xl2m(3,8),xl2s(3,8),xlpg(3,8),p23(3),p31(3),
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
      nelems=nelemloadloc
      if((ipkon(nelems).ge.0).or.(lakon(nelems)(1:1).ne.'S')) then
        write(*,*) '*WARNING in mastintpoints'
        write(*,*) '         element ',nelems,' with interface loading'
        write(*,*) '         is not a shell element with the attribute'
        write(*,*) '         DLOAD'
        return
      endif
      read(sideloadloc(2:2),'(i1)') jfaces
      ifaces=10*nelems+jfaces
!     
!     get nopes
!     
      read(lakon(nelems)(2:2),'(i1)') nopes
!     
!     actual position of the nodes belonging to the
!     slave surface
!     
      do j=1,nopes
         konl(j)=kon(-2-ipkon(nelems)+j)
      enddo
!     
      do m=1,nopes
         do j=1,3
            xl2s(j,m)=co(j,m)+
     &           vold(j,konl(m))     
         enddo
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
!     calculating a middle point in the face
!     
      do j=1,3
         pmiddle(j)=0.d0
         do m=1,nopes
            pmiddle(j)=pmiddle(j)+xl2s(j,m)
         enddo
         pmiddle(j)=pmiddle(j)/nopes
      enddo
!     
!     calculate the mean normal vector xn on the slave face
!     +
!     determine the equations of the slave face
!     (mean)plane and of the planes perpendicular to its 
!     piecewise linear approximation 
!     
      if(nopes.eq.3) then
        do i=1,3
          p12(i)=xl3s(i,2)-xl3s(i,1)
          p23(i)=xl3s(i,3)-xl3s(i,2)
          p31(i)=xl3s(i,1)-xl3s(i,3)
        enddo
!     
!     normalized vector normal to the triangle: xn = p12 x p23
!     
        xn(1)=p12(2)*p23(3)-p12(3)*p23(2)
        xn(2)=p12(3)*p23(1)-p12(1)*p23(3)
        xn(3)=p12(1)*p23(2)-p12(2)*p23(1)
        dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
        do i=1,3
          xn(i)=xn(i)/dd
        enddo
        slavstraight(16)=-xn(1)*xl3s(1,1)-xn(2)*xl3s(2,1)-
     &       xn(3)*xl3s(3,1)
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
               xn(k)=xn(k)+xsj2(k)
            enddo
         enddo 
!     
!     normalizing the mean normal on the slave surface
!     
         dd=dsqrt(xn(1)**2+xn(2)**2+xn(3)**2)
         do k=1,3
            xn(k)=xn(k)/dd
         enddo           
c         call approxplane(xl3s,slavstraight,xn,nopes)
         slavstraight(nopes*4+4)=
     &        -xn(1)*pmiddle(1)-xn(2)*pmiddle(2)-xn(3)*pmiddle(3)
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
!     moving the nodes towards the middle     
!     
      do j=1,3
         do m=1,nopes
            xl2sr(j,m)=xl2s(j,m)-2*err*(xl2s(j,m)-pmiddle(j))
         enddo
      enddo
      do j=1,nopes
         call neartriangle_load(xl2sr(1,j),xn,xo,yo,zo,x,y,z,nx,ny,nz,
     &        neigh,kneigh,ncont,straight,imastop,itri)
c         node=konl(j) 
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &          xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &          xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &     xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &           xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &           xload,nload,nload_,iamload,nam,ifaces)
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
     &          xlpg,npg,nodepg,areaslav,nopem,nelemload,sideload,
     &          xload,nload,nload_,iamload,nam,ifaces)
        endif
      enddo
!     
      return
      end
