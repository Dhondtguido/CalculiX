!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
!     Calculating the slave-master triangulation
!     and generating the integration points needed for the calculation of the coupling matrices
!     for details see phd-thesis Sitzmann Appendix A 
!     
!     [out]    nintpoint   number of generated integration points 
!     [in]     koncont     (1:3,i) nodes of triagle_i (4,i) element face
!     [in,out] islavsurf   islavsurf(1,i) slaveface i islavsurf(2,i) # integration points generated before looking at face i
!     [out]    imastsurf   pointer into pmastsurf  
!     [out]    pmastsurf   field storing position and weights for integration points on master side 
!     [in]     islavnode   fields containing nodes of slace surfaces
!     [in]     nslavnode   (i) for contraint i pointer into field islavnode
!     [in]     imastop     (l,i) for edge l in triagle i neightbouring triangle
!     [out]    gapmints    stores gaps between master and slave side
!     [in,out] islavact	 (i) indicates, if slave node i is active (=-3 no-slave-node, =-2 no-LM-node, =-1 no-gap-node, =0 inactive node, =1 sticky node, =2 slipping/active node) 
!     [in]     ncont	 number of triangles in triagulation of master side
!     [in]     ipe	 (i) pointer to ime for node i 
!     [in]     ime         ... cataloging the edges with node i
!     [out]    pslavsurf   field storing position xil, etal and weight for integration point on slave side
!     [in]     i           current tie
!     [in]     l           current face
!     [in]     ntri        # triangles
!     
      subroutine slavintmortar(ntie,itietri,ipkon,kon,lakon,straight,
     &     nintpoint,koncont,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,
     &     iinc,islavsurf,imastsurf,pmastsurf,
     &     islavnode,nslavnode,imastop,gapmints,
     &     islavact,mi,ncont,ipe,ime,pslavsurf,i,l,ntri,tietol,
     &     reltime,nmethod)
!     
!     Authors: Samoela Rakotonanahary, Saskia Sitzmann
!
!     Treating slave face l in field islavsurf    
!     
      implicit none
!     
      logical debug,shrink
!     
      character*8 lakon(*)
!     
      integer ntie,nintpoint,imastop(3,*),ncont,
     &     itietri(2,ntie),ipkon(*),kon(*),koncont(4,*),node,
     &     neigh(10),iflag,kneigh,i,j,k,l,ii,itri,nx(*),ny(*),
     &     nz(*),nmethod,nelemm,jfacem,indexe,iinc,nnodelem,nope,
     &     islavsurf(2,*),islavnode(*),nslavnode(ntie+1),
     &     imastsurf(*),ifaces,nelems,jfaces,mi(*),
     &     m,nopes,konl(26),id,islavact(*),
     &     imface(8),ntria,imfacecorner(8,8),line,
     &     iactiveline(3,3*ncont),icoveredmelem(3*ncont),
     &     nactiveline,ipe(*),ime(4,*),k1,j1,ncoveredmelem,
     &     ntri,nintpfirst,nodem(8),nodem2(8),ithree,
     &     il,ifac,getlocno,ifacem,idummy,nopemm,nmp,k2,j2
!     
      real*8 straight(16,*),co(3,*),vold(0:mi(2),*),
     &     xo(*),yo(*),zo(*),x(*),y(*),z(*),
     &     pmastsurf(2,*),xl2m(3,8),xl2s(3,8),xl2m2(3,8),
     &     pmiddle(3),xl2sr(3,8),xl2sp(3,8),xl2s2(3,8),
     &     dd,xns(3,8),areaslav,al,xn(3),gapmints(*),slavstraight(36),
     &     dist,distmin,tietol(4,*),clearance,reltime,
     &     pslavsurf(3,*),err,xquad(2,8),cl2s(3,8),cl2m(3,8),
     &     xtri(2,6),xi,et,xsj2(3),xs2(3,7),shp2(7,8),anglesm
!     
      include "gauss.f"
!     
      data iflag /2/
      data ithree /3/
!     
      data xquad /-1,-1,
     &     1,-1,
     &     1,1,
     &     -1,1,
     &     0,-1,
     &     1,0,
     &     0,1,
     &     -1,0/
!     
      data xtri /0,0,
     &     1,0,
     &     0,1,
     &     0.5,0,
     &     0.5,0.5,
     &     0,0.5/
!     
      debug=.false.
!     
      kneigh=1
      err=0.1
      areaslav=0.0
      nintpfirst=nintpoint
      islavsurf(2,l)=nintpoint
!     
      if(debug)WRITE(30,*) '#SLAVINTMORTAR iinc',iinc,'face',l
      if(debug)WRITE(20,*) '#SLAVINTMORTAR iinc',iinc,'face',l      
!     
!     get clearance and shrink
!     
      clearance=tietol(3,i)
      shrink=.true.
!     
!     Research of the contact integration points
!     
      ifaces=islavsurf(1,l)
      nelems=int(ifaces/10)
      jfaces=ifaces - nelems*10
!     
!     get nope,nopes
!     
      call getnumberofnodes(nelems,jfaces,lakon,nope,nopes,idummy)
!     
!     actual position of the nodes belonging to the
!     slave surface
!     
      do j=1,nope
        konl(j)=kon(ipkon(nelems)+j)
      enddo
!
!     xl2s: deformed coordinates
!     cl2s: undeformed coordinates
!
      do m=1,nopes
        do j=1,3
          ifac=getlocno(m,jfaces,nope)
          xl2s(j,m)=co(j,konl(ifac))+
     &         vold(j,konl(ifac))
          cl2s(j,m)=co(j,konl(ifac))   
        enddo
      enddo  
!     
!     slightly reducing the size of the slave surface;
!     needed for contact search
!
!     xl2sr: reduced deformed coordinates
!     
      do j=1,3
        pmiddle(j)=0.d0
        do m=1,nopes
          pmiddle(j)=pmiddle(j)+xl2s(j,m)
        enddo
        pmiddle(j)=pmiddle(j)/nopes
      enddo
      do j=1,3
        do m=1,nopes
          xl2sr(j,m)=xl2s(j,m)-0.05*err*(xl2s(j,m)-pmiddle(j))
        enddo
      enddo
!     
!     re-sort vertices for quadratic elements
!     e.g. instead of 1-5-2-6-3-7-4-8 in ccw-direction
!     1-2-3-4-5-6-7-8
!
!     re-sorted deformed coordinates: xl2s2      
!     
      if(nopes.eq.3 .or. nopes.eq.4)then
        do j=1,nopes
          do k=1,3
            xl2s2(k,j)=xl2s(k,j)
          enddo
        enddo
      else
        do j=1,int(nopes/2.0)
          do k=1,3
            xl2s2(k,2*j-1)=xl2s(k,j)
            xl2s2(k,2*j)=xl2s(k,(int(nopes/2.0))+j)           
          enddo
        enddo
      endif
!     
!     calculate the mean normal vector on the Slave Surface
!     +
!     determine the equations of the triangle/quadrilateral
!     (mean)plane and of the planes boardering the 
!     triangle/quadrilateral
!     
      if(nopes.eq.3) then
        call straighteq3d(xl2s2,slavstraight)
        do k=1,3
          xn(k)=slavstraight(4*nopes+k)
        enddo               
      else
        do k=1,3
          xn(k)=0.d0
        enddo
!     
        do m=1,nopes
          if(nopes.eq.4 .or. nopes.eq.8)then
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
     &         +xsj2(3)*xsj2(3))
          xsj2(1)=xsj2(1)/dd
          xsj2(2)=xsj2(2)/dd
          xsj2(3)=xsj2(3)/dd
!     
          do k=1,3
            xn(k)=xn(k)
     &           +xsj2(k)
          enddo
        enddo 
!     
!     normalizing the mean normal on the Slave surface
!     
        dd=dsqrt(xn(1)**2+xn(2)**2+xn(3)**2)
        do k=1,3
          xn(k)=xn(k)/dd
        enddo           
        call approxplane(xl2s2,slavstraight,xn,nopes)
      endif
!     
!     Project slave nodes to meanplane, needed for Sutherland-Hodgman
!
!     projected re-sorted deformed coordinates in xl2sp
!
      do j=1,nopes
        al=-xn(1)*xl2s2(1,j)-xn(2)*
     &       xl2s2(2,j)-xn(3)*xl2s2(3,j)-
     &       slavstraight(nopes*4+4)
        if(nopes.eq.3)then
          do k=1,3
            xl2sp(k,j)=xl2s2(k,j)
          enddo
        else
          do k=1,3
            xl2sp(k,j)=xl2s2(k,j)+al*xn(k)
          enddo
        endif
      enddo
!     
!     determine the triangles corresponding to the corner
!     nodes
!     
      ntria=0
      do j=1,8
        imface(j)=0
        do k=1,8
          imfacecorner(j,k)=0
        enddo
      enddo
!     
!     check for nogap-nodes: these are slave nodes for which there
!     is no opposite master face
!     
      do j=1,nopes
!     
!     here unreduced coordinates are used
!     
        call neartriangle(xl2s(1,j),xn,xo,yo,zo,x,y,z,nx,ny,nz,
     &       ntri,neigh,kneigh,itietri,ntie,straight,imastop,itri,i)
        ifac=getlocno(j,jfaces,nope)
        node=konl(ifac)
        if(debug) then
          write(20,*) neigh(1),neigh(1)+itietri(1,i)-1
          write(20,*) 'itri',itri,'node',node
        endif
        call nident(islavnode(nslavnode(i)+1),node,
     &       nslavnode(i+1)-nslavnode(i),id) 
        if(itri.eq.0) then  
!     
!     no triangle found (no opposite master face; value -1)
!     node was possibly active in last increment
!     
          islavact(nslavnode(i)+id)=-1
        else
!     
!     triangle found; slave node is by default set to inactive (value 0)
!     
          if(debug)write(20,*) ' node',node,
     &         islavnode(nslavnode(i)+id),
     &         id,islavact(nslavnode(i)+id) 
          if(islavact(nslavnode(i)+id).lt.0) then  
            islavact(nslavnode(i)+id)=0
          endif
          if(debug)write(20,*) ' node',node,
     &         islavnode(nslavnode(i)+id),
     &         id,islavact(nslavnode(i)+id) 
        endif
      enddo
!
!     nogap check finished
!
!     looking for the master faces covering the slave nodes
!     and storing them in imface(1.....max-8) in ascending order
!     if for slave node j the opposite master face is in imface(k),
!     imfacecorner(j,k)=1 and imfacecorner(j,l)=0,l.ne.k. One and the
!     same master face can cover several nodes of the same slave face      
!
!     increasing the reduction from 0.05*err to 2*err
!      
      do j=1,3
        do m=1,nopes
          xl2sr(j,m)=xl2s(j,m)-2*err*(xl2s(j,m)-pmiddle(j))
        enddo
      enddo
      distmin=1.1   
      do j=1,nopes
!     
!     for initial master elements reduced coordinates are used
!     
        call neartriangle(xl2sr(1,j),xn,xo,yo,zo,x,y,z,nx,ny,nz,
     &       ntri,neigh,kneigh,itietri,ntie,straight,imastop,itri,i)
        ifac=getlocno(j,jfaces,nope)
        node=konl(ifac) 
        if(itri.eq.0) then  
          cycle
        endif
        dist=-(straight(13,itri)*xl2sr(1,j)+
     &       straight(14,itri)*xl2sr(2,j)+
     &       straight(15,itri)*xl2sr(3,j)+
     &       straight(16,itri))/
     &       (straight(13,itri)*xn(1)+
     &       straight(14,itri)*xn(2)+
     &       straight(15,itri)*xn(3))
        if(dist.lt.distmin)distmin=dist
        if(debug)write(20,*) 'j',j,'dist',dist,distmin
        ifacem=koncont(4,itri)
        if(debug)write(20,*)'noder ',node,'itri',itri,
     &       'ifacem',ifacem
!     
        call nident(imface,ifacem,ntria,id)
        if(id.gt.0) then
          if(imface(id).eq.ifacem) then
            imfacecorner(j,id)=1
            cycle
          endif
        endif
!     
!     triangle was not covered yet: add to stack
!     
!     angle criteria: the master triangle should not deviate
!     more than 45 degrees from the mean slave face
!     
        anglesm=xn(1)*straight(13,itri)
     &       +xn(2)*straight(14,itri)
     &       +xn(3)*straight(15,itri)
        if(debug)write(20,*)'cos alpa',anglesm
!     
!     check angle to avoid problem with the dual gap         
!     
        if(anglesm.lt.-0.7)then
!     
!     angle between surface normals between 135 and 225 degree;
!     this corresponds to an          
!     angle between surfaces between 0 and 45 degree
!     
          ntria=ntria+1
          do k=ntria,id+2,-1
            imface(k)=imface(k-1)
            do m=1,j-1
              imfacecorner(m,k)=imfacecorner(m,k-1)
            enddo
          enddo
          imface(id+1)=ifacem
          imfacecorner(j,id+1)=1
          do m=1,j-1
            imfacecorner(m,id+1)=0
          enddo 
        endif              
      enddo
      if(debug)then 
        write(20,*)'ifacem n1:n8 face',l,ntria
        do k=1,ntria
          write(20,*) imface (k),imfacecorner(1:nopes,k)
        enddo
      endif
!     
      nactiveline=0
!     
!     treating the corner elements first
!     
      ncoveredmelem=0
      do j=1,ntria
        if(debug)write(20,*) 'corner triangle j',j
        ifacem=imface(j)
        nelemm=int(ifacem/10.d0)
        jfacem=ifacem-10*nelemm
!     
!     add master element to covered stack
!     
        call nident(icoveredmelem,nelemm,ncoveredmelem,id)
c     if(id.ne.0 .and. icoveredmelem(id).eq.nelemm)then
        if(id.gt.0) then
!     
!     master elment was already threated
!     
          if(icoveredmelem(id).eq.nelemm) cycle
        endif
!     
!     add master element to covered elements
!     
        ncoveredmelem=ncoveredmelem+1
        do ii=ncoveredmelem,id+2,-1
          icoveredmelem(ii)=icoveredmelem(ii-1)
        enddo
        icoveredmelem(id+1)=nelemm
!        
        if(debug)write(20,*)itri,imface(j),ifacem,nelemm,jfacem
        call getnumberofnodes(nelemm,jfacem,lakon,nopemm,
     &       nnodelem,idummy)     
!     
!     determining the nodes of the face
!     
        do j1=1,nopemm
          konl(j1)=kon(ipkon(nelemm)+j1)
        enddo
!
!     deformed and undeformed coordinates of the master face
!        
        do k1=1,nnodelem
          ifac=getlocno(k1,jfacem,nopemm)
          nodem(k1)=konl(ifac)
          do j1=1,3
            xl2m(j1,k1)=co(j1,konl(ifac))+
     &           vold(j1,konl(ifac))
            cl2m(j1,k1)=co(j1,konl(ifac))
          enddo
        enddo
!        
!     guido: xn already normalized?
!        
        dd=dsqrt(xn(1)**2+xn(2)**2+xn(3)**2)
        if(debug)then
          write(20,*) 'dd',dd    
          write(20,100)(xn(k),k=1,3)
          write(20,100)(slavstraight(nopes*4+k),k=1,3)
          write(20,*) 'SIM: xl2s'
          do j1=1,nopes
            write(20,*)(xl2s(k,j1),k=1,3)
          enddo    
          write(20,*) 'SIM: xl2m'
          do j1=1,nnodelem
            write(20,*)(xl2m(k,j1),k=1,3)
          enddo
        endif                                        
 100    format('SIM: xns',3(3x,e15.8))
!     
        if(debug) write(20,*) 'TT: itri',nelemm
!     
!     divide master element into konvex subelements
!     
        if(nnodelem.eq.3 .or.nnodelem.eq.4)then
!     
!     no loop needed
!     
          nmp=nnodelem
          do k2=1,nnodelem
            nodem2(k2)=nodem(k2)
            do j2=1,3
              xl2m2(j2,k2)=xl2m(j2,k2)
            enddo
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
        elseif(nnodelem.eq.6)then
!     
!     tri6 surface is divided into 4 tri3 surfaces
!     
!     1. triangle
c     write(20,*)'*********tria 1**************'
          nmp=3
          nodem2(1)=nodem(1)
          nodem2(2)=nodem(4)
          nodem2(3)=nodem(6)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,1)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,4)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,6)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
!     2. triangle
!     
c     write(20,*)'*********tria 2**************'
          nmp=3
          nodem2(1)=nodem(4)
          nodem2(2)=nodem(2)
          nodem2(3)=nodem(5)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,4)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,2)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,5)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
!     3. triangle
!     
c     write(20,*)'*********tria 3**************'
          nmp=3
          nodem2(1)=nodem(5)
          nodem2(2)=nodem(3)
          nodem2(3)=nodem(6)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,5)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,3)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,6)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
!     4. triangle
!     
c     write(20,*)'*********tria 4**************'
          nmp=3
          nodem2(1)=nodem(4)
          nodem2(2)=nodem(5)
          nodem2(3)=nodem(6)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,4)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,5)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,6)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
        elseif(nnodelem.eq.8)then
!     
!     quad8 surface is divided into 4 tri3 + 1 quad4 surfaces
!     
!     1. triangle
          nmp=3
          nodem2(1)=nodem(1)
          nodem2(2)=nodem(5)
          nodem2(3)=nodem(8)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,1)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,5)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,8)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
!     2. triangle
!     
          nmp=3
          nodem2(1)=nodem(5)
          nodem2(2)=nodem(2)
          nodem2(3)=nodem(6)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,5)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,2)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,6)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
!     3. triangle
!     
          nmp=3
          nodem2(1)=nodem(6)
          nodem2(2)=nodem(3)
          nodem2(3)=nodem(7)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,6)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,3)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,7)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
!     4. triangle
!     
          nmp=3
          nodem2(1)=nodem(7)
          nodem2(2)=nodem(4)
          nodem2(3)=nodem(8)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,7)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,4)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,8)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
!     quad
!     
          nmp=4
          nodem2(1)=nodem(5)
          nodem2(2)=nodem(6)
          nodem2(3)=nodem(7)
          nodem2(4)=nodem(8)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,5)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,6)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,7)
          enddo
          do j2=1,3
            xl2m2(j2,4)=xl2m(j2,8)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
        endif
      enddo
!     
!     corners of the Slave surface have already been treated
!     
      do j=1,nopes
        imfacecorner(j,1)=0
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
        if(itri.gt.itietri(2,i) .or. itri.lt.itietri(1,i))then
          itri=0
        endif                            
!     
        if(itri.eq.0) then
          nactiveline=nactiveline-1
          do il=1,nactiveline
            do k=1,3
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
c     if(id .gt. 0 .and. icoveredmelem(id).eq.nelemm)then
        if(id.gt.0) then
          if(icoveredmelem(id).eq.nelemm) then
!     
!     master elment was already threated
!     
            nactiveline=nactiveline-1
            do il=1,nactiveline
              do k=1,3
                iactiveline(k,il)=iactiveline(k,il+1)
              enddo
            enddo
            cycle
          endif
        endif
!     
!     add master element to covered elements
!     
        ncoveredmelem=ncoveredmelem+1
        do ii=ncoveredmelem,id+2,-1
          icoveredmelem(ii)=icoveredmelem(ii-1)
        enddo
        icoveredmelem(id+1)=nelemm
!     
        indexe=ipkon(nelemm)
!        
        call getnumberofnodes(nelemm,jfacem,lakon,nopemm,
     &       nnodelem,idummy)
        
!     
!     determining the nodes of the face
!     
        do j1=1,nopemm
          konl(j1)=kon(ipkon(nelemm)+j1)
        enddo
        do k1=1,nnodelem
          ifac=getlocno(k1,jfacem,nopemm)
          nodem(k1)=konl(ifac)
          do j1=1,3
            xl2m(j1,k1)=co(j1,konl(ifac))+
     &           vold(j1,konl(ifac))
            cl2m(j1,k1)=co(j1,konl(ifac))
          enddo
        enddo
        if(debug)then
          write(20,*) 'dd',dd ,'ifacem',ifacem   
          write(20,100)(xn(k),k=1,3)
          write(20,*) 'SIM: xl2s'
          do j1=1,nopes
            write(20,*)(xl2s(k,j1),k=1,3)
          enddo    
          write(20,*) 'SIM: xl2m'
          do j1=1,nnodelem
            write(20,*)(xl2m(k,j1),k=1,3)
          enddo
        endif  
!     
!     divide master element into konvex subelements
!     
        if(nnodelem.eq.3 .or.nnodelem.eq.4)then
!     
!     no loop needed
!     
          nmp=nnodelem
          do k2=1,nnodelem
            nodem2(k2)=nodem(k2)
            do j2=1,3
              xl2m2(j2,k2)=xl2m(j2,k2)
            enddo
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
        elseif(nnodelem.eq.6)then
!     
!     tri6 surface is divided into 4 tri3 surfaces
!     
!     1. triangle
c     write(20,*)'*********tria 1**************'
          nmp=3
          nodem2(1)=nodem(1)
          nodem2(2)=nodem(4)
          nodem2(3)=nodem(6)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,1)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,4)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,6)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
!     2. triangle
!     
c     write(20,*)'*********tria 2**************'
          nmp=3
          nodem2(1)=nodem(4)
          nodem2(2)=nodem(2)
          nodem2(3)=nodem(5)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,4)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,2)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,5)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
!     3. triangle
!     
c     write(20,*)'*********tria 3**************'
          nmp=3
          nodem2(1)=nodem(5)
          nodem2(2)=nodem(3)
          nodem2(3)=nodem(6)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,5)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,3)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,6)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
!     4. triangle
!     
c     write(20,*)'*********tria 4**************'
          nmp=3
          nodem2(1)=nodem(4)
          nodem2(2)=nodem(5)
          nodem2(3)=nodem(6)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,4)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,5)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,6)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
        elseif(nnodelem.eq.8)then
!     
!     quad8 surface is divided into 4 tri3 + 1 quad4 surfaces
!     
!     1. triangle
          nmp=3
          nodem2(1)=nodem(1)
          nodem2(2)=nodem(5)
          nodem2(3)=nodem(8)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,1)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,5)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,8)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
!     2. triangle
!     
          nmp=3
          nodem2(1)=nodem(5)
          nodem2(2)=nodem(2)
          nodem2(3)=nodem(6)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,5)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,2)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,6)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
!     3. triangle
!     
          nmp=3
          nodem2(1)=nodem(6)
          nodem2(2)=nodem(3)
          nodem2(3)=nodem(7)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,6)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,3)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,7)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
!     4. triangle
!     
          nmp=3
          nodem2(1)=nodem(7)
          nodem2(2)=nodem(4)
          nodem2(3)=nodem(8)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,7)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,4)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,8)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
!     
!     quad
!     
          nmp=4
          nodem2(1)=nodem(5)
          nodem2(2)=nodem(6)
          nodem2(3)=nodem(7)
          nodem2(4)=nodem(8)
          do j2=1,3
            xl2m2(j2,1)=xl2m(j2,5)
          enddo
          do j2=1,3
            xl2m2(j2,2)=xl2m(j2,6)
          enddo
          do j2=1,3
            xl2m2(j2,3)=xl2m(j2,7)
          enddo
          do j2=1,3
            xl2m2(j2,4)=xl2m(j2,8)
          enddo
          call treatmasterface_mortar(
     &         nopes,slavstraight,xn,xns,xl2s,xl2sp,
     &         ipe,ime,iactiveline,nactiveline,
     &         ifacem,
     &         nintpoint,pslavsurf,imastsurf,pmastsurf,
     &         xl2m,nnodelem,xl2m2,nmp,nodem2,
     &         gapmints,l,areaslav,debug,clearance,
     &         cl2s,cl2m,shrink,reltime)
        endif
      enddo
      islavsurf(2,l+1)=nintpoint
      if(debug) then
        if(areaslav.lt.1.e-12)write(*,*)'areaslav(',l,')=',
     &       areaslav
      endif
!     
      if(debug)then
        write(20,*) 'mint2d',(nintpoint-islavsurf(2,l))
        write(20,*)'areaslav(',l,')=',areaslav
        write(20,*)'*********************************'
      endif
!     
c            do j=1,nintpoint
c              write(*,*) 'slavintmortar', j,imastsurf(j)
c              write(*,*) 'slavintmortar', (pmastsurf(j1,j),j1=1,2)
c              write(*,*) 'slavintmortar', (pslavsurf(j1,j),j1=1,3)
c              write(*,*) 'slavintmortar', gapmints(j)
c            enddo
      return
      end
