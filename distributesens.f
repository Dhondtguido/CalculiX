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
      subroutine distributesens(istartdesi,ialdesi,ipkon,lakon,
     &   ipoface,ndesi,nodedesi,nodface,kon,co,dgdx,nobject,
     &   weightformgrad,nodedesiinv,iregion,objectset,
     &   dgdxglob,nk,physcon,nobjectstart)
!
!     calculation of the relationship between facial distributions
!     and corresponding nodal values
!
!     transforming the nodal sensitivities into facially 
!     distributed sensitivities
!
      implicit none
!
      character*81 objectset(5,*)
      character*8 lakon(*)
!
      integer idesvar,j,k,kk,l,m,m1,n,istartdesi(*),
     &   ielem,ipoface(*),nodface(5,*),ndesi,nodedesi(*),nope,
     &   ialdesi(*),nopes,indexe,ifour,indexel,
     &   mint3d,mint2d,ipkon(*),konl(26),iflag,ifaceq(8,6),
     &   ifacet(6,4),ifacew1(4,5),ifacew2(8,5),kon(*),
     &   nodes1(8),ifacel,iobject,nodes2(8),indexs,nk,
     &   ithree,i,node,ifacq(2,3,20),ifact(2,3,10),ifacw(2,3,15),
     &   nopedesi,nnodes,nodedesiinv(*),nobject,iregion,iactnode,
     &   nobjectstart
!
      real*8 xi,et,weight,xl(3,9),xs(3,2),xsj(3),shp(7,9),
     &   co(3,*),xsjj,weightformgrad(ndesi),dgdx(ndesi,nobject),
     &   dgdxglob(2,nk,nobject),physcon(*)
!
!
!
!     flag for shape functions
!
      data iflag /3/
      data indexel /0/
!
      save indexel
!
      include "gauss.f"
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
!     ifacq indicates for each local node number to which
!     faces this node belongs and the internal number
!     within this face. e.g., local node number 1 belongs to
!     face 1 with internal number 4, face 3 with internal number
!     1 and face 6 with internal number 2 (first line underneath)
!
!     ifacq, ifact and ifacw can be derived in a unique way from
!     ifaceq, ifacet and ifacew2
!
      data ifacq /1,4,3,1,6,2,
     &1,3,3,2,4,1,
     &1,2,4,2,5,1,
     &1,1,5,2,6,1,
     &2,1,3,4,6,3,
     &2,2,3,3,4,4,
     &2,3,4,3,5,4,
     &2,4,5,3,6,4,
     &1,7,3,5,0,0,
     &1,6,4,5,0,0,
     &1,5,5,5,0,0,
     &1,8,6,5,0,0,
     &2,5,3,7,0,0,
     &2,6,4,7,0,0,
     &2,7,5,7,0,0,
     &2,8,6,7,0,0,
     &3,8,6,6,0,0,
     &3,6,4,8,0,0,
     &4,6,5,8,0,0,
     &5,6,6,8,0,0/
!
      data ifact /1,1,2,1,4,1,
     &1,3,2,2,3,1,
     &1,2,3,2,4,3,
     &2,3,3,3,4,2,
     &1,6,2,4,0,0,
     &1,5,3,4,0,0,
     &1,4,4,6,0,0,
     &2,6,4,4,0,0,
     &2,5,3,6,0,0,
     &3,5,4,5,0,0/
!
      data ifacw /1,1,3,1,5,2,
     &1,3,3,2,4,1,
     &1,2,4,2,5,1,
     &2,1,3,4,5,3,
     &2,2,3,3,4,4,
     &2,3,4,3,5,4,
     &1,6,3,5,0,0,
     &1,5,4,5,0,0,
     &1,4,5,5,0,0,
     &2,4,3,7,0,0,
     &2,5,4,7,0,0,
     &2,6,5,7,0,0,
     &3,8,5,6,0,0,
     &3,6,4,8,0,0,
     &4,6,5,8,0,0/
!
!     flag for shape functions
!
      ifour=4
      ithree=3
!
!
!     Loop over all designvariables
!
      do idesvar=1,ndesi
!     
         node=nodedesi(idesvar)
!     
!     Loop over all elements belonging to the designvariable
!     
         do j=istartdesi(idesvar),istartdesi(idesvar+1)-1
            ielem=ialdesi(j)
            
            indexe=ipkon(ielem)
!     
!     mint2d: # of integration points on the surface
!     nopes: # of nodes in the surface
!     nope: # of nodes in the element
!           
            if(lakon(ielem)(4:5).eq.'8R') then
               mint2d=1
               nopes=4
               nope=8
               nopedesi=3
            elseif(lakon(ielem)(4:4).eq.'8') then
               mint2d=4
               nopes=4
               nope=8
               nopedesi=3
            elseif(lakon(ielem)(4:5).eq.'20') then
               mint2d=9
               nopes=8
               nope=20
               nopedesi=5
            elseif(lakon(ielem)(4:5).eq.'10') then
               mint2d=3
               nopes=6
               nope=10
c               nopedesi=3
               nopedesi=4
            elseif(lakon(ielem)(4:4).eq.'4') then
               mint2d=1
               nopes=3
               nope=4
               nopedesi=3
            elseif(lakon(ielem)(4:4).eq.'6') then                
               mint2d=1
               nopes=3
               nope=6
               nopedesi=3
            elseif(lakon(ielem)(4:5).eq.'15') then
               mint2d=3
               nopes=6
               nope=15
c               nopedesi=3
            else
               exit
            endif
            if(iregion.eq.0) nopedesi=0
!     
!     actual position of the nodes belonging to the
!     master surface
!     
            do l=1,nope
               konl(l)=kon(indexe+l)
            enddo
!     
            do i=1,nope
               if(kon(indexe+i).eq.node) exit
            enddo
!     
!     Loop over all faces of the actual element
!     
            loop1: do m=1,3
!     
!     hexahedral elements
!     
               if((nope.eq.20).or.(nope.eq.8)) then
                  k=ifacq(1,m,i)
                  m1=ifacq(2,m,i)
                  do l=1,nopes
                     nodes1(l)=konl(ifaceq(l,k))
                     nodes2(l)=nodes1(l)
                  enddo
                  call insertsorti(nodes1,ifour)
c                  call isortii(nodes1,iaux,ifour,kflag)
!     
!     tetrahedral element
!     
               elseif((nope.eq.10).or.(nope.eq.4)) then
                  k=ifact(1,m,i)
                  m1=ifact(2,m,i)
                  do l=1,nopes
                     nodes1(l)=konl(ifacet(l,k))
                     nodes2(l)=nodes1(l)
                  enddo
                  call insertsorti(nodes1,ithree)
c                  call isortii(nodes1,iaux,ithree,kflag)
                  nodes1(4)=0
!     
!     wedge element
!     
               elseif(nope.eq.15) then
                  k=ifacw(1,m,i)
                  m1=ifacw(2,m,i)
                  if(k.le.2) then
                     nopes=6
                     nopedesi=4
                     mint3d=3
                     do l=1,nopes
                        nodes1(l)=konl(ifacew2(l,k))
                        nodes2(l)=nodes1(l)
                     enddo  
                     call insertsorti(nodes1,ithree)
c                     call isortii(nodes1,iaux,ithree,kflag)
                     nodes1(4)=0
                  else
                     nopes=8
                     nopedesi=5
                     mint3d=4
                     do l=1,nopes
                        nodes1(l)=konl(ifacew2(l,k))
                        nodes2(l)=nodes1(l)
                     enddo  
                     call insertsorti(nodes1,ithree)
c                     call isortii(nodes1,iaux,ithree,kflag)
                  endif        
               else 
                  k=ifacw(1,m,i)
                  m1=ifacw(2,m,i)
                  if(k.le.2) then
                     nopes=3
                     mint3d=1
                     do l=1,nopes
                        nodes1(l)=konl(ifacew1(l,k))
                        nodes2(l)=nodes1(l)
                     enddo  
                     call insertsorti(nodes1,ithree)
c                     call isortii(nodes1,iaux,ithree,kflag)
                     nodes1(4)=0
                  else
                     nopes=4
                     mint3d=2
                     do l=1,nopes
                        nodes1(l)=konl(ifacew1(l,k))
                        nodes2(l)=nodes1(l)
                     enddo  
                     call insertsorti(nodes1,ithree)
c                     call isortii(nodes1,iaux,ithree,kflag)
                  endif        
               endif
               if(k.eq.0) exit
!     
!              Find external element face
!     
               indexs=ipoface(nodes1(1))
               do
                  if(indexs.eq.0) cycle loop1
                  if((nodface(1,indexs).eq.nodes1(2)).and.
     &                 (nodface(2,indexs).eq.nodes1(3))) then
!     
!                    Check if sufficient designnodes on that surface
!     
                     nnodes=0
                     do n=1,nopes
                        if(nodedesiinv(nodes1(n)).eq.1) then
                           nnodes=nnodes+1
                        endif
                     enddo
!     
!                    Assign coordinates to nodes of surface
!     
c                     write(*,*) 'formgradient',nnodes,nopedesi
                     if(nnodes.ge.nopedesi) then
                        if((nope.eq.20).or.(nope.eq.8)) then
                           do l=1,nopes
                              do n=1,3
                                 ifacel=ifaceq(l,nodface(4,indexs))
                                 xl(n,l)=co(n,konl(ifacel))
                              enddo
                           enddo
                        elseif((nope.eq.10).or.(nope.eq.4)) then
                           do l=1,nopes
                              do n=1,3
                                 ifacel=ifacet(l,nodface(4,indexs))
                                 xl(n,l)=co(n,konl(ifacel))
                              enddo
                           enddo
                        elseif(nope.eq.15) then 
                           do l=1,nopes
                              do n=1,3
                                 ifacel=ifacew2(l,nodface(4,indexs))
                                 xl(n,l)=co(n,konl(ifacel))
                              enddo
                           enddo
                        elseif(nope.eq.6) then
                           do l=1,nopes
                              do n=1,3
                                 ifacel=ifacew1(l,nodface(4,indexs))
                                 xl(n,l)=co(n,konl(ifacel))
                              enddo
                           enddo
                        endif           
!     
!                       Determine weighting factors
!     
                        do kk=1,mint2d   
                           if(lakon(ielem)(4:5).eq.'20') then
                              xi=gauss2d3(1,kk)
                              et=gauss2d3(2,kk)
                              weight=weight2d3(kk)
                              call shape8q(xi,et,xl,xsj,
     &                             xs,shp,iflag)
                           elseif(lakon(ielem)(4:5).eq.'10') then
                              xi=gauss2d5(1,kk)
                              et=gauss2d5(2,kk)
                              weight=weight2d5(kk)
                              call shape6tri(xi,et,xl,xsj,
     &                             xs,shp,iflag)
                           elseif(lakon(ielem)(4:4).eq.'4') then
                              xi=gauss2d4(1,kk)
                              et=gauss2d4(2,kk)
                              weight=weight2d4(kk)
                              call shape3tri(xi,et,xl,xsj,
     &                             xs,shp,iflag)
                           elseif(lakon(ielem)(4:5).eq.'15') then
                              if(k.le.2) then
                                 xi=gauss2d5(1,kk)
                                 et=gauss2d5(2,kk)
                                 weight=weight2d5(kk)
                                 call shape6tri(xi,et,xl,xsj,
     &                                xs,shp,iflag)
                              else
                                 xi=gauss2d3(1,kk)
                                 et=gauss2d3(2,kk)
                                 weight=weight2d3(kk)
                                 call shape8q(xi,et,xl,xsj,
     &                                xs,shp,iflag)
                              endif
                           elseif(lakon(ielem)(4:5).eq.'8R') then
                              xi=gauss2d1(1,kk)
                              et=gauss2d1(2,kk)
                              weight=weight2d1(kk)
                              call shape4q(xi,et,xl,xsj,
     &                             xs,shp,iflag)
                           elseif(lakon(ielem)(4:4).eq.'8') then
                              xi=gauss2d2(1,kk)
                              et=gauss2d2(2,kk)
                              weight=weight2d2(kk)
                              call shape4q(xi,et,xl,xsj,
     &                             xs,shp,iflag)
                           elseif(lakon(ielem)(4:4).eq.'6') then
                              if(k.le.2) then
                                 xi=gauss2d4(1,kk)
                                 et=gauss2d4(2,kk)
                                 weight=weight2d4(kk)
                                 call shape3tri(xi,et,xl,xsj,
     &                                xs,shp,iflag)
                               else
                                 xi=gauss2d2(1,kk)
                                 et=gauss2d2(2,kk)
                                 weight=weight2d2(kk)
                                 call shape4q(xi,et,xl,xsj,
     &                                xs,shp,iflag)
                              endif               
                           endif
!     
!                          Calculate Jacobian determinant
!     
                           xsjj=dsqrt(xsj(1)**2+xsj(2)**2+
     &                          xsj(3)**2)
!     
!                          Evaluate Shape functions for weightmatrix
!     
                           weightformgrad(idesvar)=weightformgrad
     &                          (idesvar)+weight*shp(4,m1)*xsjj 
!                           write(*,*) 'formgradient ',idesvar,j,m,kk,
!     &                             weightformgrad(idesvar)
                        enddo
                     endif
                     exit
                  endif
                  indexs=nodface(5,indexs)
               enddo            ! find external element faces
!     
            enddo loop1         ! Loop over all faces of the actual element
!     
         enddo                  ! Loop over all elements belonging to the designvariable
!         
      enddo                     ! Loop over all designvariables
!     
!     Scaling of sensitivities
!     
      do idesvar=1,ndesi
         do iobject=1+nobjectstart,nobject
            if(objectset(1,iobject)(4:13).eq.'MEMBERSIZE') cycle
            if(objectset(1,iobject)(1:9).eq.'MAXGROWTH') cycle
            if(objectset(1,iobject)(1:12).eq.'MAXSHRINKAGE') cycle
            if(objectset(1,iobject)(1:9).eq.'PACKAGING') cycle
            if(weightformgrad(idesvar).gt.0.d0) then
               dgdx(idesvar,iobject)=dgdx(idesvar,iobject)
     &             /weightformgrad(idesvar)
               iactnode=nodedesi(idesvar)
               dgdxglob(1,iactnode,iobject)=dgdx(idesvar,iobject)
            else
               iactnode=nodedesi(idesvar)
               dgdxglob(1,iactnode,iobject)=dgdx(idesvar,iobject)
            endif
         enddo
      enddo
!      
      return
      end
      
