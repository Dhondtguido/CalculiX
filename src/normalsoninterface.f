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
      subroutine normalsoninterface(istartset,iendset,
     &  ialset,imast,ipkon,kon,lakon,imastnode,nmastnode,
     &  xmastnor,co)
!
      logical exist
!
      character*8 lakon(*)
!
      integer istartset(*),iendset(*),ialset(*),imast,j,ifacem,
     &  nelemm,jfacem,indexe,ipkon(*),kon(*),nopem,node,
     &  ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),
     &  nmastnode,id,indexnode(9),konl(26),imastnode(*)
!
      real*8 xmastnor(3,*),xsj2(3),shp2(7,9),xs2(3,2),xi,et,dd,
     &  xquad(2,9),xtri(2,7),xl2m(3,9),co(3,*)
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
!     new added data for the local coodinates for nodes
!
      data xquad /-1.d0,-1.d0,
     &             1.d0,-1.d0,
     &             1.d0,1.d0,
     &            -1.d0,1.d0,
     &             0.d0,-1.d0,
     &             1.d0,0.d0,
     &             0.d0,1.d0,
     &            -1.d0,0.d0,
     &             0.d0,0.d0/
!
      data xtri /0.d0,0.d0,
     &           1.d0,0.d0,
     &           0.d0,1.d0,
     &           .5d0,0.d0,
     &           .5d0,.5d0,
     &           0.d0,.5d0,
     &           0.333333333333333d0,0.333333333333333d0/
!
      data iflag /2/
!
!     creates the normals for the interface condition A.n=0 
!     for electromagnetic calculations
!
      nmastnode=0
!
      do j=istartset(imast),iendset(imast)
!     
!        create imastnode, and nmastnode (number of nodes in the
!        master surface
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
     &           (lakon(nelemm)(4:4).eq.'8')) then
               node=kon(indexe+ifaceq(l,jfacem))
            elseif((lakon(nelemm)(4:4).eq.'4').or.
     &                    (lakon(nelemm)(4:5).eq.'10')) then
               node=kon(indexe+ifacet(l,jfacem))
            elseif(lakon(nelemm)(4:4).eq.'6') then
               node=kon(indexe+ifacew1(l,jfacem))
            elseif(lakon(nelemm)(4:5).eq.'15') then
               node=kon(indexe+ifacew2(l,jfacem))
            endif
            call nident(imastnode,node,
     &           nmastnode,id)
            exist=.FALSE.
            if(id.gt.0) then
               if(imastnode(id).eq.node) then
                  exist=.TRUE.
               endif
            endif
            if(exist) cycle
            nmastnode=nmastnode+1
            do k=nmastnode,id+2,-1
               imastnode(k)=imastnode(k-1)
            enddo
            imastnode(id+1)=node
         enddo
!     
      enddo
!     
!     calculate the mean normal in the master nodes
!     
      do j=istartset(imast),iendset(imast)
!     
         ifacem=ialset(j)
         nelemm=int(ifacem/10)
         jfacem=ifacem - nelemm*10
         indexe=ipkon(nelemm)
!     
!     nopem: # of nodes in the master face
!     nope: # of nodes in the element
!     
         if(lakon(nelemm)(4:5).eq.'8R') then
            nopem=4
            nope=8
         elseif(lakon(nelemm)(4:4).eq.'8') then
            nopem=4
            nope=8
         elseif(lakon(nelemm)(4:5).eq.'20') then
            nopem=8
            nope=20
         elseif(lakon(nelemm)(4:5).eq.'10') then
            nopem=6
            nope=10
         elseif(lakon(nelemm)(4:4).eq.'4') then
            nopem=3
            nope=4
!     
!     treatment of wedge faces
!     
         elseif(lakon(nelemm)(4:4).eq.'6') then
            nope=6
            if(jfacem.le.2) then
               nopem=3
            else
               nopem=4
            endif
         elseif(lakon(nelemm)(4:5).eq.'15') then
            nope=15
            if(jfacem.le.2) then
               nopem=6
            else
               nopem=8
            endif
         endif
!     
!     actual position of the nodes belonging to the
!     master surface
!     
         do k=1,nope
            konl(k)=kon(ipkon(nelemm)+k)
         enddo
!     
         if((nope.eq.20).or.(nope.eq.8)) then
            do m=1,nopem
               do k=1,3
                  xl2m(k,m)=co(k,konl(ifaceq(m,jfacem)))
               enddo
            enddo
         elseif((nope.eq.10).or.(nope.eq.4)) 
     &           then
            do m=1,nopem
               do k=1,3
                  xl2m(k,m)=co(k,konl(ifacet(m,jfacem)))
               enddo
            enddo
         elseif(nope.eq.15) then
            do m=1,nopem
               do k=1,3
                  xl2m(k,m)=co(k,konl(ifacew2(m,jfacem)))
               enddo
            enddo
         else
            do m=1,nopem
               do k=1,3
                  xl2m(k,m)=co(k,konl(ifacew1(m,jfacem)))
               enddo
            enddo
         endif
         
!     calculate the normal vector in the nodes belonging to the master surface
!     
         if(nopem.eq.8) then
            do m=1,nopem
               xi=xquad(1,m)
               et=xquad(2,m)
               call shape8q(xi,et,xl2m,xsj2,xs2,shp2,iflag)
               dd=dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2)
     &              + xsj2(3)*xsj2(3))
               xsj2(1)=xsj2(1)/dd
               xsj2(2)=xsj2(2)/dd
               xsj2(3)=xsj2(3)/dd
!     
               if(nope.eq.20) then
                  node=konl(ifaceq(m,jfacem))
               elseif(nope.eq.15) then
                  node=konl(ifacew2(m,jfacem))
               endif
!     
               call nident(imastnode,node,
     &              nmastnode,id)
               index1=id
               indexnode(m)=index1
               xmastnor(1,index1)=xmastnor(1,index1)
     &              +xsj2(1)
               xmastnor(2,index1)=xmastnor(2,index1)
     &              +xsj2(2)
               xmastnor(3,index1)=xmastnor(3,index1)
     &              +xsj2(3)
            enddo
         elseif(nopem.eq.4) then
            do m=1,nopem
               xi=xquad(1,m)
               et=xquad(2,m)
               call shape4q(xi,et,xl2m,xsj2,xs2,shp2,iflag)
               dd=dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2) 
     &              + xsj2(3)*xsj2(3))
               xsj2(1)=xsj2(1)/dd
               xsj2(2)=xsj2(2)/dd
               xsj2(3)=xsj2(3)/dd
!     
               if(nope.eq.8) then
                  node=konl(ifaceq(m,jfacem))
               elseif(nope.eq.6) then
                  node=konl(ifacew1(m,jfacem))
               endif
!     
               call nident(imastnode,node,
     &              nmastnode,id)
!     
               index1=id
               indexnode(m)=index1
               xmastnor(1,index1)=xmastnor(1,index1)
     &              +xsj2(1)
               xmastnor(2,index1)=xmastnor(2,index1)
     &              +xsj2(2)
               xmastnor(3,index1)=xmastnor(3,index1)
     &              +xsj2(3)
            enddo
         elseif(nopem.eq.6) then
            do m=1,nopem
               xi=xtri(1,m)
               et=xtri(2,m)
               call shape6tri(xi,et,xl2m,xsj2,xs2,shp2,iflag)
               dd=dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2) 
     &              + xsj2(3)*xsj2(3))
               xsj2(1)=xsj2(1)/dd
               xsj2(2)=xsj2(2)/dd
               xsj2(3)=xsj2(3)/dd
!     
               if(nope.eq.10) then
                  node=konl(ifacet(m,jfacem))
               elseif(nope.eq.15) then
                  node=konl(ifacew2(m,jfacem))
               endif
!     
               call nident(imastnode,node,
     &              nmastnode,id)
               index1=id
               indexnode(m)=index1
               xmastnor(1,index1)=xmastnor(1,index1)
     &              +xsj2(1)
               xmastnor(2,index1)=xmastnor(2,index1)
     &              +xsj2(2)
               xmastnor(3,index1)=xmastnor(3,index1)
     &              +xsj2(3)
            enddo
         else
            do m=1,nopem
               xi=xtri(1,m)
               et=xtri(2,m)
               call shape3tri(xi,et,xl2m,xsj2,xs2,shp2,iflag)
               dd=dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2) 
     &              + xsj2(3)*xsj2(3))
               xsj2(1)=xsj2(1)/dd
               xsj2(2)=xsj2(2)/dd
               xsj2(3)=xsj2(3)/dd
!     
               if(nope.eq.6) then
                  node=konl(ifacew1(m,jfacem))
               elseif(nope.eq.4) then
                  node=konl(ifacet(m,jfacem))
               endif
!     
               call nident(imastnode,node,
     &              nmastnode,id)
               index1=id
               indexnode(m)=index1
               xmastnor(1,id)=xmastnor(1,index1)
     &              +xsj2(1)
               xmastnor(2,id)=xmastnor(2,index1)
     &              +xsj2(2)
               xmastnor(3,id)=xmastnor(3,index1)
     &              +xsj2(3)
            enddo
         endif
      enddo
!     
!     normalizing the normals
!     
      do l=1,nmastnode
         dd=dsqrt(xmastnor(1,l)**2+xmastnor(2,l)**2+
     &        xmastnor(3,l)**2)
         do m=1,3
            xmastnor(m,l)=xmastnor(m,l)/dd
         enddo
      enddo    
!     
      return
      end
