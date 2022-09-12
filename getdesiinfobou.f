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
      subroutine getdesiinfobou(ndesibou,nodedesibou,nodedesiinv,
     &  lakon,ipkon,kon,ipoface,nodface,nodedesiinvbou,ndesi,
     &  nodedesi,nk)
!
!     storing the design variables on the boundary in nodedesibou
!
!     create the field nodedesiinvbou which contains
!     the following information:
!     entry > 0: node contains to the set ndesibou
!     number of entry: specifies the position of the node 
!            in the set nodedesiinvbou
!
      implicit none
!
      character*8 lakon(*)
!
      character*81 setname
!
      integer ndesibou,node,nodedesibou(*),i,k,m,index,
     &  nodedesiinv(*),nelem,nope,nopedesi,ipkon(*),nnodes,
     &  kon(*),konl(26),iaux,kflag,ipoface(*),nodface(5,*),jfacem,
     &  nopesurf(9),ifaceq(8,6),ifacet(6,4),ifacew1(4,5),
     &  ifacew2(8,5),nopem,nodedesiinvbou(*),ndesi,nodedesi(*),nk
!
      setname(1:1)=' '
      ndesibou=0
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
!     A design node is considered to be at the border if the external
!     does not contain as many designvariables as nodes on the surface
! 
      do i=1,nk
         node=i
         if(ipoface(node).eq.0) cycle
         index=ipoface(node)
         do
            nelem=nodface(3,index)
            jfacem=nodface(4,index)
!     
            if(lakon(nelem)(4:4).eq.'8') then
               nope=8
               nopedesi=4
               nopem=4
            elseif(lakon(nelem)(4:5).eq.'20') then
               nope=20
               nopedesi=8
               nopem=8
            elseif(lakon(nelem)(4:5).eq.'10') then
               nope=10
               nopedesi=4
               nopem=6
            elseif(lakon(nelem)(4:4).eq.'4') then
               nope=4
               nopedesi=3
               nopem=3
            elseif(lakon(nelem)(4:4).eq.'6') then
               nope=6
               if(jfacem.le.2) then
                  nopem=3
                  nopedesi=3
               else
                  nopem=4
                  nopedesi=4
               endif
            elseif(lakon(nelem)(4:5).eq.'15') then
               nope=15
               if(jfacem.le.2) then
                  nopem=6
                  nopedesi=6
               else
                  nopem=8
                  nopedesi=8
               endif
            endif
!     
!     actual position of the nodes belonging to the
!     master surface
!     
            do k=1,nope
               konl(k)=kon(ipkon(nelem)+k)
            enddo
!     
            if((nope.eq.20).or.(nope.eq.8)) then
               do m=1,nopem
                  nopesurf(m)=konl(ifaceq(m,jfacem))
               enddo
            elseif((nope.eq.10).or.(nope.eq.4)) then
               do m=1,nopem
                  nopesurf(m)=konl(ifacet(m,jfacem))
               enddo
            elseif(nope.eq.15) then
               do m=1,nopem
                  nopesurf(m)=konl(ifacew2(m,jfacem))
               enddo
            else
               do m=1,nopem
                  nopesurf(m)=konl(ifacew1(m,jfacem))
               enddo
            endif
!    
!     sum up how many designvariables are on that surface
!
            nnodes=0
            do m=1,nopem
               if(nodedesiinv(nopesurf(m)).ne.0) then
                  nnodes=nnodes+1
               endif
            enddo
!
            if(nnodes.lt.nopedesi) then
               do m=1,nopem
                  if(nodedesiinv(nopesurf(m)).eq.1) then
                     nodedesiinv(nopesurf(m))=-1
                     ndesibou=ndesibou+1
                     nodedesibou(ndesibou)=nopesurf(m)
                     nodedesiinvbou(nopesurf(m))=1
                  endif
               enddo
            endif
            index=nodface(5,index)
            if(index.eq.0) exit      
         enddo
      enddo
!    
!     sort the boundary nodes
!      
      kflag=1
      call isortii(nodedesibou,iaux,ndesibou,kflag)
!
      do i=1,ndesibou
         index=nodedesibou(i)
         nodedesiinvbou(index)=i
         if(nodedesiinv(index).eq.-1) then
            nodedesiinv(index)=1
         endif 
      enddo
!
      return
      end

