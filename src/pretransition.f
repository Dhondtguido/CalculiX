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
      subroutine pretransition(ipkon,kon,lakon,co,nk,ipoface,
     &      nodface,nodedesiinv,xo,yo,zo,x,y,z,nx,ny,nz,ifree)
!
      implicit none
!
      character*8 lakon(*)
!
      integer j,nelem,jface,indexe,ipkon(*),kon(*),nopem,
     &  ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),
     &  konl(26),ipoface(*),nodface(5,*),nodedesiinv(*),
     &  nopesurf(9),ifree,nope,m,k,nk,nx(*),ny(*),nz(*),
     &  actnode,kflag,ndesinode
!
      real*8 co(3,*),xo(*),yo(*),zo(*),x(*),y(*),z(*)
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
!
!     all external nodes which are not designvariables are put in a set
!    
      ifree=0 
      do j=1,nk
!        
         if(ipoface(j).eq.0) cycle
         indexe=ipoface(j)
!        
         do
!
            nelem=nodface(3,indexe)
            jface=nodface(4,indexe)
!
!           nopem: # of nodes in the surface
!           nope: # of nodes in the element
!     
            if(lakon(nelem)(4:4).eq.'8') then
               nopem=4
               nope=8
            elseif(lakon(nelem)(4:5).eq.'20') then
               nopem=8
               nope=20
            elseif(lakon(nelem)(4:5).eq.'10') then
               nopem=6
               nope=10
            elseif(lakon(nelem)(4:4).eq.'4') then
               nopem=3
               nope=4
            elseif(lakon(nelem)(4:4).eq.'6') then
               nope=6
               if(jface.le.2) then
                  nopem=3
               else
                  nopem=4
               endif
            elseif(lakon(nelem)(4:5).eq.'15') then
               nope=15
               if(jface.le.2) then
                  nopem=6
               else
                  nopem=8
               endif
            endif
!     
!     get node numbers of actual surface
!     
            do k=1,nope
               konl(k)=kon(ipkon(nelem)+k)
            enddo
!     
            ndesinode=0
            if((nope.eq.20).or.(nope.eq.8)) then
               do m=1,nopem
                  actnode=konl(ifaceq(m,jface))
                  if(nodedesiinv(actnode).ne.1) then
                     ndesinode=ndesinode+1
                     nopesurf(ndesinode)=actnode         
                  endif
               enddo
            elseif((nope.eq.10).or.(nope.eq.4)) then
               do m=1,nopem
                  actnode=konl(ifacet(m,jface))
                  if(nodedesiinv(actnode).ne.1) then
                     ndesinode=ndesinode+1
                     nopesurf(ndesinode)=actnode         
                  endif
               enddo
            elseif(nope.eq.15) then
               do m=1,nopem
                  actnode=konl(ifacew2(m,jface))
                  if(nodedesiinv(actnode).ne.1) then
                     ndesinode=ndesinode+1
                     nopesurf(ndesinode)=actnode         
                  endif
               enddo
            else
               do m=1,nopem
                  actnode=konl(ifacew1(m,jface))
                  if(nodedesiinv(actnode).ne.1) then
                     ndesinode=ndesinode+1
                     nopesurf(ndesinode)=actnode         
                  endif
               enddo
            endif
!    
!     creation of the sets
!     
            if((ndesinode.gt.0).and.(ndesinode.lt.nopem)) then 
               do m=1,ndesinode
                  actnode=nopesurf(m)
                  if(nodedesiinv(actnode).eq.-1) cycle
                  ifree=ifree+1
                  xo(ifree)=co(1,actnode)
                  x(ifree)=xo(ifree)
                  nx(ifree)=ifree
                  yo(ifree)=co(2,actnode)
                  y(ifree)=yo(ifree)          
                  ny(ifree)=ifree
                  zo(ifree)=co(3,actnode)
                  z(ifree)=zo(ifree)
                  nz(ifree)=ifree
                  nodedesiinv(actnode)=-1
               enddo
            endif          
            indexe=nodface(5,indexe)
            if(indexe.eq.0) exit           
         enddo
      enddo
!
!     Correction of nodedesiinv
!
      do m=1,nk
         if(nodedesiinv(m).eq.-1) then
            nodedesiinv(m)=0
         endif
      enddo
!
!     Sorting of node set w.r.t. coordinates
!      
      kflag=2
      call dsort(x,nx,ifree,kflag)
      call dsort(y,ny,ifree,kflag)
      call dsort(z,nz,ifree,kflag)
!
      return
      end
