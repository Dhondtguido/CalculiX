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
      subroutine smoothingvertexnodes(inn,iponn,nktet,
     &     iexternnode,netet_,kontet,cotet,ipoeln,ieln,h,quality,
     &     jfix)
!
!     smoothing the mesh (cf. Triangulation de Delaunay et maillage,
!     P-L George and H. Borouchaki, section 8.4.1)      
!
      implicit none
!
      integer netet_,kontet(4,*),ielem,i,j,k,iter,index,indexe,node,
     &     iexternnode(*),iponn(*),inn(2,*),ipoeln(*),ieln(2,*),
     &     indexeold,ielemold,nktet,jfix(*)
!
      real*8 cotet(3,*),quality(*),sumalpha,relaxw,qmax,h(*),
     &     cpycotet(3),sumalphap(3),haverage,alphaj
!
!     relaxation factor
!
      relaxw=0.5d0
!
!     number of smoothing iterations
!
      iter=10
!      
      do k=1,iter
!
!     subsurface nodes
!        
        loop1: do i=1,nktet
          if(iponn(i).eq.0) cycle
          if(jfix(i).eq.1) cycle
          if(iexternnode(i).ne.0) cycle
          indexe=ipoeln(i)
!     
!     calculate the quality of the ball (all elements containing node i)
!     
          qmax=0.d0
          do
            if(indexe.eq.0) exit
            ielem=ieln(1,indexe)
            if(kontet(1,ielem).eq.0) then
              indexe=ieln(2,indexe)
              cycle
            endif
            if(quality(ielem).gt.qmax) qmax=quality(ielem)
            indexe=ieln(2,indexe)
          enddo
!     
!     initializing numerator (sumalphap) and denominator (sumalpha)
!     of the new position
!     
          sumalpha=0.d0
          do j=1,3
            sumalphap(j)=0.d0
          enddo
!     
!     calculating the new position: loop over all neighboring nodes
!     
          index=iponn(i)
          do
            if(index.eq.0) exit
!     
!     neighboring node
!     
            node=inn(1,index)
            haverage=(h(i)+h(node))/2.d0
!     
!     weight = inverse of the square of the desired edge length
!     
            alphaj=1.d0/(haverage*haverage)
            sumalpha=sumalpha+alphaj
            do j=1,3
              sumalphap(j)=sumalphap(j)+alphaj*cotet(j,node)
            enddo
            index=inn(2,index)
          enddo
!     
!     storing old position of node i
!     determining new position
!     
          do j=1,3
            cpycotet(j)=cotet(j,i)
            cotet(j,i)=(1.d0-relaxw)*cotet(j,i)+
     &           relaxw*sumalphap(j)/sumalpha
          enddo
!     
!     check the quality of the tetrahedral elements in the ball for
!     the new position of node i; as soon as a tetrahedron is detected
!     with a quality exceeding qmax the coordinates of node i are reverted
!     to the ones before smoothing
!     
          indexe=ipoeln(i)
          do
            if(indexe.eq.0) exit
            ielem=ieln(1,indexe)
            if(kontet(1,ielem).eq.0) then
              indexe=ieln(2,indexe)
              cycle
            endif
            call meshquality(netet_,kontet,cotet,quality,ielem)
            if(quality(ielem).gt.qmax) then
!     
!     1. restore the original coordinates
!     2. recalculate the quality of the tetrahedrons in the ball
!     covered so far in the current loop
!     
              do j=1,3
                cotet(j,i)=cpycotet(j)
              enddo
              indexeold=ipoeln(i)
              do
                ielemold=ieln(1,indexeold)
                if(kontet(1,ielemold).eq.0) then
                  indexeold=ieln(2,indexeold)
                  cycle
                endif
                call meshquality(netet_,kontet,cotet,quality,ielemold)
                if(ielemold.eq.ielem) exit
                indexeold=ieln(2,indexeold)
              enddo
              cycle loop1
            endif
            indexe=ieln(2,indexe)
          enddo
        enddo loop1
!
!     surface nodes
!        
        loop2: do i=1,nktet
          if(iponn(i).eq.0) cycle
          if(iexternnode(i).eq.0) cycle
          indexe=ipoeln(i)
!     
!     calculate the quality of the ball (all elements containing node i)
!     
          qmax=0.d0
          do
            if(indexe.eq.0) exit
            ielem=ieln(1,indexe)
            if(kontet(1,ielem).eq.0) then
              indexe=ieln(2,indexe)
              cycle
            endif
            if(quality(ielem).gt.qmax) qmax=quality(ielem)
            indexe=ieln(2,indexe)
          enddo
!     
!     initializing numerator (sumalphap) and denominator (sumalpha)
!     of the new position
!     
          sumalpha=0.d0
          do j=1,3
            sumalphap(j)=0.d0
          enddo
!     
!     calculating the new position: loop over all neighboring nodes
!     
          index=iponn(i)
          do
            if(index.eq.0) exit
!     
!     neighboring node
!     
            node=inn(1,index)
            haverage=(h(i)+h(node))/2.d0
!     
!     weight = inverse of the square of the desired edge length
!     
            alphaj=1.d0/(haverage*haverage)
            sumalpha=sumalpha+alphaj
            do j=1,3
              sumalphap(j)=sumalphap(j)+alphaj*cotet(j,node)
            enddo
            index=inn(2,index)
          enddo
!     
!     storing old position of node i
!     determining new position
!     
          do j=1,3
            cpycotet(j)=cotet(j,i)
            cotet(j,i)=(1.d0-relaxw)*cotet(j,i)+
     &           relaxw*sumalphap(j)/sumalpha
          enddo
!     
!     check the quality of the tetrahedral elements in the ball for
!     the new position of node i; as soon as a tetrahedron is detected
!     with a quality exceeding qmax the coordinates of node i are reverted
!     to the ones before smoothing
!     
          indexe=ipoeln(i)
          do
            if(indexe.eq.0) exit
            ielem=ieln(1,indexe)
            if(kontet(1,ielem).eq.0) then
              indexe=ieln(2,indexe)
              cycle
            endif
            call meshquality(netet_,kontet,cotet,quality,ielem)
            if(quality(ielem).gt.qmax) then
!     
!     1. restore the original coordinates
!     2. recalculate the quality of the tetrahedrons in the ball
!     covered so far in the current loop
!     
              do j=1,3
                cotet(j,i)=cpycotet(j)
              enddo
              indexeold=ipoeln(i)
              do
                ielemold=ieln(1,indexeold)
                if(kontet(1,ielemold).eq.0) then
                  indexeold=ieln(2,indexeold)
                  cycle
                endif
                call meshquality(netet_,kontet,cotet,quality,ielemold)
                if(ielemold.eq.ielem) exit
                indexeold=ieln(2,indexeold)
              enddo
              cycle loop2
            endif
            indexe=ieln(2,indexe)
          enddo
        enddo loop2
!        
      enddo
!            
      return
      end
