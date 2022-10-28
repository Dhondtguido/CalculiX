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
      subroutine smoothingmidnodes(cotet,ipoed,kontet,iedtet,iedgmid,
     &     ipoeled,ieled,qualityjac,iponn,inn,h,iexternedg,netet_,
     &     nktet_)
!
!     smoothing a quadratic tetrahedral mesh by moving the midnodes      
!
      implicit none
!
      integer iedge,ipoed(*),iedtet(6,*),kontet(4,*),iedgmid(*),
     &     ipoeled(*),ieled(2,*),iter,iponn(*),inn(2,*),indexe,ielemold,
     &     indexeold,i,j,k,node,ielem,index,iexternedg(*),neigh,netet_,
     &     nktet_
!
      real*8 cotet(3,*),qualityjac(*),relaxw,qmax,sumalpha,cpycotet(3),
     &     haverage,alphaj,h(*),sumalphap(3)
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
        loop1: do i=1,nktet_
!
          iedge=ipoed(i)
!
          if(iedge.eq.0) cycle
          if(iexternedg(iedge).ne.0) cycle
!
          node=iedgmid(iedge)
!
          if(iponn(node).eq.0) cycle
!
!         calculate the quality of the shell (all elements containing
!         the node), it is the worst quality of all elements belonging
!         to the shell (a shell is the set of all elements containing
!         a specific edge
!
          qmax=0.d0
          indexe=ipoeled(iedge)
!
          do
            if(indexe.eq.0) exit
            ielem=ieled(1,indexe)
!
            if(qualityjac(ielem).gt.qmax) qmax=qualityjac(ielem)
            indexe=ieled(2,indexe)
          enddo
!
!         initializing the numerator (sumalphap) and denominator
!         (sumalpha) of the new position
!
          sumalpha=0.d0
          do j=1,3
            sumalphap(j)=0.d0
          enddo
!
!         calculating the new position: loop over all neighboring nodes
!
          index=iponn(node)
          do
            if(index.eq.0) exit
!
!           neighboring node
!
            neigh=inn(1,index)
            haverage=(h(node)+h(neigh))/2.d0
!
!           weight = inverse of the square of the desired edge length
!
            alphaj=1.d0/(haverage*haverage)
            sumalpha=sumalpha+alphaj
            do j=1,3
              sumalphap(j)=sumalphap(j)+alphaj*cotet(j,neigh)
            enddo
            index=inn(2,index)
          enddo
!
!         storing the old position of the node; determining the new
!         position
!
          do j=1,3
            cpycotet(j)=cotet(j,node)
            cotet(j,node)=(1.d0-relaxw)*cotet(j,node)+
     &           relaxw*sumalphap(j)/sumalpha
          enddo
!
!         check the quality of the tetrahedral elements in the shell for
!         the new position of node; as soon as a tetrahedron is detected
!         with a quality exceeding qmax the coordinates of node are reverted
!         to the ones before smoothing
!
          indexe=ipoeled(iedge)
          do
            if(indexe.eq.0) exit
            ielem=ieled(1,indexe)
!
            call quadmeshquality(netet_,cotet,kontet,iedtet,
     &           iedgmid,qualityjac,ielem)
!
            if(qualityjac(ielem).gt.qmax) then
!
!             1. restore the original coordinates
!             2. recalculate the quality of the tetrahedra in the shell 
!                covered so far in the current loop
!
              do j=1,3
                cotet(j,node)=cpycotet(j)
              enddo
              indexeold=ipoeled(iedge)
              do
                ielemold=ieled(1,indexeold)
!
                call quadmeshquality(netet_,cotet,kontet,iedtet,
     &               iedgmid,qualityjac,ielemold)
                if(ielemold.eq.ielem) exit
                indexeold=ieled(2,indexeold)
c                write(*,*) '*INFO in smoothingmidnodes: reverting to'
c                write(*,*) '      old coordinate for node ',node
c                write(*,*) '      and element ',ielem
!
              enddo
              cycle loop1
            endif
            indexe=ieled(2,indexe)
          enddo
        enddo loop1
!        
      enddo
!            
      return
      end
