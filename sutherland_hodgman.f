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
!     Sutherland-Hodgman-Algo for polygon clipping combined with 
!     active line search
!
      subroutine sutherland_hodgman(nopes,xn,xl2sp,xl2mp,
     &  nodem,ipe,ime,iactiveline,nactiveline,ifreeintersec,
     &  nelemm,nnodelem,nvertex,pvertex)
!     
      implicit none 
!     
      logical invert,oldactive,altered,border
!
      integer nvertex,nopes,ipe(*),ime(4,*),iactiveline(3,*),
     &  nactiveline,ifreeintersec,itri,nelemm,i,ii,j,k,nnodelem,id,
     &  nodem(*),ncvertex,node1,node2,modf,node,indexl,ithree, 
     &  insertl(3),ninsertl
!
      real*8 pvertex(3,*),xn(3),xl2sp(3,*),pa(3),pb(3),xinters(3),
     &  xcp(3),diff,dd,xl2mp(3,*),c_pvertex(3,13),t,cedge(3),xtest(3),
     &     eplane,area,areax,areay,areaz,p1(3),p2(3),areaface
!
!
!     
      data ithree /3/
!     
      nvertex=0
      ninsertl=0
      border=.false.
!     
!     Initialize Polygon
!     
      do j=1,nopes
         nvertex=nvertex+1
         do k=1,3
            pvertex(k,nvertex)=xl2sp(k,j)
         enddo
      enddo
        areaface=0.d0
        do k=1,nvertex-2
         p1(1)=pvertex(1,k+1)-pvertex(1,1)
         p1(2)=pvertex(2,k+1)-pvertex(2,1)
         p1(3)=pvertex(3,k+1)-pvertex(3,1)
         p2(1)=pvertex(1,k+2)-pvertex(1,1)
         p2(2)=pvertex(2,k+2)-pvertex(2,1)
         p2(3)=pvertex(3,k+2)-pvertex(3,1)
         areax=((p1(2)*p2(3))-(p2(2)*p1(3)))**2
         areay=(-(p1(1)*p2(3))+(p2(1)*p1(3)))**2
         areaz=((p1(1)*p2(2))-(p2(1)*p1(2)))**2
         areaface=areaface+dsqrt(areax+areay+areaz)/2.d0
       enddo
!     
!     loop over clipping edges
!     
      do i=1,(nnodelem)
         ncvertex=0
         altered=.false.
!     
!     generate clipping plane
!     
         node1=nodem(modf(nnodelem,i))
         node2=nodem(modf(nnodelem,i+1))
         invert=.false.
         if(node2.lt.node1)then
            node=node1
            node1=node2
            node2=node
            invert=.true.
         endif
         indexl=ipe(node1)
         do
            if(ime(1,indexl).eq.node2) exit
            indexl=ime(4,indexl)
            if(indexl.eq.0) then
               write(*,*) 
     &        '*ERROR in SH:line was not properly catalogued'
               write(*,*) 'itri',itri,'node1',node1,'node2',node2
               call exit(201)
            endif 
         enddo
         do k=1,3
            cedge(k)=xl2mp(k,modf(nnodelem,i+1))
     &           -xl2mp(k,modf(nnodelem,i))
         enddo
         xcp(1)=xn(2)*cedge(3)-xn(3)*cedge(2)
         xcp(2)=xn(3)*cedge(1)-xn(1)*cedge(3)
         xcp(3)=xn(1)*cedge(2)-xn(2)*cedge(1)
         dd=dsqrt(xcp(1)**2+xcp(2)**2+xcp(3)**2)
         do k=1,3
            xcp(k)=xcp(k)/dd
         enddo
         t=-eplane(xl2mp(1:3,modf(nnodelem,i)),xcp,0.d0)
!     
!     inside-outside-test 
!     
         do k=1,3 
            xtest(k)=xl2mp(k,modf(nnodelem,i+2))
         enddo
         if (eplane(xtest,xcp,t).gt.0)then
            t=-t
            do k=1,3
               xcp(k)=-xcp(k)
            enddo
         endif
         oldactive=.false.
         call nidentk(iactiveline,indexl, nactiveline,id,ithree)
         if(id.gt.0.and.iactiveline(1,id).eq.indexl)then
            oldactive=.true.
         endif    
         if(oldactive)then
            nactiveline=nactiveline-1
            do ii=id,nactiveline
               do k=1,3
                  iactiveline(k,ii)=iactiveline(k,ii+1)
               enddo
            enddo 
         endif
!     
         if(nvertex.lt.3) cycle
!     
!     loop over polygon vertices
!     
         do j=0, (nvertex-1)
            do k=1,3
               pa(k)=pvertex(k,modf(nvertex,j))
               pb(k)=pvertex(k,modf(nvertex,j+1))
            enddo 
            if(eplane(pa,xcp,t).le.1.0d-12) then
               if(eplane(pb,xcp,t).le.1.0d-12)then
                  ncvertex=ncvertex+1
                  do k=1,3
                     c_pvertex(k,ncvertex)=pb(k)
                  enddo 
               else
                  if(abs(eplane(pa,xcp,t)).gt.1.0d-10)then
                     call intersectionpoint(pa,pb,xcp,t,xinters)
                     diff=(xinters(1)-pa(1))**2+(xinters(2)-pa(2))**2+
     &                    (xinters(3)-pa(3))**2
                     diff=dsqrt(diff)
                     if(diff.gt.1.d-11)then
                        ncvertex=ncvertex+1
                        do k=1,3
                           c_pvertex(k,ncvertex)=xinters(k)
                        enddo
                     endif
                  endif 
!     
                  if((.not.oldactive).and.(.not.altered))then
                     altered=.true.
                     nactiveline=nactiveline+1
                     ifreeintersec=ifreeintersec+1
                     do ii=nactiveline,id+2,-1
                        do k=1,3
                           iactiveline(k,ii)=iactiveline(k,ii-1)
                        enddo
                     enddo
                     iactiveline(1,id+1)=indexl
                     iactiveline(2,id+1)=nelemm
                     iactiveline(3,id+1)=ifreeintersec
                     ninsertl=ninsertl+1
                     insertl(ninsertl)=indexl
                  endif
               endif
            else
               if(eplane(pb,xcp,t).le.1.0d-12)then
                  if(abs(eplane(pb,xcp,t)).lt.1.0d-10)then
                     do ii=1,3
                        xinters(ii)=pb(ii)
                     enddo
                     ncvertex=ncvertex+1
                     do k=1,3
                        c_pvertex(k,ncvertex)=pb(k)
                     enddo
                  else
                     call intersectionpoint(pa,pb,xcp,t,xinters)
                     ncvertex=ncvertex+2
                     do k=1,3
                        c_pvertex(k,ncvertex-1)=xinters(k)
                        c_pvertex(k,ncvertex)=pb(k)
                     enddo
                  endif       
                  if((.not.oldactive).and.(.not.altered))then
                     if(eplane(pb,xcp,t).lt.0.d0.and.nvertex.gt.2)then
                        altered=.true.
                        nactiveline=nactiveline+1
                        ifreeintersec=ifreeintersec+1
                        do ii=nactiveline,id+2,-1
                           do k=1,3
                              iactiveline(k,ii)=iactiveline(k,ii-1)
                           enddo
                        enddo
                        iactiveline(1,id+1)=indexl
                        iactiveline(2,id+1)=nelemm
                        iactiveline(3,id+1)=ifreeintersec
                        ninsertl=ninsertl+1
                        insertl(ninsertl)=indexl
                     endif
                  endif
               else
                  if((.not.oldactive).and.(.not.altered))then
                     altered=.true.
                     nactiveline=nactiveline+1
                     ifreeintersec=ifreeintersec+1
                     do ii=nactiveline,id+2,-1
                        do k=1,3
                           iactiveline(k,ii)=iactiveline(k,ii-1)
                        enddo
                     enddo
                     iactiveline(1,id+1)=indexl
                     iactiveline(2,id+1)=nelemm
                     iactiveline(3,id+1)=ifreeintersec
                     ninsertl=ninsertl+1
                     insertl(ninsertl)=indexl
                  endif    
               endif
            endif 
!     
!     end loop over polygon vertices
!     
         enddo
         do j=1,ncvertex
            do k=1,3
               pvertex(k,j)=c_pvertex(k,j)
            enddo
         enddo
         nvertex=ncvertex
!     
!     end loop over clipping edges
!     
      enddo
!     
!     remove inserted active lines,if polygon is degenerated
!     
      if(nvertex.ge.3)then
         area=0
        do k=1,nvertex-2
         p1(1)=pvertex(1,k+1)-pvertex(1,1)
         p1(2)=pvertex(2,k+1)-pvertex(2,1)
         p1(3)=pvertex(3,k+1)-pvertex(3,1)
         p2(1)=pvertex(1,k+2)-pvertex(1,1)
         p2(2)=pvertex(2,k+2)-pvertex(2,1)
         p2(3)=pvertex(3,k+2)-pvertex(3,1)
         areax=((p1(2)*p2(3))-(p2(2)*p1(3)))**2
         areay=(-(p1(1)*p2(3))+(p2(1)*p1(3)))**2
         areaz=((p1(1)*p2(2))-(p2(1)*p1(2)))**2
         area=area+dsqrt(areax+areay+areaz)/2.d0
       enddo
         if(border)write(20,*)'border reached'
      endif
      if(nvertex.lt.3.or.border)then
         do i=1,ninsertl
            oldactive=.false.
            indexl=insertl(i)
            call nidentk(iactiveline,indexl, nactiveline,id,ithree)
            if(id.gt.0)then
               if(iactiveline(1,id).eq.indexl) oldactive=.true.
            endif
            
            if(oldactive)then
               nactiveline=nactiveline-1
               do ii=id,nactiveline
                  do k=1,3
                     iactiveline(k,ii)=iactiveline(k,ii+1)
                  enddo
               enddo  
            endif
         enddo
      endif 
!
      return
      end  
