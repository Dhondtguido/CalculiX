!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine sutherland_hodgman(nopes,xn,xl3sp,xlpgp,
     &     nodepg,ipe,ime,iactiveline,nactiveline,
     &     ifacem,npg,nvertex,pvertex)
!     
      implicit none 
!     
      logical invert,oldactive,altered
!     
      integer nvertex,nopes,ipe(*),ime(4,*),iactiveline(2,*),
     &     nactiveline,itri,ifacem,i,ii,j,k,npg,id,
     &     nodepg(*),ncvertex,node1,node2,modf,node,indexline,itwo, 
     &     insertl(3),ninsertl,node1loc,node2loc,node3loc
!     
      real*8 pvertex(3,*),xn(3),xl3sp(3,*),pa(3),pb(3),xinters(3),
     &     xcp(3),diff,dd,xlpgp(3,*),c_pvertex(3,13),t,cedge(3),
     &     xtest(3),eplane
!     
      data itwo /2/
!     
      nvertex=0
      ninsertl=0
!     
!     Initialize slave polygon
!     
      do j=1,nopes
        nvertex=nvertex+1
        do k=1,3
          pvertex(k,nvertex)=xl3sp(k,j)
        enddo
      enddo
!     
!     loop over clipping edges
!     
      do i=1,npg
        ncvertex=0
        altered=.false.
!     
!     generate clipping plane
!
        node1loc=modf(npg,i)
        node2loc=modf(npg,i+1)
        node3loc=modf(npg,i+2)
!        
        node1=nodepg(node1loc)
        node2=nodepg(node2loc)
        invert=.false.
        if(node2.lt.node1) then
          node=node1
          node1=node2
          node2=node
          invert=.true.
        endif
        indexline=ipe(node1)
        do
          if(ime(1,indexline).eq.node2) exit
          indexline=ime(4,indexline)
          if(indexline.eq.0) then
            write(*,*) '*ERROR in sutherland_hodgman:'
            write(*,*) '       line was not properly catalogued'
            write(*,*) '       itri',itri,'node1',node1,'node2',node2
            call exit(201)
          endif 
        enddo
        do k=1,3
          cedge(k)=xlpgp(k,node2loc)-xlpgp(k,node1loc)
        enddo
        xcp(1)=xn(2)*cedge(3)-xn(3)*cedge(2)
        xcp(2)=xn(3)*cedge(1)-xn(1)*cedge(3)
        xcp(3)=xn(1)*cedge(2)-xn(2)*cedge(1)
        dd=dsqrt(xcp(1)**2+xcp(2)**2+xcp(3)**2)
        do k=1,3
          xcp(k)=xcp(k)/dd
        enddo
        t=-eplane(xlpgp(1,node1loc),xcp,0.d0)
!     
!     inside-outside-test 
!     
        do k=1,3 
          xtest(k)=xlpgp(k,node3loc)
        enddo
        if(eplane(xtest,xcp,t).gt.0) then
          t=-t
          do k=1,3
            xcp(k)=-xcp(k)
          enddo
        endif
        oldactive=.false.
        call nidentk(iactiveline,indexline,nactiveline,id,itwo)
        if((id.gt.0).and.(iactiveline(1,id).eq.indexline)) then
          oldactive=.true.
        endif    
        if(oldactive) then
          nactiveline=nactiveline-1
          do ii=id,nactiveline
            do k=1,2
              iactiveline(k,ii)=iactiveline(k,ii+1)
            enddo
          enddo 
        endif
!     
        if(nvertex.lt.3) cycle
!     
!     loop over slave face vertices
!     
        do j=0,nvertex-1
          do k=1,3
            pa(k)=pvertex(k,modf(nvertex,j))
            pb(k)=pvertex(k,modf(nvertex,j+1))
          enddo 
          if(eplane(pa,xcp,t).le.1.d-12) then
            if(eplane(pb,xcp,t).le.1.d-12) then
              ncvertex=ncvertex+1
              do k=1,3
                c_pvertex(k,ncvertex)=pb(k)
              enddo 
            else
              if(abs(eplane(pa,xcp,t)).gt.1.d-10) then
                call intersectionpoint(pa,pb,xcp,t,xinters)
                diff=(xinters(1)-pa(1))**2+(xinters(2)-pa(2))**2+
     &               (xinters(3)-pa(3))**2
                diff=dsqrt(diff)
                if(diff.gt.1.d-11) then
                  ncvertex=ncvertex+1
                  do k=1,3
                    c_pvertex(k,ncvertex)=xinters(k)
                  enddo
                endif
              endif 
!     
              if((.not.oldactive).and.(.not.altered)) then
                altered=.true.
                nactiveline=nactiveline+1
                do ii=nactiveline,id+2,-1
                  do k=1,2
                    iactiveline(k,ii)=iactiveline(k,ii-1)
                  enddo
                enddo
                iactiveline(1,id+1)=indexline
                iactiveline(2,id+1)=ifacem
                ninsertl=ninsertl+1
                insertl(ninsertl)=indexline
              endif
            endif
          else
            if(eplane(pb,xcp,t).le.1.d-12) then
              if(abs(eplane(pb,xcp,t)).lt.1.d-10) then
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
              if((.not.oldactive).and.(.not.altered)) then
                if((eplane(pb,xcp,t).lt.0.d0).and.(nvertex.gt.2)) then
                  altered=.true.
                  nactiveline=nactiveline+1
                  do ii=nactiveline,id+2,-1
                    do k=1,2
                      iactiveline(k,ii)=iactiveline(k,ii-1)
                    enddo
                  enddo
                  iactiveline(1,id+1)=indexline
                  iactiveline(2,id+1)=ifacem
                  ninsertl=ninsertl+1
                  insertl(ninsertl)=indexline
                endif
              endif
            else
              if((.not.oldactive).and.(.not.altered)) then
                altered=.true.
                nactiveline=nactiveline+1
                do ii=nactiveline,id+2,-1
                  do k=1,2
                    iactiveline(k,ii)=iactiveline(k,ii-1)
                  enddo
                enddo
                iactiveline(1,id+1)=indexline
                iactiveline(2,id+1)=ifacem
                ninsertl=ninsertl+1
                insertl(ninsertl)=indexline
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
      if(nvertex.lt.3) then
        do i=1,ninsertl
          oldactive=.false.
          indexline=insertl(i)
          call nidentk(iactiveline,indexline, nactiveline,id,itwo)
          if(id.gt.0) then
            if(iactiveline(1,id).eq.indexline) oldactive=.true.
          endif
          
          if(oldactive) then
            nactiveline=nactiveline-1
            do ii=id,nactiveline
              do k=1,2
                iactiveline(k,ii)=iactiveline(k,ii+1)
              enddo
            enddo  
          endif
        enddo
      endif 
!     
      return
      end  
