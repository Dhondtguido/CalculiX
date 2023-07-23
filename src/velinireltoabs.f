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
      subroutine velinireltoabs(ibody,xbody,cbody,nbody,set,
     &     istartset,iendset,ialset,nset,veold,mi,ipkon,kon,lakon,co,
     &     itreated)
!
!     assigns the body forces to the elements by use of field ipobody
!
      implicit none
!
      character*8 lakon(*)
      character*81 cbody(*),elset,set(*)
!
      integer ibody(3,*),i,j,k,l,m,istartset(*),nbody,kon(*),id,
     &     iendset(*),ialset(*),nset,istat,indexe,mi(*),ipkon(*),
     &     node,nope,itreated(*)
!
      real*8 xbody(7,*),veold(0:mi(2),*),co(3,*),om,a(3),xn(3),p(3),dd
!
      do k=1,nbody
         elset=cbody(k)
!
!        check for centrifugal loads
!
         if(ibody(1,k).eq.1) then
            ibody(1,k)=0
         else
            cycle
         endif
!
!        determine
!          om: rotational speed
!          a: point on axis
!          xn: unit vector on axis         
!
         om=dsqrt(xbody(1,k))
         a(1)=xbody(2,k)
         a(2)=xbody(3,k)
         a(3)=xbody(4,k)
         xn(1)=xbody(5,k)
         xn(2)=xbody(6,k)
         xn(3)=xbody(7,k)
!
!     check whether element number or set name
!
         read(elset,'(i21)',iostat=istat) l
         if(istat.eq.0) then
!
!           treat element l
!
            indexe=ipkon(l)
!
!           nope is the number of nodes belonging to the element
!
            if(lakon(l)(4:5).eq.'20') then
               nope=20
            elseif(lakon(l)(4:4).eq.'8') then
               nope=8
            elseif(lakon(l)(4:5).eq.'10') then
               nope=10
            elseif(lakon(l)(4:4).eq.'4') then
               nope=4
            elseif(lakon(l)(4:5).eq.'15') then
               nope=15
            elseif(lakon(l)(4:4).eq.'6') then
               nope=6
            elseif(lakon(l)(1:2).eq.'ES') then
!
!              contact elements need not be treated
!
               if(lakon(l)(7:7).ne.'C') then
                  nope=ichar(lakon(l)(8:8))-47
               endif
            elseif(lakon(l)(1:4).eq.'MASS') then
               nope=1
            endif
!
            do i=1,nope
               node=kon(indexe+i)
               if(itreated(node).eq.1) cycle
               do m=1,3
                  p(m)=co(m,node)-a(m)
               enddo
               dd=p(1)*xn(1)+p(2)*xn(2)+p(3)*xn(3)
               do m=1,3
                  p(m)=(p(m)-dd*xn(m))*om
               enddo
               veold(1,node)=veold(1,node)+
     &                       xn(2)*p(3)-xn(3)*p(2)
               veold(2,node)=veold(2,node)+
     &                       xn(3)*p(1)-xn(1)*p(3)
               veold(3,node)=veold(3,node)+
     &              xn(1)*p(2)-xn(2)*p(1)
               itreated(node)=1
            enddo
         else
!
!     set name
!
c            do i=1,nset
c               if(set(i).eq.elset) exit
c            enddo
           call cident81(set,elset,nset,id)
           i=nset+1
           if(id.gt.0) then
             if(elset.eq.set(id)) then
               i=id
             endif
           endif
!     
            do j=istartset(i),iendset(i)
               if(ialset(j).gt.0) then
                  l=ialset(j)
!
!                 treat element l
!
                  indexe=ipkon(l)
!
!                 nope is the number of nodes belonging to the element
!
                  if(lakon(l)(4:5).eq.'20') then
                     nope=20
                  elseif(lakon(l)(4:4).eq.'8') then
                     nope=8
                  elseif(lakon(l)(4:5).eq.'10') then
                     nope=10
                  elseif(lakon(l)(4:4).eq.'4') then
                     nope=4
                  elseif(lakon(l)(4:5).eq.'15') then
                     nope=15
                  elseif(lakon(l)(4:4).eq.'6') then
                     nope=6
                  elseif(lakon(l)(1:2).eq.'ES') then
!     
!                    contact elements need not be treated
!     
                     if(lakon(l)(7:7).ne.'C') then
                        nope=ichar(lakon(l)(8:8))-47
                     endif
                  elseif(lakon(l)(1:4).eq.'MASS') then
                     nope=1
                  endif
!     
                  do i=1,nope
                     node=kon(indexe+i)
                     if(itreated(node).eq.1) cycle
                     do m=1,3
                        p(m)=co(m,node)-a(m)
                     enddo
                     dd=p(1)*xn(1)+p(2)*xn(2)+p(3)*xn(3)
                     do m=1,3
                        p(m)=(p(m)-dd*xn(m))*om
                     enddo
                     veold(1,node)=veold(1,node)+
     &                    xn(2)*p(3)-xn(3)*p(2)
                     veold(2,node)=veold(2,node)+
     &                    xn(3)*p(1)-xn(1)*p(3)
                     veold(3,node)=veold(3,node)+
     &                    xn(1)*p(2)-xn(2)*p(1)
                     itreated(node)=1
                  enddo
               else
                  l=ialset(j-2)
                  do
                     l=l-ialset(j)
                     if(l.ge.ialset(j-1)) exit
!     
!                    treat element l
!     
                     indexe=ipkon(l)
!     
!                    nope is the number of nodes belonging to the element
!     
                     if(lakon(l)(4:5).eq.'20') then
                        nope=20
                     elseif(lakon(l)(4:4).eq.'8') then
                        nope=8
                     elseif(lakon(l)(4:5).eq.'10') then
                        nope=10
                     elseif(lakon(l)(4:4).eq.'4') then
                        nope=4
                     elseif(lakon(l)(4:5).eq.'15') then
                        nope=15
                     elseif(lakon(l)(4:4).eq.'6') then
                        nope=6
                     elseif(lakon(l)(1:2).eq.'ES') then
!     
!                       contact elements need not be treated
!     
                        if(lakon(l)(7:7).ne.'C') then
                           nope=ichar(lakon(l)(8:8))-47
                        endif
                     elseif(lakon(l)(1:4).eq.'MASS') then
                        nope=1
                     endif
!     
                     do i=1,nope
                        node=kon(indexe+i)
                        if(itreated(node).eq.1) cycle
                        do m=1,3
                           p(m)=co(m,node)-a(m)
                        enddo
                        dd=p(1)*xn(1)+p(2)*xn(2)+p(3)*xn(3)
                        do m=1,3
                           p(m)=(p(m)-dd*xn(m))*om
                        enddo
                        veold(1,node)=veold(1,node)+
     &                       xn(2)*p(3)-xn(3)*p(2)
                        veold(2,node)=veold(2,node)+
     &                       xn(3)*p(1)-xn(1)*p(3)
                        veold(3,node)=veold(3,node)+
     &                       xn(1)*p(2)-xn(2)*p(1)
                        itreated(node)=1
                     enddo
                  enddo
               endif
            enddo
         endif
      enddo
!     
      return
      end

