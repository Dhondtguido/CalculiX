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
      subroutine bodyadd(cbody,ibody,xbody,nbody,nbody_,set,label,
     &  iamplitude,xmagnitude,p1,p2,bodyf,xbodyold,lc,idefbody)
!
!     adds a volumetric dload condition to the data base
!
      implicit none
!
      character*20 label
      character*81 set,cbody(*)
!
      integer ibody(3,*),nbody,nbody_,id,iamplitude,ilabel,i,j,id1,lc,
     &  idefbody(*)
!
      real*8 xbody(7,*),p1(3),p2(3),bodyf(3),xmagnitude,xbodyold(7,*),
     &  dd,p(3)
!
!     assigning a number to the load type (stored in ibody(1,*))
!
      if(label(1:7).eq.'CENTRIF') then
         ilabel=1
      elseif(label(1:4).eq.'GRAV') then
         ilabel=2
      elseif(label(1:6).eq.'NEWTON') then
         ilabel=3
      elseif(label(1:5).eq.'CORIO') then
         ilabel=4
      endif
!
!     normalizing the direction for gravity forces
!
      if(ilabel.eq.2) then
         dd=dsqrt(bodyf(1)*bodyf(1)+bodyf(2)*bodyf(2)+bodyf(3)*bodyf(3))
         do i=1,3
            bodyf(i)=bodyf(i)/dd
         enddo
      endif
!
!     checking whether a similar load type was already assigned to the
!     same set
!
      call cident(cbody,set,nbody,id)
!
      if(id.ne.0) then
         do
            if(id.eq.0) exit
            if(cbody(id).eq.set) then
               if(ibody(1,id).eq.ilabel) then
!
!                 for gravity forces the gravity direction is
!                 checked; if the direction is different,it is 
!                 a new loading
!
                  if(ilabel.eq.2) then
                     if(dabs(bodyf(1)*xbody(2,id)+bodyf(2)*xbody(3,id)+
     &                  bodyf(3)*xbody(4,id)-1.d0).gt.1.d-10) then
                        id=id-1
                        cycle
                     endif
                  endif
!
!                 for centrifugal loads the centrifugal axis is 
!                 checked
!
                  if(ilabel.eq.1) then
                     if(dabs(p2(1)*xbody(5,id)+p2(2)*xbody(6,id)+
     &                       p2(3)*xbody(7,id)-1.d0).gt.1.d-10) then
                        id=id-1
                        cycle
                     endif
                     do i=1,3
                        p(i)=xbody(1+i,id)-p1(i)
                     enddo
                     dd=dsqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3))
                     if(dd.gt.1.d-10) then
                        do i=1,3
                           p(i)=p(i)/dd
                        enddo
                        if(dabs(p(1)*xbody(5,id)+p(2)*xbody(6,id)+
     &                          p(3)*xbody(7,id)-1.d0).gt.1.d-10) then
                           id=id-1
                           cycle
                        endif
                     endif
                  endif
!
!                 check for the same loadcase
!
                  if(ibody(3,id).ne.lc) then
                     id=id-1
                     cycle
                  endif
!
                  ibody(2,id)=iamplitude
                  ibody(3,id)=lc
                  if(ilabel.eq.1) then
                     if(idefbody(id).eq.0) then
                        xbody(1,id)=xmagnitude
                        idefbody(id)=1
                     else
                        if(ibody(2,id).ne.iamplitude) then
                           write(*,*) '*ERROR in bodyadd:'
                           write(*,*) '       it is not allowed to add'
                           write(*,*)'       two centrifugal loads with'
                           write(*,*) '       different amplitudes'
                           call exit(201)
                        endif
                        xbody(1,id)=xbody(1,id)+xmagnitude
                     endif
                     xbody(2,id)=p1(1)
                     xbody(3,id)=p1(2)
                     xbody(4,id)=p1(3)
                     xbody(5,id)=p2(1)
                     xbody(6,id)=p2(2)
                     xbody(7,id)=p2(3)
                  elseif(ilabel.eq.2) then
                     if(idefbody(id).eq.0) then
                        xbody(1,id)=xmagnitude
                        idefbody(id)=1
                     else
                        if(ibody(2,id).ne.iamplitude) then
                           write(*,*) '*ERROR in bodyadd:'
                           write(*,*) '       it is not allowed to add'
                           write(*,*) '       two gravity loads with'
                           write(*,*) '       different amplitudes'
                           call exit(201)
                        endif
                        xbody(1,id)=xbody(1,id)+xmagnitude
                     endif
                     xbody(2,id)=bodyf(1)
                     xbody(3,id)=bodyf(2)
                     xbody(4,id)=bodyf(3)
                  endif
                  return
               endif
               id=id-1
            else
               exit
            endif
         enddo
      endif
!
!     new set/loadtype combination
!
      nbody=nbody+1
      if(nbody.gt.nbody_) then
         write(*,*) '*ERROR in bodyadd: increase nbody_'
         call exit(201)
      endif
!
!     reordering the arrays
!
      do i=nbody,id+2,-1
         cbody(i)=cbody(i-1)
         idefbody(i)=idefbody(i-1)
         do j=1,3
            ibody(j,i)=ibody(j,i-1)
         enddo
         do j=1,7
            xbody(j,i)=xbody(j,i-1)
            xbodyold(j,i)=xbodyold(j,i-1)
         enddo
      enddo
!
!     inserting the new values
!
      id1=id+1
!
      cbody(id1)=set
      idefbody(id1)=1
      ibody(1,id1)=ilabel
      ibody(2,id1)=iamplitude
      ibody(3,id1)=lc
      if(ilabel.eq.1) then
         xbody(1,id1)=xmagnitude
         xbody(2,id1)=p1(1)
         xbody(3,id1)=p1(2)
         xbody(4,id1)=p1(3)
         xbody(5,id1)=p2(1)
         xbody(6,id1)=p2(2)
         xbody(7,id1)=p2(3)
      elseif(ilabel.eq.2) then
         xbody(1,id1)=xmagnitude
         xbody(2,id1)=bodyf(1)
         xbody(3,id1)=bodyf(2)
         xbody(4,id1)=bodyf(3)
      endif
!
      return
      end

