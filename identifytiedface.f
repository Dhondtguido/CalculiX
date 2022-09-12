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
      subroutine identifytiedface(tieset,ntie,set,nset,ifaceslave,kind)
!
!     identifies slave nodes in tied slave faces
!
      implicit none
!
      character*1 kind
      character*81 tieset(3,*),slavset,set(*)
!
      integer ifaceslave(*),i,j,nset,ipos,ntie,id
!
!     nodes per face for tet elements
!
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.kind) cycle
         slavset=tieset(2,i)
         ipos=index(slavset,' ')
         slavset(ipos:ipos)='T'
c         do j=1,nset
c            if(set(j).eq.slavset) exit
c         enddo
         call cident81(set,slavset,nset,id)
         j=nset+1
         if(id.gt.0) then
           if(slavset.eq.set(id)) then
             j=id
           endif
         endif
         if(j.gt.nset) then
            slavset(ipos:ipos)='S'
c            do j=1,nset
c               if(set(j).eq.slavset) then
c                  exit
c               endif
c           enddo
            call cident81(set,slavset,nset,id)
            j=nset+1
            if(id.gt.0) then
              if(slavset.eq.set(id)) then
                j=id
              endif
            endif
            if(j.gt.nset) then
               write(*,*) 
     &           '*ERROR in identifytiedface: ',
     &           'tied contact nodal slave surface ',
     &              slavset
               write(*,*) '       does not exist'
               call exit(201)
            else
               tieset(2,i)(ipos:ipos)='S'
               ifaceslave(i)=0
            endif
         else
            tieset(2,i)(ipos:ipos)='T'
            ifaceslave(i)=1
         endif
      enddo
      return
      end

