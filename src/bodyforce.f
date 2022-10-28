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
      subroutine bodyforce(cbody,ibody,ipobody,nbody,set,istartset,
     &  iendset,ialset,inewton,nset,ifreebody,k)
!
!     assigns the body forces to the elements by use of field ipobody
!
      implicit none
!
      character*81 cbody(*),elset,set(*)
!
      integer ibody(3,*),ipobody(2,*),i,j,l,istartset(*),nbody,
     &  iendset(*),ialset(*),kindofbodyforce,inewton,nset,istat,
     &  ifreebody,k,index,id
!
      elset=cbody(k)
      kindofbodyforce=ibody(1,k)
      if(kindofbodyforce.eq.3) inewton=1
!
!     check whether element number or set name
!
      read(elset,'(i21)',iostat=istat) l
      if(istat.eq.0) then
         if(ipobody(1,l).eq.0) then
            ipobody(1,l)=k
         else
!
            index=l
            do
               if(ipobody(1,index).eq.k) exit
               if(ipobody(2,index).eq.0) then
                  ipobody(2,index)=ifreebody
                  ipobody(1,ifreebody)=k
                  ipobody(2,ifreebody)=0
                  ifreebody=ifreebody+1
                  exit
               endif
               index=ipobody(2,index)
            enddo
         endif
         return
      endif
!
!     set name
!
c      do i=1,nset
c         if(set(i).eq.elset) exit
c      enddo
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
            if(ipobody(1,l).eq.0) then
               ipobody(1,l)=k
            else
!
               index=l
               do
                  if(ipobody(1,index).eq.k) exit
                  if(ipobody(2,index).eq.0) then
                     ipobody(2,index)=ifreebody
                     ipobody(1,ifreebody)=k
                     ipobody(2,ifreebody)=0
                     ifreebody=ifreebody+1
                     exit
                  endif
                  index=ipobody(2,index)
               enddo
            endif
         else
            l=ialset(j-2)
            do
               l=l-ialset(j)
               if(l.ge.ialset(j-1)) exit
               if(ipobody(1,l).eq.0) then
                  ipobody(1,l)=k
               else
!
                  index=l
                  do
                     if(ipobody(1,index).eq.k) exit
                     if(ipobody(2,index).eq.0) then
                        ipobody(2,index)=ifreebody
                        ipobody(1,ifreebody)=k
                        ipobody(2,ifreebody)=0
                        ifreebody=ifreebody+1
                        exit
                     endif
                     index=ipobody(2,index)
                  enddo
               endif
            enddo
         endif
      enddo
!     
      return
      end

