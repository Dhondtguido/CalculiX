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
      subroutine headings(inpc,textpart,istat,n,iline,ipol,inl,ipoinp,
     &  inp,ipoinpc,heading,istep,irstrt,ier)
!
!     reading the input deck: *HEADING
!
      implicit none
!
      character*1 inpc(*)
      character*66 heading(*)
      character*132 textpart(16)
      character*1320 text
!
      integer istat,n,key,iline,ipol,inl,ipoinp(2,*),inp(3,*),
     &  ipoinpc(0:*),i,j,nentries,nheading,istep,irstrt(*),ier
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *HEADING: *HEADING should be placed'
         write(*,*) '       before all step definitions'
         ier=1
         return
      endif
!
      do i=2,n
         write(*,*) 
     &        '*WARNING in headings: parameter not recognized:'
         write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
         call inputwarning(inpc,ipoinpc,iline,
     &"*HEADING%")
      enddo
!
      nheading=0
      do
!
!     reading a new line
!
         if(iline.eq.inp(2,inl)) then
            if(inp(3,inl).eq.0) then
               do
                  ipol=ipol+1
                  if(ipol.gt.nentries) then
                     istat=-1
                     return
                  elseif(ipoinp(1,ipol).ne.0) then
                     exit
                  endif
               enddo
               inl=ipoinp(1,ipol)
               iline=inp(1,inl)
            else
               inl=inp(3,inl)
               iline=inp(1,inl)
            endif
         else
            iline=iline+1
         endif
         j=0
         do i=ipoinpc(iline-1)+1,ipoinpc(iline)
            j=j+1
            text(j:j)=inpc(i)
         enddo
         text(j+1:j+1)=' '
!     
         istat=0
         key=0
!     
!     only free format is supported
!     
         if((text(1:1).eq.'*').and.(text(2:2).ne.'*')) then
            key=1
         endif
!     
         call splitline(text,textpart,n)
         if((key.ne.0).or.(istat.lt.0)) exit
!
         do i=j+2,66
            text(i:i)=' '
         enddo
!     
         nheading=nheading+1
         heading(nheading)=text(1:66)
      enddo
!
      return
      end

