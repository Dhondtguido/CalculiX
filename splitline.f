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
      subroutine splitline(text,textpart,n)
!
      implicit none
!
!     splits an input line (text) in n comma separated fields (textpart)
!
!     n = # comma's +1,
!
      integer n,i,j,k,ierror
!
      character*1 ctext
      character*132 textpart(16)
      character*1320 text
!
      n=1
      j=0
      do i=1,1320
         ctext=text(i:i)
         if(ctext.ne.',') then
            if(ctext.eq.' ') then
c               cycle
               exit
            else
c               if((ichar(ctext).ge.97).and.(ichar(ctext).le.122))
c     &             ctext=char(ichar(ctext)-32)
            endif
            j=j+1
            if(j.le.132) textpart(n)(j:j)=ctext
         else
            do k=j+1,132
               textpart(n)(k:k)=' '
            enddo
            j=0
            n=n+1
            if(n.gt.16) then
               ierror=0
               do k=i+1,1320
                  if(text(k:k).eq.',') cycle
                  if(text(k:k).eq.' ') then
                     if(ierror.eq.0) then
                        exit
                     else
                        write(*,*) 
     &                     '*ERROR in splitline: there should not'
                        write(*,*)'       be more than 16 entries in a '
                        write(*,*) '       line; '
                        write(*,'(a)') text(1:k-1)
c                        write(*,'(a1320)') text(1:k-1)
                        call exit(201)
                     endif
                  endif
                  ierror=1
               enddo
               exit
            endif
         endif
      enddo
      if(j.eq.0) then
         n=n-1
      else
         do k=j+1,132
            textpart(n)(k:k)=' '
         enddo
      endif
!
!     clearing all textpart fields not used
!
      do i=n+1,16
         textpart(i)='
     &
     &                '
      enddo
!
      return
      end



