!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
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
      subroutine writeinput(inpc,ipoinp,inp,nline,ninp,ipoinpc)
!
      implicit none
!
      integer nentries
      parameter(nentries=18)
!
      character*1 inpc(*)
      character*20 nameref(nentries)
!
      integer nline,i,j,ninp,ipoinp(2,nentries),inp(3,ninp),
     &  ipoinpc(0:*)
!
      data nameref /'RESTART,READ','NODE','USERELEMENT','ELEMENT',
     &              'NSET',
     &              'ELSET','SURFACE','TRANSFORM','MATERIAL',
     &              'DISTRIBUTION',
     &              'ORIENTATION','TIE','INTERACTION',
     &              'INITIALCONDITIONS','AMPLITUDE',
     &              'CONTACTPAIR','COUPLING','REST'/
!
      open(16,file='input.inpc',status='unknown',err=161)
      do i=1,nline
         write(16,'(1x,i6,1x,1320a1)') i,
     &       (inpc(j),j=ipoinpc(i-1)+1,ipoinpc(i))
      enddo
      close(16)
!
      open(16,file='input.ipoinp',status='unknown',err=162)
      do i=1,nentries
         write(16,'(1x,a20,1x,i6,1x,i6)') nameref(i),(ipoinp(j,i),j=1,2)
      enddo
      close(16)
!
      open(16,file='input.inp',status='unknown',err=163)
      do i=1,ninp
         write(16,'(1x,i3,1x,i6,1x,i6,1x,i6)') i,(inp(j,i),j=1,3)
      enddo
      close(16)
!
      return
!
 161  write(*,*) '*ERROR in writeinput: could not open file input.inpc'
      call exit(201)
!
 162  write(*,*) 
     &    '*ERROR in writeinput: could not open file input.ipoinp'
      call exit(201)
!
 163  write(*,*) '*ERROR in writeinput: could not open file input.inp'
      call exit(201)
      end
