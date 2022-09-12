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
      subroutine writesen(g0,dgdx,ndesi,nobject,nodedesi,jobnamef)
!
!     writes "raw" sensitivities to file jobname.sen
!
      implicit none
!
      character*132 jobnamef(*),cfile
!
      integer ndesi,nobject,nodedesi(*),i,j
!
      real*8 g0(*),dgdx(ndesi,*)
!
!
!     storing the objectives
!
      do i=1,132
         cfile(i:i)=' '
      enddo
      do i=1,132
         if(jobnamef(1)(i:i).eq.' ') exit
         cfile(i:i)=jobnamef(1)(i:i)
      enddo
      cfile(i:i+4)='.sen0'
      open(27,file=cfile,status='unknown')
!
      write(27,*) (g0(j),j=1,nobject)
!
      close(27)
!
!     storing the sensitivity of the objectives
!
      cfile(i+4:i+4)='1'
      open(27,file=cfile,status='unknown')
!
      do i=1,ndesi
         write(27,*) nodedesi(i),(dgdx(i,j),j=1,nobject)
      enddo
!
      close(27)
!
      return
      end

