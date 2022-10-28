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
      subroutine writeobj(objectset,iobject,g0,dgdxglob,nobject,
     &   ndesi,nodedesi,nk,nobjectstart)
!
!     writes the results design repsonse information in the .dat file
!
      implicit none
!
      character*81 objectset(5,*)
!      
      integer iobject,i,j,nobject,ndesi,nodedesi(*),nk,nobjectstart
!
      real*8 g0(*),dgdxglob(2,nk,*),dd,inode
!          
!
!     write header in .dat file
!
      if(iobject.eq.nobjectstart) then
         write(5,*)
         write(5,'(a113)') '   #########################################
     &##################################################################
     &######'
         write(5,*) '  D E S I G N   R E S P O N S E   
     &I N F O R M A T I O N'
         write(5,*)
         write(5,'(3x,16a,3x,18a,3x,18a,3x,80a)') 'FUNCTION        ',
     &   'FUNCTION VALUE    ','LENGTH GRADIENT   ','NAME                                    
     &                                        '
         write(5,'(a113)') '   #########################################
     &##################################################################
     &######'
         write(5,*)
      endif
!
!     write design repsonse in .dat file
!
      i=iobject+1
      dd=0.d0
      do j=1,ndesi
         inode=nodedesi(j)
         dd=dd+dgdxglob(1,inode,i)**2
      enddo
      dd=dsqrt(dd)
      write(5,'(3x,a16,e14.7,3x,e16.7,3x,a80)') objectset(1,i),g0(i),
     &   dd,objectset(5,i)
!
 101  format(3x,16a,3x,a16,a11,3x,a11,3x,a11,3x,a8,3x,a10,3x,a10)
 102  format(3x,i2,8x,3x,a16,a4,3x,e14.7,3x,e14.7,3x,e14.7,3x,a8,3x,a80)
!
      return
      end

