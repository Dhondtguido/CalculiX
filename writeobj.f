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
      subroutine writeobj(objectset,iobject,g0)
!
!     writes the results of the objective function in the .dat file
!
      implicit none
!
      character*81 objectset(5,*)
      integer iobject,i
      real*8 g0(*)
!          
      i=iobject+1
!
      if(i.eq.1) then
         write(5,*)
         write(5,*) '  #################################################
     &###############'
         write(5,*) '  D E S I G N   R E S P O N S E   I N F O R M A T I
     & O N'
         write(5,*)
         write(5,'(3x,a16,a14,3x,a80)') 'FUNCTION        ',
     &   'VALUE         ','NAME                                    
     &                                        '
         write(5,*) '  #################################################
     &###############'
         write(5,*)
      endif
!
      if(objectset(1,i)(1:11).ne.'PROJECTGRAD') then
         write(5,'(3x,a16,e14.7,3x,a80)') objectset(1,i),
     &   g0(i),objectset(5,i)
      endif
!      
      return
      end

