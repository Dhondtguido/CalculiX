!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine checkdispoutonly(prlab,nprint,nlabel,filab,
     &     idispfrdonly)
!     
      implicit none
!     
!     checks for massless explicit dynamics if anything else but    
!     displacements have to be stored in the .dat or .frd file    
!     
      character*6 prlab(*)
      character*87 filab(*)
!     
      integer i,nprint,nlabel,idispfrdonly
!     
      idispfrdonly=-1
!     
!     storage in the .dat file     
!
      if(nprint.gt.0) idispfrdonly=0
!     
!     storage in the .frd file        
!     
      do i=1,nlabel
        if(filab(i)(1:2).eq.'U ') then
          if(idispfrdonly.lt.0) then
            idispfrdonly=1
          endif
        elseif(filab(i)(1:4).ne.'    ') then
          idispfrdonly=0
        endif
      enddo
!
      if(idispfrdonly.lt.0) idispfrdonly=1
!     
      return
      end
