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
      subroutine checktempload(iamload,nload,sideload,ibody,nbody,
     &     masslesslinear,nloadrhs,nbodyrhs,nam)
!     
!     checks for linear massless explicit dynamic calculations
!     whether the distributed loading changes in the step     
!     (yes -> masslesslinear=1; no -> masslesslinear=2). If no,
!     rhsmain has to calculate the corresponding entries in fext     
!     only once at the start of the step (in the first increment)
!     
      implicit none
!
      character*20 sideload(*)
!
      integer nload,nbody,masslesslinear,nloadrhs,nbodyrhs,iamload(2,*),
     &     ibody(3,*),i,nam
!
!     default: no change
!
      masslesslinear=2
      nloadrhs=0
      nbodyrhs=0
!
!     checking facial distributed loading for amplitudes 
!
      if(nam.gt.0) then
        do i=1,nload
          if(iamload(1,i).ne.0) then
            masslesslinear=1
            nloadrhs=nload
            nbodyrhs=nbody
            return
          endif
        enddo
      endif
!
!     checking facial distributed loading for user subroutines
!
      do i=1,nload
        if(sideload(i)(3:4).eq.'NU') then
          masslesslinear=1
          nloadrhs=nload
          nbodyrhs=nbody
          return
        endif
      enddo
!
!     checking body distributed loading for amplitudes
!
      do i=1,nbody
        if(ibody(2,i).ne.0) then
          masslesslinear=1
          nloadrhs=nload
          nbodyrhs=nbody
          return
        endif
      enddo
!     
      return
      end
