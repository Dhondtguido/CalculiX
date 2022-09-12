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
      subroutine identifytransform(nelement,label,nelemload,sideload,
     &  nload,loadid)
!
!     checks whether a transformation was applied to a given face,
!     and if so, which one (only for CFD)
!
      implicit none
!
      character*20 label,sideload(*)
!
      integer nelemload(2,*),nelement,nload,id,loadid
!
      loadid=0
!
      call nident2(nelemload,nelement,nload,id)
      if(id.gt.0) then
!
!        it is possible that several *DLOAD, *FILM or
!        *RADIATE boundary conditions are applied to one
!        and the same element
!
         if(nelemload(1,id).eq.nelement) then
            do
               if (sideload(id).eq.label) then
                  loadid=id
               elseif(sideload(id).lt.label) then
                  exit
               endif
               id=id-1
               if((id.eq.0).or.(nelemload(1,id).ne.nelement)) then
                  exit
               endif
            enddo
         endif
      endif
!
      return
      end

