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
      subroutine checktruecontact(ntie,tieset,tietol,elcon,itruecontact,
     &  ncmat_,ntmat_)
!
!     check whether for face-to-face penalty contact the
!     surface behavior definition is such that true contact is 
!     defined and not just tied contact
!
      implicit none
!
      character*81 tieset(3,*)
!
      integer itruecontact,i,ntie,imat,ncmat_,ntmat_
!
      real*8 tietol(4,*),elcon(0:ncmat_,ntmat_,*)
!
!     if at least one tied contact is present, itruecontact
!     is set to zero and no check is performed whether tension
!     occurs in the contact areas
!
      itruecontact=1
      do i=1,ntie
         if(tieset(1,i)(81:81).eq.'C') then
            imat=int(tietol(2,i))
            if(int(elcon(3,1,imat)).eq.4) then
               itruecontact=0
               exit
            endif
         endif
      enddo
!
      return
      end
