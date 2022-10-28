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
      subroutine materialdata_cp(imat,ntmat_,t1l,shcon,nshcon,cp)
!
      implicit none
!
!     determines the specific heat
!
      integer imat,ntmat_,id,nshcon(*),four
!
      real*8 t1l,shcon(0:3,ntmat_,*),cp
!
      four=4
!     
!     calculating the specific heat
!
      call ident2(shcon(0,1,imat),t1l,nshcon(imat),four,id)
      if(nshcon(imat).eq.0) then
         continue
      elseif(nshcon(imat).eq.1) then
         cp=shcon(1,1,imat)
      elseif(id.eq.0) then
         cp=shcon(1,1,imat)
      elseif(id.eq.nshcon(imat)) then
         cp=shcon(1,id,imat)
      else
         cp=shcon(1,id,imat)+
     &        (shcon(1,id+1,imat)-shcon(1,id,imat))*
     &        (t1l-shcon(0,id,imat))/
     &        (shcon(0,id+1,imat)-shcon(0,id,imat))
      endif
!
      return
      end







