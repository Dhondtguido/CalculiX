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
      subroutine materialdata_cp_sec(imat,ntmat_,t1l,shcon,nshcon,cp,
     &  physcon)
!
      implicit none
!
!     determines the secant specific heat at constant pressure cp 
!
!     the difference with materialdata_cp is that the specific heat at
!     constant pressure cp as returned from the present routine
!     is the secant value and not the differential value. 
!     For the differential value we have:
!            dh=cp*dT
!     and consequently
!            h=int_from_0_to_T cp*dT
!     For the secant value one has:
!            h=cp_secant*T
!
      integer imat,ntmat_,id,nshcon(*),four,i
!
      real*8 t1l,shcon(0:3,ntmat_,*),cp,physcon(*)
!
      four=4
!     
!     calculating the tangent specific heat
!
      call ident2(shcon(0,1,imat),t1l,nshcon(imat),four,id)
      if(nshcon(imat).eq.0) then
         continue
      elseif(nshcon(imat).eq.1) then
         cp=shcon(1,1,imat)
      elseif(id.eq.0) then
         cp=shcon(1,1,imat)
      elseif(id.eq.nshcon(imat)) then
         cp=(shcon(0,1,imat)-physcon(1))*shcon(1,1,imat)
         do i=2,nshcon(imat)
            cp=cp+(shcon(0,i,imat)-shcon(0,i-1,imat))*
     &              (shcon(1,i,imat)+shcon(1,i-1,imat))/2.d0
         enddo
         cp=cp+(t1l-shcon(0,nshcon(imat),imat))*
     &           (shcon(1,nshcon(imat),imat))
         cp=cp/(t1l-physcon(1))
      else
         cp=shcon(1,id,imat)+
     &        (shcon(1,id+1,imat)-shcon(1,id,imat))*
     &        (t1l-shcon(0,id,imat))/
     &        (shcon(0,id+1,imat)-shcon(0,id,imat))
         cp=(t1l-shcon(0,id,imat))*(cp+shcon(1,id,imat))/2.d0
         do i=2,id
            cp=cp+(shcon(0,i,imat)-shcon(0,i-1,imat))*
     &              (shcon(1,i,imat)+shcon(1,i-1,imat))/2.d0
         enddo
         cp=cp+(shcon(0,1,imat)-physcon(1))*shcon(1,1,imat)
         cp=cp/(t1l-physcon(1))
      endif
!
      return
      end







