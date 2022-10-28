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
      subroutine materialdata_rho(rhcon,nrhcon,imat,rho,
     &  t1l,ntmat_,ithermal)
!
      implicit none
!
!     determines the density of the material
!
      integer nrhcon(*),imat,two,ntmat_,id,ithermal(*)
!
      real*8 rhcon(0:1,ntmat_,*),rho,t1l
!
      two=2
!
      if(ithermal(1).eq.0) then
         rho=rhcon(1,1,imat)
      else
         call ident2(rhcon(0,1,imat),t1l,nrhcon(imat),two,id)
         if(nrhcon(imat).eq.0) then
            continue
         elseif(nrhcon(imat).eq.1) then
            rho=rhcon(1,1,imat)
         elseif(id.eq.0) then
            rho=rhcon(1,1,imat)
         elseif(id.eq.nrhcon(imat)) then
            rho=rhcon(1,id,imat)
         else
            rho=rhcon(1,id,imat)+
     &           (rhcon(1,id+1,imat)-rhcon(1,id,imat))*
     &           (t1l-rhcon(0,id,imat))/
     &           (rhcon(0,id+1,imat)-rhcon(0,id,imat))
         endif
      endif
!
      return
      end
!     
