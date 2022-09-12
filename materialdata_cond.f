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
      subroutine materialdata_cond(imat,ntmat_,t1l,cocon,ncocon,cond)
!
      implicit none
!
!     determines the thermal conductivity
!
      integer imat,ntmat_,id,ncocon(2,*),ncoconst,seven
!
      real*8 t1l,cocon(0:6,ntmat_,*),cond
!
      seven=7
!
!     calculating the conductivity coefficients
!
      ncoconst=ncocon(1,imat)
      if(ncoconst.eq.0) then
         write(*,*) '*ERROR in materialdata_cond'
         write(*,*) 
     &        '       fluid conductivity is lacking'
         call exit(201)
      elseif(ncoconst.gt.1) then
         write(*,*) '*ERROR in materialdata_cond'
         write(*,*) 
     &        '       conductivity for fluids must be isotropic'
         call exit(201)
      endif
!     
      call ident2(cocon(0,1,imat),t1l,ncocon(2,imat),seven,id)
      if(ncocon(2,imat).eq.0) then
         cond=0.d0
         continue
      elseif(ncocon(2,imat).eq.1) then
         cond=cocon(1,1,imat)
      elseif(id.eq.0) then
         cond=cocon(1,1,imat)
      elseif(id.eq.ncocon(2,imat)) then
         cond=cocon(1,id,imat)
      else
         cond=(cocon(1,id,imat)+
     &        (cocon(1,id+1,imat)-cocon(1,id,imat))*
     &        (t1l-cocon(0,id,imat))/
     &        (cocon(0,id+1,imat)-cocon(0,id,imat)))
     &        
      endif
!     
      return
      end







