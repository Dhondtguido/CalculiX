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
      subroutine materialdata_crack(crcon,ncrconst,ncrtem,t1l,
     &     crconloc)
!     
!     calculate the crack propagation increment
!     
      implicit none
!     
      integer ncrconst,ncrtem,k,id
!     
      real*8 crcon(0:ncrconst,*),t1l,crconloc(*)
!     
        call ident2(crcon,t1l,ncrtem,ncrconst+1,id)
        if(ncrtem.eq.0) then
          continue
        elseif(ncrtem.eq.1) then
          do k=1,ncrconst
            crconloc(k)=crcon(k,1)
          enddo
        elseif(id.eq.0) then
          do k=1,ncrconst
            crconloc(k)=crcon(k,1)
          enddo
        elseif(id.eq.ncrtem) then
          do k=1,ncrconst
            crconloc(k)=crcon(k,id)
          enddo
        else
          do k=1,ncrconst
            crconloc(k)=crcon(k,id)+
     &           (crcon(k,id+1)-crcon(k,id))*
     &           (t1l-crcon(0,id))/
     &           (crcon(0,id+1)-crcon(0,id))
          enddo
        endif
!     
      return
      end
      
