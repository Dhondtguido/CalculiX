!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine mafillv2rhs(kon,ipkon,lakon,b2,v,nea,neb,mi,dtimef,
     &     ipvar,var,nk)
!     
!     filling the rhs b2 of the velocity equations (step 3)
!     
      implicit none
!     
      character*8 lakon(*)
!     
      integer kon(*),ipkon(*),nea,neb,mi(*),ipvar(*),i,j,k,
     &     node,indexe,nope,nk
!     
      real*8 b2(nk,3),v(nk,0:mi(2)),bb(3,8),dtimef,var(*)
!     
      do i=nea,neb
!     
        if(lakon(i)(1:1).ne.'F') cycle
        indexe=ipkon(i)
        if(lakon(i)(4:4).eq.'8') then
          nope=8
        elseif(lakon(i)(4:4).eq.'4') then
          nope=4
        elseif(lakon(i)(4:4).eq.'6') then
          nope=6
        else
          cycle
        endif
!     
        call e_c3d_v2rhs(kon(indexe+1),lakon(i),bb,i,v,dtimef,mi,
     &       ipvar,var,nk)
!     
        do j=1,nope
          node=kon(indexe+j)
          do k=1,3
            b2(node,k)=b2(node,k)+bb(k,j)
          enddo
        enddo
      enddo
c      write(*,*) 'mafillv2rhs '
c      do i=1,nk
c        write(*,*) i,(b2(i,j),j=1,3)
c      enddo
!     
      return
      end
