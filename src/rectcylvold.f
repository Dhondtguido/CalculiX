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
      subroutine rectcylvold(co,vold,cs,icntrl,
     &  mi,iznode,nznode,nsectors,nk)
!
!     special version of routine rectcyl for use in expand.c
!     transforms the reference displacements
!
!     icntrl=2:  rectangular to cylindrical coordinates for field
!                vold
!     icntrl=-2: cylindrical to rectangular coordinates for field
!                vold
!
!     nk: number of nodes in one segment
!     nkt: number of nodes in 360°
!
      implicit none
!
      integer i,j,icntrl,mi(*),iznode(*),nznode,nsectors,nk,
     &  ii,jj,node
!
      real*8 co(3,*),vold(0:mi(2),*),a(3,3),xr,xt,xz,cs(17,*),csab(7)
!
      do i=1,7
         csab(i)=cs(5+i,1)
      enddo
!
      if(icntrl.eq.2) then
         do ii=1,nznode
            i=iznode(ii)
            j=i
            call transformatrix(csab,co(1,i),a)
!
            xr=vold(1,j)*a(1,1)+vold(2,j)*a(2,1)+vold(3,j)*a(3,1)
            xt=vold(1,j)*a(1,2)+vold(2,j)*a(2,2)+vold(3,j)*a(3,2)
            xz=vold(1,j)*a(1,3)+vold(2,j)*a(2,3)+vold(3,j)*a(3,3)
            vold(1,j)=xr
            vold(2,j)=xt
            vold(3,j)=xz
!
         enddo
      elseif(icntrl.eq.-2) then
         do ii=1,nznode
            node=iznode(ii)
            do jj=1,nsectors
               i=node+(jj-1)*nk
               j=i
               call transformatrix(csab,co(1,i),a)
!     
               xr=vold(1,j)*a(1,1)+vold(2,j)*a(1,2)+vold(3,j)*a(1,3)
               xt=vold(1,j)*a(2,1)+vold(2,j)*a(2,2)+vold(3,j)*a(2,3)
               xz=vold(1,j)*a(3,1)+vold(2,j)*a(3,2)+vold(3,j)*a(3,3)
               vold(1,j)=xr
               vold(2,j)=xt
               vold(3,j)=xz
!     
            enddo
         enddo
      endif
!
      return
      end













