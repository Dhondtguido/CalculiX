!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2025 Guido Dhondt
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
      subroutine extrapolate_ub31(yi,yn,ipkon,inum,kon,lakon,nfield,nk,
     &  ne,mi,ndim,orab,ielorien,co,iorienloc,cflag,
     &  vold,iforce,ielmat,thicke,ielprop,prop,i)
!
!     extrapolates field values at the integration points to the 
!     nodes for user element i of type UB31
!
      implicit none
!
      character*1 cflag
      character*8 lakon(*)
!
      integer ipkon(*),inum(*),kon(*),mi(*),ne,nfield,nk,i,ndim,
     &  iorienloc,ielorien(mi(3),*),ielmat(mi(3),*),ielprop(*),iforce
!
      real*8 yi(ndim,mi(1),*),yn(nfield,*),orab(7,*),co(3,*),prop(*),
     &  vold(0:mi(2),*),thicke(mi(3),*)
!
      integer indexe,j,k,node,ncomp
!
      ncomp=min(nfield,ndim)
      if(ncomp.ge.1) then
         indexe=ipkon(i)
         do j=1,2
            node=kon(indexe+j)
            do k=1,ncomp
               yn(k,node)=yn(k,node)+yi(k,1,i)
            enddo
            inum(node)=inum(node)+1
         enddo
      endif
!
!
      return
      end

