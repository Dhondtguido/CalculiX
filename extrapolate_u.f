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
      subroutine extrapolate_u(yi,yn,ipkon,inum,kon,lakon,nfield,nk,
     &  ne,mi,ndim,orab,ielorien,co,iorienloc,cflag,
     &  vold,force,ielmat,thicke,ielprop,prop,i)
!
!     extrapolates nfield values at the integration points to the 
!     nodes for user element i
!
      implicit none
!
      logical force
!
      character*1 cflag
      character*8 lakon(*)
!
      integer ipkon(*),inum(*),kon(*),mi(*),ne,nfield,nk,i,ndim,
     &  iorienloc,ielorien(mi(3),*),ielmat(mi(3),*),ielprop(*)
!
      real*8 yi(ndim,mi(1),*),yn(nfield,*),orab(7,*),co(3,*),prop(*),
     &  vold(0:mi(2),*),thicke(mi(3),*)
!
      if(lakon(i)(2:3).eq.'1 ') then
         call extrapolate_u1(yi,yn,ipkon,inum,kon,lakon,nfield,nk,
     &        ne,mi,ndim,orab,ielorien,co,iorienloc,cflag,
     &        vold,force,ielmat,thicke,ielprop,prop,i)
      elseif(lakon(i)(2:3).eq.'S3') then
         call extrapolate_us3(yi,yn,ipkon,inum,kon,lakon,nfield,nk,
     &        ne,mi,ndim,orab,ielorien,co,iorienloc,cflag,
     &        vold,force,ielmat,thicke,ielprop,prop,i)
      elseif(lakon(i)(2:4).eq.'S45') then
         call extrapolate_us45(yi,yn,ipkon,inum,kon,lakon,nfield,nk,
     &        ne,mi,ndim,orab,ielorien,co,iorienloc,cflag,
     &        vold,force,ielmat,thicke,ielprop,prop,i)
      endif
!
      return
      end
