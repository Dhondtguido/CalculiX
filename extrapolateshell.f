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
      subroutine extrapolateshell(yi,yn,ipkon,inum,kon,lakon,nfield,nk,
     &  ne,mi,ndim,orab,ielorien,co,iorienloc,cflag,
     &  ielmat,thicke,ielprop,prop,iflag)
!
!     extrapolates field values at the integration points to the 
!     nodes for user-defined shell elements
!
!     iflag=-1: NEG-position
!     iflag=0: MID-position
!     iflag=1: POS-position
!
      implicit none
!
      character*1 cflag
      character*8 lakon(*)
!
      integer ipkon(*),inum(*),kon(*),mi(*),ne,iorienloc,nfield,nk,i,j,
     &     ndim,ielorien(mi(3),*),ielmat(mi(3),*),ielprop(*),iflag
!
      real*8 yi(ndim,mi(1),*),yn(nfield,*),orab(7,*),co(3,*),prop(*),
     &  thicke(mi(3),*)
!
      do i=1,nk
         inum(i)=0
      enddo
!
      do i=1,nk
         do j=1,nfield
            yn(j,i)=0.d0
         enddo
      enddo
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
!
         if(lakon(i)(1:4).eq.'US45') then
            call extrapolateshell_us45(yi,yn,ipkon,inum,kon,lakon,
     &           nfield,nk,ne,mi,ndim,orab,ielorien,co,iorienloc,cflag,
     &           ielmat,thicke,ielprop,prop,i,iflag)
         elseif(lakon(i)(1:3).eq.'US3') then
            call extrapolateshell_us3(yi,yn,ipkon,inum,kon,lakon,
     &           nfield,nk,ne,mi,ndim,orab,ielorien,co,iorienloc,cflag,
     &           ielmat,thicke,ielprop,prop,i,iflag)
         else
            cycle
         endif
!
      enddo
!
!     taking the mean of nodal contributions coming from different
!     elements having the node in common
!
      do i=1,nk
         if(inum(i).gt.0) then
            do j=1,nfield
               yn(j,i)=yn(j,i)/inum(i)
            enddo
         endif
      enddo
!
      return
      end
