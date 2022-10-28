!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998 Guido Dhondt
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
      subroutine reinit_refine(kontet,ifac,ieln,netet_,newsize,
     &  ifatet,itetfa,iedg,ieled);
!
      implicit none
!
      integer kontet(4,*),ifac(4,*),ieln(2,*),netet_,newsize,i,j,
     &  ifatet(4,*),itetfa(2,*),iedg(3,*),ieled(2,*)
!
!     reinitialization of kontet
!
      kontet(4,netet_)=netet_+1
      do i=netet_+1,newsize
         do j=1,3
            kontet(j,i)=0
         enddo
         kontet(4,i)=i+1
      enddo
      kontet(4,newsize)=0
!
!     reinitialization of ifatet
!
      do i=netet_+1,newsize
         do j=1,4
            ifatet(j,i)=0
         enddo
      enddo
!
!     reinitialization of ipofa and ifac
!
      ifac(4,4*netet_)=4*netet_+1
      do i=4*netet_+1,4*newsize
         do j=1,3
            ifac(j,i)=0
         enddo
         ifac(4,i)=i+1
      enddo
      ifac(4,4*newsize)=0
!
!     reinitialization of itetfa
!
      do i=4*netet_+1,4*newsize
         do j=1,2
            itetfa(j,i)=0
         enddo
      enddo
!
!     reinitialization of ipoeln and ieln
!
      ieln(2,4*netet_)=4*netet_+1
      do i=4*netet_+1,4*newsize
         ieln(1,i)=0
         ieln(2,i)=i+1
      enddo
      ieln(2,4*newsize)=0
!
!     reinitialization of noded
!
      iedg(3,6*netet_)=6*netet_+1
      do i=6*netet_+1,6*newsize
         do j=1,2
            iedg(j,i)=0
         enddo
         iedg(3,i)=i+1
      enddo
      iedg(3,6*newsize)=0
!
!     reinitialization of ieled
!
      ieled(2,6*netet_)=6*netet_+1
      do i=6*netet_+1,6*newsize
         ieled(1,i)=0
         ieled(2,i)=i+1
      enddo
      ieled(2,6*newsize)=0
!
      return
      end
