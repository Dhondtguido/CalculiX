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
      subroutine rearrange(au,irow,icol,ndim,neq)
!
!     modifies the sparse storage mode for the iterative solver of
!     Ernst Rank (pcgsolver)
!
      implicit none
!
      integer irow(*),icol(*),ndim,i,neq,k,icr,istart,idiag,kflag
      real*8 au(*)
!
      kflag=2
!
      call isortiid(irow,icol,au,ndim,kflag)
!
      istart=1
      k=irow(1)
      icr=0
      idiag=0
!
      do i=1,ndim
         if(irow(i).eq.k) then
            icr=icr+1
            cycle
         else
            call isortid(icol(istart),au(istart),icr,kflag)
            icr=1
            istart=i
            k=irow(i)
            idiag=idiag+1
            irow(idiag)=i-1
         endif
      enddo
!
!     last row
!
      call isortid(icol(istart),au(istart),icr,kflag)
      idiag=idiag+1
      irow(idiag)=ndim
!
      return
      end
