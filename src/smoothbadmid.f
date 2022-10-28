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
      subroutine smoothbadmid(cotet,kontet,ipoeln,ieln,nbadnodes,
     &     ibadnodes,iexternedg,ipoeled,ieled,iedgmid,iedtet)
!     
!     optimizing the position of bad midnodes by means of fminsi;
!     bad midnodes are midnodes on the free surface which
!     were not successfully projected in projectmidnodes.f
!     
      implicit none
!     
      integer ibadnodes(*),nbadnodes,i,k,n,
     &     neigh,ier,kontet(4,*),ipoeln(*),ieln(2,*),
     &     iedge,ipoeled(*),ieled(2,*),iedgmid(*),iedtet(6,*),
     &     iexternedg(*)
!     
      real*8 cotet(3,*),cpycotet(3),x(3),fumid,eps(3),fmin,ref
!     
      external fumid
!     
      do i=1,nbadnodes
        iedge=ibadnodes(i)
!     
!     only subsurface neighbors are optimized
!     
        if(iexternedg(iedge).ne.0) cycle
!     
!     saving the original coordinates
!     
        neigh=iedgmid(iedge)
!     
        do k=1,3
          cpycotet(k)=cotet(k,neigh)
          x(k)=cotet(k,neigh)
          eps(k)=0.d0
        enddo
!     
!     starting function value (not really necessary, just in
!     case one wants to print this value)
!     
        n=3
        fmin=fumid(n,x,cotet,kontet,ipoeln,ieln,neigh,iedge,
     &     ipoeled,ieled,iedgmid,iedtet)
!     
!     calling the optimizer
!
        ref=fmin
        do
          ier=0
          call fminsirefine(n,x,fumid,eps,fmin,ier,cotet,
     &         kontet,ipoeln,ieln,neigh,iedge,
     &         ipoeled,ieled,iedgmid,iedtet)
          if(ier.ne.0) then
            exit
          elseif(fmin.ge.ref) then
            exit
          else
            ref=fmin
          endif
        enddo
!     
!     restoring the original coordinates in case of error
!     
        if(ier.ne.0) then
          do k=1,3
            cotet(k,neigh)=cpycotet(k)
          enddo
        else
          do k=1,3
            cotet(k,neigh)=x(k)
          enddo
        endif
!     
      enddo
!     
      return
      end
