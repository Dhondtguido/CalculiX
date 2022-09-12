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
      subroutine smoothbadvertex(cotet,kontet,ipoeln,ieln,nbadnodes,
     &     ibadnodes,iponn,inn,iexternnode,ipoeled,ieled,iedgmid,iedtet)
!     
!     optimizing the position of bad vertex nodes by means of fminsi;
!     bad vertex nodes are vertex nodes on the free surface which
!     were not successfully projected in projectvertexnodes.f
!     
      implicit none
!     
      integer ibadnodes(*),nbadnodes,i,node,iexternnode(*),k,index,n,
     &     iponn(*),neigh,inn(2,*),ier,kontet(4,*),ipoeln(*),ieln(2,*),
     &     ipoeled(*),ieled(2,*),iedgmid(*),iedtet(6,*),iedge
!     
      real*8 cotet(3,*),cpycotet(3),x(3),fuvertex,eps(3),fmin,ref
!     
      external fuvertex
!     
      do i=1,nbadnodes
        neigh=ibadnodes(i)
!     
!     only subsurface neighbors are optimized
!     
        if(iexternnode(neigh).ne.0) cycle
!     
!     saving the original coordinates
!     
        do k=1,3
          cpycotet(k)=cotet(k,neigh)
          x(k)=cotet(k,neigh)
          eps(k)=0.d0
        enddo
!
!     starting function value (not really necessary, just in
!     case one want to print this value)
!
        n=3
        fmin=fuvertex(n,x,cotet,kontet,ipoeln,ieln,neigh,iedge,
     &     ipoeled,ieled,iedgmid,iedtet)
!     
!     calling the optimizer
!
        ref=fmin
        do
          ier=0
          call fminsirefine(n,x,fuvertex,eps,fmin,ier,cotet,
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
