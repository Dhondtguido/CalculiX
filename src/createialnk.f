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
      subroutine createialnk(nk,iponoel,inoel,istartnk,ialnk,ipkon)
!
      implicit none
!
      integer nk,iponoel(*),inoel(2,*),ipkon(*),ielem,
     &   istartnk(*),ialnk(*),ifree,index,i
!
!     determining the elements belonging to a node i.
!     They are stored in ialnk(istartnk(i))..
!     ...up to..... ialnk(istartnk(i+1)-1)
!
      ifree=1
      do i=1,nk
         istartnk(i)=ifree
         index=iponoel(i)
         do
            if(index.eq.0) exit
            ielem=inoel(1,index)
            if(ipkon(ielem).ge.0) then
               ialnk(ifree)=ielem
               ifree=ifree+1
            endif
            index=inoel(2,index)
         enddo
      enddo
      istartnk(nk+1)=ifree
!
      return
      end
