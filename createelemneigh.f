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
      subroutine createelemneigh(nk,iponoel,inoel,istartnneigh,
     &   ialnneigh,icheckelems,istarteneigh,ialeneigh)
!
      implicit none
!
      integer nk,iponoel(*),inoel(2,*),istartnneigh(*),ialnneigh(*),
     &   ifree,index,i,j,ipos,na,nb,node,istarteneigh(*),ialeneigh(*),
     &   icheckelems(*)
!
!     determining all the elements to which the objective
!     nodes of the neighboring elements of a node i belong
!     They are stored in ialeneigh(istarteneigh(i))..
!     ...up to..... ialeneigh(istarteneigh(i+1)-1)
!
      ifree=1
      do i=1,nk
!
         istarteneigh(i)=ifree 
         index=iponoel(i)
         if(index.eq.0) cycle 
         na=istartnneigh(i)
         nb=istartnneigh(i+1)-1
!   
         do j=na,nb
!   
            node=ialnneigh(j)
            index=iponoel(node)
!      
            do
               if(index.eq.0) exit
               ipos=inoel(1,index)
               if(icheckelems(ipos).ne.i) then
                  ialeneigh(ifree)=inoel(1,index)
                  ifree=ifree+1
                  icheckelems(ipos)=i
               endif
               index=inoel(2,index)
            enddo
         enddo
      enddo   
      istarteneigh(nk+1)=ifree
!
      return
      end
