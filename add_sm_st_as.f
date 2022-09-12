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
      subroutine add_sm_st_as(au,ad,jq,irow,i,j,value,i0,i1,nzs)
!
!     stores the stiffness coefficient (i,j) with value "value"
!     in the stiffness matrix stored in spare matrix format
!
      implicit none
!
      integer jq(*),irow(*),i,j,ii,jj,ipointer,id,i0,i1,nzs(3),ioffset
      real*8 ad(*),au(*),value
!
      if(i.eq.j) then
         if(i0.eq.i1) then
            ad(i)=ad(i)+value
         else
            ad(i)=ad(i)+2.d0*value
         endif
         return
      elseif(i.gt.j) then
         ioffset=0
         ii=i
         jj=j
      else
         ioffset=nzs(3)
         ii=j
         jj=i
      endif
!
      call nident(irow(jq(jj)),ii,jq(jj+1)-jq(jj),id)
!
      ipointer=jq(jj)+id-1
!
      if(irow(ipointer).ne.ii) then
         write(*,*) '*ERROR in add_sm_st_as: coefficient should be 0'
         call exit(201)
      else
         ipointer=ipointer+ioffset
         au(ipointer)=au(ipointer)+value
      endif
!
      return
      end













