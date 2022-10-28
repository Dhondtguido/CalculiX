!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine add_sm_fl(aub,adb,jq,irow,i,j,value,
     &  i0,i1)
!
!     stores the coefficient (i,j) with value "value" in the
!     fluid matrix
!
      implicit none
!
      integer jq(*),irow(*),i,j,ii,jj,ipointer,id,i0,i1
      real*8 adb(*),aub(*),value
!
      if(i.eq.j) then
         if(i0.eq.i1) then
            adb(i)=adb(i)+value
         else
            adb(i)=adb(i)+2.d0*value
         endif
         return
      elseif(i.gt.j) then
         ii=i
         jj=j
      else
         ii=j
         jj=i
      endif
!
      call nident(irow(jq(jj)),ii,jq(jj+1)-jq(jj),id)
!
      ipointer=jq(jj)+id-1
!
      if(irow(ipointer).ne.ii) then
         write(*,*) '*ERROR in add_sm_ei: coefficient should be 0'
c         write(*,*) i,j,ii,jj,ipointer,irow(ipointer)
      else
         aub(ipointer)=aub(ipointer)+value
      endif
!
      return
      end













