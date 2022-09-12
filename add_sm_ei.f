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
      subroutine add_sm_ei(au,ad,aub,adb,jq,irow,i,j,value,valuem,
     &  i0,i1)
!
!     stores the stiffness coefficient (i,j) with value "value"
!     in the stiffness matrix stored in spare matrix format and 
!     the mass coefficient (i,j) with value "valuem" in the lumped 
!     mass matrix
!
      implicit none
!
      integer jq(*),irow(*),i,j,ii,jj,ipointer,id,i0,i1
!
      real*8 ad(*),au(*),adb(*),aub(*),value,valuem
!
!
!
      if(i.eq.j) then
         if(i0.eq.i1) then
            ad(i)=ad(i)+value
            adb(i)=adb(i)+valuem
         else
            ad(i)=ad(i)+2.d0*value
            adb(i)=adb(i)+2.d0*valuem
         endif
         return
      elseif(i.gt.j) then
         ii=i
         jj=j
      else
         ii=j
         jj=i
      endif
c      write(*,*) ii,jj,value,valuem
!
      call nident(irow(jq(jj)),ii,jq(jj+1)-jq(jj),id)
!
      ipointer=jq(jj)+id-1
!
      if(irow(ipointer).ne.ii) then
         write(*,*) '*ERROR in add_sm_ei: coefficient should be 0'
c         write(*,*) i,j,ii,jj,ipointer,irow(ipointer)
c         call exit(201)
      else
         au(ipointer)=au(ipointer)+value
         aub(ipointer)=aub(ipointer)+valuem
      endif
!
      return
      end













