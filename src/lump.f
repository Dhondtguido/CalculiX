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
      subroutine lump(adb,aub,adl,irow,jq,neq)
!
!     lumping the matrix stored in adb,aub and storing the result
!     in adl
!
      implicit none
!
      integer irow(*),jq(*),neq,i,j,k
!
      real*8 adb(*),aub(*),adl(*)
!
      do i=1,neq
         adl(i)=adb(i)
      enddo
!
      do j=1,neq
         do k=jq(j),jq(j+1)-1
            i=irow(k)
            adl(i)=adl(i)+aub(k)
            adl(j)=adl(j)+aub(k)
         enddo
      enddo
!
!     change of meaning of adb and adl
!     first adb is replaced by adb-adl
!     then, adl is replaced by 1./adl
!
      do i=1,neq
         adb(i)=adb(i)-adl(i)
         adl(i)=1.d0/adl(i)
      enddo
!
      return
      end
