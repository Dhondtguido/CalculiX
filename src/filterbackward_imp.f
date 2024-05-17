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
      subroutine filterbackward_imp(ndesi,au,ad,aub,adb,jq,
     &   objectset)
!
!     Assembly of implicit filter matrix: 
!     au=filterrad*au+aub
!
      implicit none
!
      integer idesvar,ndesi,jq(*),jj,istat
!
      character*81 objectset(5,*)
!
      real*8 au(*),ad(*),aub(*),adb(*),filterrad 
!
!     assigning the filterradius
!
      read(objectset(2,1)(21:40),'(f20.0)',iostat=istat) filterrad     
!
      do idesvar=1,ndesi
         ad(idesvar)=filterrad**2*ad(idesvar)+adb(idesvar)
         do jj=jq(idesvar),jq(idesvar+1)-1
            au(jj)=filterrad**2*au(jj)+aub(jj)
         enddo
      enddo
!
      return
      end
      
