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
      subroutine mafillpbc(nef,au,ad,jq,irow,
     &  b,iatleastonepressurebc,nzs)
!
!     filling the lhs and rhs to calculate p
!
      implicit none
!
      integer i,nef,irow(*),jq(*),iatleastonepressurebc,nzs
!
      real*8 ad(*),au(*),b(*)
!     
!     at least one pressure bc is needed. If none is applied,
!     the last dof is set to 0
!     
!     a pressure bc is only recognized if not all velocity degrees of
!     freedom are prescribed on the same face
!     
c      write(*,*) 'mafillpbc', iatleastonepressurebc
      if(iatleastonepressurebc.eq.0) then
         ad(nef)=1.d0
         b(nef)=0.d0
         do i=2,nef
            if(jq(i)-1>0) then
               if(irow(jq(i)-1).eq.nef) then
                  au(jq(i)-1)=0.d0
               endif
            endif
         enddo
      endif
!     
c      do i=1,nzs
c         write(*,*) 'mafillp irow,au',i,au(i)
c      enddo
c      do i=1,nef
c         write(*,*) 'mafillp ad b',i,ad(i),b(i)
c      enddo
!     
      return
      end
