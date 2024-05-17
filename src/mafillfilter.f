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
      subroutine mafillfilter(adf,auf,jqf,irowf,ndesi,
     &   nodedesi,filterrad,co,weighting,objectset,xdesi,
     &   area)           
!     
!     calculates the filtervalues of the filter matrix
!     
      implicit none
!
      character*81 objectset(5,*)
!
      integer jqf(*),irowf(*),ndesi,nodedesi(*),kk,jj,
     &   inode1,inode2,ipos,actdir
!     
      real*8 auf(*),co(3,*),filterrad,dist,dx,dy,dz,distmin,
     &   weighting(*),scalar,xdesi(3,*),adf(*),filterval,
     &   area(*)
! 
!     Check if direction weighting is turned on
!
      if(objectset(2,1)(14:16).eq.'DIR') then
         actdir=1
      else
         actdir=0
      endif
! 
!     calculation of distmin
      if(actdir.eq.1) then
         distmin=dsqrt(xdesi(1,1)**2+xdesi(2,1)**2+xdesi(3,1)**2)
      endif
!
!     loop over all columns
      do kk=1,ndesi
         inode1=nodedesi(kk)
! 
!        entry on main diagonal
         weighting(kk)=weighting(kk)+area(kk)*1.0d0
         adf(kk)=1.0d0
!
!        loop over all rows of sub-diagonal
         do jj=jqf(kk),jqf(kk+1)-1
            ipos=irowf(jj)
            inode2=nodedesi(ipos)
!
            dx=co(1,inode1)-co(1,inode2)
            dy=co(2,inode1)-co(2,inode2)
            dz=co(3,inode1)-co(3,inode2)
            dist=dsqrt(dx**2+dy**2+dz**2)
            if(actdir.eq.1) then
               scalar=(xdesi(1,kk)*xdesi(1,ipos)
     &                +xdesi(2,kk)*xdesi(2,ipos)
     &                +xdesi(3,kk)*xdesi(3,ipos))/(distmin**2)
               if(scalar.lt.0.d0) then
                  scalar=0.d0
               endif
            else
               scalar=1.d0
            endif         
!  
!           Linear filter function
            filterval=max(0.d0,(filterrad-dist)/filterrad)
!
!           Entry caused by sub-diagonal
            weighting(kk)=weighting(kk)+area(ipos)*filterval
!
!           Entry caused by top-diagonal as filter matrix is symmetric
            weighting(ipos)=weighting(ipos)+area(kk)*filterval
!
!           Sparse matrix entry
            auf(jj)=filterval*scalar
!
         enddo
      enddo
!  
      return        
      end
