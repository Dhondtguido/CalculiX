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
      subroutine filtermatrix(au,jq,irow,icol,ndesi,nodedesi,
     &     filterrad,co,nk,denominator,objectset,filterval,
     &     xdesi,distmin)           
!     
!     calculates the filtervalues of the filter matrix
!     
      implicit none
!
      character*81 objectset(5,*)
!
      integer jq(*),irow(*),icol(*),ndesi,nodedesi(*),i,ii,ll,
     &   kk,jj,inode1,inode2,nk,ipos,icolumn,irowstart,actdir
!     
      real*8 au(*),co(3,*),pi,filterrad,dist,dx,dy,dz,distmin,
     &   denominator(*),filterval(*),sigma,scalar,xdesi(3,*)
! 
!     Check if direction weighting is turned on
!
      if(objectset(2,1)(14:16).eq.'DIR') then
         actdir=1
      else
         actdir=0
      endif
! 
!     loop over all columns
      do kk=1,ndesi
         inode1=nodedesi(kk)
! 
!        loop over all rows
         do jj=jq(kk),jq(kk+1)-1
            ipos=irow(jj)
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
!           calculate filter filtervalue 
!
            if(objectset(2,1)(1:3).eq.'LIN') then
              filterval(ipos)=(1-dist/filterrad)*filterrad            
            elseif(objectset(2,1)(1:4).eq.'QUAD') then
              filterval(ipos)=(-1*(1+dist/filterrad)         
     &           *(-1+dist/filterrad))*filterrad
            elseif(objectset(2,1)(1:5).eq.'CUBIC') then
              filterval(ipos)=(2*(dist/filterrad)**3
     &           -3*(dist/filterrad)**2+1)*filterrad
            elseif(objectset(2,1)(1:5).eq.'GAUSS') then
               pi=4.d0*datan(1.d0)
               sigma=filterrad/3
               filterval(ipos)=1/(dsqrt(2*pi)*sigma)*exp(-(dist**2)
     &            /(2*sigma**2))
            endif
!
            denominator(kk)=denominator(kk)+filterval(ipos)
            filterval(ipos)=filterval(ipos)*scalar
         enddo
!
!        loop over all rows in column
         do jj=jq(kk),jq(kk+1)-1
            ipos=irow(jj)
            if(denominator(kk).ne.0d0) then
               au(jj)=filterval(ipos)/denominator(kk) 
            endif
         enddo    
      enddo
!  
!      do kk=1,ndesi
!         do jj=jq(kk),jq(kk+1)-1
!      write(5,*) kk,jj,au(jj)
!         enddo
!      enddo
      return        
      end
