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
      subroutine prefilter(co,nodedesi,ndesi,xo,yo,zo,x,y,z,
     &   nx,ny,nz,objectset,filterrad)               
!
!     Filtering of sensitivities and assigning the filterradius      
!
      implicit none
!
      character*81 objectset(5,*)
!
      integer nodedesi(*),ndesi,m,nx(ndesi),istat,ny(ndesi),
     &   nz(ndesi),kflag
!
      real*8 xo(ndesi),yo(ndesi),zo(ndesi),filterrad,x(ndesi),
     &   y(ndesi),z(ndesi),co(3,*)
!   
!     Create set of designnodes and perform the sorting
!     needed for near3d_se
!
      do m=1,ndesi
         xo(m)=co(1,nodedesi(m))
         x(m)=xo(m)
         nx(m)=m
         yo(m)=co(2,nodedesi(m))
         y(m)=yo(m)
         ny(m)=m
         zo(m)=co(3,nodedesi(m))
         z(m)=zo(m)
         nz(m)=m
      enddo
      kflag=2
      call dsort(x,nx,ndesi,kflag)
      call dsort(y,ny,ndesi,kflag)
      call dsort(z,nz,ndesi,kflag)
!
!     assinging the filterradius
!
      read(objectset(2,1)(21:40),'(f20.0)',iostat=istat) filterrad     
!
!     For the GAUSS filter search in the 3sigma distance
! 
      if(objectset(2,1)(1:5).eq.'GAUSS') then
         filterrad=3*filterrad
      endif
!
      return        
      end
