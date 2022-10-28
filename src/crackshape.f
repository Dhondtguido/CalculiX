      
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
      subroutine crackshape(nnfront,ifront,istartfront,iendfront,
     &     isubsurffront,angle,posfront,shape)
!
!     User Subroutine
!      
!     determine the shape factor for each crack front node
!
!     INPUT:
!
!     nnfront            number of crack fronts
!     ifront(i)          node number of front node i; the fronts are
!                        stored consecutively in ifront; within each
!                        front the nodes or stored by adjacency      
!     istartfront(j)     start of front j in field ifront
!     iendfront(j)       end of front j in field ifront
!     isubsurffront(j)   0: front j is a front belonging to a surface
!                           crack
!                        1: front j belongs to a subsurface crack
!     angle(j)           angle between tangents at start and end of      
!                        front j
!     posfront(i)        relative position of node ifront(i) within the
!                        the front it belongs to; 0<=posfront(i)<=1
!
!
!     OUTPUT (general):
!
!     shape(k,i)         shape factor for mode k (1<=k<=3) at front
!                        node ifront(i)
!     
      implicit none
!     
      integer i,j,k,nnfront,isubsurffront(*),istartfront(*),
     &     iendfront(*),ifront(*)
!     
      real*8 pi,shape(3,*),angle(*),posfront(*),s,twodpi,
     &     shape0,shapepi,shapeangle
!
      pi=4.d0*datan(1.d0)
      twodpi=2.d0/pi
!
      do i=1,nnfront
        if(isubsurffront(i).eq.1) then
          do j=istartfront(i),iendfront(i)
            do k=1,3
              shape(k,j)=twodpi
            enddo
          enddo
        else
          if(angle(i).lt.pi) then
            do j=istartfront(i),iendfront(i)
              s=2.d0*dabs(posfront(j)-0.5d0)
              shape0=1.12d0*(1.d0-s*s*0.02d0)
              shapepi=twodpi*(1.04+s*s*1.1d0)
              shapeangle=shape0*(1.d0-angle(i)/pi)+
     &             shapepi*angle(i)/pi
              do k=1,3
                shape(k,j)=shapeangle
              enddo
            enddo
          else
            do j=istartfront(i),iendfront(i)
              s=2.d0*dabs(posfront(j)-0.5d0)
              shapepi=twodpi*(1.04+s*s*1.1d0)
              do k=1,3
                shape(k,j)=shapepi
              enddo
            enddo
          endif
        endif
      enddo
!     
      return
      end

