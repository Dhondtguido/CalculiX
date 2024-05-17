
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
      subroutine extrapol2dto3d(dgdxglob,nod2nd3rd,ndesi,nodedesi,
     &     nobject,nk)
!
!     for 2D models extrapolate the results of the midnodes 
!     to the 2 symmetry planes
!
      implicit none
!
      integer nod2nd3rd(2,*),ndesi,nodedesi(*),nobject,node,nk,
     &     idesvar,i,j,l
!
      real*8 dgdxglob(2,nk,*)
!
!     Loop over all designvariables
!
      do i=1,ndesi
         idesvar=nodedesi(i)
         do j=1,2
            node=nod2nd3rd(j,idesvar)
            if(node.eq.0) cycle
!
!           Loop over all objectives/constraints
!         
            do l=1,nobject   
               dgdxglob(1,node,l)=dgdxglob(1,idesvar,l)
               dgdxglob(2,node,l)=dgdxglob(2,idesvar,l)
            enddo
         enddo
      enddo
!
      return
      end        
