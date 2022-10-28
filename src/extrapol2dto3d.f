
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
      subroutine extrapol2dto3d(dgdxglob,iponod2dto3d,ndesi,
     &   nodedesi,nobject,nk,xinterpol,nnodes,ipkon,lakon,kon,
     &   nodedesimpc,ndesimpc,ne,iponoel,inoel)
!
!     for 2D models extrapolate the results of the midnodes 
!     to the 2 symmetry planes
!
      implicit none
!
      character*8 lakon(*)
!
      integer iponod2dto3d(2,*),ndesi,nodedesi(*),nobject,
     &  node,nk,idesvar,ipkon(*),konl(26),ifaceq(2,20),
     &  ifacet(2,10),ifacew(2,15),kon(*),nnodes(nk),nope,ielem,
     &  start,nodedesimpc(*),ndesimpc,ne,indexe,i,j,l,ii,
     &  iponoel(*),inoel(2,*),nodecor1,nodecor2
!
      real*8 dgdxglob(2,nk,nobject),xinterpol(2,nk,nobject)
!
!     cornernodes next to the midnode for quadratic hex element
!
      data ifaceq /0,0,
     &             0,0,
     &             0,0,
     &             0,0,
     &             0,0,
     &             0,0,
     &             0,0,
     &             0,0,
     &             1,2,
     &             2,3,
     &             3,4,
     &             1,4,
     &             5,6,
     &             6,7,
     &             7,8,
     &             5,8,
     &             1,5,
     &             2,6,
     &             3,7,
     &             4,8/
!
!     cornernodes next to the midnode for quadratic wedge elements
!
      data ifacew /0,0,
     &             0,0,
     &             0,0,
     &             0,0,
     &             0,0,
     &             0,0,
     &             1,2,
     &             2,3,
     &             1,3,
     &             4,5,
     &             5,6,
     &             4,6,
     &             1,4,
     &             2,5,
     &             3,6/
!
!     Loop over all designvariables
!
      do i=1,ndesi
         idesvar=nodedesi(i)
         do j=1,2
            node=iponod2dto3d(j,idesvar)
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
