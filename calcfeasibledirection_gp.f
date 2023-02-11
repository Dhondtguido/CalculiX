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
      subroutine calcfeasibledirection_gp(ndesi,nodedesi,dgdxglob,
     &   nactive,nobject,nk,gradproj)         
!
!     calculates the projected gradient
!
      implicit none
!
      integer ndesi,nodedesi(*),nactive,nobject,nk,node,i 
!
      real*8 dgdxglob(2,nk,*),gradproj(3,*),dd1,dd2
!
      dd1=0.d0
      dd2=0.d0
!
!     calc projected gradient if nactive greater than 0
!     else, take sensitivity of objective function directly
!
      if(nactive.gt.0) then
         do i=1,ndesi
            node=nodedesi(i)
            gradproj(2,node)=dgdxglob(2,node,1)-gradproj(2,node)  
         enddo  
         write(*,*)
         write(*,*) '*INFO: at least 1 constraint active,'
         write(*,*) '       projected gradient calculated'   
         write(*,*)
      else
         do i=1,ndesi
            node=nodedesi(i)
            gradproj(1,node)=0.d0
            gradproj(2,node)=dgdxglob(2,node,1)  
         enddo
         write(*,*)
         write(*,*) '*INFO: no constraint active, sensitivity of'
         write(*,*) '      objective function taken as projected'
         write(*,*) '       gradient' 
         write(*,*)
      endif               
!
!     calc inifinity norms of the correction vector and 
!     gradient projection vector and normalize gradient fields
!
      do i=1,ndesi
         node=nodedesi(i)
         dd1=max(dd1,abs(gradproj(1,node)))
         dd2=max(dd2,abs(gradproj(2,node)))
      enddo
      if(dd1.le.0.d0) then
         dd1=1.d0
      endif
      if(dd2.le.0.d0) then
         dd2=1.d0
      endif
!
      do i=1,ndesi
         node=nodedesi(i)
         gradproj(1,node)=gradproj(1,node)/dd1
         gradproj(2,node)=gradproj(2,node)/dd2
         gradproj(3,node)=gradproj(2,node)
      enddo
!
      return        
      end

