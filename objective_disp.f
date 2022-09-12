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
      subroutine objective_disp(nodeset,istartset,iendset,ialset,
     &  nk,idesvarc,iobject,mi,g0,nobject,vold,objectset)
!
!     calculates the sum of the square of the displacements of a node
!     set
!
      implicit none
!
      character*81 objectset(5,*)
!
      integer nk,istartset(*),iendset(*),ialset(*),nodeset,idir,
     &  idesvarc,iobject,mi(*),j,k,nobject,idesvar
!
      real*8 g0(nobject),vold(0:mi(2),*)
!
!
!
      idesvar=idesvarc+1
!
      g0(iobject)=0.d0
!
!     check for the existence of a set, else take the complete mesh
!
      if(nodeset.eq.0) then
         if(objectset(1,iobject)(1:8).eq.'ALL-DISP') then
            do j=1,nk           
               do idir=1,3
                  g0(iobject)=g0(iobject)+vold(idir,j)**2
               enddo      
            enddo
         elseif(objectset(1,iobject)(1:6).eq.'X-DISP') then
            do j=1,nk           
               g0(iobject)=g0(iobject)+vold(1,j)**2      
            enddo
         elseif(objectset(1,iobject)(1:6).eq.'Y-DISP') then
            do j=1,nk           
               g0(iobject)=g0(iobject)+vold(2,j)**2      
            enddo
         elseif(objectset(1,iobject)(1:6).eq.'Z-DISP') then
            do j=1,nk           
               g0(iobject)=g0(iobject)+vold(3,j)**2      
            enddo
         endif
      else
         do j=istartset(nodeset),iendset(nodeset)
            if(ialset(j).gt.0) then
               if(objectset(1,iobject)(1:8).eq.'ALL-DISP') then
                  do idir=1,3
                     g0(iobject)=g0(iobject)+vold(idir,ialset(j))**2
                  enddo    
               elseif(objectset(1,iobject)(1:6).eq.'X-DISP') then  
                  g0(iobject)=g0(iobject)+vold(1,ialset(j))**2
               elseif(objectset(1,iobject)(1:6).eq.'Y-DISP') then  
                  g0(iobject)=g0(iobject)+vold(2,ialset(j))**2
               elseif(objectset(1,iobject)(1:6).eq.'Z-DISP') then  
                  g0(iobject)=g0(iobject)+vold(3,ialset(j))**2
               endif
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  if(objectset(1,iobject)(1:8).eq.'ALL-DISP') then
                     do idir=1,3
                       g0(iobject)=g0(iobject)+vold(idir,k)**2
                     enddo     
                  elseif(objectset(1,iobject)(1:6).eq.'X-DISP') then     
                     g0(iobject)=g0(iobject)+vold(1,k)**2
                 elseif(objectset(1,iobject)(1:6).eq.'Y-DISP') then     
                     g0(iobject)=g0(iobject)+vold(2,k)**2
                 elseif(objectset(1,iobject)(1:6).eq.'Z-DISP') then     
                     g0(iobject)=g0(iobject)+vold(3,k)**2
                 endif
               enddo
            endif
         enddo
      endif
!
!     Euclidian length
!
      g0(iobject)=dsqrt(g0(iobject))
!
      return
      end
      
