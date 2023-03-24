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
!
!     center of gravity of the projection of the vertices for
!     visibility purposes
!     exact integration for one triangle: routine cubtri
!     if the surfaces are far enough away, one-point integration
!     is used
! 
      subroutine postview(ntr,sideload,nelemload,kontri,ntri,nloadtr,
     &     tenv,adview,auview,area,fenv,jqrad,irowrad,nzsrad)
!     
      implicit none
!
!     change following line if nlabel is increased
!     
      character*20 sideload(*)
!     
      integer ntr,iptri,irowrad(*),jqrad(*),nzsrad,nloadtr(*),
     &   nelemload(2,*),i,j,ntri,kontri(4,*)
!
      real*8 area(*),fenv(*),auview(*),adview(*),tenv(*),pi,
     &  totarea
!     
      pi=4.d0*datan(1.d0)
!     
!     division through total area and through pi
!     
      iptri=1
      do i=1,ntr
         totarea=0.d0
         do
            if(iptri.gt.ntri) exit
            if(kontri(4,iptri).ne.i) exit
            totarea=totarea+area(iptri)
            iptri=iptri+1
         enddo
c            
c         if(i.lt.ntr) then
c            do j=iptri(i),iptri(i+1)-1
c               totarea=totarea+area(j)
c            enddo
c         else
c            do j=iptri(i),ntri
c               totarea=totarea+area(j)
c            enddo
c         endif
         totarea=totarea*pi
!     
!     auxiliary field
!     
         fenv(i)=totarea
      enddo
!     
!     diagonal entries
!     
      do i=1,ntr
         adview(i)=adview(i)/fenv(i)
      enddo
!     
!     lower triangular entries
!     
      do i=1,nzsrad
         auview(i)=auview(i)/fenv(irowrad(i))
      enddo
!     
!     upper triangular entries
!     
      do i=1,ntr
         do j=nzsrad+jqrad(i),nzsrad+jqrad(i+1)-1
            auview(j)=auview(j)/fenv(i)
         enddo
      enddo
!     
!     checking whether the sum of the viewfactors does not
!     exceed 1 (summing the rows of the viewfactor matrix)
!     
!     diagonal entries
!     
      do i=1,ntr
         fenv(i)=adview(i)
      enddo
!     
!     lower triangular entries
      
      do i=1,nzsrad
         fenv(irowrad(i))=fenv(irowrad(i))+auview(i)
      enddo
!     
!     upper triangular entries
!     
      do i=1,ntr
         do j=nzsrad+jqrad(i),nzsrad+jqrad(i+1)-1
            fenv(i)=fenv(i)+auview(j)
         enddo
      enddo
!     
!     if fenv exceeds 1 or if the user explicitly asked
!     for a scaling to 1 (negative environmental temperature)
!     -> scale to 1 (only if fenv exceeds 0)
!     
!     diagonal entries
!     
      do i=1,ntr
         if((fenv(i).gt.0.d0).and.
     &        ((fenv(i).gt.1.d0).or.(tenv(i).lt.0))) 
     &        adview(i)=adview(i)/fenv(i)
      enddo
!     
!     lower triangular entries
!     
      do i=1,nzsrad
         if((fenv(irowrad(i)).gt.0.d0).and.
     &        ((fenv(irowrad(i)).gt.1.d0).or.(tenv(irowrad(i)).lt.0))) 
     &        auview(i)=auview(i)/fenv(irowrad(i))
      enddo
!     
!     upper triangular entries
!     
      do i=1,ntr
         if((fenv(i).gt.0.d0).and.
     &        ((fenv(i).gt.1.d0).or.(tenv(i).lt.0))) then
            do j=nzsrad+jqrad(i),nzsrad+jqrad(i+1)-1
               auview(j)=auview(j)/fenv(i)
            enddo
         endif
      enddo
!     
!     updating fenv accordingly
!     
      do i=1,ntr
         if((fenv(i).gt.1.d0).or.(tenv(i).lt.0)) then
            if(fenv(i).gt.0.d0) then
               fenv(i)=1.d0
            else
               write(*,*) '*WARNING in radmatrix: viewfactors'
               write(*,*) '         for 3D-face''',
     &              sideload(nloadtr(i)),''''
               write(*,*) '         of element',
     &              nelemload(1,nloadtr(i))
               write(*,*) '         cannot be scaled since they are'
               write(*,*) '         all zero'
               write(*,*)
            endif
         endif
c         write(*,*) 'postview ',i,nelemload(1,nloadtr(i)),
c     &        sideload(nloadtr(i)),fenv(i)
         fenv(i)=1.d0-fenv(i)
      enddo
!
      return
      end
