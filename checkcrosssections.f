     
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
      subroutine checkcrosssections(co,doubleglob,integerglob,stress,
     &     nnfront,ifront,ifrontrel,costruc,temp,nstep,istartfront,
     &     iendfront)
!     
!     adjusting the projection on the structure of the external front
!     nodes 
!     
      implicit none
!     
      integer ifront(*),ifrontrel(*),nnfront,
     &     i,j,k,l,m,istartfront(*),iendfront(*),
     &     nstep,nktet,netet,integerglob(*),ne,nkon,
     &     nfaces,nfield,i1,nodeext,nodeextrel,nodeint,iselect(7),
     &     nselect,istartset(1),iendset(1),imastset,ialset(1),
     &     nterms,konl(20),nelem,loopa,iconstant
!     
      real*8 stress(6,nstep,*),costruc(3,*),
     &     temp(nstep,*),co(3,*),doubleglob(*),v1(3),v2(3),alambda,
     &     dist,cosang,coords(3),value(7),ratio(20)
!     
      nktet=integerglob(1)
      netet=integerglob(2)
      ne=integerglob(3)
      nkon=integerglob(4)
      nfaces=integerglob(5)
      nfield=13
!
!     select the stresses for interpolation
!
      nselect=7
      iselect(1)=1
      do k=2,7
         iselect(k)=k+3
      enddo
!
      imastset=0
      loopa=8
!
!     check the cross section of the crack fronts with the external
!     surface of the structure: first node in the fronts
!
      do i=1,nnfront
        i1=istartfront(i)
        nodeext=ifront(i1)
        nodeextrel=ifrontrel(i1)
        nodeint=ifront(i1+1)
        do k=1,3
          v1(k)=co(k,nodeint)-co(k,nodeext)
          v2(k)=costruc(k,nodeextrel)-co(k,nodeext)
        enddo
        alambda=dsqrt(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))
!
!       internal and external node are very close together
!
        if(alambda.lt.1.d-10) cycle
!
!       unit vector connecting external with internal node on front
!
        do k=1,3
          v1(k)=v1(k)/alambda
        enddo
!
!       distance between external node and its projection
!
        dist=dsqrt(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))
!
!       external node is very close to its projection
!
        if(dist.lt.1.d-10) cycle
!
!       cosine of angle between front direction and projection
!       direction
!
        cosang=(v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))/dist
!
!       do nothing if angle is small enough (roughly smaller than
!       25 degrees)
!
        if(cosang.gt.0.9d0) cycle
!
!       approach external node from internal node in steps of 
!       one tenth of the distance between both
!
        alambda=alambda/10.d0
!     
        do m=1,10
!     
          do k=1,3
            coords(k)=co(k,nodeint)-m*alambda*v1(k)
          enddo
!     
          call basis(doubleglob(1),doubleglob(netet+1),
     &         doubleglob(2*netet+1),
     &         doubleglob(3*netet+1),doubleglob(4*netet+1),
     &         doubleglob(5*netet+1),integerglob(6),
     &         integerglob(netet+6),
     &         integerglob(2*netet+6),doubleglob(6*netet+1),
     &         integerglob(3*netet+6),nktet,netet,
     &         doubleglob(4*nfaces+6*netet+1),nfield,
     &         doubleglob(nstep*13*nktet+4*nfaces+6*netet+1),
     &         integerglob(7*netet+6),integerglob(ne+7*netet+6),
     &         integerglob(2*ne+7*netet+6),
     &         integerglob(nkon+2*ne+7*netet+6),
     &         coords(1),coords(2),coords(3),value,ratio,iselect,
     &         nselect,
     &         istartset,iendset,ialset,imastset,
     &         integerglob(nkon+2*ne+8*netet+6),nterms,konl,nelem,loopa,
     &         dist)
!
!         exit as soon as the node is external
!
          if(dist.gt.1.d-6) exit
!
        enddo
!
!     store the coordinates of the node; these may not be the same
!     as in co due to the projection of external nodes onto the
!     structure        
!
        do j=1,3
          costruc(j,nodeextrel)=coords(j)
        enddo
!
        temp(1,nodeextrel)=value(1)
        do j=1,6
          stress(j,1,nodeextrel)=value(j+1)
        enddo
!
!       interpolating the values of the subsequent steps
!
        do l=2,nstep
          iconstant=13*(l-1)*nktet+4*nfaces+6*netet
          do k=1,nselect
            m=iselect(k)
            value(k)=0.d0
            do j=1,nterms
              value(k)=value(k)+ratio(j)*
     &             doubleglob(iconstant+(konl(j)-1)*13+m)
            enddo
          enddo
          temp(l,nodeextrel)=value(1)
          do j=1,6
            stress(j,l,nodeextrel)=value(j+1)
          enddo
        enddo
      enddo
!
!     check the cross section of the crack fronts with the external
!     surface of the structure: last node in the fronts
!
      do i=1,nnfront
        i1=iendfront(i)
        nodeext=ifront(i1)
        nodeextrel=ifrontrel(i1)
        nodeint=ifront(i1-1)
        do k=1,3
          v1(k)=co(k,nodeint)-co(k,nodeext)
          v2(k)=costruc(k,nodeextrel)-co(k,nodeext)
        enddo
        alambda=dsqrt(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))
!
!       internal and external node are very close together
!
        if(alambda.lt.1.d-10) cycle
!
!       unit vector connecting external with internal node on front
!
        do k=1,3
          v1(k)=v1(k)/alambda
        enddo
!
!       distance between external node and its projection
!
        dist=dsqrt(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))
!
!       external node is very close to its projection
!
        if(dist.lt.1.d-10) cycle
!
!       cosine of angle between front direction and projection
!       direction
!
        cosang=(v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))/dist
!
!       do nothing if angle is small enough (roughly smaller than
!       25 degrees)
!
        if(cosang.gt.0.9d0) cycle
!
!       approach external node from internal node in steps of 
!       one tenth of the distance between both
!
        alambda=alambda/10.d0
!     
        do m=1,10
!     
          do k=1,3
            coords(k)=co(k,nodeint)-m*alambda*v1(k)
          enddo
!     
          call basis(doubleglob(1),doubleglob(netet+1),
     &         doubleglob(2*netet+1),
     &         doubleglob(3*netet+1),doubleglob(4*netet+1),
     &         doubleglob(5*netet+1),integerglob(6),
     &         integerglob(netet+6),
     &         integerglob(2*netet+6),doubleglob(6*netet+1),
     &         integerglob(3*netet+6),nktet,netet,
     &         doubleglob(4*nfaces+6*netet+1),nfield,
     &         doubleglob(nstep*13*nktet+4*nfaces+6*netet+1),
     &         integerglob(7*netet+6),integerglob(ne+7*netet+6),
     &         integerglob(2*ne+7*netet+6),
     &         integerglob(nkon+2*ne+7*netet+6),
     &         coords(1),coords(2),coords(3),value,ratio,iselect,
     &         nselect,
     &         istartset,iendset,ialset,imastset,
     &         integerglob(nkon+2*ne+8*netet+6),nterms,konl,nelem,loopa,
     &         dist)
!
!         exit as soon as the node is external
!
          if(dist.gt.1.d-6) exit
!
        enddo
!
!     store the coordinates of the node; these may not be the same
!     als in co due to the projection of external nodes onto the
!     structure        
!
        do j=1,3
          costruc(j,nodeextrel)=coords(j)
        enddo
!
        temp(1,nodeextrel)=value(1)
        do j=1,6
          stress(j,1,nodeextrel)=value(j+1)
        enddo
!
!       interpolating the values of the subsequent steps
!
        do l=2,nstep
          iconstant=13*(l-1)*nktet+4*nfaces+6*netet
          do k=1,nselect
            m=iselect(k)
            value(k)=0.d0
            do j=1,nterms
              value(k)=value(k)+ratio(j)*
     &             doubleglob(iconstant+(konl(j)-1)*13+m)
            enddo
          enddo
          temp(l,nodeextrel)=value(1)
          do j=1,6
            stress(j,l,nodeextrel)=value(j+1)
          enddo
        enddo
      enddo
      return
      end

