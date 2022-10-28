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
      subroutine genratio(co,doubleglob,integerglob,nkold,nk,iprfn,
     &     konrfn,ratiorfn)
!     
!     create MPC's connecting the nodes in the refined tet mesh to
!     the nodes in the initial tet mesh; this is needed for the
!     temperature boundary conditions in the new mesh (t0 and t1)
!     
      implicit none
!     
      integer integerglob(*),nktet,netet,ne,nkon,
     &     nfaces,nfield,nselect,imastset,iselect(1),nterms,
     &     nelem,ialset(1),iprfn(*),konrfn(*),nkold,nk,
     &     iendset(1),istartset(1),konl(20),loopa,i,j
!     
      real*8 co(3,*),doubleglob(*),coords(3),ratio(20),value,
     &     dist,ratiorfn(*)
!     
      nktet=integerglob(1)
      netet=integerglob(2)
      ne=integerglob(3)
      nkon=integerglob(4)
      nfaces=integerglob(5)
      nfield=0
      nselect=0
      imastset=0
      loopa=8
!
      iprfn(1)=0
!     
      loop:do i=1,nk-nkold
!     
        do j=1,3
          coords(j)=co(j,nkold+i)
        enddo
!     
        call basis(doubleglob(1),doubleglob(netet+1),
     &       doubleglob(2*netet+1),
     &       doubleglob(3*netet+1),doubleglob(4*netet+1),
     &       doubleglob(5*netet+1),integerglob(6),integerglob(netet+6),
     &       integerglob(2*netet+6),doubleglob(6*netet+1),
     &       integerglob(3*netet+6),nktet,netet,
     &       doubleglob(4*nfaces+6*netet+1),nfield,
     &       doubleglob(13*nktet+4*nfaces+6*netet+1),
     &       integerglob(7*netet+6),integerglob(ne+7*netet+6),
     &       integerglob(2*ne+7*netet+6),
     &       integerglob(nkon+2*ne+7*netet+6),
     &       coords(1),coords(2),coords(3),value,ratio,iselect,nselect,
     &       istartset,iendset,ialset,imastset,
     &       integerglob(nkon+2*ne+8*netet+6),nterms,konl,nelem,loopa,
     &       dist)
!
!       coefficients of the equation
!
        do j=1,nterms
          ratiorfn(iprfn(i)+j)=ratio(j)
          konrfn(iprfn(i)+j)=konl(j)
        enddo
        iprfn(i+1)=iprfn(i)+nterms
!     
      enddo loop
!     
      return
      end

