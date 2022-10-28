     
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
      subroutine interpolextnodes(ibounnod,nbounnod,co,doubleglob,
     &     integerglob,stress,ifront,nfront,ifrontrel,costruc,
     &     temp,nstep)
!     
!     interpolation of the master file stress results to the 
!     boundary nodes of the crack shape(s)
!     
      implicit none
!     
      integer ibounnod(*),nbounnod,integerglob(*),nktet,netet,ne,
     &     nkon,nfaces,nfield,nselect,imastset,iselect(7),nterms,
     &     nelem,ialset(1),ifront(*),nfront,ifrontrel(*),kflag,
     &     iendset(1),istartset(1),konl(20),loopa,nstep,
     &     node,i,j,k,l,m,iconstant
!     
      real*8 co(3,*),doubleglob(*),coords(3),ratio(20),costruc(3,*),
     &     stress(6,nstep,*),dist,temp(nstep,*),value(7)
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
      nfront=0
      loop:do i=1,nbounnod
        node =ibounnod(i)
c        write(*,*) 'interpolextnodes ',node
!     
        do j=1,3
          coords(j)=co(j,node)
        enddo
!     
        call basis(doubleglob(1),doubleglob(netet+1),
     &       doubleglob(2*netet+1),
     &       doubleglob(3*netet+1),doubleglob(4*netet+1),
     &       doubleglob(5*netet+1),integerglob(6),integerglob(netet+6),
     &       integerglob(2*netet+6),doubleglob(6*netet+1),
     &       integerglob(3*netet+6),nktet,netet,
     &       doubleglob(4*nfaces+6*netet+1),nfield,
     &       doubleglob(nstep*13*nktet+4*nfaces+6*netet+1),
     &       integerglob(7*netet+6),integerglob(ne+7*netet+6),
     &       integerglob(2*ne+7*netet+6),
     &       integerglob(nkon+2*ne+7*netet+6),
     &       coords(1),coords(2),coords(3),value,ratio,iselect,
     &       nselect,
     &       istartset,iendset,ialset,imastset,
     &       integerglob(nkon+2*ne+8*netet+6),nterms,konl,nelem,loopa,
     &       dist)
!
!     store the coordinates of the node; these may not be the same
!     als in co due to the projection of external nodes onto the
!     structure        
!
        do j=1,3
          costruc(j,i)=coords(j)
        enddo
!
        temp(1,i)=value(1)
        do j=1,6
          stress(j,1,i)=value(j+1)
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
          temp(l,i)=value(1)
          do j=1,6
            stress(j,l,i)=value(j+1)
          enddo
        enddo
!
!       check whether inside structure
!
        if(dist.lt.1.d-6) then
          nfront=nfront+1
          ifront(nfront)=node
          ifrontrel(nfront)=i
        endif
!     
      enddo loop
!     
!     sorting the front nodes in ascending order
!     
      kflag=2
      call isortii(ifront,ifrontrel,nfront,kflag)
!     
      return
      end

