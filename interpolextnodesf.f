     
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
      subroutine interpolextnodesf(ibounnod,nbounnod,co,doubleglobf,
     &     integerglobf,stressf,tempf,nstepf,mode)
!     
!     interpolation of the hcf master file stress results for mode "mode"
!     to the boundary nodes of the crack shape(s)
!     
      implicit none
!     
      integer ibounnod(*),nbounnod,integerglobf(*),nktet,netet,ne,
     &     nkon,nfaces,nfield,nselect,imastset,iselect(7),nterms,
     &     nelem,ialset(1),iendset(1),istartset(1),konl(20),loopa,
     &     nstepf,mode,node,i,j,k
!     
      real*8 co(3,*),doubleglobf(*),coords(3),ratio(20),
     &     stressf(6,*),dist,tempf(*),value(7)
!     
      nktet=integerglobf(1)
      netet=integerglobf(2)
      ne=integerglobf(3)
      nkon=integerglobf(4)
      nfaces=integerglobf(5)
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
      loop:do i=1,nbounnod
        node =ibounnod(i)
c        write(*,*) 'interpolextnodes ',node
!     
        do j=1,3
          coords(j)=co(j,node)
        enddo
!     
        call basis(doubleglobf(1),doubleglobf(netet+1),
     &       doubleglobf(2*netet+1),
     &       doubleglobf(3*netet+1),doubleglobf(4*netet+1),
     &       doubleglobf(5*netet+1),integerglobf(6),
     &       integerglobf(netet+6),
     &       integerglobf(2*netet+6),doubleglobf(6*netet+1),
     &       integerglobf(3*netet+6),nktet,netet,
     &       doubleglobf(4*nfaces+6*netet+13*(mode-1)*nktet+1),nfield,
     &       doubleglobf(nstepf*13*nktet+4*nfaces+6*netet+1),
     &       integerglobf(7*netet+6),integerglobf(ne+7*netet+6),
     &       integerglobf(2*ne+7*netet+6),
     &       integerglobf(nkon+2*ne+7*netet+6),
     &       coords(1),coords(2),coords(3),value,ratio,iselect,
     &       nselect,
     &       istartset,iendset,ialset,imastset,
     &       integerglobf(nkon+2*ne+8*netet+6),nterms,konl,nelem,loopa,
     &       dist)
!
        tempf(i)=value(1)
        do j=1,6
          stressf(j,i)=value(j+1)
        enddo
!     
      enddo loop
!     
      return
      end

