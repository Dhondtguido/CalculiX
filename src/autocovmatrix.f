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
      subroutine autocovmatrix(co,ad,au,jqs,irows,ndesi,nodedesi,
     &     corrlen,randomval,irobustdesign)         
!     
!     calculates the values of the autocovariance matrix
!     
      implicit none
!     
      integer jqs(*),irows(*),ndesi,nodedesi(*),idof,j,jdof,node1,
     &     node2,irobustdesign(2)
!     
      real*8 co(3,*),ad(*),au(*),dist,corrlen,
     &     sigma1,sigma2,randomval(2,*)
!     
      if(irobustdesign(2).eq.1) then
!     case of homogeneous random field
        sigma1=randomval(2,1)         
        do idof=1,ndesi
          ad(idof)=sigma1*sigma1
          do j=jqs(idof),jqs(idof+1)-1
            jdof=irows(j)
            node1=nodedesi(idof)
            node2=nodedesi(jdof)
            dist=dsqrt((co(1,node1)-co(1,node2))**2+
     &           (co(2,node1)-co(2,node2))**2+
     &           (co(3,node1)-co(3,node2))**2)
!     
!     assign the value to the autocovariance matrix
!     
            au(j)=ad(idof)*dexp(-(dist/corrlen)**2)
          enddo
        enddo
      else
!        
!     case of inhomogeneous random field
!        
        do idof=1,ndesi
          node1=nodedesi(idof)
          sigma1=randomval(2,node1)
          ad(idof)=sigma1*sigma1
          do j=jqs(idof),jqs(idof+1)-1
            jdof=irows(j)              
            node2=nodedesi(jdof)
            sigma2=randomval(2,node2)
            dist=dsqrt((co(1,node1)-co(1,node2))**2+
     &           (co(2,node1)-co(2,node2))**2+
     &           (co(3,node1)-co(3,node2))**2)
!     
!     assign the value to the autocovariance matrix
!     
            au(j)=sigma1*sigma2*dexp(-(dist/corrlen)**2)
          enddo
        enddo
      endif
!     
      return        
      end
