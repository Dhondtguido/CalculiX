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
      subroutine objective_modalstress(ndesi,neq,b,daldx,bfix,jqs,
     &     irows,df,iev,nev,z,dgduz,d,iobject,dgdx,dfm)
!
      implicit none
!
      integer i,j,idesvar,ndesi,neq(*),irows(*),iev,nev,iobject,jqs(*)
!      
      real*8 b(*),daldx(*),bfix(*),df(*),sensi,sum,z(neq(2),*),
     &     d(*),dgdx(ndesi,*),dgduz(*),dfm(*)
!
!     calculation of the modal stress sensitivity
!
      do idesvar=1,ndesi
!
!       dlambda_i/ds*M*U_i
!
        do j=1,neq(2)
          b(j)=daldx(idesvar)*bfix(j)
        enddo
!
!       dF/ds := (-dK/ds+lambda_i*dM/ds+dlambda_i/ds*M)*U_i
!
        do j=jqs(idesvar),jqs(idesvar+1)-1
          b(irows(j))=b(irows(j))+df(j)
        enddo
!
        sensi=0.d0
        do i=1,nev
          if(i.eq.iev+1) then
!
!           calculation of U_i^T*(dM/ds*U_i) = -2*betas
!
            sum=0.d0
            do j=jqs(idesvar),jqs(idesvar+1)-1
              sum=sum+z(irows(j),i)*dfm(j)
            enddo
            sensi=sensi-sum*dgduz(i)/2.d0
            cycle
          endif
!
!         calculation of U_i^T . dF/ds
!
          sum=0.d0
          do j=1,neq(2)
            sum=sum+z(j,i)*b(j)
          enddo
          sensi=sensi+sum*dgduz(i)/(d(iev+1)-d(i))
        enddo
        dgdx(idesvar,iobject)=dgdx(idesvar,iobject)+sensi
      enddo
!
      return
      end
