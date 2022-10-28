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
      subroutine mafillfreq_em(ad,au,adb,aub,irow,jq,neq,
     &  adfreq,aufreq,irowfreq,iaux,jqfreq,icolfreq,neqfreq,
     &  nzsfreq,om,symmetryflag,inputformat,b,bfreq)
!
!     expanding K and M into a new matrix needed for steady state
!     harmonic electromagnetic field calculations. The matrix satisfies:
!
!     __                __
!     |                   |
!     | K        -Omega*M |
!     |                   |
!     | Omega*M         K |
!     |_                 _|
!
      implicit none
!
      integer irow(*),jq(*),neq,irowfreq(*),iaux(*),jqfreq(*),
     &  icolfreq(*),neqfreq,nzsfreq,i,j,k,kflag,symmetryflag,
     &  inputformat
!
      real*8 ad(*),au(*),adb(*),aub(*),adfreq(*),aufreq(*),om,b(*),
     &  bfreq(*)
!
!     non-symmetric matrix
!     diagonal terms are stored in adfreq, off-diagonal terms are
!     stored column by column in aufreq
!
      symmetryflag=2
      inputformat=3
!
!     row numbers are stored in irowfreq
!     column numbers are stored in iaux
!     values are stored in aufreq
!
      k=0
      do i=1,neq
!
!        diagonal K-values
!
         adfreq(i)=ad(i)
         adfreq(neq+i)=ad(i)
!
!        off-diagonal K-values
!
         do j=jq(i),jq(i+1)-1
!
!           upper left matrix, lower triangle
!
            k=k+1
            irowfreq(k)=irow(j)
            iaux(k)=i
            aufreq(k)=au(j)
!
!           upper left matrix, upper triangle
!
            k=k+1
            irowfreq(k)=i
            iaux(k)=irow(j)
            aufreq(k)=au(j)
!
!           lower right matrix, lower triangle
!
            k=k+1
            irowfreq(k)=neq+irow(j)
            iaux(k)=neq+i
            aufreq(k)=au(j)
!
!           lower right matrix, upper triangle
!
            k=k+1
            irowfreq(k)=neq+i
            iaux(k)=neq+irow(j)
            aufreq(k)=au(j)
         enddo
!
!        diagonal M-values
!
         k=k+1
         irowfreq(k)=neq+i
         iaux(k)=i
         aufreq(k)=om*adb(i)
         k=k+1
         irowfreq(k)=i
         iaux(k)=neq+i
         aufreq(k)=-om*adb(i)
!
!        off-diagonal M-values
!
         do j=jq(i),jq(i+1)-1
!
!           lower left matrix, lower triangle
!
            k=k+1
            irowfreq(k)=neq+irow(j)
            iaux(k)=i
            aufreq(k)=om*aub(j)
!
!           lower left matrix, upper triangle
!
            k=k+1
            irowfreq(k)=neq+i
            iaux(k)=irow(j)
            aufreq(k)=om*aub(j)
!
!           upper right matrix, lower triangle
!
            k=k+1
            irowfreq(k)=irow(j)
            iaux(k)=neq+i
            aufreq(k)=-om*aub(j)
!
!           upper right matrix, upper triangle
!
            k=k+1
            irowfreq(k)=i
            iaux(k)=neq+irow(j)
            aufreq(k)=-om*aub(j)
         enddo
         bfreq(i)=b(i)
      enddo
!
!     sorting the column numbers
!
      kflag=2
      call isortiid(iaux,irowfreq,aufreq,nzsfreq,kflag)
!
!     determining the location of the start of each new column      
!
      k=0
      do i=1,nzsfreq
         if(iaux(i).gt.k) then
            k=k+1
            jqfreq(k)=i
            cycle
         endif
      enddo
      jqfreq(neqfreq+1)=nzsfreq+1
!
!     determining the length of each column
!
      do i=1,neqfreq
         icolfreq(i)=jqfreq(i+1)-jqfreq(i)
      enddo
!
!     sorting each column
!
      do i=1,neqfreq
        call isortid(irowfreq(jqfreq(i)),aufreq(jqfreq(i)),
     &       jqfreq(i+1)-jqfreq(i),kflag)
      enddo
!
      return
      end
