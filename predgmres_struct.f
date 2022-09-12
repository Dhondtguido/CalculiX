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
!     y=A*x for real sparse matrices (symmetric and non-symmetric)
!
!     storage of the matrix in a:
!        - first the lower triangular terms
!        - then, if the matrix is non-symmetric, the upper triangular terms
!        - finally the diagonal terms
!
      subroutine predgmres_struct(n,b,x,nelt,ia,ja,a,isym,itol,tol,
     &  itmax,iter,err,ierr,iunit,sb,sx,rgwk,lrgw,igwk,ligw,rwork,iwork)
!
      implicit none
!
      integer n,nelt,ia(*),ja(*),isym,itol,itmax,iter,ierr,
     &  iunit,lrgw,igwk(*),ligw,iwork(*)
!
      real*8 b(*),x(*),a(*),tol,err,sb(*),sx(*),rgwk(*),
     &  rwork(*)
!
      external matvec_struct,msolve_struct
!
      itol=0
      tol=1.e-6
      itmax=0
      iunit=0
!
      igwk(1)=10
      igwk(2)=10
      igwk(3)=0
      igwk(4)=1
      igwk(5)=10
      ligw=20
!
      call dgmres(n,b,x,nelt,ia,ja,a,isym,matvec_struct,
     &  msolve_struct,itol,tol,itmax,
     &  iter,err,ierr,iunit,sb,sx,rgwk,lrgw,igwk,ligw,rwork,iwork)
!
      return
      end
