!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2021 Guido Dhondt
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
      subroutine auglag_inclusion(gmatrix,cvec,iacti,nacti,
     &     fric,atol,rtol,alglob,kitermax,
     &     auw,jqw,iroww,nslavs,al,alnew,r,omega)
!     
      implicit none
!
      include 'omp_lib.h'
!
      logical iscvg
!
      character*1 uplo
!     
      integer i,j,nacti,kitermax,iacti(*),incx,incy,
     &     icont,inorm,irow,jqw(*),iroww(*),nslavs
!     
!     al is in local contact coordinates!
!     alglob is in global coordinates
!     
      real*8 fric(*),atol,rtol,alglob(*),alsize,al(*),err,
     &     alnew(*),r(*),cvec(*),gmatrix(nacti,nacti),
     &     omega,value,auw(*),alpha,beta,altan,altanmax,ratio
!
      iscvg=.false.
!     
!     arguments for dsymv
!     
      uplo='U'
      alpha=1.d0
      incx=1
      beta=0.d0
      incy=1
!
      err=1.d30
      icont=0
!
!     determine the relaxation parameter
!
      call relaxval_al(r,gmatrix,nacti)
!
      do while((icont.le.kitermax).and.(.not.(iscvg)))
!        
!     G*lam via BLAS symmetric matrix vector multiplication
!     BLAS is part of the ARPACK library
!
        call dsymv(uplo,nacti,alpha,gmatrix,
     &       nacti,al,incx,beta,alnew,incy)
!
!     al-omega*r*(alnew+c)
!
        do i=1,nacti
          alnew(i)=al(i)-omega*r(i)*(alnew(i)+cvec(i))
        enddo
!
!     Projection operations for normal and tangential contact (friction)
!
        do i=1,nslavs
          if(iacti(3*i).eq.0) cycle
          inorm=iacti(3*i)-2
!     
!         F_normal
!     
          alnew(inorm)=max(alnew(inorm),0.d0)
!     
!         F_tangential
!     
          altan=sqrt(alnew(inorm+1)**2+alnew(inorm+2)**2 )
!     
!         F_normal * fric_coefficient
!     
          altanmax=fric(i)*alnew(inorm)
!     
!         Eval. stick or slip
!     
          if(altan.gt.altanmax) then
!     
!           if slip, adjust velocities according to Coulomb model
!     
            ratio=altanmax/altan
            alnew(inorm+1)=ratio*alnew(inorm+1)
            alnew(inorm+2)=ratio*alnew(inorm+2)
          endif
        enddo
!
!       determining the change in solution
!
        err=0.d0
        alsize=0.d0
!        
        do i=1,nacti
          err=err+(alnew(i)-al(i))**2
          alsize=alsize+al(i)**2
        enddo
!        
        err=dsqrt(err)
        alsize=dsqrt(alsize)
!
        do i=1,nacti
          al(i)=alnew(i)
        enddo
!
!       check for convergence          
!          
        if(err.le.(alsize*rtol+atol)) then
          iscvg=.true.
        endif
!
        icont=icont+1
      enddo
!
      if(icont.gt.kitermax)then
        write(*,*) '*WARNING!!: maximum iterations for massless'
        write(*,*) ' contact solution reached:' ,kitermax
        write(*,*) ' with error norm:', err
        ! call exit(201) ! TODO CMT should it be error or not? 
c      else
c        write(*,*) 'AL converged! NactiveDOF,it,err:',nacti,icont,err
      endif
!
!     Expansion of pk to global cys: alglob = Wb*al
!
      do i=1,3*nslavs
        if(iacti(i).ne.0)then
            do j=jqw(i),jqw(i+1)-1 ! each row 
              value=auw(j)
              irow=iroww(j)
              alglob(irow)=alglob(irow)+value*al(iacti(i)) !
            enddo
        endif
      enddo
!
      return
      end
