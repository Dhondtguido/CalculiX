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
      subroutine subspace(d,aa,bb,cc,alpham,betam,nev,xini,
     &  cd,cv,time,rwork,lrw,m,jout,rpar,bj,iwork,liw,iddebdf,bjp)
!
!     solves the linear dynamic equations mapped on the subspace
!     of the eigenvectors (only if there are dashpots in the
!     model)
!
      implicit none
!
      save time0
!
      integer nev,nev2,info(15),idid,lrw,iwork(*),liw,jout,id,
     &  iab,iaa,ibb,i,j,m,iddebdf
!
      real*8 d(*),aa(nev),bb(nev),cc(nev,*),alpham,betam,
     &  xini(*),cd(*),cv(*),time,time0,rtol,atol,rwork(*),rpar(*),
     &  bj(*),bjp(*)
!
      external df,djac
!
!
      nev2=2*nev
!
!     transferring fields into global field rpar
!     (needed for subroutine fd)
!     rpar contains (field, size): m+0.5, 1
!                                  alpham, 1
!                                  betam, 1
!                                  cc, nev**2
!                                  d, nev
!                                  time
!                                  aa(1)..aa(nev), nev
!                                  bb(1)..bb(nev), nev
!
      if(m.eq.1) then
         rpar(2)=alpham
         rpar(3)=betam
         do j=1,nev
            do i=1,nev
               rpar(3+(j-1)*nev+i)=cc(i,j)
            enddo
         enddo
         id=3+nev*nev
         do i=1,nev
            rpar(id+i)=d(i)
         enddo
!
!        copying the initial conditions for the system of first order
!        differential equations
!
         do i=1,nev
            xini(i)=cd(i)
            xini(nev+i)=cv(i)
         enddo
       endif
!
      iaa=3+nev*(1+nev)+1
      rpar(iaa)=time
      ibb=iaa+nev
      do i=1,nev
         rpar(iaa+i)=aa(i)
         rpar(ibb+i)=bb(i)
      enddo
!
      do i=1,3
         info(i)=0
      enddo
      info(4)=1
      info(5)=1
      info(6)=0
      rwork(1)=time
!
!     absolute and relative tolerance for dderkf
!
      rtol=1.d-5
      atol=1.d-3
!
      if(iddebdf.eq.0) then
         call ddeabm(df,nev2,time0,xini,time,info,rtol,atol,idid,rwork,
     &        lrw,iwork,liw,rpar,nev)
!
         if((idid.ne.2).and.(idid.ne.3)) then
            write(*,*) 
     &         '*WARNING in subspace: ddeabm did not converge properly'
            write(*,*) '         idid= ',idid
            write(*,*) '         switch to routine ddebdf'
            iddebdf=2
            return
         endif
      else
         call ddebdf(df,nev2,time0,xini,time,info,rtol,atol,idid,rwork,
     &        lrw,iwork,liw,rpar,nev,djac)
         if((idid.ne.2).and.(idid.ne.3)) then
            write(*,*) 
     &           '*ERROR in subspace: ddebdf did not converge properly'
            write(*,*) '       idid= ',idid
            call exit(201)
         endif
      endif
!
!     copying the solution into field bj
!
      do i=1,nev
         bj(i)=xini(i)
         bjp(i)=xini(nev+i)
      enddo
!
      return
      end
!
!     subroutine df expressing the first order derivative as a function
!     of time and the function itself
!
      subroutine df(x,u,uprime,rpar,nev)
!
      implicit none
!
      integer nev,i,j,id,iab,iaa,ibb,nev2p1,m
!
      real*8 rpar(*),x,u(*),uprime(*)
!
      iaa=4+nev*(nev+1)
      ibb=iaa+nev
      id=3+nev*nev
!
      do i=1,nev
         uprime(i)=u(nev+i)
         uprime(nev+i)=rpar(iaa+i)+x*rpar(ibb+i)
     &             -rpar(id+i)*rpar(id+i)*u(i)
     &             -(rpar(2)+rpar(3)*rpar(id+i)*rpar(id+i))*u(nev+i)
!
!        contribution of the dashpots
!   
         do j=1,nev
            uprime(nev+i)=uprime(nev+i)-rpar(3+(j-1)*nev+i)*u(nev+j)
         enddo
      enddo
!
      return
      end
!
!     subroutine djac 
!
      subroutine djac(x,u,pd,nrowpd,rpar,nev)
!
      implicit none
!
      integer nrowpd,nev,id,i,j
!
      real*8 rpar(*),x,u(*),pd(nrowpd,*)
!
      id=3+nev*nev
!
      do i=1,nev
         pd(i,nev+i)=1.d0
         pd(nev+i,i)=-rpar(id+i)*rpar(id+i)
         pd(nev+i,nev+i)=-(rpar(2)+rpar(3)*rpar(id+i)*rpar(id+i))
!
!        contribution of the dashpots
!   
         do j=1,nev
            pd(nev+i,nev+j)=pd(nev+i,nev+j)-rpar(3+(j-1)*nev+i)
         enddo
      enddo
!
      return
      end
