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
      subroutine coriolissolve(cc,nev,a,b,x,eiga,eigb,eigcorio,iter,d,
     &  temp,nevcomplex)
!
!     solves for the complex eigenfrequencies due to Coriolis 
!     forces
!
      implicit none
!
      logical wantx
!
      integer nev,iter(*),n,i,j,k,l,m,iii,ksmall,nincrement,
     &  nevcomplex,ninc,id,ipos,nevcomplexmax
!
      real*8 cc(nev,*),d(*),size,dmax
!
      complex*16 a(nev,*),b(nev,*),x(nev,*),eiga(*),eigb(*),eigcorio(*),
     &  temp(nev,*),deig,eig,eignew
!
      n=nev
      nincrement=10
      nevcomplexmax=nevcomplex
c      write(*,*) 'coriolissolve ',nevcomplexmax
      nevcomplex=0
      ninc=2
      ipos=0
!
!     determining the largest eigenvalue square
!
      dmax=0.d0
      do i=1,nev
         if(dabs(d(i)).gt.dmax) dmax=dabs(d(i))
      enddo
!
!     loop over all modes
!
      loop2: do i=1,nev
!
!        negative eigenvalues square are not taken into account as
!        initial conditions
!
         if(d(i).lt.0.d0) cycle
!
!        eigenvalues which are very small are not taken into account
!        as initial conditions
!
         if(dabs(d(i)).lt.dmax*1.d-6) cycle
!
!        initial values are taken at regular intervals between
!        the eigenvalues without Coriolis
!
         loop1: do m=1,nincrement
!
!           initial frequency
!     
            if(ipos.eq.0) then
!
!              for the lowest initial values: bias away from zero
!
               eig=dsqrt(d(i))*(1.d0,0.d0)*
     &              (1.d0-((nincrement-m)/(1.d0*nincrement))**2)
c               eig=dsqrt(d(i))*(1.d0,0.d0)*m/nincrement
c         write(*,*) 'coriolissolve ',i,m,eig
            else
               eig=(dsqrt(d(i-1))+
     &          (dsqrt(d(i))-dsqrt(d(i-1)))*m/nincrement)*(1.d0,0.d0)
            endif
!     
            iii=0
            do
               do k=1,nev
                  do l=1,nev
                     a(k,l)=eig*cc(k,l)*(0.d0,1.d0)
                     b(k,l)=-cc(k,l)*(0.d0,1.d0)
                  enddo
                  a(k,k)=a(k,k)-eig*eig+d(k)
                  b(k,k)=b(k,k)+2.d0*eig
               enddo
!     
!     solving for the complex eigenvalues
!     
               wantx=.false.
               call dlzhes(n,a,n,b,n,temp,n,wantx)
               call dlzit(n,a,n,b,n,temp,n,wantx,iter,eiga,eigb)
!     
               if(iter(1).eq.-1) then
                  write(*,*) '*ERROR in coriolissolve: fatal error'
                  write(*,*) '       in dlzit'
                  call exit(201)
               elseif(cdabs(eigb(1)).lt.1.d-10) then
                  write(*,*) '*WARNING in coriolissolve: eigenvalue'
                  write(*,*) '         out of bounds'
c                  cycle loop1
                  deig=0.d0
               else
!
!                 looking of the smallest change in eigenvalue
!
                  size=1.d30
                  do k=1,nev
                     deig=eiga(k)/eigb(k)
                     if(cdabs(deig).lt.size) then
                        ksmall=k
                        size=cdabs(deig)
                     endif
                  enddo
                  deig=eiga(ksmall)/eigb(ksmall)
               endif
               iii=iii+1
!     
               if((cdabs(deig).lt.1.d-3*cdabs(eig)).or.
     &              (cdabs(deig).lt.1.d-3).or.(iii.eq.100)) then
!
!                 eigenvalue found
!
                  id=0
                  eignew=eig+deig
!                  
!                 check whether eigenvalue is identical to one which
!                 was already found
!
c                  if((i.ne.1).or.(m.ne.1)) then
                  if(nevcomplex.gt.0) then
                     call ident2(eigcorio,eignew,nevcomplex,ninc,id)
                     if(id.ne.0) then
                        if((cdabs(eignew-eigcorio(id)).lt.
     &                       2.d-3*cdabs(eignew)).or.
     &                     (cdabs(eignew).lt.1.d-3)) cycle loop1
                     endif
                     if(id.ne.nevcomplex) then
                        if((cdabs(eignew-eigcorio(id+1)).lt.
     &                       2.d-3*cdabs(eignew)).or.
     &                     (cdabs(eignew).lt.1.d-3)) cycle loop1
                     endif
                  endif
!     
!     solve for the eigenvalues AND eigenvectors
!     
                  do k=1,nev
                     do l=1,nev
                        a(k,l)=eig*cc(k,l)*(0.d0,1.d0)
                        b(k,l)=-cc(k,l)*(0.d0,1.d0)
                     enddo
                     a(k,k)=a(k,k)-eig*eig+d(k)
                     b(k,k)=b(k,k)+2.d0*eig
                  enddo
                  wantx=.true.
                  call dlzhes(n,a,n,b,n,temp,n,wantx)
                  call dlzit(n,a,n,b,n,temp,n,wantx,iter,eiga,eigb)
!     
!     copying the eigenvector into the right place
!     
                  nevcomplex=nevcomplex+1
c                  write(*,*) 'coriolissolve ',i,m,nevcomplex,eignew
!
                  do k=nevcomplex,id+2,-1
                     eigcorio(k)=eigcorio(k-1)
                     do j=1,nev
                        x(j,k)=x(j,k-1)
                     enddo
                  enddo
!
                  do j=1,nev
                     x(j,id+1)=temp(j,ksmall)
                  enddo
                  eigcorio(id+1)=eignew
!
!                 check whether enough eigenvalues were found
!
                  if(nevcomplex.eq.nevcomplexmax) exit loop2
                  exit
               else
                  eig=eig+deig
               endif
!     
            enddo
         enddo loop1
         ipos=1
      enddo loop2
!
      if(nevcomplex.lt.nevcomplexmax) then
         write(*,*) '*WARNING in coriolissolve: not all requested'
         write(*,*) '         eigenvalues found'
         write(*,*) '         number of requested eigenvalues:',
     &       nevcomplexmax
         write(*,*) '         number of eigenvalues found:',nevcomplex
         write(*,*)
      endif
!     
      return
      end
      
