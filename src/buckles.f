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
      subroutine buckles(inpc,textpart,nmethod,mei,fei,
     &  nforc,nload,ithermal,iprestr,nbody,t0,t1,nk,iperturb,
     &  istep,istat,n,iline,ipol,inl,ipoinp,inp,isolver,ipoinpc,
     &  ier)
!
!     reading the input deck: *BUCKLE
!
      implicit none
!
      character*1 inpc(*)
      character*20 solver
      character*132 textpart(16)
!
      integer nmethod,mei(4),istep,istat,n,key,ncv,mxiter,
     &  nforc,nload,ithermal(*),iprestr,i,nk,iperturb(*),iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),nev,isolver,nbody,ipoinpc(0:*),ier
!
      real*8 fei(3),t0(*),t1(*),tol
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *BUCKLE: *BUCKLE can only be used'
         write(*,*) '  within a STEP'
         ier=1
         return
      endif
!
!     no heat transfer analysis
!
      if(ithermal(1).gt.1) then
         ithermal(1)=1
      endif
!
!     default solver
!
      solver='                    '
      if(isolver.eq.0) then
         solver(1:7)='SPOOLES'
      elseif(isolver.eq.2) then
         solver(1:16)='ITERATIVESCALING'
      elseif(isolver.eq.3) then
         solver(1:17)='ITERATIVECHOLESKY'
      elseif(isolver.eq.4) then
         solver(1:3)='SGI'
      elseif(isolver.eq.5) then
         solver(1:5)='TAUCS'
      elseif(isolver.eq.7) then
         solver(1:7)='PARDISO'
      elseif(isolver.eq.8) then
         solver(1:6)='PASTIX'
      endif
!
      do i=2,n
         if(textpart(i)(1:7).eq.'SOLVER=') then
            read(textpart(i)(8:27),'(a20)') solver
         else
            write(*,*) 
     &        '*WARNING reading *BUCKLE: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*BUCKLE%")
         endif
      enddo
!
      if(solver(1:7).eq.'SPOOLES') then
         isolver=0
      elseif(solver(1:16).eq.'ITERATIVESCALING') then
         write(*,*) '*WARNING in frequencies: the iterative scaling'
         write(*,*) '         procedure is not available for buckling'
         write(*,*) '         calculations; the default solver is used'
      elseif(solver(1:17).eq.'ITERATIVECHOLESKY') then
         write(*,*) '*WARNING in frequencies: the iterative scaling'
         write(*,*) '         procedure is not available for buckling'
         write(*,*) '         calculations; the default solver is used'
      elseif(solver(1:3).eq.'SGI') then
         isolver=4
      elseif(solver(1:5).eq.'TAUCS') then
         isolver=5
      elseif(solver(1:13).eq.'MATRIXSTORAGE') then
         isolver=isolver+10
      elseif(solver(1:7).eq.'PARDISO') then
         isolver=7
      elseif(solver(1:6).eq.'PASTIX') then
         isolver=8
      else
         write(*,*) '*WARNING reading *BUCKLE: unknown solver;'
         write(*,*) '         the default solver is used'
      endif
!
      if((isolver.eq.2).or.(isolver.eq.3)) then
         write(*,*) '*ERROR reading *BUCKLE: the default solver ',
     & solver
         write(*,*) '       cannot be used for buckling calculations '
         ier=1
         return
      endif
!
      nmethod=3
      if(iperturb(1).gt.1) iperturb(1)=0
      iperturb(2)=0
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*) '*ERROR reading *BUCKLE: definition not complete'
         write(*,*) '  '
         call inputerror(inpc,ipoinpc,iline,
     &        "*BUCKLE%",ier)
         return
      endif
      read(textpart(1)(1:10),'(i10)',iostat=istat) nev
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*BUCKLE%",ier)
         return
      endif
      if(nev.le.0) then
         write(*,*) '*ERROR reading *BUCKLE: less than 1 eigenvalue req
     &uested'
         ier=1
         return
      endif
      read(textpart(2)(1:20),'(f20.0)',iostat=istat) tol
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*BUCKLE%",ier)
         return
      endif
      if(tol.le.0.d0) then
         tol=1.d-2
      endif
      read(textpart(3)(1:10),'(i10)',iostat=istat) ncv
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*BUCKLE%",ier)
         return
      endif
      if(ncv.le.0) then
         ncv=4*nev
      endif
      ncv=ncv+nev
      read(textpart(4)(1:10),'(i10)',iostat=istat) mxiter
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*BUCKLE%",ier)
         return
      endif
      if(mxiter.le.0) then
         mxiter=1000
      endif
!
!     removing the natural boundary conditions
!
      nforc=0
      nload=0
      nbody=0
      iprestr=0
      if(ithermal(1).eq.1) then
         do i=1,nk
            t1(i)=t0(i)
         enddo
      endif
!
      mei(1)=nev
      mei(2)=ncv
      mei(3)=mxiter
      fei(1)=tol
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end


