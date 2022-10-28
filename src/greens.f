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
      subroutine greens(inpc,textpart,nmethod,
     &  mei,iperturb,istep,istat,n,iline,ipol,inl,
     &  ipoinp,inp,ithermal,isolver,xboun,nboun,ipoinpc,
     &  ier)
!
!     reading the input deck: *GREEN
!
      implicit none
!
      character*1 inpc(*)
      character*20 solver
      character*132 textpart(16)
!
      integer nmethod,mei(4),istep,istat,iperturb(*),i,nboun,ier,
     &  n,key,iline,ipol,inl,ipoinp(2,*),inp(3,*),ithermal(*),isolver,
     &  ipoinpc(0:*)
!
      real*8 xboun(*)
!
      mei(4)=0
!
      if(istep.lt.1) then
         write(*,*) 
     &      '*ERROR reading *GREEN: *GREEN can only be used'
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
         elseif(textpart(i)(1:11).eq.'STORAGE=YES') then
            mei(4)=1
         else
            write(*,*) 
     &        '*WARNING reading *GREEN: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*GREEN%")
         endif
      enddo
!
      if(solver(1:7).eq.'SPOOLES') then
         isolver=0
      elseif(solver(1:16).eq.'ITERATIVESCALING') then
         write(*,*) '*WARNING reading *GREEN: the iterative scaling'
         write(*,*) '         procedure is not available for green'
         write(*,*) '         calculations; the default solver is used'
      elseif(solver(1:17).eq.'ITERATIVECHOLESKY') then
         write(*,*) '*WARNING reading *GREEN: the iterative scaling'
         write(*,*) '         procedure is not available for green'
         write(*,*) '         calculations; the default solver is used'
      elseif(solver(1:3).eq.'SGI') then
         isolver=4
      elseif(solver(1:5).eq.'TAUCS') then
         isolver=5
      elseif(solver(1:7).eq.'PARDISO') then
         isolver=7
      elseif(solver(1:6).eq.'PASTIX') then
         isolver=8
      else
         write(*,*) '*WARNING reading *GREEN: unknown solver;'
         write(*,*) '         the default solver is used'
      endif
!
      if((isolver.eq.2).or.(isolver.eq.3)) then
         write(*,*) '*ERROR reading *GREEN: the default solver ',
     & solver
         write(*,*) '       cannot be used for green calculations '
         ier=1
         return
      endif
!
      nmethod=13
      if(iperturb(1).gt.1) iperturb(1)=0
!
!     removing nonzero boundary conditions
!
      do i=1,nboun
         xboun(i)=0.d0
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

