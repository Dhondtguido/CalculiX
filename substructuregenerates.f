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
      subroutine substructuregenerates(inpc,textpart,nmethod,iperturb,
     &  isolver,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &  ithermal,ipoinpc,filab,ier)
!
!     reading the input deck: *SUBSTRUCTURE GENERATE
!
!     isolver=0: SPOOLES
!             7: pardiso
!             8: pastix
!
      implicit none
!
      character*1 inpc(*)
      character*20 solver
      character*87 filab(*)
      character*132 textpart(16)
!
      integer nmethod,iperturb(*),isolver,istep,istat,n,key,i,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),ithermal(*),ipoinpc(0:*),
     &  ier
!
      if((iperturb(1).eq.1).and.(istep.ge.1)) then
         write(*,*) '*ERROR reading *SUBSTRUCTURE GENERATE:'
         write(*,*) '       perturbation analysis is'
         write(*,*) '       not provided in a *SUBSTRUCTURE'
         write(*,*) '       GENERATE step.'
         ier=1
         return
      endif
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *SUBSTRUCTURE GENERATE:'
         write(*,*) '       *SUBSTRUCTURE GENERATE can only be used'
         write(*,*) '       within a STEP'
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
            write(*,*) '*WARNING reading *SUBSTRUCTURE GENERATE:'
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*SUBSTRUCTURE GENERATE%")
         endif
      enddo
!
      if(solver(1:7).eq.'SPOOLES') then
         isolver=0
      elseif(solver(1:7).eq.'PARDISO') then
         isolver=7
      elseif(solver(1:6).eq.'PASTIX') then
         isolver=8
      else
         write(*,*) '*ERROR reading *SUBSTRUCTURE GENERATE:'
         write(*,*) '       solver:',solver,'is not allowed.'
         write(*,*) '       please specify SPOOLES or PARDISO'
         ier=1
         return
      endif
!
      nmethod=11
      filab(5)(1:4)='RF  '
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

