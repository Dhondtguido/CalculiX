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
      subroutine viscos(inpc,textpart,nmethod,iperturb,isolver,istep,
     &  istat,n,tinc,tper,tmin,tmax,idrct,iline,ipol,inl,ipoinp,
     &  inp,ipoinpc,nelcon,nmat,ncmat_,ctrl,ier)
!
!     reading the input deck: *VISCO (provided for compatibility
!     reasons with ABAQUS)
!
!     isolver=0: SPOOLES
!             2: iterative solver with diagonal scaling
!             3: iterative solver with Cholesky preconditioning
!             4: sgi solver
!             5: TAUCS
!             7: pardiso
!             8: pastix
!
      implicit none
!
      character*1 inpc(*)
      character*20 solver
      character*132 textpart(16)
!
      integer nmethod,iperturb(*),isolver,istep,istat,n,key,i,idrct,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),ipoinpc(0:*),nelcon(2,*),
     &  nmat,ncmat_,ier
!
      real*8 tinc,tper,tmin,tmax,ctrl(*)
!
      idrct=0
      tmin=0.d0
      tmax=0.d0
!
      if((iperturb(1).eq.1).and.(istep.gt.1)) then
         write(*,*) '*ERROR reading *VISCO: perturbation analysis is'
         write(*,*) '       not provided in a *VISCO step. Perform'
         write(*,*) '       a genuine nonlinear geometric calculation'
         write(*,*) '       instead (parameter NLGEOM)'
         ier=1
         return
      endif
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *VISCO: *VISCO can only be used'
         write(*,*) '  within a STEP'
         ier=1
         return
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
         elseif((textpart(i)(1:6).eq.'DIRECT').and.
     &          (textpart(i)(1:9).ne.'DIRECT=NO')) then
            idrct=1
         elseif(textpart(i)(1:6).eq.'CETOL=') then
            read(textpart(i)(7:26),'(f20.0)',iostat=istat) ctrl(40)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*VISCO%",ier)
               return
            endif
         else
            write(*,*) 
     &        '*WARNING reading *VISCO: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*VISCO%")
         endif
      enddo
!
      if(ctrl(40).le.0.d0) then
         write(*,*) '*ERROR reading *VISCO:'
         write(*,*) '       no strictly positive value'
         write(*,*) '       for the parameter CETOL given'
         ier=1
         return
      endif
      write(*,*) 'INFO reading *VISCO:'
      write(*,*) '     although the specification of CETOL for a'
      write(*,*) '     *VISCO procedure is mandatory, it is only'
      write(*,*) '     used for creep in materials with a linear'
      write(*,*) '     elastic behavior which is ISOTROPIC'
      write(*,*)
!
      if(solver(1:7).eq.'SPOOLES') then
         isolver=0
      elseif(solver(1:16).eq.'ITERATIVESCALING') then
         isolver=2
      elseif(solver(1:17).eq.'ITERATIVECHOLESKY') then
         isolver=3
      elseif(solver(1:3).eq.'SGI') then
         isolver=4
      elseif(solver(1:5).eq.'TAUCS') then
         isolver=5
      elseif(solver(1:7).eq.'PARDISO') then
         isolver=7
      elseif(solver(1:6).eq.'PASTIX') then
         isolver=8
      else
         write(*,*) '*WARNING reading *VISCO: unknown solver;'
         write(*,*) '         the default solver is used'
      endif
!
      nmethod=-1
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         if(iperturb(1).ge.2) then
            write(*,*) '*WARNING reading *VISCO: a nonlinear analysis is
     & requested'
            write(*,*) '         but no time increment nor step is speci
     &fied'
            write(*,*) '         the defaults (1,1) are used'
            tinc=1.d0
            tper=1.d0
            tmin=1.d-5
            tmax=1.d+30
         endif
         return
      endif
!
      read(textpart(1)(1:20),'(f20.0)',iostat=istat) tinc
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*VISCO%",ier)
         return
      endif
      read(textpart(2)(1:20),'(f20.0)',iostat=istat) tper
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*VISCO%",ier)
         return
      endif
      read(textpart(3)(1:20),'(f20.0)',iostat=istat) tmin
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*VISCO%",ier)
         return
      endif
      read(textpart(4)(1:20),'(f20.0)',iostat=istat) tmax
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*VISCO%",ier)
         return
      endif
!
      if(tinc.le.0.d0) then
         write(*,*) '*ERROR reading *VISCO: initial increment size is ne
     &gative'
      endif
      if(tper.le.0.d0) then
         write(*,*) '*ERROR reading *VISCO: step size is negative'
      endif
      if(tinc.gt.tper) then
         write(*,*) '*ERROR reading *VISCO: initial increment size excee
     &ds step size'
      endif
!      
      if(idrct.ne.1) then
         if(dabs(tmin).lt.1.d-6*tper) then
            tmin=min(tinc,1.d-6*tper)
         endif
         if(dabs(tmax).lt.1.d-10) then
            tmax=1.d+30
         if(tinc.gt.dabs(tmax)) then
            write(*,*) '*WARNING reading *VISCO:'
            write(*,*) '         the initial increment ',tinc
            write(*,*) '         exceeds the maximum increment ',
     &          tmax
            write(*,*) '         the initial increment is reduced'
            write(*,*) '         to the maximum value'
            tinc=dabs(tmax)
         endif
         endif
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

