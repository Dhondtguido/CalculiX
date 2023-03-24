!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine dynamics(inpc,textpart,nmethod,iperturb,tinc,tper,
     &     tmin,tmax,idrct,alpha,iexpl,isolver,istep,istat,n,iline,
     &     ipol,inl,ipoinp,inp,ithermal,ipoinpc,nef,ctrl,tincf,nener,
     &     ier,mortar)
!     
!     reading the input deck: *DYNAMIC
!     
!     isolver=0: SPOOLES
!     2: iterative solver with diagonal scaling
!     3: iterative solver with Cholesky preconditioning
!     4: sgi solver
!     5: TAUCS
!     7: pardiso
!     8: pastix
!     
!     iexpl=0:  structure:implicit, fluid:incompressible
!     iexpl=2:  structure:explicit, fluid:incompressible
!     
      implicit none
!     
      character*1 inpc(*)
      character*20 solver
      character*132 textpart(16)
!     
      integer nmethod,istep,istat,n,key,i,iperturb(*),idrct,iexpl,
     &     isolver,iline,ipol,inl,ipoinp(2,*),inp(3,*),ithermal(*),
     &     ipoinpc(0:*),nef,nener,ier,mortar
!     
      real*8 tinc,tper,tmin,tmax,alpha(*),ctrl(*),tincf
!     
      if(istep.lt.1) then
        write(*,*) '*ERROR reading *DYNAMIC: *DYNAMIC can only'
        write(*,*) '  be used within a STEP'
        ier=1
        return
      endif
!     
!     default is implicit
!     
      iexpl=0
!     
!     no heat transfer analysis
!     
      if(ithermal(1).gt.1) then
        ithermal(1)=1
      endif
!     
!     only nonlinear analysis allowed for this procedure
!     
      if(iperturb(1).lt.2) iperturb(1)=2
!     
!     default values
!     
      idrct=0
      alpha(1)=-0.05d0
      alpha(2)=0.5d0
      tmin=0.d0
      tmax=0.d0
      tincf=1.d-2
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
        if(textpart(i)(1:6).eq.'ALPHA=') then
          read(textpart(i)(7:26),'(f20.0)',iostat=istat) alpha(1)
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*DYNAMIC%",ier)
            return
          endif
          if(alpha(1).lt.-1.d0/3.d0) then
            write(*,*) '*WARNING reading *DYNAMIC: alpha is smaller'
            write(*,*) '  than -1/3 and is reset to -1/3'
            alpha(1)=-1.d0/3.d0
          elseif(alpha(1).gt.0.d0) then
            write(*,*) '*WARNING reading *DYNAMIC: alpha is greater'
            write(*,*) '  than 0 and is reset to 0'
            alpha(1)=0.d0
          endif
        elseif(textpart(i)(1:8).eq.'EXPLICIT') then
          iexpl=2
        elseif((textpart(i)(1:6).eq.'DIRECT').and.
     &         (textpart(i)(1:9).ne.'DIRECT=NO')) then
          idrct=1
        elseif(textpart(i)(1:7).eq.'SOLVER=') then
          read(textpart(i)(8:27),'(a20)') solver
        elseif(textpart(i)(1:18).eq.'RELATIVETOABSOLUTE') then
          alpha(2)=1.5
        else
          write(*,*) 
     &         '*WARNING reading *DYNAMIC: parameter not recognized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*DYNAMIC%")
        endif
      enddo
!     
      if((iexpl.eq.0).and.(mortar.eq.-1)) then
        write(*,*) '*ERROR reading *DYNAMIC; massless contact is not'
        write(*,*) '       allowed in an implicit dynamic calculation' 
        ier=1
        return
      endif
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
        write(*,*) '*WARNING reading *DYNAMIC: unknown solver;'
        write(*,*) '         the default solver is used'
      endif
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
        if((iperturb(1).ge.2).or.(nef.gt.0)) then
          write(*,*)'*WARNING reading *DYNAMIC: a nonlinear analysis is
     &requested'
          write(*,*) '         but no time increment nor step is specifi
     &ed'
          write(*,*) '         the defaults (1,1) are used'
          tinc=1.d0
          tper=1.d0
c          tmin=1.d-5
          tmin=1.d-10
          tmax=1.d+30
          tincf=1.d-2
        endif
        nmethod=4
        return
      endif
!     
      read(textpart(1)(1:20),'(f20.0)',iostat=istat) tinc
      if(istat.gt.0) then
        call inputerror(inpc,ipoinpc,iline,
     &       "*DYNAMIC%",ier)
        return
      endif
      read(textpart(2)(1:20),'(f20.0)',iostat=istat) tper
      if(istat.gt.0) then
        call inputerror(inpc,ipoinpc,iline,
     &       "*DYNAMIC%",ier)
        return
      endif
      read(textpart(3)(1:20),'(f20.0)',iostat=istat) tmin
      if(istat.gt.0) then
        call inputerror(inpc,ipoinpc,iline,
     &       "*DYNAMIC%",ier)
        return
      endif
      read(textpart(4)(1:20),'(f20.0)',iostat=istat) tmax
      if(istat.gt.0) then
        call inputerror(inpc,ipoinpc,iline,
     &       "*DYNAMIC%",ier)
        return
      endif
      read(textpart(4)(1:20),'(f20.0)',iostat=istat) tincf
      if(istat.gt.0) then
        call inputerror(inpc,ipoinpc,iline,
     &       "*DYNAMIC%",ier)
        return
      endif
!     
      if(tper.le.0.d0) then
        write(*,*) '*ERROR reading *DYNAMIC: step size is negative'
        ier=1
        return
      endif
!
      if(idrct.ne.1) then
        if(tinc.gt.0.d0) then
          write(*,*) '*WARNING reading *DYNAMIC: the initial time'
          write(*,*) '         increment defined by the user will not'
          write(*,*) '         be used since the time increment is'
          write(*,*) '         determined automatically by the dynamic'
          write(*,*) '         procedure based on stability'
          write(*,*) '         considerations'
        endif
      else
        if(tinc.le.0.d0) then
          write(*,*)'*ERROR reading *DYNAMIC: initial increment size is 
     &negative'
          ier=1
          return
        endif
        if(tinc.gt.tper) then
          write(*,*)'*ERROR reading *DYNAMIC: initial increment size exc
     &eeds step size'
          ier=1
          return
        endif
      endif
!      
      if((nef.gt.0).and.(tincf.le.0.d0)) then
        write(*,*) '*WARNING reading *DYNAMIC: initial CFD increment siz
     &e is zero or negative; the default of 0.01 is taken'
        tincf=1.d-2
      endif
!     
      if(idrct.ne.1) then
        if(dabs(tmin).lt.1.d-10*tper) then
          tmin=min(tinc,1.e-10*tper)
        endif
        if(dabs(tmax).lt.1.d-10) then
          tmax=1.d+30
        endif
        if(dabs(tmax).gt.tper) then
          tmax=tper
        endif
        if(tinc.gt.dabs(tmax)) then
          write(*,*) '*WARNING reading *DYNAMIC:'
          write(*,*) '         the initial increment ',tinc
          write(*,*) '         exceeds the maximum increment ',
     &         tmax
          write(*,*) '         the initial increment is reduced'
          write(*,*) '         to the maximum value'
          tinc=dabs(tmax)
        endif
      else
        tmin=tinc
        tmax=tinc
      endif
!     
!     10 cutbacks allowed for dynamics (because of contact)
!     
      ctrl(8)=10.5d0
!     
      nmethod=4
!
!     for implicit dynamic calculations the energy balance is always
!     calculated
!     for explicit dynamic calculations it is only calculated if ENER
!     output was requested on a *EL PRINT, *EL FILE OR *ELEMENT OUTPUT
!     card
!
      if(iexpl.le.1) nener=1
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      return
      end

