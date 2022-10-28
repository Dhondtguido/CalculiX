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
      subroutine uncouptempdisps(inpc,textpart,
     &  nmethod,iperturb,isolver,
     &  istep,istat,n,tinc,tper,tmin,tmax,idrct,ithermal,iline,ipol,
     &  inl,ipoinp,inp,ipoinpc,alpha,ctrl,ttime,nener,ier)
!
!     reading the input deck: *UNCOUPLED TEMPERATURE-DISPLACEMENT
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
      logical timereset
!
      character*1 inpc(*)
      character*20 solver
      character*132 textpart(16)
!
      integer nmethod,iperturb(*),isolver,istep,istat,n,key,i,idrct,
     &  ithermal(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),ipoinpc(0:*),
     &  nener,ier
!
      real*8 tinc,tper,tmin,tmax,alpha(*),ctrl(*),ttime
!
      idrct=0
      alpha(1)=-0.05d0
      tmin=0.d0
      tmax=0.d0
      nmethod=4
      timereset=.false.
!
      if(iperturb(1).eq.0) then
         iperturb(1)=2
      elseif((iperturb(1).eq.1).and.(istep.gt.1)) then
         write(*,*) 
     &       '*ERROR reading *UNCOUPLED TEMPERATURE-DISPLACEMENT:'
         write(*,*) '       perturbation analysis is not provided in a '
         write(*,*) '       *UNCOUPLED TEMPERATURE-DISPLACEMENT step.'
         ier=1
         return
      endif
!
      if(istep.lt.1) then
         write(*,*) 
     &       '*ERROR reading *UNCOUPLED TEMPERATURE-DISPLACEMENT:'
         write(*,*) '       *UNCOUPLED TEMPERATURE-DISPLACMENT '
         write(*,*) '       can only be used within a STEP'
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
         if(textpart(i)(1:6).eq.'ALPHA=') then
            read(textpart(i)(7:26),'(f20.0)',iostat=istat) alpha(1)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &            "*UNCOUPLED TEMPERATURE-DISPLACEMENT%",ier)
               return
            endif
            if(alpha(1).lt.-1.d0/3.d0) then
               write(*,*) '*WARNING in dynamics: alpha is smaller'
               write(*,*) '  than -1/3 and is reset to -1/3'
               alpha(1)=-1.d0/3.d0
            elseif(alpha(1).gt.0.d0) then
               write(*,*) '*WARNING in dynamics: alpha is greater'
               write(*,*) '  than 0 and is reset to 0'
               alpha(1)=0.d0
            endif
         elseif(textpart(i)(1:7).eq.'SOLVER=') then
            read(textpart(i)(8:27),'(a20)') solver
         elseif((textpart(i)(1:6).eq.'DIRECT').and.
     &          (textpart(i)(1:9).ne.'DIRECT=NO')) then
            idrct=1
         elseif(textpart(i)(1:11).eq.'STEADYSTATE') then
            nmethod=1
         elseif(textpart(i)(1:7).eq.'DELTMX=') then
            read(textpart(i)(8:27),'(f20.0)',iostat=istat) ctrl(27)
         elseif(textpart(i)(1:9).eq.'TIMERESET') then
            timereset=.true.
         elseif(textpart(i)(1:17).eq.'TOTALTIMEATSTART=') then
            read(textpart(i)(18:37),'(f20.0)',iostat=istat) ttime
         else
            write(*,*) 
     &        '*WARNING in uncouptempdisps: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*UNCOUPLED TEMPERATURE-DISPLACEMENT%")
         endif
      enddo
      if(nmethod.eq.1) ctrl(27)=1.d30
!
      if((ithermal(1).eq.0).and.(nmethod.ne.1).and.
     &   (nmethod.ne.2).and.(iperturb(1).ne.0)) then
         write(*,*) 
     &        '*ERROR reading *UNCOUPLED TEMPERATURE-DISPLACEMENT:'
         write(*,*) '       please define initial '
         write(*,*) '       conditions for the temperature'
         ier=1
         return
      else
         ithermal(1)=4
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
         write(*,*) 
     &       '*WARNING reading *UNCOUPLED TEMPERATURE-DISPLACEMENT:'
         write(*,*) '         unknown solver;'
         write(*,*) '         the default solver is used'
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         if(iperturb(1).ge.2) then
            write(*,*) 
     &       '*WARNING reading *UNCOUPLED TEMPERATURE-DISPLACEMENT:'
            write(*,*) '         a nonlinear analysis is requested'
            write(*,*) '         but no time increment nor step is speci
     &fied'
            write(*,*) '         the defaults (1,1) are used'
            tinc=1.d0
            tper=1.d0
            tmin=1.d-5
            tmax=1.d+30
         endif
         if(timereset)ttime=ttime-tper
         if(nmethod.eq.4) nener=1
         return
      endif
!
      read(textpart(1)(1:20),'(f20.0)',iostat=istat) tinc
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*UNCOUPLED TEMPERATURE-DISPLACEMENT%",ier)
         return
      endif
      read(textpart(2)(1:20),'(f20.0)',iostat=istat) tper
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*UNCOUPLED TEMPERATURE-DISPLACEMENT%",ier)
         return
      endif
      read(textpart(3)(1:20),'(f20.0)',iostat=istat) tmin
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*UNCOUPLED TEMPERATURE-DISPLACEMENT%",ier)
         return
      endif
      read(textpart(4)(1:20),'(f20.0)',iostat=istat) tmax
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*UNCOUPLED TEMPERATURE-DISPLACEMENT%",ier)
         return
      endif
!
      if(tinc.le.0.d0) then
         write(*,*) 
     &       '*ERROR reading *UNCOUPLED TEMPERATURE-DISPLACEMENT:'
         write(*,*) '       initial increment size is 
     &negative'
      endif
      if(tper.le.0.d0) then
         write(*,*) 
     &      '*ERROR reading *UNCOUPLED TEMPERATURE-DISPLACEMENT:'
         write(*,*) '       step size is negative'
      endif
      if(tinc.gt.tper) then
         write(*,*) 
     &     '*ERROR reading *UNCOUPLED TEMPERATURE-DISPLACEMENT:'
         write(*,*) '       initial increment size exc
     &eeds step size'
      endif
!      
      if(idrct.ne.1) then
         if(dabs(tmin).lt.1.d-6*tper) then
            tmin=min(tinc,1.d-6*tper)
         endif
         if(dabs(tmax).lt.1.d-10) then
            tmax=1.d+30
         endif
         if(tinc.gt.dabs(tmax)) then
            write(*,*) 
     &       '*WARNING reading *UNCOUPLED TEMPERATURE-DISPLACEMENT'
            write(*,*) '         the initial increment ',tinc
            write(*,*) '         exceeds the maximum increment ',
     &          tmax
            write(*,*) '         the initial increment is reduced'
            write(*,*) '         to the maximum value'
            tinc=dabs(tmax)
         endif
      endif
!
      if(timereset)ttime=ttime-tper
!
      if(nmethod.eq.4) nener=1
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end








