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
      subroutine statics(inpc,textpart,nmethod,iperturb,isolver,istep,
     &  istat,n,tinc,tper,tmin,tmax,idrct,iline,ipol,inl,ipoinp,inp,
     &  ithermal,cs,ics,tieset,istartset,
     &  iendset,ialset,ipompc,nodempc,coefmpc,nmpc,nmpc_,ikmpc,
     &  ilmpc,mpcfree,mcs,set,nset,labmpc,ipoinpc,iexpl,nef,ttime,
     &  iaxial,nelcon,nmat,tincf,ier)
!
!     reading the input deck: *STATIC
!
!     isolver=0: SPOOLES
!             2: iterative solver with diagonal scaling
!             3: iterative solver with Cholesky preconditioning
!             4: sgi solver
!             5: TAUCS
!             7: pardiso
!             8: pastix
!
!      iexpl==0:  structure:implicit, fluid:incompressible
!
      implicit none
!
      logical timereset
!
      character*1 inpc(*)
      character*20 labmpc(*),solver
      character*81 set(*),tieset(3,*)
      character*132 textpart(16)
!
      integer nmethod,iperturb(*),isolver,istep,istat,n,key,i,idrct,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),ithermal(*),ics(*),iexpl,
     &  istartset(*),iendset(*),ialset(*),ipompc(*),nodempc(3,*),
     &  nmpc,nmpc_,ikmpc(*),ilmpc(*),mpcfree,nset,mcs,ipoinpc(0:*),
     &  nef,iaxial,nelcon(2,*),nmat,ier
!
      real*8 tinc,tper,tmin,tmax,cs(17,*),coefmpc(*),ttime,tincf
!
      idrct=0
      tmin=0.d0
      tmax=0.d0
      timereset=.false.
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *STATIC: *STATIC can only be used'
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
         elseif((textpart(i)(1:6).eq.'DIRECT').and.
     &          (textpart(i)(1:9).ne.'DIRECT=NO')) then
            idrct=1
         elseif(textpart(i)(1:9).eq.'TIMERESET') then
            timereset=.true.
         elseif(textpart(i)(1:17).eq.'TOTALTIMEATSTART=') then
            read(textpart(i)(18:37),'(f20.0)',iostat=istat) ttime
         else
            write(*,*) 
     &        '*WARNING reading *STATIC: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*STATIC%")
         endif
      enddo
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
         write(*,*) '*WARNING reading *STATIC: unknown solver;'
         write(*,*) '         the default solver is used'
         write(*,*)
      endif
!
      nmethod=1
!
!     check for nodes on a cyclic symmetry axis
!
      if((mcs.eq.0).or.(iaxial.eq.180)) then
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
      else
         n=3
         textpart(2)='NMIN=0
     &
     &                '
         textpart(3)='NMAX=0
     &
     &                '
         nmethod=2
         call selectcyclicsymmetrymodess(inpc,textpart,cs,ics,tieset,
     &        istartset,
     &        iendset,ialset,ipompc,nodempc,coefmpc,nmpc,nmpc_,ikmpc,
     &        ilmpc,mpcfree,mcs,set,nset,labmpc,istep,istat,n,iline,
     &        ipol,inl,ipoinp,inp,nmethod,key,ipoinpc)
         nmethod=1
         do i=1,mcs
            cs(2,i)=-0.5d0
            cs(3,i)=-0.5d0
         enddo
      endif
!
      if((istat.lt.0).or.(key.eq.1)) then
         if((iperturb(1).ge.2).or.(nef.gt.0)) then
            write(*,*) '*WARNING reading *STATIC: a nonlinear analysis i
     &s requested'
            write(*,*) '         but no time increment nor step is speci
     &fied'
            write(*,*) '         the defaults (1,1) are used'
            write(*,*)
            tinc=1.d0
            tper=1.d0
c            tmin=1.d-5
            tmin=1.d-6
            tmax=1.d+30
            tincf=-1.d0
         else
            tper=1.d0
         endif
         if(timereset)ttime=ttime-tper
         return
      endif
!
      read(textpart(1)(1:20),'(f20.0)',iostat=istat) tinc
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*STATIC%",ier)
         return
      endif
      read(textpart(2)(1:20),'(f20.0)',iostat=istat) tper
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*STATIC%",ier)
         return
      endif
      read(textpart(3)(1:20),'(f20.0)',iostat=istat) tmin
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*STATIC%",ier)
         return
      endif
      read(textpart(4)(1:20),'(f20.0)',iostat=istat) tmax
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*STATIC%",ier)
         return
      endif
      read(textpart(5)(1:20),'(f20.0)',iostat=istat) tincf
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*STATIC%",ier)
         return
      endif
!
      if(tper.lt.0.d0) then
         write(*,*) '*ERROR reading *STATIC: step size is negative'
         ier=1
         return
      elseif(tper.le.0.d0) then
         tper=1.d0
      endif
      if(tinc.lt.0.d0) then
         write(*,*) '*ERROR reading *STATIC: initial increment size is n
     &egative'
         ier=1
         return
      elseif(tinc.le.0.d0) then
         tinc=tper
      endif
      if(tinc.gt.tper) then
         write(*,*) '*ERROR reading *STATIC: initial increment size exce
     &eds step size'
         ier=1
         return
      endif
!      
      if(idrct.ne.1) then
         if(dabs(tmin).lt.1.d-6*tper) then
            write(*,*) '*WARNING reading *STATIC:'
            write(*,*) '         the minimum increment ',tmin
            write(*,*) '         is smaller then 1.e-6 times the '
            write(*,*) '         step time;'
            write(*,*) '         the minimum increment is changed'
            write(*,*) '         to ',min(tinc,1.d-6*tper)
            write(*,*) '         which is the minimum of the initial'
            write(*,*) 
     &         '         increment time and 1.e-6 times the step time'
            write(*,*)
            tmin=min(tinc,1.d-6*tper)
         endif
         if(dabs(tmax).lt.1.d-10) then
            tmax=1.d+30
         endif
         if(tinc.gt.dabs(tmax)) then
            write(*,*) '*WARNING reading *STATIC:'
            write(*,*) '         the initial increment ',tinc
            write(*,*) '         exceeds the maximum increment ',
     &          tmax
            write(*,*) '         the initial increment is reduced'
            write(*,*) '         to the maximum value'
            write(*,*)
            tinc=dabs(tmax)
         endif
      endif
!
      if(timereset)ttime=ttime-tper
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

