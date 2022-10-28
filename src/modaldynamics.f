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
      subroutine modaldynamics(inpc,textpart,nmethod,tinc,tper,iexpl,
     &  istep,istat,n,iline,ipol,inl,ipoinp,inp,iperturb,isolver,
     &  cs,mcs,ipoinpc,idrct,ctrl,tmin,tmax,nforc,nload,nbody,iprestr,
     &  t0,t1,ithermal,nk,vold,veold,xmodal,set,nset,mi,cyclicsymmetry,
     &  ier)
!
!     reading the input deck: *MODAL DYNAMIC
!
      implicit none
!
      logical steadystate,nodalset
!
      character*1 inpc(*)
      character*20 solver
      character*81 set(*),noset
      character*132 textpart(16)
!
      integer nmethod,istep,istat,n,key,iexpl,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),iperturb(*),isolver,i,mcs,ipoinpc(0:*),
     &  idrct,nforc,nload,nbody,iprestr,ithermal(*),j,nk,ipos,nset,
     &  mi(*),cyclicsymmetry,ier
!
      real*8 tinc,tper,cs(17,*),ctrl(*),tmin,tmax,t0(*),t1(*),
     &  vold(0:mi(2),*),veold(0:mi(2),*),xmodal(*)
!
      iexpl=0
c      iperturb(1)=0
      iperturb(2)=0
      idrct=1
      tmin=0.d0
      tmax=0.d0
      steadystate=.false.
      if((mcs.ne.0).and.(cs(2,1).ge.0.d0)) then
         cyclicsymmetry=1
      endif
      nodalset=.false.
!
      if(istep.lt.1) then
         write(*,*) 
     &      '*ERROR reading *MODAL DYNAMIC: *MODAL DYNAMIC can only'
         write(*,*) '  be used within a STEP'
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
         elseif(textpart(i)(1:9).eq.'DIRECT=NO') then
            idrct=0
         elseif(textpart(i)(1:7).eq.'DELTMX=') then
            read(textpart(i)(8:27),'(f20.0)',iostat=istat) ctrl(27)
         elseif(textpart(i)(1:11).eq.'STEADYSTATE') then
            steadystate=.true.
         else
            write(*,*) 
     &      '*WARNING reading *MODAL DYNAMIC: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*MODAL DYNAMIC%")
         endif
      enddo
!
      if(solver(1:7).eq.'SPOOLES') then
         isolver=0
      elseif(solver(1:16).eq.'ITERATIVESCALING') then
         write(*,*) 
     &      '*WARNING reading *MODAL DYNAMIC: the iterative scaling'
         write(*,*) '         procedure is not available for modal'
         write(*,*) '         dynamic calculations; the default solver'
         write(*,*) '         is used'
      elseif(solver(1:17).eq.'ITERATIVECHOLESKY') then
         write(*,*) 
     &       '*WARNING reading *MODAL DYNAMIC: the iterative scaling'
         write(*,*) '         procedure is not available for modal'
         write(*,*) '         dynamic calculations; the default solver'
         write(*,*) '         is used'
      elseif(solver(1:3).eq.'SGI') then
         isolver=4
      elseif(solver(1:5).eq.'TAUCS') then
         isolver=5
      elseif(solver(1:7).eq.'PARDISO') then
         isolver=7
      elseif(solver(1:6).eq.'PARDISO') then
         isolver=8
      else
         write(*,*) '*WARNING reading *MODAL DYNAMIC: unknown solver;'
         write(*,*) '         the default solver is used'
      endif
!
c      if((isolver.eq.2).or.(isolver.eq.3)) then
c        write(*,*) '*ERROR reading *MODAL DYNAMIC: the default solver ',
c     & solver
c         write(*,*) '       cannot be used for modal dynamic'
c         write(*,*) '       calculations '
c         ier=1
c         return
c      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*) 
     &       '*ERROR reading *MODAL DYNAMIC: definition not complete'
         write(*,*) '       '
         call inputerror(inpc,ipoinpc,iline,
     &        "*MODAL DYNAMIC%",ier)
         return
      endif
      read(textpart(1)(1:20),'(f20.0)',iostat=istat)tinc
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*MODAL DYNAMIC%",ier)
         return
      endif
      read(textpart(2)(1:20),'(f20.0)',iostat=istat)tper
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*MODAL DYNAMIC%",ier)
         return
      endif
      read(textpart(3)(1:20),'(f20.0)',iostat=istat) tmin
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*MODAL DYNAMIC%",ier)
         return
      endif
      read(textpart(4)(1:20),'(f20.0)',iostat=istat) tmax
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*MODAL DYNAMIC%",ier)
         return
      endif
!
      if(steadystate) then
!
!        modal dynamics calculation till steady state
!
         if(tper.le.0.d0) then
            write(*,*) '*ERROR reading *MODAL DYNAMIC: relative error'
            write(*,*) '       is nonpositive'
            ier=1
            return
         endif
         tper=-tper
         if(tinc.le.0.d0) then
            write(*,*) 
     &         '*ERROR reading *MODAL DYNAMIC: initial increment'
            write(*,*) '       size is nonpositive'
            ier=1
            return
         endif
         if(tmin.lt.0.d0) then
            tmin=1.d-10
         endif
         if(tmax.lt.1.d-10) then
            tmax=1.d+30
         endif
      else
!
!        transient modal dynamics calculation
!
         if(tper.lt.0.d0) then
            write(*,*) 
     &        '*ERROR reading *MODAL DYNAMIC: step size is negative'
            ier=1
            return
         elseif(tper.le.0.d0) then
            tper=1.d0
         endif
         if(tinc.lt.0.d0) then
            write(*,*) 
     &            '*ERROR reading *MODAL DYNAMIC: initial increment size 
     &is negative'
            ier=1
            return
         elseif(tinc.le.0.d0) then
            tinc=tper
         endif
         if(tinc.gt.tper) then
            write(*,*) 
     &            '*ERROR reading *MODAL DYNAMIC: initial increment size 
     &exceeds step size'
            ier=1
            return
         endif
!      
         if(idrct.ne.1) then
            if(tmin.lt.1.d-10*tper) then
               tmin=min(tinc,1.d-10*tper)
            endif
            if(tmax.lt.1.d-10) then
               tmax=1.d+30
            endif
         endif
      endif
!
!     removing the present loading
!
c      nforc=0
c      nload=0
c      nbody=0
c      iprestr=0
c      if((ithermal.eq.1).or.(ithermal.eq.3)) then
c         do j=1,nk
c            t1(j)=t0(j)
c         enddo
c      endif
!
!     resetting fields vold and veold after a frequency or
!     buckling step
!
      if((nmethod.eq.2).or.(nmethod.eq.3)) then
         do i=1,nk
            do j=1,3
               vold(j,i)=0.d0
               veold(j,i)=0.d0
            enddo
         enddo
      endif
!
      nmethod=4
!
!     correction for cyclic symmetric structures:
!     if the present step was not preceded by a frequency step
!     no nodal diameter has been selected. To make sure that
!     mastructcs is called instead of mastruct a fictitious
!     minimum nodal diameter is stored
!
      if((cyclicsymmetry.eq.1).and.(mcs.ne.0).and.(cs(2,1).lt.0.d0)) 
     &       cs(2,1)=0.d0
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end
