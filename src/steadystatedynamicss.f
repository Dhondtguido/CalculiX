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
      subroutine steadystatedynamicss(inpc,textpart,nmethod,
     &  iexpl,istep,istat,n,iline,ipol,inl,ipoinp,inp,iperturb,isolver,
     &  xmodal,cs,mcs,ipoinpc,nforc,nload,nbody,iprestr,t0,t1,ithermal,
     &  nk,set,nset,cyclicsymmetry,ibody,ier)
!
!     reading the input deck: *STEADY STATE DYNAMICS
!
      implicit none
!
      logical nodalset
!
      character*1 inpc(*)
      character*3 harmonic
      character*20 solver
      character*81 set(*),noset
      character*132 textpart(16)
!
      integer nmethod,istep,istat,n,key,iexpl,iline,ipol,inl,nset,
     &  ipoinp(2,*),inp(3,*),iperturb(*),isolver,i,ndata,nfour,mcs,
     &  ipoinpc(0:*),nforc,nload,nbody,iprestr,ithermal(*),j,nk,ipos,
     &  cyclicsymmetry,ier,ibody(3,*)
!
      real*8 fmin,fmax,bias,tmin,tmax,xmodal(*),cs(17,*),t0(*),t1(*)
!
      iexpl=0
      iperturb(2)=0
      harmonic='YES'
      if((mcs.ne.0).and.(cs(2,1).ge.0.d0)) then
         cyclicsymmetry=1
      endif
      nodalset=.false.
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *STEADY STATE DYNAMICS:'
         write(*,*) '               *STEADY STATE DYNAMICS'
         write(*,*) '  can only be used within a STEP'
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
         elseif(textpart(i)(1:9).eq.'HARMONIC=') then
            read(textpart(i)(10:12),'(a3)') harmonic
         else
            write(*,*) '*WARNING reading *STEADY STATE DYNAMICS:' 
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*STEADY STATE DYNAMICS%")
         endif
      enddo
!
      if(solver(1:7).eq.'SPOOLES') then
         isolver=0
      elseif(solver(1:16).eq.'ITERATIVESCALING') then
         write(*,*) '*WARNING reading *STEADY STATE DYNAMICS:'
         write(*,*) '         the iterative scaling'
         write(*,*) '         procedure is not available for modal'
         write(*,*) '         dynamic calculations; the default solver'
         write(*,*) '         is used'
      elseif(solver(1:17).eq.'ITERATIVECHOLESKY') then
         write(*,*) '*WARNING reading *STEADY STATE DYNAMICS:'
         write(*,*) '         the iterative scaling'
         write(*,*) '         procedure is not available for modal'
         write(*,*) '         dynamic calculations; the default solver'
         write(*,*) '         is used'
      elseif(solver(1:3).eq.'SGI') then
         isolver=4
      elseif(solver(1:5).eq.'TAUCS') then
         isolver=5
      elseif(solver(1:7).eq.'PARDISO') then
         isolver=7
      elseif(solver(1:6).eq.'PASTIX') then
         isolver=8
      else
         write(*,*) '*WARNING reading *STEADY STATE DYNAMICS:'
         write(*,*) '         unknown solver;'
         write(*,*) '         the default solver is used'
      endif
!
c      if((isolver.eq.2).or.(isolver.eq.3)) then
c         write(*,*) '*ERROR reading *STEADY STATE DYNAMICS:'
c         write(*,*) '       the default solver ',solver
c         write(*,*) '       cannot be used for modal dynamic'
c         write(*,*) '       calculations '
c         ier=1
c         return
c      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*) '*ERROR reading *STEADY STATE DYNAMICS:'
         write(*,*) '       definition not complete'
         write(*,*) '       '
         call inputerror(inpc,ipoinpc,iline,
     &        "*STEADY STATE DYNAMICS%",ier)
         return
      endif
      read(textpart(1)(1:20),'(f20.0)',iostat=istat) fmin
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*STEADY STATE DYNAMICS%",ier)
         return
      endif
      xmodal(3)=fmin
      read(textpart(2)(1:20),'(f20.0)',iostat=istat) fmax
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*STEADY STATE DYNAMICS%",ier)
         return
      endif
      xmodal(4)=fmax
      read(textpart(3)(1:20),'(i10)',iostat=istat) ndata
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*STEADY STATE DYNAMICS%",ier)
         return
      endif
      if(textpart(3)(1:1).eq.' ') then
         ndata=20
      elseif(ndata.lt.2) then
         write(*,*) '*WARNING reading *STEADY STATE DYNAMICS:'
         write(*,*) '         number of data points:',ndata
         write(*,*) '         is smaller than 2;'
         write(*,*) '         2 data points are taken instead'
         ndata=2
      endif
      xmodal(5)=ndata+0.5d0
      read(textpart(4)(1:20),'(f20.0)',iostat=istat) bias
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*STEADY STATE DYNAMICS%",ier)
         return
      endif
      if(bias.lt.1.) bias=3.d0
      xmodal(6)=bias
!
      if(harmonic.eq.'YES') then
         xmodal(7)=-0.5d0
      else
         read(textpart(5)(1:10),'(i10)',iostat=istat) nfour
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*STEADY STATE DYNAMICS%",ier)
            return
         endif
         if(nfour.le.0) nfour=20
         if(n.ge.6) then
            read(textpart(6)(1:20),'(f20.0)',iostat=istat) tmin
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*STEADY STATE DYNAMICS%",ier)
               return
            endif
         else
            tmin=0.d0
         endif
         if(n.ge.7) then
            read(textpart(7)(1:20),'(f20.0)',iostat=istat) tmax
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*STEADY STATE DYNAMICS%",ier)
               return
            endif
         else
            tmax=1.d0
         endif
         xmodal(7)=nfour+0.5d0
         xmodal(8)=tmin
         xmodal(9)=tmax
      endif
!
!     removing the present loading except centrigufal loading
!
      nforc=0
      nload=0
      do i=1,nbody
         if(ibody(1,i).ne.1) ibody(1,i)=0
      enddo
c      nbody=0
      iprestr=0
      if((ithermal(1).eq.1).or.(ithermal(1).eq.3)) then
         do j=1,nk
            t1(j)=t0(j)
         enddo
      endif
!
      nmethod=5
!
!     correction for cyclic symmetric structures:
!       if the present step was not preceded by a frequency step
!       no nodal diameter has been selected. To make sure that
!       mastructcs is called instead of mastruct a fictitious
!       minimum nodal diameter is stored
!
      if((cyclicsymmetry.eq.1).and.(mcs.ne.0).and.(cs(2,1)<0.d0)) 
     &         cs(2,1)=0.d0
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

