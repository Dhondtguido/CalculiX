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
      subroutine frequencys(inpc,textpart,nmethod,
     &  mei,fei,iperturb,istep,istat,n,iline,ipol,inl,
     &  ipoinp,inp,ithermal,isolver,xboun,nboun,ipoinpc,
     &  ipompc,labmpc,fmpc,ikmpc,ilmpc,nmpc,ier)
!
!     reading the input deck: *FREQUENCY
!
      implicit none
!
      logical global,cycmpcactive
!
      character*1 inpc(*)
      character*20 solver,labmpc(*)
      character*132 textpart(16)
!
      integer nmethod,mei(4),ncv,mxiter,istep,istat,iperturb(*),i,
     &  nboun,ier,isolver,
     &  n,key,iline,ipol,inl,ipoinp(2,*),inp(3,*),nev,ithermal(*),
     &  ipoinpc(0:*),nmpcred,kflag,ipompc(*),ikmpc(*),ilmpc(*),nmpc
!
      real*8 fei(3),pi,fmin,fmax,tol,xboun(*),fmpc(*)
!
      global=.true.
      cycmpcactive=.true.
!
      pi=4.d0*datan(1.d0)
      mei(4)=0
!
!     defaults for fmin and fmax
!
      fmin=-1.d0
      fmax=-1.d0
!
      if(istep.lt.1) then
         write(*,*) 
     &      '*ERROR reading *FREQUENCY: *FREQUENCY can only be used'
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
         elseif(textpart(i)(1:9).eq.'GLOBAL=NO') then
            global=.false.
         elseif(textpart(i)(1:15).eq.'CYCMPC=INACTIVE') then
            cycmpcactive=.false.
         else
            write(*,*) 
     &        '*WARNING reading *FREQUENCY: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*FREQUENCY%")
         endif
      enddo
!
      if(solver(1:7).eq.'SPOOLES') then
         isolver=0
      elseif(solver(1:16).eq.'ITERATIVESCALING') then
         write(*,*) '*WARNING reading *FREQUENCY: the iterative scaling'
         write(*,*) '         procedure is not available for frequency'
         write(*,*) '         calculations; the default solver is used'
      elseif(solver(1:17).eq.'ITERATIVECHOLESKY') then
         write(*,*) '*WARNING reading *FREQUENCY: the iterative scaling'
         write(*,*) '         procedure is not available for frequency'
         write(*,*) '         calculations; the default solver is used'
      elseif(solver(1:3).eq.'SGI') then
         isolver=4
      elseif(solver(1:5).eq.'TAUCS') then
         isolver=5
      elseif(solver(1:13).eq.'MATRIXSTORAGE') then
         isolver=6
      elseif(solver(1:7).eq.'PARDISO') then
         isolver=7
      elseif(solver(1:6).eq.'PASTIX') then
         isolver=8
      else
         write(*,*) '*WARNING reading *FREQUENCY: unknown solver;'
         write(*,*) '         the default solver is used'
      endif
!
      if((isolver.eq.2).or.(isolver.eq.3)) then
         write(*,*) '*ERROR reading *FREQUENCY: the default solver ',
     & solver
         write(*,*) '       cannot be used for frequency calculations '
         ier=1
         return
      endif
!
      nmethod=2
      if(iperturb(1).gt.1) iperturb(1)=0
c      iperturb(2)=0
!
      if(isolver.ne.6) then
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) then
            write(*,*) 
     &       '*ERROR reading *FREQUENCY: definition not complete'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &           "*FREQUENCY%",ier)
            return
         endif
         read(textpart(1)(1:10),'(i10)',iostat=istat) nev
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*FREQUENCY%",ier)
            return
         endif
         if(nev.le.0) then
            write(*,*) '*ERROR reading *FREQUENCY: less than 1 eigenvalu
     &e requested'
            ier=1
            return
         endif
         tol=1.d-2
         ncv=4*nev
         ncv=ncv+nev
         mxiter=1000
         if(textpart(2)(1:1).ne.' ') then
            read(textpart(2)(1:20),'(f20.0)',iostat=istat) fmin
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*FREQUENCY%",ier)
               return
            endif
         endif
         if(textpart(3)(1:1).ne.' ') then
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) fmax
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*FREQUENCY%",ier)
               return
            endif
         endif
!     
         mei(1)=nev
         mei(2)=ncv
         mei(3)=mxiter
         fei(1)=tol
         fei(2)=fmin
         fei(3)=fmax
!     
      else
!
!        matrix storage
!
         if(global) then
            mei(1)=1
         else
            mei(1)=0
         endif
!
         if(.not.cycmpcactive) then
!
!           remove the CYCLIC and SUBCYCLIC MPC's
!
            nmpcred=0
!
            kflag=2
            call isortii(ilmpc,ikmpc,nmpc,kflag)
!
            do i=1,nmpc
               if((labmpc(i)(1:6).eq.'CYCLIC').or.
     &            (labmpc(i)(1:9).eq.'SUBCYCLIC')) cycle
               nmpcred=nmpcred+1
               ipompc(nmpcred)=ipompc(i)
               labmpc(nmpcred)=labmpc(i)
               fmpc(nmpcred)=fmpc(i)
               ikmpc(nmpcred)=ikmpc(i)
               ilmpc(nmpcred)=nmpcred
            enddo
!
            nmpc=nmpcred
            call isortii(ikmpc,ilmpc,nmpc,kflag)
         endif
      endif
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

