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
      subroutine complexfrequencys(inpc,textpart,nmethod,
     &  mei,iperturb,istep,istat,n,iline,ipol,inl,
     &  ipoinp,inp,ithermal,xboun,nboun,ipoinpc,mcs,cs,cyclicsymmetry,
     &  ier)
!
!     reading the input deck: *COMPLEX FREQUENCY
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nmethod,mei(4),istep,istat,iperturb(*),i,nboun,
     &  n,key,iline,ipol,inl,ipoinp(2,*),inp(3,*),nev,ithermal(*),
     &  ipoinpc(0:*),mcs,cyclicsymmetry,ier
!
      real*8 xboun(*),cs(17,*)
!
      mei(4)=0
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *COMPLEX FREQUENCY:'
         write(*,*) '       *COMPLEX FREQUENCY can only be used'
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
!     check for cyclic symmetry
!
      if((mcs.ne.0).and.(cs(2,1).ge.0.d0)) then
         cyclicsymmetry=1
      endif
!
      nmethod=0
      do i=2,n
         if(textpart(i)(1:8).eq.'CORIOLIS') then
            nmethod=6
         elseif(textpart(i)(1:7).eq.'FLUTTER') then
            nmethod=7
         elseif(textpart(i)(1:11).eq.'STORAGE=YES') then
            write(*,*) '*WARNING reading *COMPLEX FREQUENCY:'
            write(*,*) '         for this keyword'
            write(*,*) '         STORAGE=YES is deactivated'
            write(*,*) '         in the CalculiX code'

c            mei(4)=1
         else
            write(*,*) 
     &           '*WARNING reading *COMPLEX FREQUENCY:'
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',
     &           textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*COMPLEX FREQUENCY%")
         endif
      enddo
      if(nmethod.eq.0) then
         write(*,*) 
     &        '*ERROR reading *COMPLEX FREQUENCY:'
         write(*,*) '       either parameter CORIOLIS'
         write(*,*) '       or parameter FLUTTER is required'
         call inputerror(inpc,ipoinpc,iline,
     &        "*COMPLEX FREQUENCY%",ier)
         return
      endif
!
      if(iperturb(1).gt.1) iperturb(1)=0
      iperturb(2)=0
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*) '*ERROR reading *COMPLEX FREQUENCY:'
         write(*,*) '       definition not complete'
         write(*,*) '  '
         call inputerror(inpc,ipoinpc,iline,
     &        "*COMPLEX FREQUENCY%",ier)
         return
      endif
      read(textpart(1)(1:10),'(i10)',iostat=istat) nev
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*COMPLEX FREQUENCY%",ier)
         return
      endif
      if(nev.le.0) then
         write(*,*) '*ERROR reading *COMPLEX FREQUENCY:'
         write(*,*) '       less than 1 eigenvalue requested'
         ier=1
         return
      endif
!
      mei(1)=nev
!
!     removing nonzero boundary conditions
!
      do i=1,nboun
         xboun(i)=0.d0
      enddo
!
!     correction for cyclic symmetric structures:
!     if the present step was not preceded by a frequency step
!     no nodal diameter has been selected. To make sure that
!     mastructcs is called instead of mastruct a fictitious
!     minimum nodal diameter is stored
!
      if((cyclicsymmetry.eq.1).and.(mcs.ne.0).and.(cs(2,1)<0.d0)) 
     &       cs(2,1)=0.d0
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

