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
      subroutine robustdesigns(inpc,textpart,nmethod,
     &  istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,tieset,ipoinpc,ntie,tinc,tper,tmin,tmax,tincf,isens,
     &  ier,physcon,irobustdesign)
!
!     reading the input deck: *ROBUST DESIGN
!
      implicit none
!
      logical iread,iwrite
!
      character*1 inpc(*)
      character*81 tieset(3,*)
      character*132 textpart(16)
!
      integer nmethod,istep,istat,n,key,i,isens,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),irobustdesign(*),
     &  ipoinpc(0:*),ntie,ier
!
      real*8 tinc,tper,tmin,tmax,tincf,reliability,physcon(*)
!
!     Read in *ROBUST DESIGN
!
      if(isens.eq.1) then
         write(*,*) '*ERROR reading *ROBUST DESIGN:'
         write(*,*) '      no more than one *ROBUST DESIGN'
         write(*,*) '      is allowed per input deck'
         ier=1
         return
      endif
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *ROBUST DESIGN: *ROBUST DESIGN can
     &         only be used within a STEP'     
         ier=1
         return
      endif
!
c      if(istep.lt.2) then
c         write(*,*) '*ERROR reading *ROBUST DESIGN: *ROBUST DESIGN'
c         write(*,*) '      requires a previous *STATIC step'
c         ier=1
c         return
c      endif
!
      tinc=0.d0
      tper=0.d0
      tmin=0.d0
      tmax=0.d0
      tincf=0.d0
!
      iwrite=.false.
      iread=.false.
!
      nmethod=14
!
!     check whether design variables were defined
!
c      do i=1,ntie
c         if(tieset(1,i)(81:81).eq.'D') exit
c      enddo
c      if(i.gt.ntie) then
c         write(*,*) '*ERROR reading *ROBUST DESIGN'
c         write(*,*) '      no design variables were defined'
c         call inputerror(inpc,ipoinpc,iline,
c     &        "*ROBUST DESIGN%",ier)
c         return
c      endif
!
!     check what information is requested by the user
!     irobustdesign(1)=1 --> the full stochastic perturbation method 
!                         is performed (default)
!     irobustdesign(1)=2 --> only the eigenvectors of the randomfield is
!                         calculated 
!
c      if(n.eq.2) then
c         if(textpart(2)(1:15).eq.'RANDOMFIELDONLY') then
            irobustdesign(1)=2
c         else
c            irobustdesign(1)=1
c            write(*,*) '*WARNING Keyword in *ROBUST DESIGN'
c            write(*,*) '   not known, keyword ignored'
c         endif
c      endif
!
!     Read in the reliability of the random field
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
!
      read(textpart(1)(1:20),'(f20.0)',iostat=istat) reliability
      if(istat.gt.0) then
         write(*,*) '*ERROR in *ROBUST DESIGN reliability of'
         write(*,*) '       the random field not specified'
         call inputerror(inpc,ipoinpc,iline,
     &        "*ROBUST DESIGN%",ier)
         return
      endif
      if((reliability.le.0.d0).or.(reliability.ge.1.d0)) then
         write(*,*) '*ERROR reading *ROBUST DESIGN'
         write(*,*) '       Reliability of the random field'
         write(*,*) '       has to be in the range'
         write(*,*) '       between 0 and 1'
         write(*,*) 
         call inputerror(inpc,ipoinpc,iline,
     &          "*ROBUST DESIGN%",ier)
         return
      endif
      physcon(11)=reliability
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
!
      return
      end
