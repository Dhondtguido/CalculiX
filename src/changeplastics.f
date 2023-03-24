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
      subroutine changeplastics(inpc,textpart,imat,ntmat_,npmat_,
     &        plicon,nplicon,plkcon,nplkcon,istep,istat,n,iline,ipol,
     &        inl,ipoinp,inp,ipoinpc,nelcon,ier)
!
!     reading the input deck: *CHANGE PLASTIC
!
      implicit none
!
      logical iso
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer imat,ntmat_,ntmat,npmat_,npmat,istep,nelcon(2,*),
     &  n,key,i,nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),istat,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),ipoinpc(0:*),ier
!
      real*8 plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     & temperature
!
      iso=.true.
!
      ntmat=0
      npmat=0
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *CHANGE PLASTIC: *CHANGE PLASTIC'
         write(*,*) '       should only be used within a STEP'
         ier=1
         return
      endif
!
      if((nelcon(1,imat).ne.-51).and.
     &   (nelcon(1,imat).ne.-52)) then
         write(*,*) '*ERROR reading *CHANGE PLASTIC: *CHANGE PLASTIC'
         write(*,*) '       can only be used to change the plastic'
         write(*,*) '       definition of an elastically isotropic'
         write(*,*) '       material with *PLASTIC data'
         ier=1
         return
      endif
!
      do i=2,n
         if(textpart(i)(1:10).eq.'HARDENING=') then
            if(textpart(i)(11:19).eq.'KINEMATIC') then
               iso=.false.
            elseif(textpart(i)(11:18).eq.'COMBINED') then
               write(*,*) '*ERROR reading *CHANGE PLASTIC'
               write(*,*) '       combined hardening is not allowed'
               ier=1
               return
            elseif(textpart(i)(11:14).eq.'USER') then
               write(*,*) '*ERROR reading *CHANGE PLASTIC'
               write(*,*) '       parameter USER is not allowed'
               ier=1
               return
            endif
            exit
         else
            write(*,*) 
     &     '*WARNING reading *CHANGE PLASTIC: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*CHANGE PLASTIC%")
         endif
      enddo
!
      if(iso) then
!
!        isotropic hardening coefficients
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) temperature
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*CHANGE PLASTIC%",ier)
               return
            endif
!
!           first temperature
!
            if(ntmat.eq.0) then
               npmat=0
               ntmat=ntmat+1
               if(ntmat.gt.ntmat_) then
                  write(*,*) '*ERROR reading *CHANGE PLASTIC:'
                  write(*,*) '       more temperature data points'
                  write(*,*) '       than underneath the *PLASTIC card'
                  ier=1
                  return
               endif
               nplicon(0,imat)=ntmat
               plicon(0,ntmat,imat)=temperature
!
!           new temperature
!
            elseif(plicon(0,ntmat,imat).ne.temperature) then
               npmat=0
               ntmat=ntmat+1
               if(ntmat.gt.ntmat_) then
                  write(*,*) '*ERROR reading *CHANGE PLASTIC:'
                  write(*,*) '       more temperature data points'
                  write(*,*) '       than underneath the *PLASTIC card'
                  ier=1
                  return
               endif
               nplicon(0,imat)=ntmat
               plicon(0,ntmat,imat)=temperature
            endif
            do i=1,2
               read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &              plicon(2*npmat+i,ntmat,imat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*CHANGE PLASTIC%",ier)
                  return
               endif
            enddo
            npmat=npmat+1
            if(npmat.gt.npmat_) then
                  write(*,*) '*ERROR reading *CHANGE PLASTIC:'
                  write(*,*) '       more stress versus equivalent'
                  write(*,*) '       plastic strain data points'
                  write(*,*) '       than underneath the *PLASTIC card'
               ier=1
               return
            endif
            nplicon(ntmat,imat)=npmat
         enddo
      else
!
!        kinematic hardening coefficients
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) temperature
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*CHANGE PLASTIC%",ier)
               return
            endif
!
!           first temperature
!
            if(ntmat.eq.0) then
               npmat=0
               ntmat=ntmat+1
               if(ntmat.gt.ntmat_) then
                  write(*,*) '*ERROR reading *CHANGE PLASTIC:'
                  write(*,*) '       more temperature data points'
                  write(*,*) '       than underneath the *PLASTIC card'
                  ier=1
                  return
               endif
               nplkcon(0,imat)=ntmat
               plkcon(0,ntmat,imat)=temperature
!
!           new temperature
!
            elseif(plkcon(0,ntmat,imat).ne.temperature) then
               npmat=0
               ntmat=ntmat+1
               if(ntmat.gt.ntmat_) then
                  write(*,*) '*ERROR reading *CHANGE PLASTIC:'
                  write(*,*) '       more temperature data points'
                  write(*,*) '       than underneath the *PLASTIC card'
                  ier=1
                  return
               endif
               nplkcon(0,imat)=ntmat
               plkcon(0,ntmat,imat)=temperature
            endif
            do i=1,2
               read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &              plkcon(2*npmat+i,ntmat,imat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*CHANGE PLASTIC%",ier)
                  return
               endif
            enddo
            npmat=npmat+1
            if(npmat.gt.npmat_) then
               write(*,*) '*ERROR reading *CHANGE PLASTIC:'
               write(*,*) '       more stress versus equivalent'
               write(*,*) '       plastic strain data points'
               write(*,*) '       than underneath the *PLASTIC card'
               ier=1
               return
            endif
            nplkcon(ntmat,imat)=npmat
         enddo
      endif
!
      if(ntmat.eq.0) then
         write(*,*) '*ERROR reading *CHANGE PLASTIC:'
         write(*,*) '       *CHANGE PLASTIC card without data'
         ier=1
         return
      endif
!
      return
      end

