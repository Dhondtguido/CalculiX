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
      subroutine plastics(inpc,textpart,nelcon,nmat,ntmat_,npmat_,
     &     plicon,nplicon,plkcon,nplkcon,iplas,iperturb,nstate_,
     &     ncmat_,elcon,matname,irstrt,istep,istat,n,iline,ipol,
     &     inl,ipoinp,inp,ipoinpc,ier)
!     
!     reading the input deck: *PLASTIC
!     
      implicit none
!     
      logical iso,johnsoncook
!     
      character*1 inpc(*)
      character*80 matname(*)
      character*132 textpart(16)
!     
      integer nelcon(2,*),nmat,ntmat_,ntmat,npmat_,npmat,istep,
     &     n,key,i,nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),ncmat_,
     &     iplas,iperturb(*),istat,nstate_,
     &     irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),ipoinpc(0:*),
     &     ier
!     
      real*8 plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &     temperature,elcon(0:ncmat_,ntmat_,*)
!     
      iso=.true.
      johnsoncook=.false.
!     
      ntmat=0
      npmat=0
!     
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
        write(*,*) '*ERROR reading *PLASTIC: *PLASTIC should be placed'
        write(*,*) '  before all step definitions'
        ier=1
        return
      endif
!     
      if(nmat.eq.0) then
        write(*,*) 
     &       '*ERROR reading *PLASTIC: *PLASTIC should be preceded'
        write(*,*) '  by a *MATERIAL card'
        ier=1
        return
      endif
!     
      if((nelcon(1,nmat).ne.2).and.(nelcon(1,nmat).ne.9)) then
        write(*,*) 
     &       '*ERROR reading *PLASTIC: *PLASTIC should be preceded'
        write(*,*) '       by an *ELASTIC,TYPE=ISO card or'
        write(*,*) '       by an *ELASTIC,TYPE=ORTHO card'
        ier=1
        return
      endif
!     
      iperturb(1)=3
!     
      do i=2,n
        if(textpart(i)(1:10).eq.'HARDENING=') then
          if(textpart(i)(11:19).eq.'KINEMATIC') then
            iso=.false.
          elseif(textpart(i)(11:18).eq.'COMBINED') then
            iso=.false.
          elseif(textpart(i)(11:14).eq.'USER') then
            if(nelcon(1,nmat).ne.2) then
              write(*,*) '*ERROR reading *PLASTIC: user defined '
              write(*,*) '       hardening is not allowed for '
              write(*,*) '       elastically anisotropic materials'
              ier=1
              return
            endif
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            return
          elseif(textpart(i)(11:21).eq.'JOHNSONCOOK') then
            johnsoncook=.true.
          endif
c          exit
c        elseif(textpart(i)(1:21).eq.'HARDENING=JOHNSONCOOK') then
c          johnsoncook=.true.
        else
          write(*,*) 
     &         '*WARNING reading *PLASTIC: parameter not recognized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*PLASTIC%")
        endif
      enddo
!     
!     special case: Johnson-Cook
!     user material; npmat=0; ntmat=1;
!     
      if(johnsoncook) then
!     
        if(matname(nmat)(1:11).ne.'JOHNSONCOOK') then
          if(matname(nmat)(70:80).ne.'           ') then
            write(*,*)
     &           '*ERROR reading *PLASTIC: the material name'
            write(*,*) '       for a Johnson-Cook material must'
            write(*,*) '       not exceed 69 characters'
            ier=1
            return
          else
            do i=80,12,-1
              matname(nmat)(i:i)=matname(nmat)(i-11:i-11)
            enddo
            matname(nmat)(1:11)='JOHNSONCOOK'
          endif
c          write(*,*) '*ERROR reading *PLASTIC'
c          write(*,*) '       the name of a Johnson Cook material'
c          write(*,*) '       must start with JOHNSONCOOK'
c          write(*,*) '       (blanks are allowed at any location'
c          write(*,*) '        and the string is not case sensitive)'
c          call inputerror(inpc,ipoinpc,iline,
c     &         "*PLASTIC%",ier)
c          return
        endif
!     
        if(nelcon(1,nmat).ne.2) then
          write(*,*) 
     &         '*ERROR reading *PLASTIC: for a Johnson Cook material'
          write(*,*) '       *PLASTIC should be preceded'
          write(*,*) '       by an *ELASTIC,TYPE=ISO card'
          ier=1
          return
        endif
!     
        nelcon(1,nmat)=-111
        nstate_=max(nstate_,9)
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
        if((istat.lt.0).or.(key.eq.1).or.(n.lt.4)) then
          write(*,*) '*ERROR reading *PLASTIC'
          write(*,*) '       for the Johnson-Cook model at least'
          write(*,*) '       A, B, n and m should be given'
          call inputerror(inpc,ipoinpc,iline,
     &         "*PLASTIC%",ier)
          return
        endif
        ntmat=1
        if(ntmat.gt.ntmat_) then
          write(*,*) 
     &         '*ERROR reading *PLASTIC: increase ntmat_'
          ier=1
          return
        endif
!     
!     reading A, B and n (elcon(4,5,6))
!     
        do i=1,3
          read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &         elcon(3+i,ntmat,nmat)
        enddo
        if(elcon(6,ntmat,nmat).le.0.d0) then
          write(*,*) '*ERROR reading *PLASTIC'
          write(*,*) '       n must be strictly positive'
          call inputerror(inpc,ipoinpc,iline,
     &         "*PLASTIC%",ier)
          return
        endif
!     
!     reading m (elcon(11))
!     
        read(textpart(4)(1:20),'(f20.0)',iostat=istat)
     &       elcon(11,ntmat,nmat)
        if(elcon(11,ntmat,nmat).le.0.d0) then
          write(*,*) '*ERROR reading *PLASTIC'
          write(*,*) '       m must be strictly positive'
          call inputerror(inpc,ipoinpc,iline,
     &         "*PLASTIC%",ier)
          return
        endif
!     
!     reading, if given, Tmelt and Ttransition (elcon(9,10))
!     
        if(n.gt.4) then
          do i=5,min(n,6)
            read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &           elcon(4+i,ntmat,nmat)
          enddo
        endif
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
        return
      endif
!
      iplas=1
      if(nelcon(1,nmat).eq.2) then
        nelcon(1,nmat)=-51
        nstate_=max(nstate_,13)
      else
        nelcon(1,nmat)=-54
        nstate_=max(nstate_,20)
      endif
!     
      if(iso) then
!     
!     isotropic hardening coefficients
!     
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          read(textpart(3)(1:20),'(f20.0)',iostat=istat) temperature
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*PLASTIC%",ier)
            return
          endif
!     
!     first temperature
!     
          if(ntmat.eq.0) then
            npmat=0
            ntmat=ntmat+1
            if(ntmat.gt.ntmat_) then
              write(*,*) '*ERROR reading *PLASTIC: increase ntmat_'
              ier=1
              return
            endif
            nplicon(0,nmat)=ntmat
            plicon(0,ntmat,nmat)=temperature
!     
!     new temperature
!     
          elseif(plicon(0,ntmat,nmat).ne.temperature) then
            npmat=0
            ntmat=ntmat+1
            if(ntmat.gt.ntmat_) then
              write(*,*) '*ERROR reading *PLASTIC: increase ntmat_'
              ier=1
              return
            endif
            nplicon(0,nmat)=ntmat
            plicon(0,ntmat,nmat)=temperature
          endif
          do i=1,2
            read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &           plicon(2*npmat+i,ntmat,nmat)
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*PLASTIC%",ier)
              return
            endif
          enddo
          npmat=npmat+1
          if(npmat.gt.npmat_) then
            write(*,*) '*ERROR reading *PLASTIC: increase npmat_'
            ier=1
            return
          endif
          nplicon(ntmat,nmat)=npmat
        enddo
      else
!     
!     kinematic hardening coefficients
!     
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          read(textpart(3)(1:20),'(f20.0)',iostat=istat) temperature
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*PLASTIC%",ier)
            return
          endif
!     
!     first temperature
!     
          if(ntmat.eq.0) then
            npmat=0
            ntmat=ntmat+1
            if(ntmat.gt.ntmat_) then
              write(*,*) '*ERROR reading *PLASTIC: increase ntmat_'
              ier=1
              return
            endif
            nplkcon(0,nmat)=ntmat
            plkcon(0,ntmat,nmat)=temperature
!     
!     new temperature
!     
          elseif(plkcon(0,ntmat,nmat).ne.temperature) then
            npmat=0
            ntmat=ntmat+1
            if(ntmat.gt.ntmat_) then
              write(*,*) '*ERROR reading *PLASTIC: increase ntmat_'
              ier=1
              return
            endif
            nplkcon(0,nmat)=ntmat
            plkcon(0,ntmat,nmat)=temperature
          endif
          do i=1,2
            read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &           plkcon(2*npmat+i,ntmat,nmat)
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*PLASTIC%",ier)
              return
            endif
          enddo
          npmat=npmat+1
          if(npmat.gt.npmat_) then
            write(*,*) '*ERROR reading *PLASTIC: increase npmat_'
            ier=1
            return
          endif
          nplkcon(ntmat,nmat)=npmat
        enddo
      endif
!     
      if(ntmat.eq.0) then
        write(*,*) 
     &       '*ERROR reading *PLASTIC: *PLASTIC card without data'
        ier=1
        return
      endif
!     
      return
      end

