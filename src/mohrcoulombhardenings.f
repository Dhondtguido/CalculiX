!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2023 Guido Dhondt
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
      subroutine mohrcoulombhardenings(inpc,textpart,nelcon,nmat,
     &        plicon,nplicon,plkcon,nplkcon,iplas,iperturb,nstate_,
     &        ncmat_,elcon,matname,irstrt,istep,istat,n,iline,ipol,
     &        inl,ipoinp,inp,ipoinpc,ianisoplas,ier,ntmat_,npmat_)
!
!     reading the input deck: *MOHR COULOMB HARDENING
!
      implicit none
!
      character*1 inpc(*)
      character*80 matname(*)
      character*132 textpart(16)
!     
      integer nelcon(2,*),nmat,ntmat_,ntmat,npmat_,npmat,istep,
     &     n,key,i,nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),ncmat_,
     &     iplas,iperturb(*),istat,nstate_,kin,itemp,ndata,ndatamax,id,
     &     irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),ipoinpc(0:*),
     &     ianisoplas,ier
!     
      real*8 plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &     temperature,plconloc(802),t1l,elcon(0:ncmat_,ntmat_,*)
!     
      ntmat=0
      npmat=0
!     
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
        write(*,*) '*ERROR reading *MOHR COULOMB HARDENING:'
        write(*,*) '       *MOHR COULOMB HARDENING should be placed'
        write(*,*) '       before all step definitions'
        ier=1
        return
      endif
!     
      if(nmat.eq.0) then
        write(*,*) '*ERROR reading *MOHR COULOMB HARDENING:'
        write(*,*) '       *MOHR COULOMB HARDENING should be preceded'
        write(*,*) '       by a *MATERIAL card'
        ier=1
        return
      endif
!     
      iperturb(1)=3
!     
      do i=2,n
          write(*,*) '*WARNING reading *MOHR COULOMB HARDENING:'
          write(*,*) '         parameter not recognized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*MOHR COULOMB HARDENING%")
      enddo
!     
      iplas=1
      nelcon(1,nmat)=-53
      nstate_=max(nstate_,7)
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
     &           "*MOHR COULOMB HARDENING%",ier)
            return
          endif
!     
!     first temperature
!     
          if(ntmat.eq.0) then
            npmat=0
            ntmat=ntmat+1
            if(ntmat.gt.ntmat_) then
              write(*,*) '*ERROR reading *MOHR COULOMB HARDENING:'
              write(*,*) '       increase ntmat_'
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
              write(*,*) '*ERROR reading *MOHR COULOMB HARDENING:'
              write(*,*) '       increase ntmat_'
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
     &             "*MOHR COULOMB HARDENING%",ier)
              return
            endif
          enddo
          npmat=npmat+1
          if(npmat.gt.npmat_) then
            write(*,*) '*ERROR reading *MOHR COULOMB HARDENING:'
            write(*,*) '       increase npmat_'
            ier=1
            return
          endif
          nplicon(ntmat,nmat)=npmat
        enddo
!     
      if(ntmat.eq.0) then
        write(*,*) '*ERROR reading *MOHR COULOMB HARDENING:'
        write(*,*) '       *MOHR COULOMB HARDENING card without data'
        ier=1
        return
      endif
!
      nplkcon(0,nmat)=0
!     
      return
      end

