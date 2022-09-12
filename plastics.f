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
     &        plicon,nplicon,plkcon,nplkcon,iplas,iperturb,nstate_,
     &        ncmat_,elcon,matname,irstrt,istep,istat,n,iline,ipol,
     &        inl,ipoinp,inp,ipoinpc,ianisoplas,ier)
!
!     reading the input deck: *PLASTIC
!
      implicit none
!
      logical iso
!
      character*1 inpc(*)
      character*80 matname(*)
      character*132 textpart(16)
!
      integer nelcon(2,*),nmat,ntmat_,ntmat,npmat_,npmat,istep,
     &  n,key,i,nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),ncmat_,
     &  iplas,iperturb(*),istat,nstate_,kin,itemp,ndata,ndatamax,id,
     &  irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),ipoinpc(0:*),
     &  ianisoplas,ier
!
      real*8 plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     & temperature,plconloc(802),t1l,elcon(0:ncmat_,ntmat_,*)
!
      iso=.true.
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
     &      '*ERROR reading *PLASTIC: *PLASTIC should be preceded'
         write(*,*) '  by a *MATERIAL card'
         ier=1
         return
      endif
!
      if((nelcon(1,nmat).ne.2).and.(nelcon(1,nmat).ne.9)) then
         write(*,*) 
     &        '*ERROR reading *PLASTIC: *PLASTIC should be preceded'
         write(*,*) '  by an *ELASTIC,TYPE=ISO card or'
         write(*,*) '  by an *ELASTIC,TYPE=ORTHO card'
         ier=1
         return
      endif
!
      iperturb(1)=3
c      iperturb(2)=1
c      write(*,*) '*INFO reading *PLASTIC: nonlinear geometric'
c      write(*,*) '      effects are turned on'
c      write(*,*)
!
      if(nelcon(1,nmat).eq.2) then
         iplas=1
         nelcon(1,nmat)=-51
         nstate_=max(nstate_,13)
      else
         ianisoplas=1
         nelcon(1,nmat)=-114
         nstate_=max(nstate_,14)
      endif
!
      do i=2,n
         if(textpart(i)(1:10).eq.'HARDENING=') then
            if(textpart(i)(11:19).eq.'KINEMATIC') then
               iso=.false.
            elseif(textpart(i)(11:18).eq.'COMBINED') then
               iso=.false.
            elseif(textpart(i)(11:14).eq.'USER') then
               if(nelcon(1,nmat).eq.-114) then
                  write(*,*) '*ERROR reading *PLASTIC: user defined '
                  write(*,*) '       hardening is not allowed for '
                  write(*,*) '       elastically anisotropic materials'
                  ier=1
                  return
               endif
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               return
            endif
            exit
         else
            write(*,*) 
     &        '*WARNING reading *PLASTIC: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*PLASTIC%")
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
     &              "*PLASTIC%",ier)
               return
            endif
!
!           first temperature
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
!           new temperature
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
     &              plicon(2*npmat+i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*PLASTIC%",ier)
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
!        kinematic hardening coefficients
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) temperature
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*PLASTIC%",ier)
               return
            endif
!
!           first temperature
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
!           new temperature
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
     &              plkcon(2*npmat+i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*PLASTIC%",ier)
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
!     elastically anisotropic materials: recasting the input data
!     in a format conform to the user routine umat_aniso_plas.f
!
      if(nelcon(1,nmat).eq.-114) then
         if(matname(nmat)(71:80).ne.'          ') then
            write(*,*) 
     &          '*ERROR reading *PLASTIC: the material name for an'
            write(*,*) '       elastically anisotropic material with'
            write(*,*) '       isotropic plasticity must not exceed 70'
            write(*,*) '       characters'
            ier=1
            return
         else
            do i=80,11,-1
               matname(nmat)(i:i)=matname(nmat)(i-10:i-10)
            enddo
c            matname(nmat)(11:80)=matname(nmat)(1:70)
            matname(nmat)(1:10)='ANISO_PLAS'
         endif
!
         if(iso) then
!
!           isotropic hardening
!
!           interpolating the plastic data at the elastic temperature
!           data points
!
            ndatamax=0
            do i=1,nelcon(2,nmat)
               t1l=elcon(0,i,nmat)
c               plconloc(1)=0.d0
c               plconloc(2)=0.d0
c               plconloc(3)=0.d0
c               plconloc(81)=nplicon(1,nmat)+0.5d0
!
               if(nplicon(0,nmat).eq.1) then
                  id=-1
               else
                  call ident2(plicon(0,1,nmat),t1l,nplicon(0,nmat),
     &                        2*npmat_+1,id)
               endif
!
               if(nplicon(0,nmat).eq.0) then
                  continue
               elseif((nplicon(0,nmat).eq.1).or.(id.eq.0).or.
     &                 (id.eq.nplicon(0,nmat))) then
                  if(id.le.0) then
                     itemp=1
                  else
                     itemp=id
                  endif
                  kin=0
                  call plcopy(plicon,nplicon,plconloc,npmat_,ntmat_,
     &                 nmat,itemp,i,kin)
                  if((id.eq.0).or.(id.eq.nplicon(0,nmat))) then
                  endif
               else
                  kin=0
                  call plmix(plicon,nplicon,plconloc,npmat_,ntmat_,
     &                 nmat,id+1,t1l,i,kin)
               endif
!
               ndata=int(plconloc(801))
               if(ndata.eq.1) then
                  elcon(10,i,nmat)=plconloc(2)
                  elcon(11,i,nmat)=0.d0
                  elcon(12,i,nmat)=0.d0
                  elcon(13,i,nmat)=-1.d0
                  elcon(14,i,nmat)=1.d0
               else
                  elcon(10,i,nmat)=plconloc(2)
                  elcon(11,i,nmat)=(plconloc(4)-plconloc(2))/
     &                             (plconloc(3)-plconloc(1))
                  elcon(12,i,nmat)=0.d0
                  elcon(13,i,nmat)=-1.d0
                  elcon(14,i,nmat)=1.d0
               endif
               ndatamax=max(ndata,ndatamax)
            enddo
            if(ndatamax.gt.2) then
               write(*,*) 
     &            '*WARNING reading *PLASTIC: isotropic hardening'
               write(*,*) '         curve is possibly nonlinear for'
               write(*,*) '         the elastically anisotropic'
               write(*,*) '         material ',matname(nmat)(11:80)
            endif
         else
!
!           kinematic hardening
!
!           interpolating the plastic data at the elastic temperature
!           data points
!
            ndatamax=0
            do i=1,nelcon(2,nmat)
               t1l=elcon(0,i,nmat)
c               plconloc(1)=0.d0
c               plconloc(2)=0.d0
c               plconloc(3)=0.d0
c               plconloc(82)=nplkcon(1,nmat)+0.5d0
!
               if(nplkcon(0,nmat).eq.1) then
                  id=-1
               else
                  call ident2(plkcon(0,1,nmat),t1l,nplkcon(0,nmat),
     &                        2*npmat_+1,id)
               endif
!
               if(nplkcon(0,nmat).eq.0) then
                  continue
               elseif((nplkcon(0,nmat).eq.1).or.(id.eq.0).or.
     &                 (id.eq.nplkcon(0,nmat))) then
                  if(id.le.0) then
                     itemp=1
                  else
                     itemp=id
                  endif
                  kin=1
                  call plcopy(plkcon,nplkcon,plconloc,npmat_,ntmat_,
     &                 nmat,itemp,i,kin)
                  if((id.eq.0).or.(id.eq.nplkcon(0,nmat))) then
                  endif
               else
                  kin=1
                  call plmix(plkcon,nplkcon,plconloc,npmat_,ntmat_,
     &                 nmat,id+1,t1l,i,kin)
               endif
!
               ndata=int(plconloc(802))
               if(ndata.eq.1) then
                  elcon(10,i,nmat)=plconloc(402)
                  elcon(11,i,nmat)=0.d0
                  elcon(12,i,nmat)=0.d0
                  elcon(13,i,nmat)=-1.d0
                  elcon(14,i,nmat)=1.d0
               else
                  elcon(10,i,nmat)=plconloc(402)
                  elcon(11,i,nmat)=0.d0
                  elcon(12,i,nmat)=(plconloc(404)-plconloc(402))/
     &                             (plconloc(403)-plconloc(401))
                  elcon(13,i,nmat)=-1.d0
                  elcon(14,i,nmat)=1.d0
               endif
               ndatamax=max(ndata,ndatamax)
            enddo
            if(ndatamax.gt.2) then
               write(*,*) 
     &             '*WARNING reading *PLASTIC: kinematic hardening'
               write(*,*) '         curve is possibly nonlinear for'
               write(*,*) '         the elastically anisotropic'
               write(*,*) '         material ',matname(nmat)(11:80)
            endif
         endif
      endif
!
c      if(nelcon(1,nmat).eq.-114) then
c         write(*,*) 'anisotropic elasticity+viscoplasticity'
c         do i=1,nelcon(2,nmat)
c            write(*,*) (elcon(j,i,nmat),j=0,14)
c         enddo
c      endif
!
      return
      end

