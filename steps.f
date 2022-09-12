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
      subroutine steps(inpc,textpart,iperturb,iprestr,nbody,
     &  nforc,nload,ithermal,t0,t1,nk,irstrt,istep,istat,n,jmax,ctrl,
     &  iline,ipol,inl,ipoinp,inp,newstep,ipoinpc,network,
     &  iamplitudedefault,amname,nam,nam_,namta,amta,namtot,
     &  nstam,ier,namtot_,physcon)
!
!     reading the input deck: *STEP
!
      implicit none
!
      character*1 inpc(*)
      character*80 amname(*)
      character*132 textpart(16)
!
      integer iperturb(*),nforc,nload,ithermal(*),nk,istep,istat,n,key,
     &  i,j,iprestr,jmax(2),irstrt(*),iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),namtot_,
     &  newstep,nbody,ipoinpc(0:*),network,iamplitudedefault,nam,
     &  nam_,namta(3,*),namtot,nstam,ier
!
      real*8 t0(*),t1(*),ctrl(*),amta(2,*),physcon(*)
!
      if(newstep.eq.1) then
         write(*,*) '*ERROR reading *STEP: *STEP statement detected'
         write(*,*) '       within step ',istep
         ier=1
         return
      else
         newstep=1
      endif
!
      if(iperturb(1).lt.2) iperturb(1)=0
      if(irstrt(1).lt.0) irstrt(1)=0
      istep=istep+1
      jmax(1)=100
      jmax(2)=10000
      physcon(14)=0.d0
!
!     adding a ramp and a step amplitude if AMPLITUDE is active in
!     any step
!
      if((istep.eq.1).and.(nstam.eq.1)) then
!
!        adding a ramp amplitude
!
         nam=nam+1
         if(nam.gt.nam_) then
            write(*,*) 
     &           '*ERROR reading *STEP: increase nam_'
            ier=1
            return
         endif
         namta(3,nam)=nam
         amname(nam)(1:15)='RAMP12357111317'
         do j=16,80
            amname(nam)(j:j)=' '
         enddo
         namta(1,nam)=namtot+1
         namtot=namtot+1
         if(namtot.gt.namtot_) then
            write(*,*) 
     &           '*ERROR reading *STEP: increase namtot_'
            ier=1
            return
         endif
         amta(1,namtot)=0.d0
         amta(2,namtot)=0.d0
         namtot=namtot+1
         if(namtot.gt.namtot_) then
            write(*,*) 
     &           '*ERROR reading *STEP: increase namtot_'
            ier=1
            return
         endif
         amta(1,namtot)=1.d20
         amta(2,namtot)=1.d20
         namta(2,nam)=namtot
!     
!        adding a step amplitude
!     
         nam=nam+1
         if(nam.gt.nam_) then
            write(*,*) 
     &           '*ERROR reading *STEP: increase nam_'
            ier=1
            return
         endif
         namta(3,nam)=nam
         amname(nam)(1:15)='STEP12357111317'
         do j=16,80
            amname(nam)(j:j)=' '
         enddo
         namta(1,nam)=namtot+1
         namtot=namtot+1
         if(namtot.gt.namtot_) then
            write(*,*) 
     &           '*ERROR reading *STEP: increase namtot_'
            ier=1
            return
         endif
         amta(1,namtot)=0.d0
         amta(2,namtot)=1.d0
         namtot=namtot+1
         if(namtot.gt.namtot_) then
            write(*,*) 
     &           '*ERROR reading *STEP: increase namtot_'
            ier=1
            return
         endif
         amta(1,namtot)=1.d20
         amta(2,namtot)=1.d0
         namta(2,nam)=namtot
      endif
!     
      do i=2,n
         if(textpart(i)(1:12).eq.'PERTURBATION') then
            iperturb(1)=1
            write(*,*) '*INFO reading *STEP: nonlinear geometric'
            write(*,*) '      effects are turned on'
            write(*,*)
!
!           removing the present loading (check!!)
!
            nforc=0
            iprestr=0
            if((ithermal(1).eq.1).or.(ithermal(1).eq.3)) then
               do j=1,nk
                  t1(j)=t0(j)
               enddo
            endif
!
         elseif((textpart(i)(1:6).eq.'NLGEOM').and.
     &          (textpart(i)(7:9).ne.'=NO')) then
!
!           geometrically nonlinear calculations
!
            iperturb(2)=1
            write(*,*) '*INFO reading *STEP: nonlinear geometric'
            write(*,*) '      effects are turned on'
            write(*,*)
            if(iperturb(1).eq.0) then
               iperturb(1)=2
            elseif(iperturb(1).eq.1) then
               write(*,*) 
     &            '*ERROR reading *STEP: PERTURBATION and NLGEOM'
               write(*,*) '       are mutually exclusive; '
               call inputerror(inpc,ipoinpc,iline,
     &              "*STEP%",ier)
               return
            endif
!
         elseif(textpart(i)(1:9).eq.'NLGEOM=NO') then
!
!           geometrically linear calculations
!           iperturb(1) is not changed, i.e. iterations will
!           occur (e.g. they may be needed due to contact,
!           nonlinear material laws....)
!
            iperturb(2)=0
            write(*,*) '*INFO reading *STEP: nonlinear geometric'
            write(*,*) '      effects are turned off'
            write(*,*)
!
         elseif(textpart(i)(1:4).eq.'INC=') then
!
!           maximum number of increments
!
            read(textpart(i)(5:14),'(i10)',iostat=istat) jmax(1)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*STEP%",ier)
               return
            endif
!
         elseif(textpart(i)(1:5).eq.'INCF=') then
!
!           maximum number of fluid increments
!
            read(textpart(i)(6:15),'(i10)',iostat=istat) jmax(2)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*STEP%",ier)
               return
            endif
         elseif(textpart(i)(1:14).eq.'THERMALNETWORK') then
            if(istep.ne.1) then
               write(*,*) '*ERROR reading *STEP'
               write(*,*) '       THERMAL NETWORK can only be used'
               write(*,*) '       on the first step'
               ier=1
               return
            endif
!
!           purely thermal network, i.e. Dx-elements are defined
!           for thermal purposes only (calculation of heat transfer
!           coefficient based on Reynolds number et...)
!
            network=1
         elseif(textpart(i)(1:9).eq.'AMPLITUDE') then
            if(textpart(i)(1:14).eq.'AMPLITUDE=RAMP') then
               do j=1,nam
                  if(amname(j)(1:15).eq.'RAMP12357111317') then
                     iamplitudedefault=j
                     exit
                  endif
               enddo
            elseif(textpart(i)(1:14).eq.'AMPLITUDE=STEP') then
               do j=1,nam
                  if(amname(j)(1:15).eq.'STEP12357111317') then
                     iamplitudedefault=j
                     exit
                  endif
               enddo
            endif
         elseif(textpart(i)(1:15).eq.'SHOCKSMOOTHING=') then
!
!           reading the shock smoothing parameter for compressible
!           cfd-calculations
!
            read(textpart(i)(16:35),'(f20.0)',iostat=istat) physcon(14)
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &              "*STEP%",ier)
         else
            write(*,*) 
     &          '*WARNING reading *STEP: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*STEP%")
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end




