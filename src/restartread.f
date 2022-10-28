!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine restartread(istep,nset,nload,nforc, nboun,nk,ne,
     &     nmpc,nalset,nmat,ntmat_,npmat_,norien,nam,nprint,mi,
     &     ntrans,ncs_,namtot,ncmat_,mpcfree,maxlenmpc,
     &     ne1d,ne2d,nflow,nlabel,iplas,
     &     nkon,ithermal,nmethod,iperturb,nstate_,nener,set,istartset,
     &     iendset,ialset,co,kon,ipkon,lakon,nodeboun,ndirboun,iamboun,
     &     xboun,ikboun,ilboun,ipompc,nodempc,coefmpc,labmpc,ikmpc,
     &     ilmpc,nodeforc,ndirforc,iamforc,xforc,ikforc,ilforc,
     &     nelemload,iamload,sideload,xload,elcon,nelcon,rhcon,nrhcon,
     &     alcon,nalcon,alzero,plicon,nplicon,plkcon,nplkcon,orname,
     &     orab,ielorien,trab,inotr,amname,amta,namta,t0,t1,iamt1,veold,
     &     ielmat,matname,prlab,prset,filab,vold,nodebounold,
     &     ndirbounold,xbounold,xforcold,xloadold,t1old,eme,
     &     iponor,xnor,knor,thicke,offset,iponoel,inoel,rig,
     &     shcon,nshcon,cocon,ncocon,ics,sti,
     &     ener,xstate,jobnamec,infree,irestartstep,prestr,iprestr,
     &     cbody,ibody,xbody,nbody,xbodyold,ttime,qaold,cs,mcs,
     &     output,physcon,ctrl,typeboun,fmpc,tieset,ntie,tietol,nslavs,
     &     t0g,t1g,nprop,ielprop,prop,mortar,nintpoint,ifacecount,
     &     islavsurf,pslavsurf,clearini,irstrt,vel,nef,velo,veloo,
     &     ne2boun,heading,network,nfc,ndc,coeffc,ikdc,edc)
!     
      implicit none
!     
      character*1 typeboun(*)
      character*4 output
      character*6 prlab(*)
      character*8 lakon(*)
      character*20 labmpc(*),sideload(*)
      character*66 heading(*)
      character*80 orname(*),amname(*),matname(*),version
      character*81 set(*),prset(*),tieset(*),cbody(*)
      character*87 filab(*)
      character*132 fnrstrt,jobnamec(*)
!     
      integer istep,nset,nload,nforc,nboun,nk,ne,nmpc,nalset,nmat,
     &     ntmat_,npmat_,norien,nam,nprint,mi(*),ntrans,ncs_,nheading_,
     &     namtot,ncmat_,mpcfree,ne1d,ne2d,nflow,nlabel,iplas,nkon,
     &     ithermal(*),nmethod,iperturb(*),nstate_,istartset(*),
     &     ialset(*),kon(*),ipkon(*),nodeboun(*),ndirboun(*),iamboun(*),
     &     ikboun(*),ilboun(*),ipompc(*),nodempc(*),ikmpc(*),ilmpc(*),
     &     nodeforc(*),ndirforc(*),iamforc(*),ikforc(*),ilforc(*),
     &     nelemload(*),iamload(*),nelcon(*),mt,nprop,ielprop(*),
     &     nrhcon(*),nalcon(*),nplicon(*),nplkcon(*),ielorien(*),
     &     inotr(*),mortar,nintpoint,ifacecount,islavsurf(*),iendset(*),
     &     namta(*),iamt1(*),ielmat(*),nodebounold(*),ndirbounold(*),
     &     iponor(*),knor(*),iponoel(*),inoel(*),rig(*),nfc,ndc,
     &     nshcon(*),ncocon(*),ics(*),infree(*),i,ipos,ikdc(*),
     &     nener,irestartstep,istat,iprestr,irstrt(*),
     &     maxlenmpc,mcs,mpcend,ntie,ibody(*),nbody,nslavs,nef,
     &     ne2boun(*),memmpc_,network
!     
      real*8 co(*),xboun(*),coefmpc(*),xforc(*),xload(*),elcon(*),
     &     rhcon(*),alcon(*),alzero(*),plicon(*),plkcon(*),orab(*),
     &     trab(*),amta(*),t0(*),t1(*),veold(*),tietol(*),pslavsurf(*),
     &     vold(*),xbounold(*),xforcold(*),xloadold(*),t1old(*),eme(*),
     &     xnor(*),thicke(*),offset(*),t0g(*),t1g(*),clearini(*),
     &     shcon(*),cocon(*),sti(*),ener(*),xstate(*),prestr(*),ttime,
     &     qaold(2),physcon(*),ctrl(*),cs(*),fmpc(*),xbody(*),coeffc(*),
     &     xbodyold(*),prop(*),vel(*),velo(*),veloo(*),edc(*)
!     
      ipos=index(jobnamec(1),char(0))
      fnrstrt(1:ipos-1)=jobnamec(1)(1:ipos-1)
      fnrstrt(ipos:ipos+3)=".rin"
      do i=ipos+4,132
        fnrstrt(i:i)=' '
      enddo
!     
      open(15,file=fnrstrt,ACCESS='SEQUENTIAL',FORM='UNFORMATTED',
     &     err=15)
!     
      do
!     
        read(15,iostat=istat) version
!     
        if(istat.lt.0) then
          if(irestartstep.eq.0) then
!     
!     reading the last step
!     
            irestartstep=istep
            close(15)
            open(15,file=fnrstrt,ACCESS='SEQUENTIAL',
     &           FORM='UNFORMATTED',err=15)
            read(15) version
          else
            write(*,*) '*ERROR reading *RESTART,READ: requested step'
            write(*,*) '       is not in the restart file'
            call exit(201)
          endif
        endif
        read(15)istep
!     
!     set size
!     
        read(15)nset
        read(15)nalset
!     
!     load size
!     
        read(15)nload
        read(15)nbody
        read(15)nforc
        read(15)nboun
        read(15)nflow
!     
!     mesh size
!     
        read(15)nk
        read(15)ne
        read(15)nef
        read(15)nkon
        read(15)(mi(i),i=1,3)
        mt=mi(2)+1
!     
!     constraint size
!     
        read(15)nmpc
        read(15)mpcend
        read(15)maxlenmpc
        read(15)memmpc_
!     
!     coupling, distributing
!     
        read(15)nfc
        read(15)ndc
!     
!     material size
!     
        read(15)nmat
        read(15)ntmat_
        read(15)npmat_
        read(15)ncmat_
!     
!     property info
!     
        read(15)nprop
!     
!     transformation size
!     
        read(15)norien
        read(15)ntrans
!     
!     amplitude size
!     
        read(15)nam
        read(15)namtot
!     
!     print size
!     
        read(15)nprint
        read(15)nlabel
!     
!     tie size
!     
        read(15)ntie
!     
!     cyclic symmetry size
!     
        read(15)ncs_
        read(15)mcs
!     
!     1d and 2d element size
!     
        read(15)ne1d 
        read(15)ne2d 
        read(15)(infree(i),i=1,4)
!     
!     procedure info
!     
        read(15)nmethod
        read(15)(iperturb(i),i=1,2)
        read(15)nener
        read(15)iplas
        read(15)ithermal(1)
        read(15)nstate_
        read(15)nslavs
        read(15)iprestr
        read(15)mortar
        if(mortar.eq.1) then
          read(15)ifacecount
          read(15)nintpoint
        endif
        read(15)nheading_
        read(15)network
!     
        if(istep.eq.irestartstep) exit
!     
!     skipping the next entries until the requested step is found
!     
        call skip(nset,nalset,nload,nbody,nforc,nboun,nk,ne,nkon,
     &       mi(1),nmpc,mpcend,nmat,ntmat_,npmat_,ncmat_,norien,ntrans,
     &       nam,nprint,nlabel,ncs_,ne1d,ne2d,infree,nmethod,
     &       iperturb,nener,ithermal,nstate_,iprestr,mcs,ntie,
     &       nslavs,nprop,mortar,ifacecount,nintpoint,nef,nheading_,
     &       nfc,ndc)
!     
      enddo
!     
!     sets
!     
      read(15)(set(i),i=1,nset)
      read(15)(istartset(i),i=1,nset)
      read(15)(iendset(i),i=1,nset)
      do i=1,nalset
        read(15)ialset(i)
      enddo
!     
!     header lines
!     
      read(15)(heading(i),i=1,nheading_)
!     
!     mesh
!     
      read(15)(co(i),i=1,3*nk)
      read(15)(kon(i),i=1,nkon)
      read(15)(ipkon(i),i=1,ne)
      read(15)(lakon(i),i=1,ne)
!     
!     single point constraints
!     
      read(15)(nodeboun(i),i=1,nboun)
      read(15)(ndirboun(i),i=1,nboun)
      read(15)(typeboun(i),i=1,nboun)
      read(15)(xboun(i),i=1,nboun)
      read(15)(ikboun(i),i=1,nboun)
      read(15)(ilboun(i),i=1,nboun)
      if(nam.gt.0) read(15)(iamboun(i),i=1,nboun)
      read(15)(nodebounold(i),i=1,nboun)
      read(15)(ndirbounold(i),i=1,nboun)
      read(15)(xbounold(i),i=1,nboun)
!     
!     multiple point constraints
!     
      read(15)(ipompc(i),i=1,nmpc)
      read(15)(labmpc(i),i=1,nmpc)
      read(15)(ikmpc(i),i=1,nmpc)
      read(15)(ilmpc(i),i=1,nmpc)
      read(15)(fmpc(i),i=1,nmpc)
      read(15)(nodempc(i),i=1,3*mpcend)
      read(15)(coefmpc(i),i=1,mpcend)
      mpcfree=mpcend+1
!     
!     force constraints
!     
      if(nfc.gt.0) then
        read(15)(coeffc(i),i=1,7*nfc)
        read(15)(ikdc(i),i=1,ndc)
        read(15)(edc(i),i=1,12*ndc)
      endif
!     
!     point forces
!     
      read(15)(nodeforc(i),i=1,2*nforc)
      read(15)(ndirforc(i),i=1,nforc)
      read(15)(xforc(i),i=1,nforc)
      read(15)(ikforc(i),i=1,nforc)
      read(15)(ilforc(i),i=1,nforc)
      if(nam.gt.0) read(15)(iamforc(i),i=1,nforc)
      read(15)(xforcold(i),i=1,nforc)
!     
!     distributed loads
!     
      read(15)(nelemload(i),i=1,2*nload)
      read(15)(sideload(i),i=1,nload)
      read(15)(xload(i),i=1,2*nload)
      if(nam.gt.0) read(15)(iamload(i),i=1,2*nload)
      read(15)(xloadold(i),i=1,2*nload)
      read(15)(cbody(i),i=1,nbody)
      read(15)(ibody(i),i=1,3*nbody)
      read(15)(xbody(i),i=1,7*nbody)
      read(15)(xbodyold(i),i=1,7*nbody)
!     
!     prestress
!     
      if(iprestr.gt.0) read(15) (prestr(i),i=1,6*mi(1)*ne)
!     
!     labels
!     
      read(15)(prlab(i),i=1,nprint)
      read(15)(prset(i),i=1,nprint)
      read(15)(filab(i),i=1,nlabel)
!     
!     elastic constants
!     
      read(15)(elcon(i),i=1,(ncmat_+1)*ntmat_*nmat)
      read(15)(nelcon(i),i=1,2*nmat)
!     
!     density
!     
      read(15)(rhcon(i),i=1,2*ntmat_*nmat)
      read(15)(nrhcon(i),i=1,nmat)
!     
!     specific heat
!     
      read(15)(shcon(i),i=1,4*ntmat_*nmat)
      read(15)(nshcon(i),i=1,nmat)
!     
!     conductivity
!     
      read(15)(cocon(i),i=1,7*ntmat_*nmat)
      read(15)(ncocon(i),i=1,2*nmat)
!     
!     expansion coefficients
!     
      read(15)(alcon(i),i=1,7*ntmat_*nmat)
      read(15)(nalcon(i),i=1,2*nmat)
      read(15)(alzero(i),i=1,nmat)
!     
!     physical constants
!     
      read(15)(physcon(i),i=1,14)
!     
!     plastic data
!     
      if(npmat_.ne.0)then
        read(15)(plicon(i),i=1,(2*npmat_+1)*ntmat_*nmat)
        read(15)(nplicon(i),i=1,(ntmat_+1)*nmat)
        read(15)(plkcon(i),i=1,(2*npmat_+1)*ntmat_*nmat)
        read(15)(nplkcon(i),i=1,(ntmat_+1)*nmat)
      endif
!     
!     material orientation
!     
      if(norien.ne.0)then
        read(15)(orname(i),i=1,norien)
        read(15)(orab(i),i=1,7*norien)
        read(15)(ielorien(i),i=1,mi(3)*ne)
      endif
!     
!     fluid section properties
!     
      if(nprop.ne.0) then
        read(15)(ielprop(i),i=1,ne)
        read(15)(prop(i),i=1,nprop)
      endif
!     
!     transformations
!     
      if(ntrans.ne.0)then
        read(15)(trab(i),i=1,7*ntrans)
        read(15)(inotr(i),i=1,2*nk)
      endif
!     
!     amplitudes
!     
      if(nam.gt.0)then
        read(15)(amname(i),i=1,nam)
        read(15)(namta(i),i=1,3*nam-1)
        read(15) namta(3*nam)
        read(15)(amta(i),i=1,2*namta(3*nam-1))
      endif
!     
!     temperatures
!     
      if(ithermal(1).gt.0)then
        read(15)(t0(i),i=1,nk)
        read(15)(t1(i),i=1,nk)
        if((ne1d.gt.0).or.(ne2d.gt.0))then
          read(15)(t0g(i),i=1,2*nk)
          read(15)(t1g(i),i=1,2*nk)
        endif
        if(nam.gt.0) read(15)(iamt1(i),i=1,nk)
        read(15)(t1old(i),i=1,nk)
      endif
!     
!     materials
!     
      read(15)(matname(i),i=1,nmat)
      read(15)(ielmat(i),i=1,mi(3)*ne)
!     
!     temperature, displacement, static pressure, velocity and acceleration
!     
      read(15)(vold(i),i=1,mt*nk)
      if((nmethod.eq.4).or.((nmethod.eq.1).and.(iperturb(1).ge.2))) then
        read(15)(veold(i),i=1,mt*nk)
      endif
!     
!     CFD results at the element centers
!     
      if(nef.gt.0) then
        read(15)(vel(i),i=1,8*nef)
        read(15)(velo(i),i=1,8*nef)
        read(15)(veloo(i),i=1,8*nef)
      endif
!     
!     1d and 2d elements
!     
      if((ne1d.gt.0).or.(ne2d.gt.0))then
        read(15)(iponor(i),i=1,2*nkon)
        read(15)(xnor(i),i=1,infree(1))
        read(15)(knor(i),i=1,infree(2))
        read(15)(thicke(i),i=1,mi(3)*nkon)
        read(15)(offset(i),i=1,2*ne)
        read(15)(iponoel(i),i=1,infree(4))
        read(15)(inoel(i),i=1,3*(infree(3)-1))
        read(15)(rig(i),i=1,infree(4))
        read(15)(ne2boun(i),i=1,2*infree(4))
      endif
!     
!     tie constraints
!     
      if(ntie.gt.0) then
        read(15)(tieset(i),i=1,3*ntie)
        read(15)(tietol(i),i=1,4*ntie)
      endif
!     
!     cyclic symmetry
!     
      if(ncs_.gt.0)then
        read(15)(ics(i),i=1,ncs_)
      endif
      if(mcs.gt.0) then
        read(15)(cs(i),i=1,17*mcs)
      endif
!     
!     integration point variables
!     
      read(15)(sti(i),i=1,6*mi(1)*ne)
      read(15)(eme(i),i=1,6*mi(1)*ne)
      if(nener.eq.1) then
        if(mortar.ne.1) then
          read(15)(ener(i),i=1,2*mi(1)*(ne+nslavs))
        else
          read(15)(ener(i),i=1,2*mi(1)*ne)
        endif
      endif
      if(nstate_.gt.0)then
        if(mortar.ne.1) then
          read(15)(xstate(i),i=1,nstate_*mi(1)*(ne+nslavs))
        else
          read(15)(xstate(i),i=1,nstate_*mi(1)*(ne+nintpoint))
        endif
      endif
!     
!     face-to-face penalty contact variables
!     
      if(mortar.eq.1) then
        read(15) (islavsurf(i),i=1,2*ifacecount+2)
        read(15) (pslavsurf(i),i=1,3*nintpoint)
        read(15) (clearini(i),i=1,3*9*ifacecount)
      endif
!     
!     control parameters
!     
      read(15) (ctrl(i),i=1,52)
      read(15) (qaold(i),i=1,2)
      read(15) output
      read(15) ttime
!     
!     restart parameters
!     
      read(15) (irstrt(i),i=1,2)
!     
      close(15)
!     
      return
!     
 15   write(*,*) '*ERROR reading *RESTART,READ: could not open file ',
     &     fnrstrt
      call exit(201)
      end
