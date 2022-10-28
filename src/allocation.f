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
      subroutine allocation(nload_,nforc_,nboun_,nk_,ne_,nmpc_,
     &     nset_,nalset_,nmat_,ntmat_,npmat_,norien_,nam_,nprint_,
     &     mi,ntrans_,set,meminset,rmeminset,ncs_,
     &     namtot_,ncmat_,memmpc_,ne1d,ne2d,nflow,jobnamec,irstrt,
     &     ithermal,nener,nstate_,irestartstep,inpc,ipoinp,inp,
     &     ntie_,nbody_,nprop_,ipoinpc,nevdamp_,npt_,nslavs,nkon_,mcs,
     &     mortar,ifacecount,nintpoint,infree,nheading_,nobject_,
     &     iuel,iprestr,nstam,ndamp,nef,nbounold,nforcold,nloadold,
     &     nbodyold,mpcend,irobustdesign,nfc_,ndc_)
!     
!     calculates a conservative estimate of the size of the 
!     fields to be allocated
!     
!     meminset=total # of terms in sets
!     rmeminset=total # of reduced terms (due to use of generate) in
!     sets
!     
!     nstate_ needs only be assigned for
!     a. restart (read from file)
!     b. initial conditions (defined by *depvar)
!     
      implicit none
!     
      logical igen,lin,frequency,cyclicsymmetry,composite,
     &     tabular,massflow,beamgeneralsection,johnsoncook
!     
      character*1 selabel,sulabel,inpc(*)
      character*5 llab
      character*8 label
      character*20 mpclabel
      character*81 set(*),noset,elset,slavset,mastset,noelset,submset,
     &     surface,slavsets,slavsett,mastsets,mastsett,surfset
      character*132 jobnamec(*),textpart(16)
!     
      integer nload_,nforc_,nboun_,nk_,ne_,nmpc_,nset_,nalset_,
     &     nmat_,ntmat_,npmat_,norien_,nam_,nprint_,kode,iline,
     &     istat,n,key,meminset(*),i,js,inoset,mi(*),ii,ipol,inl,
     &     ibounstart,ibounend,ibound,ntrans_,ntmatl,npmatl,ityp,l,
     &     ielset,nope,nteller,nterm,ialset(16),ncs_,rmeminset(*),
     &     islavset,imastset,namtot_,ncmat_,nconstants,memmpc_,j,ipos,
     &     maxrmeminset,ne1d,ne2d,necper,necpsr,necaxr,nesr,
     &     neb32,nn,nflow,nradiate,irestartread,irestartstep,icntrl,
     &     irstrt(*),ithermal(*),nener,nstate_,ipoinp(2,*),inp(3,*),
     &     ntie_,nbody_,nprop_,ipoinpc(0:*),nevdamp_,npt_,nentries,
     &     iposs,iposm,nslavs,nlayer,nkon_,nopeexp,iremove,mcs,
     &     ifacecount,nintpoint,mortar,infree(4),nheading_,icfd,
     &     multslav,multmast,nobject_,numnodes,iorientation,id,
     &     nuel,iuel(4,*),number,four,
     &     iprestr,nstam,ier,ndamp,nef,nbounold,nforcold,nloadold,
     &     nbodyold,mpcend,irobustdesign(3),iflag,network,
     &     nsubmodel,nfc_,ndc_
!     
      real*8 temperature,tempact,xfreq,tpinc,tpmin,tpmax
!     
      parameter(nentries=18)
!     
!     icfd=-1: initial value
!     =0: pure mechanical analysis
!     =1: pure CFD analysis
!     =2: mixed mechanical/cfd analysis
!     
!     mi(1): # of integration points
!     mi(2): # of dofs per node
!     mi(3): # of layers in the elements
!     
      icfd=-1
!     
      ier=0
!     
!     in the presence of mechanical steps the highest number
!     of DOF is at least 3
!     
      if(ithermal(2).ne.2) mi(2)=3
!     
!     initialisation of ipoinp
!     
      do i=1,nentries
        if(ipoinp(1,i).ne.0) then
          ipol=i
          inl=ipoinp(1,i)
          iline=inp(1,inl)-1
          exit
        endif
      enddo
!     
      istat=0
!
      nsubmodel=0
      nset_=0
      maxrmeminset=0
      necper=0
      necpsr=0
      necaxr=0
      nesr=0
      neb32=0
      nradiate=0
      nkon_=0
      nuel=0
!     
      four=4
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      loop: do
      if(istat.lt.0) then
        exit
      endif
!     
      if(textpart(1)(1:10).eq.'*AMPLITUDE') then
        nam_=nam_+1
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          namtot_=namtot_+4
        enddo
      elseif(textpart(1)(1:19).eq.'*BEAMGENERALSECTION') then
        mi(3)=max(mi(3),2)
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) then
            exit
          endif
          nprop_=nprop_+8
        enddo
      elseif(textpart(1)(1:12).eq.'*BEAMSECTION') then
        mi(3)=max(mi(3),2)
        beamgeneralsection=.false.
        do i=2,n
          if((textpart(i)(1:11).eq.'SECTION=BOX').or.
     &         (textpart(i)(1:11).eq.'SECTION=PIP').or.
     &         (textpart(i)(1:11).eq.'SECTION=GEN')) then
            beamgeneralsection=.true.
            exit
          endif
        enddo
        if(beamgeneralsection) then
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
              exit
            endif
            nprop_=nprop_+8
          enddo
        else
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &         inl,ipoinp,inp,ipoinpc)
        endif
      elseif(textpart(1)(1:10).eq.'*BOUNDARYF') then
        nam_=nam_+1
        namtot_=namtot_+1
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
!     
          read(textpart(3)(1:10),'(i10)',iostat=istat) ibounstart
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*BOUNDARYF%",ier)
            exit
          endif
!     
          if(textpart(4)(1:1).eq.' ') then
            ibounend=ibounstart
          else
            read(textpart(4)(1:10),'(i10)',iostat=istat) ibounend
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*BOUNDARYF%",ier)
              exit
            endif
          endif
          ibound=ibounend-ibounstart+1
          ibound=max(1,ibound)
!     
          read(textpart(1)(1:10),'(i10)',iostat=istat) l
          if(istat.eq.0) then
            nboun_=nboun_+ibound
            if(ntrans_.gt.0) then
              nmpc_=nmpc_+ibound
              memmpc_=memmpc_+4*ibound
              nk_=nk_+1
            endif
          else
            read(textpart(1)(1:80),'(a80)',iostat=istat) elset
            elset(81:81)=' '
            ipos=index(elset,' ')
!     
!     check for element set
!     
            elset(ipos:ipos)='E'
            call cident81(set,elset,nset_,id)
            i=nset_+1
            if(id.gt.0) then
              if(set(id).eq.elset) then
                i=id
                nboun_=nboun_+ibound*meminset(i)
                if(ntrans_.gt.0)then
                  nmpc_=nmpc_+ibound*meminset(i)
                  memmpc_=memmpc_+4*ibound*meminset(i)
                  nk_=nk_+meminset(i)
                endif
              endif
            endif
            if(i.gt.nset_) then
!     
!     check for facial surface
!     
              elset(ipos:ipos)='T'
              call cident81(set,elset,nset_,i)
              if(i.gt.0) then
                if(set(i).eq.elset) then
                  nboun_=nboun_+ibound*meminset(i)
                  if(ntrans_.gt.0)then
                    nmpc_=nmpc_+ibound*meminset(i)
                    memmpc_=memmpc_+4*ibound*meminset(i)
                    nk_=nk_+meminset(i)
                  endif
                endif
              endif
            endif
          endif
        enddo
      elseif(textpart(1)(1:9).eq.'*BOUNDARY') then
        nam_=nam_+1
        namtot_=namtot_+1
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
!     
          read(textpart(2)(1:10),'(i10)',iostat=istat) ibounstart
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*BOUNDARY%",ier)
            exit
          endif
!     
          if(textpart(3)(1:1).eq.' ') then
            ibounend=ibounstart
          else
            read(textpart(3)(1:10),'(i10)',iostat=istat) ibounend
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*BOUNDARY%",ier)
              exit
            endif
          endif
          ibound=ibounend-ibounstart+1
          ibound=max(1,ibound)
!     
          read(textpart(1)(1:10),'(i10)',iostat=istat) l
          if(istat.eq.0) then
            nboun_=nboun_+ibound
            if(ntrans_.gt.0) then
              nmpc_=nmpc_+ibound
              memmpc_=memmpc_+4*ibound
              nk_=nk_+1
            endif
          else
            read(textpart(1)(1:80),'(a80)',iostat=istat) noset
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
            call cident81(set,noset,nset_,i)
            if(i.gt.0) then
              if(set(i).eq.noset) then
                nboun_=nboun_+ibound*meminset(i)
                if(ntrans_.gt.0)then
                  nmpc_=nmpc_+ibound*meminset(i)
                  memmpc_=memmpc_+4*ibound*meminset(i)
                  nk_=nk_+meminset(i)
                endif
              endif
            endif
          endif
        enddo
      elseif(textpart(1)(1:4).eq.'*CFD') then
        iflag=1
        do i=2,n
          if(textpart(i)(1:10).eq.'TURBULENCE') then
            iflag=iflag+1
          endif
        enddo
        if(iflag.eq.2) mi(2)=max(mi(2),6)
        write(*,*) 'mi(2)= ',mi(2)
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:6).eq.'*CFLUX') then
        nam_=nam_+1
        namtot_=namtot_+1
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
!     
          read(textpart(1)(1:10),'(i10)',iostat=istat) l
          if(istat.eq.0) then
            nforc_=nforc_+1
          else
            read(textpart(1)(1:80),'(a80)',iostat=istat) noset
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
            call cident81(set,noset,nset_,i)
            if(i.gt.0) then
              if(set(i).eq.noset) then
                nforc_=nforc_+meminset(i)
              endif
            endif
          endif
        enddo
      elseif(textpart(1)(1:6).eq.'*CLOAD') then
        nam_=nam_+1
        namtot_=namtot_+1
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
!     
          read(textpart(1)(1:10),'(i10)',iostat=istat) l
          if(istat.eq.0) then
            if(ntrans_.eq.0) then
              nforc_=nforc_+1
            else
              nforc_=nforc_+3
            endif
          else
            read(textpart(1)(1:80),'(a80)',iostat=istat) noset
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
            call cident81(set,noset,nset_,i)
            if(i.gt.0) then
              if(set(i).eq.noset) then
                if(ntrans_.eq.0) then
                  nforc_=nforc_+meminset(i)
                else
                  nforc_=nforc_+3*meminset(i)
                endif
              endif
            endif
          endif
        enddo
      elseif((textpart(1)(1:13).eq.'*CONDUCTIVITY').or.
     &       (textpart(1)(1:8).eq.'*DENSITY').or.
     &       (textpart(1)(1:10).eq.'*EXPANSION').or.
     &       (textpart(1)(1:15).eq.'*FLUIDCONSTANTS').or.
     &       (textpart(1)(1:13).eq.'*SPECIFICHEAT').or.
     &       (textpart(1)(1:23).eq.'*ELECTRICALCONDUCTIVITY')) then
        ntmatl=0
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          ntmatl=ntmatl+1
          ntmat_=max(ntmatl,ntmat_)
        enddo
      elseif(textpart(1)(1:11).eq.'*CONSTRAINT') then
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          nobject_=nobject_+1
        enddo    
      elseif(textpart(1)(1:15).eq.'*CONTACTDAMPING') then
        ncmat_=max(8,ncmat_)
        ntmat_=max(1,ntmat_)
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:12).eq.'*CONTACTPAIR') then
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          ntie_=ntie_+1
        enddo
      elseif(textpart(1)(1:13).eq.'*CONTACTPRINT') then
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          nprint_=nprint_+n
        enddo
      elseif(textpart(1)(1:9).eq.'*COUPLING') then
        surface(1:1)=' '
        iorientation=0
        do i=2,n
          if(textpart(i)(1:8).eq.'SURFACE=') then
            surface=textpart(i)(9:88)
            ipos=index(surface,' ')
            surface(ipos:ipos)='T'
          elseif(textpart(i)(1:12).eq.'ORIENTATION=') then
            iorientation=1
          endif
        enddo
        if(surface(1:1).ne.' ') then
          surface(ipos:ipos)='T'
          call cident81(set,surface,nset_,i)
          if(i.gt.0) then
            if(set(i).eq.surface) then
              numnodes=8*meminset(i)
            endif
          endif
          surface(ipos:ipos)='S'
          call cident81(set,surface,nset_,i)
          if(i.gt.0) then
            if(set(i).eq.surface) then
              numnodes=meminset(i)
            endif
          endif
        endif
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:6).eq.'*CREEP') then
        ntmatl=0
        npmat_=max(2,npmat_)
        if(ncmat_.le.2) then
!     elastic isotropic
          ncmat_=max(9,ncmat_)
        else
!     elastic anisotropic
          ncmat_=max(19,ncmat_)
        endif
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          ntmatl=ntmatl+1
        enddo
        ntmat_=max(ntmatl,ntmat_)
      elseif(textpart(1)(1:16).eq.'*CYCLICHARDENING') then
        ntmatl=0
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          read(textpart(3)(1:20),'(f20.0)',iostat=istat) 
     &         temperature
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*CYCLIC HARDENING%",ier)
            exit
          endif
          if(ntmatl.eq.0) then
            npmatl=0
            ntmatl=ntmatl+1
            ntmat_=max(ntmatl,ntmat_)
            tempact=temperature
          elseif(temperature.ne.tempact) then
            npmatl=0
            ntmatl=ntmatl+1
            ntmat_=max(ntmatl,ntmat_)
            tempact=temperature
          endif
          npmatl=npmatl+1
          npmat_=max(npmatl,npmat_)
        enddo
      elseif(textpart(1)(1:20).eq.'*CYCLICSYMMETRYMODEL') then
!     
!     possible MPC's: static temperature, displacements(velocities)
!     and static pressure
!     
        nk_=nk_+1
!     
!     change on 11th of Dec. 2020
!     
        nmpc_=nmpc_+(mi(2)+1)*ncs_
        memmpc_=memmpc_+25*(mi(2)+1)*ncs_
        ntrans_=ntrans_+1
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
        enddo
      elseif(textpart(1)(1:8).eq.'*DAMPING') then
        do i=2,n
          if(textpart(i)(1:11).eq.'STRUCTURAL=') then
            ndamp=1
            exit
          endif
        enddo
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:8).eq.'*DASHPOT') then
        nmat_=nmat_+1
        frequency=.false.
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
        if((istat.lt.0).or.(key.eq.1)) then
          call inputerror(inpc,ipoinpc,iline,
     &         "*DASHPOT%",ier)
          cycle
        endif
        read(textpart(2)(1:20),'(f20.0)',iostat=istat)
     &       xfreq
        if(istat.gt.0) then
          call inputerror(inpc,ipoinpc,iline,
     &         "*DASHPOT%",ier)
          cycle
        endif
        if(xfreq.gt.0.d0) frequency=.true.
        iline=iline-1
        if(.not.frequency) then
          ntmatl=0
          ncmat_=max(2,ncmat_)
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            ntmatl=ntmatl+1
            ntmat_=max(ntmatl,ntmat_)
          enddo
        else
          ntmatl=0
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) 
     &           temperature
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*DASHPOT%",ier)
              exit
            endif
            if(ntmatl.eq.0) then
              npmatl=0
              ntmatl=ntmatl+1
              ntmat_=max(ntmatl,ntmat_)
              tempact=temperature
            elseif(temperature.ne.tempact) then
              npmatl=0
              ntmatl=ntmatl+1
              ntmat_=max(ntmatl,ntmat_)
              tempact=temperature
            endif
            npmatl=npmatl+1
            npmat_=max(npmatl,npmat_)
          enddo
          if(ncmat_.ge.9) ncmat_=max(19,ncmat_)
        endif
      elseif(textpart(1)(1:22).eq.'*DEFORMATIONPLASTICITY') then
        ncmat_=max(5,ncmat_)
        ntmatl=0
!
!     the deformation plasticity model does not have internal
!     variables. Nevertheless, space is reserved here for an
!     equivalent plastic strain (used in sensitivity procedures
!     to minimize plasticity)        
!
        nstate_=max(1,nstate_)
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          ntmatl=ntmatl+1
          ntmat_=max(ntmatl,ntmat_)
        enddo
      elseif(textpart(1)(1:7).eq.'*DEPVAR') then
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          read(textpart(1)(1:10),'(i10)',iostat=istat) l
          if(istat.lt.0) exit
          nstate_=max(l,nstate_)
        enddo
      elseif(textpart(1)(1:16).eq.'*DESIGNVARIABLES') then
        ntie_=ntie_+1   
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)        
      elseif(textpart(1)(1:15).eq.'*DESIGNRESPONSE') then
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          nobject_=nobject_+1
        enddo    
      elseif(textpart(1)(1:21).eq.'*DISTRIBUTINGCOUPLING') then
        nmpc_=nmpc_+3
        memmpc_=memmpc_+3
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
!     
          read(textpart(1)(1:10),'(i10)',iostat=istat) l
          if(istat.eq.0) then
            memmpc_=memmpc_+3
          else
            read(textpart(1)(1:80),'(a80)',iostat=istat) noset
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
            call cident81(set,noset,nset_,i)
            if(i.gt.0) then
              if(set(i).eq.noset) then
                memmpc_=memmpc_+3*meminset(i)
              endif
            endif
          endif
        enddo
      elseif(textpart(1)(2:13).eq.'DISTRIBUTING') then
c        irotation=0
c        itranslation=0
!     
!     translational dofs are always active
!     
        npt_=max(npt_,numnodes)
!
!       number of force prescriptions nfc_
!       number of forces nforc_
!
        nfc_=nfc_+3*numnodes
        nforc_=nforc_+3*numnodes
        ndc_=ndc_+3
!            
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
!     
          read(textpart(1)(1:10),'(i10)',iostat=istat) ibounstart
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*BOUNDARY%",ier)
            exit
          endif
!     
          if(textpart(2)(1:1).eq.' ') then
            ibounend=ibounstart
          else
            read(textpart(2)(1:10),'(i10)',iostat=istat) ibounend
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*BOUNDARY%",ier)
              exit
            endif
          endif
          ibounstart=max(4,ibounstart)
          ibounend=min(6,ibounend)
          ibound=max(0,ibounend-ibounstart+1)
!     
c          if(itranslation.eq.0) then
c!     
c!     translational dofs 3 MPC's + a two-term MPC for each
c!     participating node
c!     
c            npt_=max(npt_,numnodes)
c!     
c            nfc_=nfc_+3*numnodes
c            nforc_=nforc_+3*numnodes
c            ndc_=ndc_+3
c            itranslation=1
c          endif
!     
!     rotational dofs
!     
          if(ibound.gt.0) then
            ndc_=ndc_+ibound
          endif
        enddo
!     
      elseif(textpart(1)(2:18).eq.'DISTRIBUTIONTABLE') then
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(2:13).eq.'DISTRIBUTION') then
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          norien_=norien_+1
        enddo
      elseif((textpart(1)(1:6).eq.'*DLOAD').or.
     &       (textpart(1)(1:7).eq.'*DSLOAD').or.
     &       (textpart(1)(1:6).eq.'*DFLUX').or.
     &       (textpart(1)(1:9).eq.'*MASSFLOW').or.
     &       (textpart(1)(1:5).eq.'*FILM')) then
        massflow=.false.
        if((textpart(1)(1:5).ne.'*FILM').and.
     &       (textpart(1)(1:9).ne.'*MASSFLOW')) then
          nam_=nam_+1
          namtot_=namtot_+1
        elseif(textpart(1)(1:9).ne.'*MASSFLOW') then
          nam_=nam_+2
          namtot_=namtot_+2
        else
          massflow=.true.
        endif
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          read(textpart(2)(1:5),'(a5)',iostat=istat) llab
          if((llab.eq.'GRAV ').or.(llab.eq.'CENTR').or.
     &         (llab.eq.'NEWTO')) then
            nbody_=nbody_+1
            cycle
          endif
          read(textpart(1)(1:10),'(i10)',iostat=istat) l
          if(istat.eq.0) then
            nload_=nload_+1
            if(massflow) then
              nmpc_=nmpc_+1
              memmpc_=memmpc_+3
            endif
          else
            read(textpart(1)(1:80),'(a80)',iostat=istat) elset
            elset(81:81)=' '
            ipos=index(elset,' ')
!     
!     check for element set
!     
            elset(ipos:ipos)='E'
            call cident81(set,elset,nset_,id)
            i=nset_+1
            if(id.gt.0) then
              if(set(id).eq.elset) then
                i=id
                nload_=nload_+meminset(i)
                if(massflow) then
                  nmpc_=nmpc_+meminset(i)
                  memmpc_=memmpc_+3*meminset(i)
                endif
              endif
            endif
            if(i.gt.nset_) then
!     
!     check for facial surface
!     
              elset(ipos:ipos)='T'
              call cident81(set,elset,nset_,i)
              if(i.gt.0) then
                if(set(i).eq.elset) then
                  nload_=nload_+meminset(i)
                  if(massflow) then
                    nmpc_=nmpc_+meminset(i)
                    memmpc_=memmpc_+3*meminset(i)
                  endif
                endif
              endif
            endif
          endif
        enddo
      elseif((textpart(1)(1:8).eq.'*DYNAMIC').or.
     &       (textpart(1)(1:32).eq.'*COUPLEDTEMPERATURE-DISPLACEMENT')
     &       .or.
     &       (textpart(1)(1:34).eq.
     &       '*UNCOUPLEDTEMPERATURE-DISPLACEMENT'))then
!     
!     change of number of integration points except for a pure
!     CFD-calculation
!     
        if(icfd.ne.1) then
          if((mi(1).eq.1).or.(mi(1).eq.8)) then
            mi(1)=27
          elseif(mi(1).eq.4) then
            mi(1)=15
          elseif(mi(1).eq.2) then
            mi(1)=9
          endif
        endif
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:8).eq.'*ELPRINT') then
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          nprint_=nprint_+n
        enddo
      elseif(textpart(1)(1:8).eq.'*ELASTIC') then
        ntmatl=0
        ityp=2
        ncmat_=max(2,ncmat_)
        do i=2,n
          if(textpart(i)(1:5).eq.'TYPE=') then
            if(textpart(i)(6:8).eq.'ISO') then
              ityp=2
              ncmat_=max(2,ncmat_)
            elseif((textpart(i)(6:10).eq.'ORTHO').or.
     &             (textpart(i)(6:10).eq.'ENGIN')) then
              ityp=9
              ncmat_=max(9,ncmat_)
            elseif(textpart(i)(6:10).eq.'ANISO') then
              ityp=21
              ncmat_=max(21,ncmat_)
            endif
            exit
          endif
        enddo
        if(ityp.eq.2) then
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            ntmatl=ntmatl+1
          enddo
          ntmat_=max(ntmatl,ntmat_)
        elseif(ityp.eq.9) then
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            ntmatl=ntmatl+1
            iline=iline+1
          enddo
          ntmat_=max(ntmatl,ntmat_)
        elseif(ityp.eq.21) then
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            ntmatl=ntmatl+1
            iline=iline+2
          enddo
          ntmat_=max(ntmatl,ntmat_)
        endif
      elseif(textpart(1)(1:17).eq.'*ELECTROMAGNETICS') then
        mi(2)=max(mi(2),5)
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif((textpart(1)(1:8).eq.'*ELEMENT').and.
     &       (textpart(1)(1:14).ne.'*ELEMENTOUTPUT')) then
        ielset=0
!     
        loop1: do i=2,n
        if(textpart(i)(1:6).eq.'ELSET=') then
          elset=textpart(i)(7:86)
          elset(81:81)=' '
          ipos=index(elset,' ')
          elset(ipos:ipos)='E'
          ielset=1
          call cident81(set,elset,nset_,id)
          js=nset_+1
          if(id.gt.0) then
            if(set(id).eq.elset) js=id
          endif
          if(js.gt.nset_) then
            nset_=nset_+1
            do j=nset_,id+2,-1
              meminset(j)=meminset(j-1)
              rmeminset(j)=rmeminset(j-1)
              set(j)=set(j-1)
            enddo
            js=id+1
            set(js)=elset
            meminset(js)=0
            rmeminset(js)=0
          endif
        elseif(textpart(i)(1:5).eq.'TYPE=') then
          read(textpart(i)(6:13),'(a8)') label
          if(label.eq.'        ') then
            write(*,*) 
     &           '*ERROR in allocation: element type is lacking'
            write(*,*) '       '
            call inputerror(inpc,ipoinpc,iline,
     &           "*ELEMENT or *ELEMENT OUTPUT%",ier)
            exit
          endif
          if((label(1:2).eq.'DC').and.(label(1:7).ne.'DCOUP3D'))
     &         then
            label(1:7)=label(2:8)
            label(8:8)=' '
          endif
!     
          nopeexp=0
!     
          if(label.eq.'C3D20   ') then
            mi(1)=max(mi(1),27)
            nope=20
            nopeexp=20
          elseif(label(1:8).eq.'C3D20R  ') then
            mi(1)=max(mi(1),8)
            nope=20
            nopeexp=20
          elseif((label.eq.'C3D8R   ').or.(label.eq.'F3D8R   '))
     &           then
            mi(1)=max(mi(1),1)
            nope=8
            nopeexp=8
          elseif((label.eq.'C3D10   ').or.
     &           (label.eq.'C3D10T  ')) then
            mi(1)=max(mi(1),4)
            nope=10
            nopeexp=10
          elseif((label.eq.'C3D4    ').or.
     &           (label.eq.'F3D4    ')) then
            mi(1)=max(mi(1),1)
            nope=4
            nopeexp=4
          elseif(label.eq.'C3D15   ') then
            mi(1)=max(mi(1),9)
            nope=15
            nopeexp=15
          elseif(label.eq.'C3D6    ') then
            mi(1)=max(mi(1),2)
            nope=6
            nopeexp=6
          elseif(label.eq.'F3D6    ') then
            mi(1)=max(mi(1),1)
            nope=6
            nopeexp=6
          elseif((label.eq.'C3D8    ').or.(label.eq.'F3D8    '))
     &           then
            mi(1)=max(mi(1),8)
            nope=8
            nopeexp=8
c     Bernhardi start
          elseif(label.eq.'C3D8I   ') then
            mi(1)=max(mi(1),8)
            nope=8
            nopeexp=11
c     Bernhardi end
          elseif((label.eq.'CPE3    ').or.
     &           (label.eq.'CPS3    ').or.
     &           (label.eq.'CAX3    ').or.
     &           (label.eq.'M3D3    ').or.
     &           (label.eq.'S3      ')) then
            mi(1)=max(mi(1),2)
            nope=3
            nopeexp=9
          elseif((label.eq.'CPE4R   ').or.
     &           (label.eq.'CPS4R   ').or.
     &           (label.eq.'CAX4R   ').or.
     &           (label.eq.'M3D4R   ').or.
     &           (label.eq.'S4R     ')) then
            mi(1)=max(mi(1),1)
            nope=4
            nopeexp=12
          elseif((label.eq.'CPE4    ').or.
     &           (label.eq.'CPS4    ').or.
     &           (label.eq.'CAX4    ').or.
     &           (label.eq.'M3D4    ')) then
            mi(1)=max(mi(1),8)
            nope=4
            nopeexp=12
          elseif(label.eq.'S4      ') then
            mi(1)=max(mi(1),8)
            nope=4
!     modified into C3D8I (11 nodes)
            nopeexp=15
          elseif((label.eq.'CPE6    ').or.
     &           (label.eq.'CPS6    ').or.
     &           (label.eq.'CAX6    ').or.
     &           (label.eq.'M3D6    ').or.
     &           (label.eq.'S6      ')) then
            mi(1)=max(mi(1),9)
            nope=6
            nopeexp=21
          elseif((label.eq.'CPE8R   ').or.
     &           (label.eq.'CPS8R   ').or.
     &           (label.eq.'CAX8R   ').or.
     &           (label.eq.'M3D8R   ').or.
     &           (label.eq.'S8R     ')) then
            mi(1)=max(mi(1),8)
            nope=8
            nopeexp=28
          elseif((label.eq.'CPE8    ').or.
     &           (label.eq.'CPS8    ').or.
     &           (label.eq.'CAX8    ').or.
     &           (label.eq.'M3D8    ').or.
     &           (label.eq.'S8      ')) then
            mi(1)=max(mi(1),27)
            nope=8
            nopeexp=28
          elseif((label.eq.'B31     ').or.
     &           (label.eq.'B21     ').or.
     &           (label.eq.'T3D2    ').or.
     &           (label.eq.'T2D2    ')) then
            mi(1)=max(mi(1),8)
            mi(3)=max(mi(3),2)
            nope=2
!     modified into C3D8I (11 nodes)
            nopeexp=13
          elseif(label.eq.'B31R    ') then
            mi(1)=max(mi(1),1)
            nope=2
            nopeexp=10
          elseif((label.eq.'B32     ').or.
     &           (label.eq.'T3D3    ')) then
            mi(1)=max(mi(1),27)
            mi(3)=max(mi(3),2)
            nope=3
            nopeexp=23
          elseif(label.eq.'B32R    ') then
            mi(1)=max(mi(1),50)
            nope=3
            nopeexp=23
          elseif(label(1:8).eq.'DASHPOTA') then
            label='EDSHPTA1'
            nope=2
            nopeexp=2
          elseif(label(1:7).eq.'DCOUP3D') then
            nope=1
            nopeexp=1
          elseif(label(1:1).eq.'D') then
            nope=3
            nopeexp=3
            mi(2)=max(3,mi(2))
          elseif(label(1:7).eq.'SPRINGA') then
            mi(1)=max(mi(1),1)
            label='ESPRNGA1'
            nope=2
            nopeexp=2
          elseif(label(1:7).eq.'SPRING1') then
            mi(1)=max(mi(1),1)
            label='ESPRNG10'
            nope=1
            nopeexp=1
            ncmat_=max(3,ncmat_)
          elseif(label(1:7).eq.'SPRING2') then
            mi(1)=max(mi(1),1)
            label='ESPRNG21'
            nope=2
            nopeexp=2
            ncmat_=max(4,ncmat_)
          elseif(label.eq.'GAPUNI  ') then
            mi(1)=max(mi(1),1)
            label='ESPGAPA1'
            nope=2
            nopeexp=2
          elseif(label(1:4).eq.'MASS') then
            nope=1
            nopeexp=1
          elseif(label(1:1).eq.'U') then
!     
!     the number uniquely characterizes the
!     element name (consisting of 4 freely
!     selectable characters in position 2..5)
!     
            number=ichar(label(2:2))*256**3+
     &           ichar(label(3:3))*256**2+
     &           ichar(label(4:4))*256+
     &           ichar(label(5:5))
            nope=-1
            call nidentk(iuel,number,nuel,id,four)
            if(id.gt.0) then
              if(iuel(1,id).eq.number) then
                mi(1)=max(mi(1),iuel(2,id))
                mi(2)=max(mi(2),iuel(3,id))
                nope=iuel(4,id)
                nopeexp=nope
              endif
            endif
            if(nope.eq.-1) then
              write(*,*) '*ERROR reading *ELEMENT'
              write(*,*) '       nonexistent element type:'
              write(*,*) '       ',label
              call inputerror(inpc,ipoinpc,iline,
     &             "*ELEMENT%",ier)
              call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &             ipoinp,inp,ipoinpc)
              cycle loop
            endif
          endif
          if(label(1:1).eq.'F') then
            mi(2)=max(mi(2),4)
            if(icfd.eq.-1) then
              icfd=1
            elseif(icfd.eq.0) then
              icfd=2
            endif
          else
            if(icfd.eq.-1) then
              icfd=0
            elseif(icfd.eq.1) then
              icfd=2
            endif
          endif
        endif
      enddo loop1
!     
      loop2:do
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) exit
      read(textpart(1)(1:10),'(i10)',iostat=istat) i
      if(istat.gt.0) then
        call inputerror(inpc,ipoinpc,iline,
     &       "*ELEMENT or *ELEMENT OUTPUT%",ier)
        exit
      endif
c     Bernhardi start
c     space for incompatible mode nodes
      if(label(1:5).eq.'C3D8I') then
        nk_=nk_+3
      endif
c     Bernhardi end
      if(label(1:2).ne.'C3') then
        if(label(1:3).eq.'CPE') then
          necper=necper+1
        elseif(label(1:2).eq.'CP') then
          necpsr=necpsr+1
        elseif(label(1:1).eq.'C') then
          necaxr=necaxr+1
        elseif((label(1:1).eq.'S').or.
     &         ((label(1:1).eq.'M').and.(label(1:4).ne.'MASS')))
     &         then
          nesr=nesr+1
        elseif((label(1:1).eq.'B').or.
     &         (label(1:1).eq.'T')) then
          neb32=neb32+1
        elseif(label(1:1).eq.'D') then
          nflow=nflow+1
        elseif(label(1:1).eq.'F') then
          nef=nef+1
        endif
      endif
      nteller=n-1
      if(nteller.lt.nope) then
        do
          call getnewline(inpc,textpart,istat,n,key,iline,
     &         ipol,inl,ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit loop2
          if(nteller+n.gt.nope) n=nope-nteller
          nteller=nteller+n
          if(nteller.eq.nope) exit
        enddo
      endif
      ne_=max(ne_,i)
      nkon_=nkon_+nopeexp
      if(ielset.eq.1) then
        meminset(js)=meminset(js)+1
        rmeminset(js)=rmeminset(js)+1
      endif
c     !
c     !              up to 8 new mpc's with 22 terms in each mpc
c     !              (21 = 7 nodes x 3 dofs + inhomogeneous term)
c     !
      enddo loop2
      elseif((textpart(1)(1:5).eq.'*NSET').or.
     &     (textpart(1)(1:6).eq.'*ELSET')) then
        if(textpart(1)(1:5).eq.'*NSET')
     &       then
          noelset=textpart(2)(6:85)
          noelset(81:81)=' '
          ipos=index(noelset,' ')
          noelset(ipos:ipos)='N'
          kode=0
        else
          noelset=textpart(2)(7:86)
          noelset(81:81)=' '
          ipos=index(noelset,' ')
          noelset(ipos:ipos)='E'
          kode=1
        endif
!     
!     check whether new set name or old one
!     
        call cident81(set,noelset,nset_,id)
        nn=nset_+1
        if(id.gt.0) then
          if(set(id).eq.noelset) nn=id
        endif
        if(nn.gt.nset_) then
          nset_=nset_+1
          do j=nset_,id+2,-1
            meminset(j)=meminset(j-1)
            rmeminset(j)=rmeminset(j-1)
            set(j)=set(j-1)
          enddo
          nn=id+1
          set(nn)=noelset
          meminset(nn)=0
          rmeminset(nn)=0
        endif
!     
        if((n.gt.2).and.(textpart(3)(1:8).eq.'GENERATE')) then
          igen=.true.
        else
          igen=.false.
        endif
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          if(igen) then
            if(textpart(2)(1:1).eq.' ')
     &           textpart(2)=textpart(1)
            if(textpart(3)(1:1).eq.' ')
     &           textpart(3)='1        '
            do i=1,3
              read(textpart(i)(1:10),'(i10)',iostat=istat) 
     &             ialset(i)
              if(istat.gt.0) then
                call inputerror(inpc,ipoinpc,iline,
     &               "*NSET or *ELSET%",ier)
                exit
              endif
            enddo
            meminset(nn)=meminset(nn)+
     &           (ialset(2)-ialset(1))/ialset(3)+1
            rmeminset(nn)=rmeminset(nn)+3
          else
            do i=1,n
              read(textpart(i)(1:10),'(i10)',iostat=istat) 
     &             ialset(i)
              if(istat.gt.0) then
                noelset=textpart(i)(1:80)
                noelset(81:81)=' '
                ipos=index(noelset,' ')
                if(kode.eq.0) then
                  noelset(ipos:ipos)='N'
                else
                  noelset(ipos:ipos)='E'
                endif
                call cident81(set,noelset,nset_,j)
                if(j.gt.0) then
                  if(set(j).eq.noelset) then
                    meminset(nn)=meminset(nn)+
     &                   meminset(j)
                    rmeminset(nn)=rmeminset(nn)+
     &                   rmeminset(j)
                  endif
                endif
              else
                meminset(nn)=meminset(nn)+1
                rmeminset(nn)=rmeminset(nn)+1
              endif
            enddo
          endif
        enddo
      elseif((textpart(1)(1:9).eq.'*EQUATION').or.
     &       (textpart(1)(1:10).eq.'*EQUATIONF')) then
        iremove=0
        do i=2,n
          if(textpart(i)(1:6).eq.'REMOVE') iremove=1
        enddo
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if(iremove.eq.1) exit
          if((istat.lt.0).or.(key.eq.1)) exit 
          read(textpart(1)(1:10),'(i10)',iostat=istat) nterm
          if(ntrans_.eq.0) then
            nmpc_=nmpc_+1
            memmpc_=memmpc_+nterm
          else
            nmpc_=nmpc_+3
            memmpc_=memmpc_+3*nterm
          endif
          ii=0
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            ii=ii+n/3
            if(ii.eq.nterm) exit
          enddo
        enddo
      elseif(textpart(1)(1:13).eq.'*FLUIDSECTION') then
        nconstants=-1
        do i=2,n
          if(textpart(i)(1:10).eq.'CONSTANTS=') then
            read(textpart(i)(11:20),'(i10)',iostat=istat) 
     &           nconstants
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*FLUID SECTION%",ier)
              exit
            endif
            nprop_=nprop_+nconstants
            exit
          endif
        enddo
        if(nconstants.lt.0) nprop_=nprop_+65
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
        enddo
      elseif(textpart(1)(1:9).eq.'*FRICTION') then
!     
!     '8' is for Mortar.
!     
        ncmat_=max(8,ncmat_)
        ntmat_=max(1,ntmat_)
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:5).eq.'*GAP ') then
        nmat_=nmat_+1
        ncmat_=max(6,ncmat_)
        ntmat_=max(1,ntmat_)
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &         inl,ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
        enddo
      elseif(textpart(1)(1:15).eq.'*GAPCONDUCTANCE') then
        ntmatl=0
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &         inl,ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          read(textpart(3)(1:20),'(f20.0)',iostat=istat) 
     &         temperature
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*GAP CONDUCTANCE%",ier)
            exit
          endif
          if(ntmatl.eq.0) then
            npmatl=0
            ntmatl=ntmatl+1
            ntmat_=max(ntmatl,ntmat_)
            tempact=temperature
          elseif(temperature.ne.tempact) then
            npmatl=0
            ntmatl=ntmatl+1
            ntmat_=max(ntmatl,ntmat_)
            tempact=temperature
          endif
          npmatl=npmatl+1
          npmat_=max(npmatl,npmat_)
        enddo
      elseif(textpart(1)(1:18).eq.'*GAPHEATGENERATION') then
        ncmat_=max(11,ncmat_)
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:20).eq.'*GEOMETRICCONSTRAINT') then
        do   
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &      inl,ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          nobject_ = nobject_ + 1
        enddo
      elseif(textpart(1)(1:8).eq.'*HEADING') then
        if(nheading_.ne.0) then
          write(*,*) '*ERROR in allocation: more than 1'
          write(*,*) '       *HEADING card in the input deck'
          call inputerror(inpc,ipoinpc,iline,
     &         "*HEADING%",ier)
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          cycle loop
        endif
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          nheading_=nheading_+1
        enddo
      elseif(textpart(1)(1:13).eq.'*HYPERELASTIC') then
        ntmatl=0
        ityp=-7
        do i=2,n
          if(textpart(i)(1:12).eq.'ARRUDA-BOYCE') then
            ityp=-1
            ncmat_=max(3,ncmat_)
          elseif(textpart(i)(1:13).eq.'MOONEY-RIVLIN') then
            ityp=-2
            ncmat_=max(3,ncmat_)
          elseif(textpart(i)(1:8).eq.'NEOHOOKE') then
            ityp=-3
            ncmat_=max(2,ncmat_)
          elseif(textpart(i)(1:5).eq.'OGDEN') then
            ityp=-4
            ncmat_=max(3,ncmat_)
          elseif(textpart(i)(1:10).eq.'POLYNOMIAL') then
            ityp=-7
            ncmat_=max(3,ncmat_)
          elseif(textpart(i)(1:17).eq.'REDUCEDPOLYNOMIAL')
     &           then
            ityp=-10
            ncmat_=max(2,ncmat_)
          elseif(textpart(i)(1:11).eq.'VANDERWAALS') then
            ityp=-13
            ncmat_=max(5,ncmat_)
          elseif(textpart(i)(1:4).eq.'YEOH') then
            ityp=-14
            ncmat_=max(6,ncmat_)
          elseif(textpart(i)(1:2).eq.'N=') then
            if(textpart(i)(3:3).eq.'1') then
            elseif(textpart(i)(3:3).eq.'2') then
              if(ityp.eq.-4) then
                ityp=-5
                ncmat_=max(6,ncmat_)
              elseif(ityp.eq.-7) then
                ityp=-8
                ncmat_=max(7,ncmat_)
              elseif(ityp.eq.-10) then
                ityp=-11
                ncmat_=max(4,ncmat_)
              endif
            elseif(textpart(i)(3:3).eq.'3') then
              if(ityp.eq.-4) then
                ityp=-6
                ncmat_=max(9,ncmat_)
              elseif(ityp.eq.-7) then
                ityp=-9
                ncmat_=max(12,ncmat_)
              elseif(ityp.eq.-10) then
                ityp=-12
                ncmat_=max(6,ncmat_)
              endif
            endif
          endif
        enddo
        if((ityp.ne.-6).and.(ityp.ne.-9)) then
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            ntmatl=ntmatl+1
            ntmat_=max(ntmatl,ntmat_)
          enddo
        else
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            ntmatl=ntmatl+1
            ntmat_=max(ntmatl,ntmat_)
            iline=iline+1
          enddo
        endif
      elseif(textpart(1)(1:10).eq.'*HYPERFOAM') then
        ntmatl=0
        ityp=-15
        ncmat_=max(3,ncmat_)
        do i=2,n
          if(textpart(i)(1:2).eq.'N=') then
            if(textpart(i)(3:3).eq.'1') then
            elseif(textpart(i)(3:3).eq.'2') then
              ityp=-16
              ncmat_=max(6,ncmat_)
            elseif(textpart(i)(3:3).eq.'3') then
              ityp=-17
              ncmat_=max(9,ncmat_)
            endif
          endif
        enddo
        if(ityp.ne.-17) then
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            ntmatl=ntmatl+1
            ntmat_=max(ntmatl,ntmat_)
          enddo
        else
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            ntmatl=ntmatl+1
            ntmat_=max(ntmatl,ntmat_)
            iline=iline+1
          enddo
        endif
      elseif(textpart(1)(2:10).eq.'KINEMATIC') then
        npt_=max(npt_,numnodes)
!     
!     connection of rotational dofs in refnode to
!     translational dofs in rotational node
!     
        nk_=nk_+1
        nmpc_=nmpc_+3
        memmpc_=memmpc_+6
!     
!     local system
!     
        if(iorientation.ne.0) then
          nk_=nk_+2*numnodes
          nmpc_=nmpc_+3*numnodes
          memmpc_=memmpc_+3*6*numnodes
          nboun_=nboun_+3*numnodes
        endif
!     
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
!     
          read(textpart(1)(1:10),'(i10)',iostat=istat) ibounstart
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*BOUNDARY%",ier)
            exit
          endif
!     
          if(textpart(2)(1:1).eq.' ') then
            ibounend=ibounstart
          else
            read(textpart(2)(1:10),'(i10)',iostat=istat) ibounend
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*BOUNDARY%",ier)
              exit
            endif
          endif
          ibound=ibounend-ibounstart+1
          ibound=max(1,ibound)
          ibound=min(3,ibound)
!     
          if(iorientation.eq.0) then
            nk_=nk_+numnodes
            nmpc_=nmpc_+ibound*numnodes
            memmpc_=memmpc_+6*ibound*numnodes
            nboun_=nboun_+ibound*numnodes
          else
            nmpc_=nmpc_+ibound*numnodes
            memmpc_=memmpc_+ibound*6*numnodes
          endif
        enddo
      elseif(textpart(1)(1:21).eq.'*MAGNETICPERMEABILITY') then
        ntmatl=0
        ityp=2
        ncmat_=max(2,ncmat_)
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &         inl,ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          ntmatl=ntmatl+1
          ntmat_=max(ntmatl,ntmat_)
        enddo
      elseif(textpart(1)(1:5).eq.'*MASS') then
        nmat_=nmat_+1
        ntmat_=max(1,ntmat_)
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &       inl,ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:9).eq.'*MATERIAL') then
        nmat_=nmat_+1
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:13).eq.'*MODALDAMPING') then
        if(textpart(2)(1:8).ne.'RAYLEIGH') then
          nevdamp_=0
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            read(textpart(1)(1:10),'(i10)',iostat=istat) i
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*MODAL DAMPING%",ier)
              exit
            endif
            nevdamp_ = max(nevdamp_,i)
            read(textpart(2)(1:10),'(i10)',iostat=istat) i
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*MODAL DAMPING%",ier)
              exit
            endif
            nevdamp_ = max(nevdamp_,i)
          enddo
        else
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
        endif
      elseif(textpart(1)(1:12).eq.'*MODELCHANGE') then
        if(iprestr.ne.2) then
          do i=2,n
            if(textpart(i)(1:14).eq.'ADD=STRAINFREE') then
              iprestr=2
              exit
            elseif(textpart(i)(1:14).eq.'ADD=WITHSTRAIN') then
            elseif(textpart(i)(1:3).eq.'ADD') then
              iprestr=2
              exit
            elseif(textpart(i)(1:20).eq.'MECHSTRAINTORESIDUAL')
     &             then
              iprestr=2
              exit
            endif
          enddo
        endif
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:4).eq.'*MPC') then
        mpclabel='                    '
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          do i=1,n
            read(textpart(i)(1:10),'(i10)',iostat=istat) ialset(i)
            if(mpclabel.eq.'                    ') then
              mpclabel=textpart(i)(1:20)
              if((mpclabel(1:8).ne.'STRAIGHT').and.
     &             (mpclabel(1:4).ne.'PLANE')) then
                nk_=nk_+1
                nmpc_=nmpc_+1
                nboun_=nboun_+1
                memmpc_=memmpc_+1
              endif
            elseif(istat.gt.0) then
              noelset=textpart(i)(1:80)
              noelset(81:81)=' '
              ipos=index(noelset,' ')
              noelset(ipos:ipos)='N'
              call cident81(set,noelset,nset_,j)
              if(j.gt.0) then
                if(set(j).eq.noelset) then
                  if(mpclabel(1:8).eq.'STRAIGHT') then
                    nk_=nk_+2*meminset(j)
                    nmpc_=nmpc_+2*meminset(j)
                    nboun_=nboun_+2*meminset(j)
                    memmpc_=memmpc_+14*meminset(j)
                  elseif(mpclabel(1:5).eq.'PLANE') then
                    nk_=nk_+meminset(j)
                    nmpc_=nmpc_+meminset(j)
                    nboun_=nboun_+meminset(j)
                    memmpc_=memmpc_+13*meminset(j)
                  elseif(mpclabel(1:4).eq.'BEAM') then
                    memmpc_=memmpc_+3*meminset(j)
                  else
                    memmpc_=memmpc_+meminset(j)
                  endif
                endif
              endif
            else
              if(mpclabel(1:8).eq.'STRAIGHT') then
                nk_=nk_+2
                nmpc_=nmpc_+2
                nboun_=nboun_+2
                memmpc_=memmpc_+14
              elseif(mpclabel(1:5).eq.'PLANE') then
                nk_=nk_+1
                nmpc_=nmpc_+1
                nboun_=nboun_+1
                memmpc_=memmpc_+13
              elseif(mpclabel(1:4).eq.'BEAM') then
                memmpc_=memmpc_+3
              else
                memmpc_=memmpc_+1
              endif
            endif
          enddo
        enddo
      elseif(textpart(1)(1:11).eq.'*NETWORKMPC') then
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit 
          read(textpart(1)(1:10),'(i10)',iostat=istat) nterm
          if(ntrans_.eq.0) then
            nmpc_=nmpc_+1
            memmpc_=memmpc_+nterm
          else
            nmpc_=nmpc_+3
            memmpc_=memmpc_+3*nterm
          endif
          ii=0
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            ii=ii+n/3
            if(ii.eq.nterm) exit
          enddo
        enddo
      elseif((textpart(1)(1:5).eq.'*NODE').and.
     &       (textpart(1)(1:10).ne.'*NODEPRINT').and.
     &       (textpart(1)(1:9).ne.'*NODEFILE').and.
     &       (textpart(1)(1:11).ne.'*NODEOUTPUT')) then
        inoset=0
        loop3: do i=2,n
        if(textpart(i)(1:5).eq.'NSET=') then
          noset=textpart(i)(6:85)
          noset(81:81)=' '
          ipos=index(noset,' ')
          noset(ipos:ipos)='N'
          inoset=1
          call cident81(set,noset,nset_,id)
          js=nset_+1
          if(id.gt.0) then
            if(set(id).eq.noset) js=id
          endif
          if(js.gt.nset_) then
            nset_=nset_+1
            do j=nset_,id+2,-1
              meminset(j)=meminset(j-1)
              rmeminset(j)=rmeminset(j-1)
              set(j)=set(j-1)
            enddo
            js=id+1
            set(js)=noset
            meminset(js)=0
            rmeminset(js)=0
          endif
        endif
      enddo loop3
!     
      do
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
        if((istat.lt.0).or.(key.eq.1)) exit
        read(textpart(1)(1:10),'(i10)',iostat=istat) i
        if(istat.gt.0) then
          call inputerror(inpc,ipoinpc,iline,
     &         "*NODE or *NODE PRINT or *NODE FILE or *NODE OUTPUT%",
     &         ier)
          exit
        endif
        nk_=max(nk_,i)
        if(inoset.eq.1) then
          meminset(js)=meminset(js)+1
          rmeminset(js)=rmeminset(js)+1
        endif
      enddo
      elseif(textpart(1)(1:10).eq.'*NODEPRINT') then
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          nprint_=nprint_+n
        enddo
      elseif(textpart(1)(1:10).eq.'*OBJECTIVE') then
        nobject_=nobject_+1
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc) 
      elseif(textpart(1)(1:12).eq.'*ORIENTATION') then
        norien_=norien_+1
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
        enddo
      elseif(textpart(1)(1:8).eq.'*PLASTIC') then
!
        johnsoncook=.false.
        do i=2,n
          if(textpart(i)(1:21).eq.'HARDENING=JOHNSONCOOK') then
            johnsoncook=.true.
            ntmat_=max(1,ntmat_)
            ncmat_=max(11,ncmat_)
          endif
        enddo
        if(johnsoncook) then
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
        else
!     
          ntmatl=0
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) 
     &           temperature
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*PLASTIC%",ier)
              exit
            endif
            if(ntmatl.eq.0) then
              npmatl=0
              ntmatl=ntmatl+1
              ntmat_=max(ntmatl,ntmat_)
              tempact=temperature
            elseif(temperature.ne.tempact) then
              npmatl=0
              ntmatl=ntmatl+1
              ntmat_=max(ntmatl,ntmat_)
              tempact=temperature
            endif
            npmatl=npmatl+1
            npmat_=max(npmatl,npmat_)
          enddo
          if(ncmat_.ge.9) ncmat_=max(19,ncmat_)
        endif
      elseif(textpart(1)(1:19).eq.'*PRE-TENSIONSECTION') then
        surface(1:1)=' '
        do i=2,n
          if(textpart(i)(1:8).eq.'SURFACE=') then
            surface=textpart(i)(9:88)
            ipos=index(surface,' ')
            surface(ipos:ipos)='T'
            exit
          elseif(textpart(i)(1:8).eq.'ELEMENT=') then
            nmpc_=nmpc_+1
            memmpc_=memmpc_+7
            exit
          endif
        enddo
        if(surface(1:1).ne.' ') then
          call cident81(set,surface,nset_,i)
          if(i.gt.0) then
            if(set(i).eq.surface) then
!     
!     worst case: 8 nodes per element face
!     
              nk_=nk_+8*meminset(i)
              npt_=npt_+8*meminset(i)
!     
!     2 MPC's per node perpendicular to tension direction
!     + 1 thermal MPC per node
!     + 1 MPC per node in tension direction (the total of
!     which is divided into one global tension MPC and the
!     rest are MPC's specifying that the distance in tension
!     direction in all nodes should be the same)
!     
              nmpc_=nmpc_+32*meminset(i)+1
!     
!     6 terms per MPC perpendicular to tension direction
!     + 2 thermal terms per MPC
!     + 6 terms * # of nodes +1 parallel to tension
!     direction
!     + 12 terms per MPC parallel to tension direction
!     
              memmpc_=memmpc_+96*meminset(i)
     &             +16*meminset(i)
     &             +48*meminset(i)+1
     &             +12*(8*meminset(i)-1)
!     
            endif
          endif
        endif
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:8).eq.'*RADIATE') then
        nam_=nam_+2
        namtot_=namtot_+2
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          read(textpart(2)(1:5),'(a5)',iostat=istat) llab
          if((llab.eq.'GRAV ').or.(llab.eq.'CENTR')) exit
          read(textpart(1)(1:10),'(i10)',iostat=istat) l
          if(istat.eq.0) then
            nload_=nload_+1
            nradiate=nradiate+1
          else
            read(textpart(1)(1:80),'(a80)',iostat=istat) elset
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
            call cident81(set,elset,nset_,i)
            if(i.gt.0) then
              if(set(i).eq.elset) then
                nload_=nload_+meminset(i)
                nradiate=nradiate+meminset(i)
              endif
            endif
          endif
        enddo
      elseif(textpart(1)(1:13).eq.'*RATEDEPENDENT') then
        ntmat_=max(1,ntmat_)
        ncmat_=max(11,ncmat_)
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
!     
      elseif(textpart(1)(1:8).eq.'*RESTART') then
        irestartread=0
        irestartstep=0
        do i=1,n
          if(textpart(i)(1:4).eq.'READ') then
            irestartread=1
          endif
          if(textpart(i)(1:5).eq.'STEP=') then
            read(textpart(i)(6:15),'(i10)',iostat=istat) 
     &           irestartstep
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*RESTART%",ier)
              exit
            endif
          endif
        enddo
        if(irestartread.eq.1) then
          icntrl=1
          call restartshort(nset_,nload_,nbody_,nforc_,nboun_,nk_,ne_,
     &         nmpc_,nalset_,nmat_,ntmat_,npmat_,norien_,nam_,nprint_,
     &         mi,ntrans_,ncs_,namtot_,ncmat_,memmpc_,
     &         ne1d,ne2d,nflow,set,meminset,rmeminset,jobnamec,
     &         irestartstep,icntrl,ithermal,nener,nstate_,ntie_,
     &         nslavs,nkon_,mcs,nprop_,mortar,ifacecount,nintpoint,
     &         infree,nef,mpcend,nheading_,network,nfc_,ndc_,iprestr)
          irstrt(1)=-1
          nbounold=nboun_
          nforcold=nforc_
          nloadold=nload_
          nbodyold=nbody_
        else
        endif
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:18).eq.'*RETAINEDNODALDOFS') then
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
!     
          read(textpart(2)(1:10),'(i10)',iostat=istat) ibounstart
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*BOUNDARY%",ier)
            exit
          endif
!     
          if(textpart(3)(1:1).eq.' ') then
            ibounend=ibounstart
          else
            read(textpart(3)(1:10),'(i10)',iostat=istat) ibounend
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*BOUNDARY%",ier)
              exit
            endif
          endif
          ibound=ibounend-ibounstart+1
          ibound=max(1,ibound)
          ibound=min(3,ibound)
!     
          read(textpart(1)(1:10),'(i10)',iostat=istat) l
          if(istat.eq.0) then
            nboun_=nboun_+ibound
          else
            read(textpart(1)(1:80),'(a80)',iostat=istat) noset
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
            call cident81(set,noset,nset_,i)
            if(i.gt.0) then
              if(set(i).eq.noset) then
                nboun_=nboun_+ibound*meminset(i)
              endif
            endif
          endif
        enddo
      elseif(textpart(1)(1:10).eq.'*RIGIDBODY') then
        noset='
     &'
        elset='
     &'
        do i=2,n
          if(textpart(i)(1:5).eq.'NSET=')
     &         then
            noset=textpart(i)(6:85)
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
            exit
          elseif(textpart(i)(1:6).eq.'ELSET=')
     &           then
            elset=textpart(i)(7:86)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
            exit
          endif
        enddo
        if(noset(1:1).ne.' ') then
          call cident81(set,noset,nset_,i)
          if(i.gt.0) then
            if(set(i).eq.noset) then
              nk_=nk_+2+meminset(i)
              nmpc_=nmpc_+3*meminset(i)
              memmpc_=memmpc_+18*meminset(i)
              nboun_=nboun_+3*meminset(i)
            endif
          endif
        elseif(elset(1:1).ne.' ') then
          call cident81(set,elset,nset_,i)
          if(i.gt.0) then
            if(set(i).eq.elset) then
              nk_=nk_+2+20*meminset(i)
              nmpc_=nmpc_+60*meminset(i)
              memmpc_=memmpc_+360*meminset(i)
              nboun_=nboun_+60*meminset(i)
            endif
          endif
        endif
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:13).eq.'*ROBUSTDESIGN') then
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
        if((istat.lt.0).or.(key.eq.1)) exit
        irobustdesign(1)=1
      elseif(textpart(1)(1:16).eq.'*SECTIONPRINT') then
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          nprint_=nprint_+n
        enddo
      elseif(textpart(1)(1:13).eq.'*SHELLSECTION') then
        composite=.false.
        do i=2,n
          if(textpart(i)(1:9).eq.'COMPOSITE') then
            composite=.true.
            nlayer=0
          elseif(textpart(i)(1:6).eq.'ELSET=') then
            elset=textpart(i)(7:86)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
            call cident81(set,elset,nset_,id)
            js=nset_+1
            if(id.gt.0) then
              if(set(id).eq.elset) then
                js=id
              endif
            endif
          endif
        enddo
        if(composite) then
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
!     
!     conservative upper limit
!     "label" is not necessary the label of the
!     composite shell element
!     
              mi(1)=max(mi(1),8*nlayer)
              mi(3)=max(mi(3),nlayer)
              if(js.le.nset_) then
                nk_=nk_+20*nlayer*meminset(js)
                nkon_=nkon_+20*nlayer*meminset(js)
              endif
              exit
            endif
            nlayer=nlayer+1
          enddo
        else
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &         inl,ipoinp,inp,ipoinpc)
        endif
      elseif(textpart(1)(1:7).eq.'*SPRING') then
        nmat_=nmat_+1
        lin=.true.
        do i=2,n
          if(textpart(i)(1:9).eq.'NONLINEAR') then
            lin=.false.
            exit
          endif
        enddo
        if(lin) then
          ntmatl=0
          ncmat_=max(2,ncmat_)
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            ntmatl=ntmatl+1
            ntmat_=max(ntmatl,ntmat_)
          enddo
        else
          ntmatl=0
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) 
     &           temperature
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*SPRING%",ier)
              exit
            endif
            if(ntmatl.eq.0) then
              npmatl=0
              ntmatl=ntmatl+1
              ntmat_=max(ntmatl,ntmat_)
              tempact=temperature
            elseif(temperature.ne.tempact) then
              npmatl=0
              ntmatl=ntmatl+1
              ntmat_=max(ntmatl,ntmat_)
              tempact=temperature
            endif
            npmatl=npmatl+1
            npmat_=max(npmatl,npmat_)
          enddo
          if(ncmat_.ge.9) ncmat_=max(19,ncmat_)
        endif
      elseif(textpart(1)(1:5).eq.'*STEP') then
        if(nstam.eq.0) then
          do i=1,n
            if((textpart(i)(1:14).eq.'AMPLITUDE=STEP').or.
     &           (textpart(i)(1:14).eq.'AMPLITUDE=RAMP')) then
              nam_=nam_+2
              namtot_=namtot_+4
              nstam=1
              exit
            endif
          enddo
        endif
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:9).eq.'*SUBMODEL') then
        nsubmodel=nsubmodel+1
        ntie_=ntie_+1
        nam_=nam_+1
        namtot_=namtot_+4
!     
!     global element set
!     
        do j=2,n
          if(textpart(j)(1:12).eq.'GLOBALELSET=')
     &         then
            mastset(1:80)=textpart(j)(13:92)
            mastset(81:81)=' '
            ipos=index(mastset,' ')
            mastset(ipos:ipos)='E'
            call cident81(set,mastset,nset_,id)
            i=nset_+1
            if(id.gt.0) then
              if(set(id).eq.mastset) then
                i=id
              endif
            endif
            if(i.le.nset_) then
              meminset(i)=meminset(i)+meminset(i)
              rmeminset(i)=rmeminset(i)+meminset(i)
            endif
          elseif(textpart(j)(1:5).eq.'TYPE=') then
            if(textpart(j)(6:12).eq.'SURFACE') then
              selabel='T'
            else
              selabel='N'
            endif
          endif
        enddo
!     
!     local node or element face set
!     
        submset(1:8)='SUBMODEL'
        if(nsubmodel.lt.10) then
          submset(9:10)='00'
          write(submset(11:11),'(i1)') nsubmodel
        elseif(nsubmodel.lt.100) then
          submset(9:9)='0'
          write(submset(10:11),'(i2)') nsubmodel
        elseif(nsubmodel.lt.1000) then
          write(submset(9:11),'(i3)') nsubmodel
        else
          write(*,*) '*ERROR reading *SUBMODEL: no more than 999'
          write(*,*) '       submodels allowed'
          ier=1
          return
        endif
        submset(12:12)=selabel
        do i=13,81
          submset(i:i)=' '
        enddo
!        
        call cident81(set,submset,nset_,id)
        nset_=nset_+1
        do j=nset_,id+2,-1
          meminset(j)=meminset(j-1)
          rmeminset(j)=rmeminset(j-1)
          set(j)=set(j-1)
        enddo
        js=id+1
        set(js)=submset
        meminset(js)=0
        rmeminset(js)=0
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          read(textpart(1)(1:10),'(i10)',iostat=istat) ialset(1)
          if(istat.gt.0) then
            noset=textpart(1)(1:80)
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)=selabel
            call cident81(set,noset,nset_,i)
            if(i.gt.0) then
              if(set(i).eq.noset) then
                meminset(js)=meminset(js)+meminset(i)
!     
!     surfaces are stored in expanded form 
!     (no equivalent to generate)
!     
                rmeminset(js)=rmeminset(js)+meminset(i)
              endif
            endif
          else
            meminset(js)=meminset(js)+1
            rmeminset(js)=rmeminset(js)+1
          endif
        enddo
      elseif(textpart(1)(1:9).eq.'*SURFACE ') then
        sulabel='T'
        do i=2,n
          if(textpart(i)(1:5).eq.'NAME=')
     &         then
            surfset=textpart(i)(6:85)
            surfset(81:81)=' '
          elseif(textpart(i)(1:9).eq.'TYPE=NODE') then
            sulabel='S'
          endif
        enddo
        ipos=index(surfset,' ')
        surfset(ipos:ipos)=sulabel
        call cident81(set,surfset,nset_,id)
        nset_=nset_+1
        do j=nset_,id+2,-1
          meminset(j)=meminset(j-1)
          rmeminset(j)=rmeminset(j-1)
          set(j)=set(j-1)
        enddo
        js=id+1
        set(js)=surfset
        meminset(js)=0
        rmeminset(js)=0
!        
        if(sulabel.eq.'S') then
          selabel='N'
        else
          selabel='E'
        endif
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          read(textpart(1)(1:10),'(i10)',iostat=istat) ialset(1)
          if(istat.gt.0) then
            noset=textpart(1)(1:80)
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)=selabel
            call cident81(set,noset,nset_,i)
            if(i.gt.0) then
              if(set(i).eq.noset) then
                meminset(js)=meminset(js)+meminset(i)
!     
!     surfaces are stored in expanded form 
!     (no equivalent to generate)
!     
                rmeminset(js)=rmeminset(js)+meminset(i)
              endif
            endif
          else
            meminset(js)=meminset(js)+1
            rmeminset(js)=rmeminset(js)+1
          endif
        enddo
!     
!     for CFD-calculations: local coordinate systems are
!     stored as distributed load
!     
        if(icfd>0) nload_=nload_+rmeminset(js)
      elseif(textpart(1)(1:16).eq.'*SURFACEBEHAVIOR') then
        ncmat_=max(4,ncmat_)
        ntmat_=max(1,ntmat_)
        tabular=.false.
        do i=1,n
          if(textpart(i)(1:38).eq.'PRESSURE-OVERCLOSURE=TABULAR') 
     &         tabular=.true.
        enddo
        if(tabular) then
          ntmatl=0
          do
            call getnewline(inpc,textpart,istat,n,key,iline,
     &           ipol,inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) 
     &           temperature
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*SURFACE BEHAVIOR%",ier)
              exit
            endif
            if(ntmatl.eq.0) then
              npmatl=0
              ntmatl=ntmatl+1
              ntmat_=max(ntmatl,ntmat_)
              tempact=temperature
            elseif(temperature.ne.tempact) then
              npmatl=0
              ntmatl=ntmatl+1
              ntmat_=max(ntmatl,ntmat_)
              tempact=temperature
            endif
            npmatl=npmatl+1
            npmat_=max(npmatl,npmat_)
          enddo
        else
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
        endif
      elseif(textpart(1)(1:19).eq.'*SURFACEINTERACTION') then
        nmat_=nmat_+1
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:12).eq.'*TEMPERATURE') then
        nam_=nam_+1
        namtot_=namtot_+1
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:4).eq.'*TIE') then
        ntie_=ntie_+1
        cyclicsymmetry=.false.
        do i=1,n
          if((textpart(i)(1:14).eq.'CYCLICSYMMETRY').or.
     &         (textpart(i)(1:10).eq.'MULTISTAGE')) then
            cyclicsymmetry=.true.
          endif
        enddo
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
        if(.not.cyclicsymmetry) cycle
        if((istat.lt.0).or.(key.eq.1)) cycle
!     
        slavset=textpart(1)(1:80)
        slavset(81:81)=' '
        iposs=index(slavset,' ')
        slavsets=slavset
        slavsets(iposs:iposs)='S'
        slavsett=slavset
        slavsett(iposs:iposs)='T'
!     
        mastset=textpart(2)(1:80)
        mastset(81:81)=' '
        iposm=index(mastset,' ')
        mastsets=mastset
        mastsets(iposm:iposm)='S'
        mastsett=mastset
        mastsett(iposm:iposm)='T'
!     
        islavset=0
        imastset=0
!     
        call cident81(set,slavsets,nset_,i)
        if(i.gt.0) then
          if(set(i).eq.slavsets) then
            islavset=i
            multslav=1
          endif
        endif
        call cident81(set,slavsett,nset_,i)
        if(i.gt.0) then
          if(set(i).eq.slavsett) then
            islavset=i
            multslav=8
          endif
        endif
        call cident81(set,mastsets,nset_,i)
        if(i.gt.0) then
          if(set(i).eq.mastsets) then
            imastset=i
            multmast=1
          endif
        endif
        call cident81(set,mastsett,nset_,i)
        if(i.gt.0) then
          if(set(i).eq.mastsett) then
            imastset=i
            multmast=8
          endif
        endif
!
        if((islavset.ne.0).and.(imastset.ne.0)) then
          ncs_=ncs_+max(multslav*meminset(islavset),
     &         multmast*meminset(imastset))
        else
          write(*,*) '*ERROR in allocation: either the slave'
          write(*,*) '       surface or the master surface in a'
          write(*,*) '       cyclic symmetry *TIE option or both'
          write(*,*) '       do not exist or are no nodal surfaces'
          write(*,*) '       slave set:',slavset(1:iposs-1)
          write(*,*) '       master set:',mastset(1:iposm-1)
          call inputerror(inpc,ipoinpc,iline,
     &         "*TIE%",ier)
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          cycle loop
        endif
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:11).eq.'*TIMEPOINTS') then
        igen=.false.
        nam_=nam_+1
        do i=2,n
          if(textpart(i)(1:8).eq.'GENERATE') then
            igen=.true.
            exit
          endif
        enddo
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          if(igen)then
            if(n.lt.3)then
              write(*,*)'*ERROR in allocation:'
              call inputerror(inpc,ipoinpc,iline,
     &             "*TIME POINTS%",ier)
              exit
            else
              read(textpart(1)(1:20),'(f20.0)',iostat=istat) 
     &             tpmin
              if(istat.gt.0) then
                call inputerror(inpc,ipoinpc,iline,
     &               "*TIME POINTS%",ier)
                exit
              endif
              read(textpart(2)(1:20),'(f20.0)',iostat=istat) 
     &             tpmax
              if(istat.gt.0) then
                call inputerror(inpc,ipoinpc,iline,
     &               "*TIME POINTS%",ier)
                exit
              endif
              read(textpart(3)(1:20),'(f20.0)',iostat=istat) 
     &             tpinc
              if(istat.gt.0) then
                call inputerror(inpc,ipoinpc,iline,
     &               "*TIME POINTS%",ier)
                exit
              endif
!     
              if((tpinc.le.0).or.(tpmin.ge.tpmax)) then
                write(*,*) '*ERROR in allocation:'
                call inputerror(inpc,ipoinpc,iline,
     &               "*TIME POINTS%",ier)
                exit
              else
                namtot_=namtot_+2+INT((tpmax-tpmin)/tpinc)
              endif

            endif
          else
            namtot_=namtot_+8
          endif
        enddo
      elseif(textpart(1)(1:10).eq.'*TRANSFORM') then
        ntrans_=ntrans_+1
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
        enddo
      elseif(textpart(1)(1:11).eq.'*TRANSFORMF') then
        ntrans_=ntrans_+1
        surface(1:1)=' '
        do i=2,n
          if(textpart(i)(1:8).eq.'SURFACE=') then
            surface=textpart(i)(9:88)
            ipos=index(surface,' ')
            surface(ipos:ipos)='T'
            exit
          endif
        enddo
        if(surface(1:1).ne.' ') then
          call cident81(set,surface,nset_,i)
          if(i.gt.0) then
            if(set(i).eq.surface) then
              nload_=nload_+meminset(i)
            endif
          endif
        endif
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
        enddo
      elseif(textpart(1)(1:12).eq.'*USERELEMENT') then
        call userelements(textpart,n,iuel,nuel,inpc,ipoinpc,iline,
     &       ier,ipoinp,inp,inl,ipol)
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      elseif(textpart(1)(1:13).eq.'*USERMATERIAL') then
        ntmatl=0
        do i=2,n
          if(textpart(i)(1:10).eq.'CONSTANTS=') then
            read(textpart(i)(11:20),'(i10)',iostat=istat) 
     &           nconstants
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*USER MATERIAL%",ier)
              exit
            endif
            ncmat_=max(nconstants,ncmat_)
            exit
          endif
        enddo
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          ntmatl=ntmatl+1
          ntmat_=max(ntmatl,ntmat_)
          do i=2,nconstants/8+1
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
          enddo
        enddo
      elseif(textpart(1)(1:12).eq.'*USERSECTION') then
        nconstants=0
        do i=2,n
          if(textpart(i)(1:10).eq.'CONSTANTS=') then
            read(textpart(i)(11:20),'(i10)',iostat=istat) 
     &           nconstants
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*USER SECTION%",ier)
              exit
            endif
            nprop_=nprop_+nconstants
            exit
          endif
        enddo
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
        enddo
      else
!     
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      endif
      enddo loop
!     
      do i=1,nset_
        nalset_=nalset_+rmeminset(i)
        maxrmeminset=max(maxrmeminset,rmeminset(i))
      enddo
!     
!     extra space needed for rearrangement in elements.f and
!     noelsets.f
!     
      nalset_=nalset_+maxrmeminset
!     
      nmpc_=nmpc_+1
      memmpc_=memmpc_+1
!     
      if(irstrt(1).eq.0) then
        ne1d=neb32
        ne2d=necper+necpsr+necaxr+nesr
      endif
!     
!     introducing a fake tie for axisymmetric elements
!     (needed for cavity radiation)
!     
      if(necaxr.gt.0) ntie_=max(1,ntie_)
!     
!     providing space for the expansion of shell and beam elements
!     to genuine volume elements (no distinction is made between
!     linear and quadratic elements. The worst case (quadratic)
!     is taken
!     
      nk_=nk_+3*8*ne2d+8*3*ne1d
      if(ne1d.gt.0) then
        nboun_=nboun_*9
        nforc_=nforc_*9
      elseif(ne2d.gt.0) then
        nboun_=nboun_*4
        nforc_=nforc_*4
      endif
!     
!     providing for rigid nodes (knots)
!     
!     number of knots: 8*ne2d+3*ne1d
!     number of expanded nodes: 3*8*ne2d+8*3*ne1d
!     
!     number of extra nodes (1 rotational node and
!     1 expansion node per knot
!     and one inhomogeneous term node per expanded node)
!     
      nk_=nk_+(2+3)*8*ne2d+(2+8)*3*ne1d
!     
!     number of equations (3 per expanded node)
!     
      nmpc_=nmpc_+3*(3*8*ne2d+8*3*ne1d)
!     
!     number of terms: 9 per equation
!     
      memmpc_=memmpc_+9*3*(3*8*ne2d+8*3*ne1d)
!     
!     number of SPC's: 1 per DOF per expanded node
!     
      nboun_=nboun_+3*(3*8*ne2d+8*3*ne1d)
!     
!     temperature DOF in knots
!     
      nmpc_=nmpc_+(3*8*ne2d+8*3*ne1d)
      memmpc_=memmpc_+2*(3*8*ne2d+8*3*ne1d)
!     
!     extra MPCs to avoid undefinid rotation of rigid body nodes
!     lying on a line
!     
      nmpc_=nmpc_+3*8*ne2d+8*3*ne1d
      memmpc_=memmpc_+3*(3*8*ne2d+8*3*ne1d)
!     
!     extra nodes for the radiation boundary conditions
!     
      nk_=nk_+nradiate
!     
!     each layer in each shell has a local orientation
!     
      norien_=norien_+nesr*mi(3)
!     
      write(*,*)
      write(*,*) ' The numbers below are estimated upper bounds'
      write(*,*)
      write(*,*) ' number of:'
      write(*,*)
      write(*,*) '  nodes: ',nk_
      write(*,*) '  elements: ',ne_
      write(*,*) '  one-dimensional elements: ',ne1d
      write(*,*) '  two-dimensional elements: ',ne2d
      write(*,*) '  integration points per element: ',mi(1)
      write(*,*) '  degrees of freedom per node: ',mi(2)
      write(*,*) '  layers per element: ',mi(3)
      write(*,*)
      write(*,*) '  distributed facial loads: ',nload_
      write(*,*) '  distributed volumetric loads: ',nbody_
      write(*,*) '  concentrated loads: ',nforc_
      write(*,*) '  single point constraints: ',nboun_
      write(*,*) '  multiple point constraints: ',nmpc_
      write(*,*) '  terms in all multiple point constraints: ',memmpc_
      write(*,*) '  tie constraints: ',ntie_
      write(*,*) '  dependent nodes tied by cyclic constraints: ',ncs_
      write(*,*) '  dependent nodes in pre-tension constraints: ',npt_
      write(*,*)
      write(*,*) '  sets: ',nset_
      write(*,*) '  terms in all sets: ',nalset_
      write(*,*)
      write(*,*) '  materials: ',nmat_
      write(*,*) '  constants per material and temperature: ',ncmat_
      write(*,*) '  temperature points per material: ',ntmat_
      write(*,*) '  plastic data points per material: ',npmat_
      write(*,*)
      write(*,*) '  orientations: ',norien_
      write(*,*) '  amplitudes: ',nam_
      write(*,*) '  data points in all amplitudes: ',namtot_
      write(*,*) '  print requests: ',nprint_
      write(*,*) '  transformations: ',ntrans_
      write(*,*) '  property cards: ',nprop_
      write(*,*)
!     
      if(ier.eq.1) then
        write(*,*) '*ERROR in allocation: at least one fatal'
        write(*,*) '       error message while reading the'
        write(*,*) '       input deck: CalculiX stops.'
        write(*,*)
        call exit(201)
      endif
!     
      return
      end
