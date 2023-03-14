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
      subroutine calinput(co,nk,kon,ipkon,lakon,nkon,
     &     ne,nodeboun,ndirboun,xboun,nboun,
     &     ipompc,nodempc,coefmpc,nmpc,nmpc_,nodeforc,ndirforc,xforc,
     &     nforc,
     &     nforc_,nelemload,sideload,xload,nload,nload_,
     &     nprint,prlab,prset,mpcfree,nboun_,
     &     mei,set,istartset,iendset,ialset,nset,nalset,elcon,nelcon,
     &     rhcon,
     &     nrhcon,alcon,nalcon,alzero,t0,t1,
     &     matname,ielmat,orname,orab,ielorien,amname,amta,namta,nam,
     &     nmethod,iamforc,iamload,iamt1,
     &     ithermal,iperturb,istat,istep,nmat,ntmat_,norien,
     &     prestr,iprestr,isolver,fei,veold,timepar,
     &     xmodal,filab,jout,nlabel,idrct,jmax,
     &     iexpl,alpha,iamboun,plicon,nplicon,plkcon,
     &     nplkcon,iplas,npmat_,mi,nk_,trab,inotr,ntrans,ikboun,
     &     ilboun,ikmpc,ilmpc,ics,dcs,ncs_,namtot_,cs,nstate_,ncmat_,
     &     mcs,labmpc,iponor,xnor,knor,thickn,thicke,ikforc,ilforc,
     &     offset,iponoel,inoel,rig,infree,nshcon,shcon,cocon,ncocon,
     &     physcon,nflow,ctrl,maxlenmpc,ne1d,
     &     ne2d,nener,vold,nodebounold,ndirbounold,xbounold,
     &     xforcold,xloadold,t1old,eme,sti,ener,xstate,jobnamec,
     &     irstrt,ttime,qaold,output,typeboun,inpc,
     &     ipoinp,inp,tieset,tietol,ntie,fmpc,cbody,ibody,xbody,
     &     nbody,nbody_,xbodyold,nam_,ielprop,nprop,nprop_,prop,itpamp,
     &     iviewfile,ipoinpc,nslavs,t0g,t1g,network,cyclicsymmetry,
     &     idefforc,idefload,idefbody,mortar,ifacecount,islavsurf,
     &     pslavsurf,clearini,heading,iaxial,nobject,objectset,nprint_,
     &     iuel,nuel_,nodempcref,coefmpcref,ikmpcref,memmpcref_,
     &     mpcfreeref,maxlenmpcref,memmpc_,isens,namtot,nstam,dacon,
     &     vel,nef,velo,veloo,ne2boun,itempuser,irobustdesign,
     &     irandomtype,randomval,nfc,nfc_,coeffc,ikdc,ndc,ndc_,edc,
     &     coini)
!     
      implicit none
!     
!     nmethod: -1: visco (=static+creep) 
!               0: no analysis 
!               1: static
!               2: frequency 
!               3: buckling 
!               4: linear dynamic
!               5: steady state dynamics
!               6: Coriolis frequency calculation
!               7: flutter frequency calculation
!               8: magnetostatics
!               9: magnetodynamics (inductive heating)
!              10: electromagnetic eigenvalue problems
!              11: superelement creation
!              12: sensitivity
!              13: Green function calculation
!              14: robustness w.r.t. geometric tolerances
!              15: crack propagation
!              16: constrained sensitivity
!      
!     iprestr: 0: no residual stresses; 1: residual stresses;
!     2; residual strains
!     iperturb: 0:no perturbation; 1:perturbation; 2: nonlinear
!     geometric analysis; 3: material and geometrical
!     nonlinearities
!     istep: step number
!     istat: error indicator (<0:EOF, >0: input error)
!     
      logical boun_flag,cload_flag,dload_flag,temp_flag,elprint_flag,
     &     nodeprint_flag,elfile_flag,nodefile_flag,contactfile_flag,
     &     dflux_flag,cflux_flag,film_flag,radiate_flag,out3d,
     &     solid,sectionprint_flag,contactprint_flag,pretension,
     &     beamgeneralsection,objective_flag,constraint_flag
!     
      character*1 typeboun(*),inpc(*)
      character*4 output
      character*6 prlab(*)
      character*8 lakon(*)
      character*20 labmpc(*),sideload(*)
      character*66 heading(*)
      character*80 matname(*),orname(*),amname(*)
      character*81 set(*),prset(*),tieset(3,*),cbody(*),objectset(4,*)
      character*87 filab(*)
      character*132 jobnamec(*),textpart(16)
!     
      integer kon(*),nodeboun(*),ndirboun(*),ipompc(*),nodempc(3,*),
     &     nodeforc(2,*),ndirforc(*),nelemload(2,*),iaxial,j,mi(*),
     &     istartset(*),iendset(*),ialset(*),ipkon(*),ics(*),nodedep,
     &     nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),nodeind,
     &     ielorien(mi(3),*),icomposite,nsubmodel,mortar,
     &     namta(3,*),iamforc(*),iamload(2,*),iamt1(*),ipoinpc(0:*),
     &     iamboun(*),inotr(2,*),ikboun(*),ilboun(*),ikmpc(*),ilmpc(*),
     &     iponor(2,*),knor(*),ikforc(*),ilforc(*),iponoel(*),
     &     inoel(3,*),infree(4),ixfree,ikfree,inoelfree,iponoelmax,
     &     rig(*),nshcon(*),ncocon(2,*),nodebounold(*),ielprop(*),nprop,
     &     nprop_,maxsectors,irestartread,
     &     ndirbounold(*),ipoinp(2,*),inp(3,*),nintpoint,ifacecount,
     &     ianisoplas,ifile_output,ichangefriction,nslavs,
     &     nalset,nalset_,nmat,nmat_,ntmat_,norien,norien_,
     &     islavsurf(2,*),
     &     nmethod,nk,ne,nboun,nmpc,nmpc_,mpcfree,i,istat,n,
     &     key,nk_,ne_,nboun_,ncs_,namtot_,nstate_,iviewfile,
     &     isolver,ithermal(*),iperturb(*),iprestr,istep,mei(4),nkon,
     &     nprint,nload,nload_,nforc,nforc_,nlabel,imat,
     &     nset,nset_,nprint_,nam,nam_,jout(2),ncmat_,itpamp,
     &     ierror,idrct,jmax(2),iexpl,iplas,npmat_,ntrans,ntrans_,
     &     M_or_SPC,nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),nflow,
     &     ne1d,ne2d,nener,irstrt(*),ii,maxlenmpc,inl,ipol,network,
     &     iline,mcs,ntie,ntie_,lprev,newstep,nbody,nbody_,ibody(3,*),
     &     cyclicsymmetry,idefforc(*),idefload(*),idefbody(*),
     &     ichangesurfacebehavior,nobject,ibasemotion,iuel(4,*),nuel_,
     &     nodempcref(3,*),ikmpcref(*),memmpcref_,mpcfreeref,
     &     maxlenmpcref,memmpc_,isens,iamplitudedefault,namtot,
     &     nstam,ier,nef,ne2boun(2,*),itempuser(*),irobustdesign(3),
     &     irandomtype(*),iparentel,nfc,nfc_,ikdc(*),ndc,ndc_
!     
      real*8 co(3,*),xboun(*),coefmpc(*),xforc(*),fmpc(*),
     &     xload(2,*),alzero(*),offset(2,*),prop(*),pslavsurf(3,*),
     &     elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),clearini(3,9,*),
     &     alcon(0:6,ntmat_,*),thicke(mi(3),*),thickn(2,*),xnor(*),
     &     t1(*),orab(7,*),prestr(6,mi(1),*),amta(2,*),dacon(*),
     &     veold(0:mi(2),*),t0(*),plicon(0:2*npmat_,ntmat_,*),
     &     plkcon(0:2*npmat_,ntmat_,*),trab(7,*),dcs(*),
     &     shcon(0:3,ntmat_,*),cocon(0:6,ntmat_,*),timepar(*),
     &     ctrl(*),vold(0:mi(2),*),xbounold(*),xforcold(*),
     &     xloadold(*),t1old(*),eme(*),sti(*),ener(*),
     &     xstate(nstate_,mi(1),*),ttime,qaold(2),cs(17,*),tietol(4,*),
     &     xbody(7,*),xbodyold(7,*),t0g(2,*),t1g(2,*),
     &     fei(3),tinc,tper,xmodal(*),tmin,tmax,tincf,
     &     alpha(*),physcon(*),coefmpcref(*),vel(nef,*),velo(*),
     &     veloo(*),randomval(2,*),coeffc(0:6,*),edc(12,*),coini(3,*)
!     
      save solid,ianisoplas,out3d,pretension
!     
      integer nentries
      parameter(nentries=18)
!
      irestartread=0
      newstep=0
      iviewfile=0
      ichangefriction=0
      ichangesurfacebehavior=0
      icomposite=0
      ibasemotion=0
      iamplitudedefault=0
      ier=0
!     
      maxsectors=1
      if(mcs.ne.0) then
        do i=1,mcs
          maxsectors=max(maxsectors,int(cs(1,i)))
        enddo
      endif
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
      ixfree=infree(1)
      ikfree=infree(2)
      inoelfree=infree(3)
      iponoelmax=infree(4)
!     
      iexpl=0
!     
!     the following flag is used to check whether any SPC's or MPC's
!     are used before transformation definitions
!     
      M_or_SPC=0
!     
!     the flags indicate whether some specific keyword cards already
!     occurred (needed to determine the effect of OP=NEW or to check
!     whether the element or nodal output selection should be reset)
!     
      boun_flag=.false.
      cload_flag=.false.
      dload_flag=.false.
      temp_flag=.false.
      elprint_flag=.false.
      nodeprint_flag=.false.
      sectionprint_flag=.false.
      contactprint_flag=.false.
      contactfile_flag=.false.
      elfile_flag=.false.
      nodefile_flag=.false.
      film_flag=.false.
      dflux_flag=.false.
      radiate_flag=.false.
      cflux_flag=.false.
      objective_flag=.false.
      constraint_flag=.false.
!     
      if(istep.eq.0) then
!     
!     initializing the maxima
!     
        ne_=ne
        nset_=nset
        nalset_=nalset
        nmat_=nmat
        norien_=norien
        ntrans_=ntrans
        ntie_=ntie
!     
        nmethod=0
!     
        ne=0
c        nset=0
        nalset=0
        nmat=0
        norien=0
        ntrans=0
        ntie=0
        nsubmodel=0
        mortar=-2
!     
        imat=0
        lprev=0
!     
        do i=1,ne_
          ipkon(i)=-1
        enddo
!     
        if((ne1d.gt.0).or.(ne2d.gt.0)) then
          do i=1,nlabel
            filab(i)='    I '
          enddo
          out3d=.false.
        else
          do i=1,nlabel
            filab(i)='      '
          enddo
          out3d=.true.
        endif
!     
        solid=.false.
        ianisoplas=0
        pretension=.false.
!     
      endif
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      loop: do
!     
      if(istat.lt.0) then
        write(*,*)
        write(*,*) 'Job finished'
        write(*,*)
        return
      endif
!     
      if(textpart(1)(1:10).eq.'*AMPLITUDE') then
        call amplitudes(inpc,textpart,amname,amta,namta,nam,
     &       nam_,namtot_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,
     &       inp,ipoinpc,namtot,ier)
!     
      elseif(textpart(1)(1:11).eq.'*BASEMOTION') then
        call basemotions(inpc,textpart,amname,nam,ibasemotion,
     &       xboun,ndirboun,iamboun,typeboun,nboun,istep,istat,n,
     &       iline,ipol,inl,ipoinp,inp,ipoinpc,iamplitudedefault,
     &       ier,xmodal,nmethod)
!     
      elseif(textpart(1)(1:19).eq.'*BEAMGENERALSECTION') then
        write(*,*) '*WARNING in calinput'
        write(*,*) '         *BEAM GENERAL SECTION is obsolete'
        write(*,*) '         please use *BEAM SECTION'
        write(*,*)
        call beamgeneralsections(inpc,textpart,set,istartset,
     &       iendset,ialset,nset,ielmat,matname,nmat,ielorien,
     &       orname,norien,thicke,ipkon,iponor,xnor,ixfree,
     &       offset,lakon,irstrt,istep,istat,n,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc,mi,ielprop,nprop,nprop_,prop,
     &       nelcon,ier)
!     
      elseif(textpart(1)(1:12).eq.'*BEAMSECTION') then
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
          call beamgeneralsections(inpc,textpart,set,istartset,
     &         iendset,ialset,nset,ielmat,matname,nmat,ielorien,
     &         orname,norien,thicke,ipkon,iponor,xnor,ixfree,
     &         offset,lakon,irstrt,istep,istat,n,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc,mi,ielprop,nprop,nprop_,prop,
     &         nelcon,ier)
        else
          call beamsections(inpc,textpart,set,istartset,iendset,
     &         ialset,nset,ielmat,matname,nmat,ielorien,orname,norien,
     &         thicke,ipkon,iponor,xnor,ixfree,
     &         offset,lakon,irstrt,istep,istat,n,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc,mi,nelcon,ier)
        endif
!     
      elseif(textpart(1)(1:10).eq.'*BOUNDARYF') then
        M_or_SPC=1
        call boundaryfs(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,nodeboun,ndirboun,xboun,nboun,nboun_,nk,
     &       iamboun,amname,nam,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &       mpcfree,trab,ntrans,ikboun,ilboun,ikmpc,ilmpc,
     &       nk_,co,labmpc,typeboun,istat,n,iline,
     &       ipol,inl,ipoinp,inp,nam_,namtot_,namta,amta,nmethod,
     &       iperturb,ipoinpc,vold,mi,xload,sideload,nload,nelemload,
     %       lakon,kon,ipkon,ne,iamplitudedefault,namtot,ier)
!     
      elseif(textpart(1)(1:9).eq.'*BOUNDARY') then
        M_or_SPC=1
        call boundarys(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,nodeboun,ndirboun,xboun,nboun,nboun_,nk,
     &       iamboun,amname,nam,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &       mpcfree,inotr,trab,ntrans,ikboun,ilboun,ikmpc,ilmpc,
     &       nk_,co,labmpc,boun_flag,typeboun,istat,n,iline,
     &       ipol,inl,ipoinp,inp,nam_,namtot_,namta,amta,nmethod,
     &       iperturb,iaxial,ipoinpc,vold,mi,iamplitudedefault,
     &       namtot,ier)
        boun_flag=.true.
!     
      elseif(textpart(1)(1:7).eq.'*BUCKLE') then
        call buckles(inpc,textpart,nmethod,mei,fei,
     &       nforc,nload,ithermal,iprestr,nbody,t0,t1,nk,iperturb,
     &       istep,istat,n,iline,ipol,inl,ipoinp,inp,isolver,ipoinpc,
     &       ier)
!     
      elseif(textpart(1)(1:4).eq.'*CFD')
     &       then
        call cfds(inpc,textpart,nmethod,
     &       iperturb,isolver,
     &       istep,istat,n,tinc,tper,tmin,tmax,idrct,ithermal,iline,
     &       ipol,inl,ipoinp,inp,ipoinpc,alpha,ctrl,iexpl,tincf,
     &       ttime,physcon,ier)
!     
      elseif(textpart(1)(1:6).eq.'*CFLUX') then
        call cfluxs(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,nodeforc,ndirforc,xforc,nforc,nforc_,iamforc,
     &       amname,nam,ntrans,trab,inotr,co,ikforc,ilforc,nk,
     &       cflux_flag,istep,istat,n,iline,ipol,inl,ipoinp,inp,nam_,
     &       namtot_,namta,amta,iaxial,ipoinpc,idefforc,ipompc,nodempc,
     &       nmpc,ikmpc,ilmpc,labmpc,iamplitudedefault,namtot,ier)
        cflux_flag=.true.
!     
      elseif(textpart(1)(1:15).eq.'*CHANGEFRICTION') then
        ichangefriction=1
        call changefrictions(inpc,textpart,matname,nmat,nmat_,
     &       irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &       ipoinpc,imat,ier)
!     
      elseif(textpart(1)(1:15).eq.'*CHANGEMATERIAL') then
        call changematerials(inpc,textpart,matname,nmat,nmat_,
     &       irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &       imat,ier)
!     
      elseif(textpart(1)(1:14).eq.'*CHANGEPLASTIC') then
        call changeplastics(inpc,textpart,imat,ntmat_,npmat_,
     &       plicon,nplicon,plkcon,nplkcon,istep,istat,n,iline,ipol,
     &       inl,ipoinp,inp,ipoinpc,nelcon,ier)
!     
      elseif(textpart(1)(1:19).eq.'*CHANGESOLIDSECTION') then
        call changesolidsections(inpc,textpart,set,istartset,
     &       iendset,ialset,nset,ielmat,matname,nmat,ielorien,orname,
     &       norien,lakon,thicke,kon,ipkon,irstrt,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,cs,mcs,iaxial,ipoinpc,mi,nelcon,ier)
!     
      elseif(textpart(1)(1:22).eq.'*CHANGESURFACEBEHAVIOR') then
        ichangesurfacebehavior=1
        call changesurfacebehaviors(inpc,textpart,matname,nmat,
     &       ntmat_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &       ipoinpc,imat,ier)
!     
      elseif(textpart(1)(1:10).eq.'*CLEARANCE') then
        call clearances(inpc,textpart,tieset,istat,n,iline,
     &       ipol,inl,ipoinp,inp,ntie,ipoinpc,istep,tietol,irstrt,
     &       ier)
!     
      elseif(textpart(1)(1:6).eq.'*CLOAD') then
        call cloads(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,nodeforc,ndirforc,xforc,nforc,nforc_,
     &       iamforc,amname,nam,ntrans,trab,inotr,co,ikforc,ilforc,
     &       nk,cload_flag,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &       nam_,namtot_,namta,amta,nmethod,iaxial,iperturb,ipoinpc,
     &       maxsectors,idefforc,ipompc,nodempc,
     &       nmpc,ikmpc,ilmpc,labmpc,iamplitudedefault,namtot,ier,
     &       edc,orab,coeffc,ikdc,ndc)
        cload_flag=.true.
c        do i=1,nforc
c          write(*,*) i,nodeforc(1,i),ndirforc(i),xforc(i)
c        enddo
c        call exit(201)
!     
      elseif(textpart(1)(1:17).eq.'*COMPLEXFREQUENCY') then
        call complexfrequencys(inpc,textpart,nmethod,
     &       mei,iperturb,istep,istat,n,iline,ipol,inl,
     &       ipoinp,inp,ithermal,xboun,nboun,ipoinpc,mcs,cs,
     &       cyclicsymmetry,ier)
!     
      elseif(textpart(1)(1:13).eq.'*CONDUCTIVITY') then
        call conductivitys(inpc,textpart,cocon,ncocon,
     &       imat,ntmat_,irstrt,istep,istat,n,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:11).eq.'*CONSTRAINT') then
        call constraints(inpc,textpart,istat,n,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc,nobject,objectset,ier,
     &       nmethod)
        constraint_flag=.false.
!     
      elseif(textpart(1)(1:15).eq.'*CONTACTDAMPING') then
        call contactdampings(inpc,textpart,elcon,nelcon,
     &       nmat,ntmat_,ncmat_,irstrt,istep,istat,n,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc,imat,ier)
      elseif((textpart(1)(1:12).eq.'*CONTACTFILE').or.
     &       (textpart(1)(1:14).eq.'*CONTACTOUTPUT')) then
        if(textpart(1)(1:12).eq.'*CONTACTFILE') then
          output='asc '
        else
          output='bin '
        endif
        ifile_output=3
        call noelfiles(inpc,textpart,jout,filab,nmethod,
     &       nodefile_flag,elfile_flag,ifile_output,nener,ithermal,
     &       istep,istat,n,iline,ipol,inl,ipoinp,inp,out3d,nlabel,
     &       amname,nam,itpamp,idrct,ipoinpc,nef,contactfile_flag,
     &       set,nset,xmodal,ier,physcon,output)
        contactfile_flag=.true.
!     
      elseif(textpart(1)(1:12).eq.'*CONTACTPAIR') then
        call contactpairs(inpc,textpart,tieset,istep,
     &       istat,n,iline,ipol,inl,ipoinp,inp,ntie,ntie_,
     &       iperturb,matname,nmat,ipoinpc,tietol,set,nset,
     &       mortar,ncmat_,ntmat_,elcon,ier)
!     
      elseif(textpart(1)(1:13).eq.'*CONTACTPRINT') then
        call contactprints(inpc,textpart,nprint,nprint_,jout,
     &       prlab,prset,contactprint_flag,ithermal,istep,istat,n,
     &       iline,ipol,inl,ipoinp,inp,amname,nam,itpamp,idrct,
     &       ipoinpc,nener,ier,ntie,tieset)
        contactprint_flag=.true.
!     
      elseif(textpart(1)(1:9).eq.'*CONTROLS') then
        call controlss(inpc,textpart,ctrl,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:18).eq.'*CORRELATIONLENGTH') then
         call correlationlengths(inpc,textpart,istep,istat,n,iline,ipol,
     &       inl,ipoinp,inp,ipoinpc,physcon,ier)
!     
      elseif(textpart(1)(1:32).eq.'*COUPLEDTEMPERATURE-DISPLACEMENT')
     &       then
        call couptempdisps(inpc,textpart,nmethod,
     &       iperturb,isolver,
     &       istep,istat,n,tinc,tper,tmin,tmax,idrct,ithermal,iline,
     &       ipol,inl,ipoinp,inp,ipoinpc,alpha,ctrl,iexpl,tincf,
     &       ttime,nener,ier)
!     
      elseif(textpart(1)(1:9).eq.'*COUPLING')
     &       then
        call couplings(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,nboun,nk,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &       mpcfree,ikboun,ikmpc,ilmpc,co,labmpc,istat,n,iline,ipol,
     &       inl,ipoinp,inp,ipoinpc,norien,orname,orab,irstrt,ipkon,
     &       kon,lakon,istep,ics,dcs,nk_,nboun_,nodeboun,ndirboun,
     &       typeboun,ilboun,ier,nfc,nfc_,coeffc,ikdc,ndc,ndc_,
     &       edc)
!     
      elseif(textpart(1)(1:17).eq.'*CRACKPROPAGATION') then
        call crackpropagations(inpc,textpart,nmethod,iperturb,isolver,
     &       istep,
     &       istat,n,tinc,tper,tmin,tmax,idrct,iline,ipol,inl,ipoinp,
     &       inp,ithermal,cs,ics,tieset,istartset,
     &       iendset,ialset,ipompc,nodempc,coefmpc,nmpc,nmpc_,ikmpc,
     &       ilmpc,mpcfree,mcs,set,nset,labmpc,ipoinpc,iexpl,nef,ttime,
     &       iaxial,nelcon,nmat,ier,jobnamec,matname)
!     
      elseif(textpart(1)(1:6).eq.'*CREEP') then
        call creeps(inpc,textpart,nelcon,imat,ntmat_,npmat_,
     &       plicon,nplicon,elcon,iplas,iperturb,nstate_,ncmat_,
     &       matname,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &       ipoinpc,ianisoplas,ier)
!     
      elseif(textpart(1)(1:16).eq.'*CYCLICHARDENING') then
        call cyclichardenings(inpc,textpart,nelcon,imat,ntmat_,
     &       npmat_,plicon,nplicon,ncmat_,elcon,matname,
     &       irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &       ier)
!     
c     ics(15*ncx_+1) was removed: can be cleaned up.
c     
      elseif(textpart(1)(1:20).eq.'*CYCLICSYMMETRYMODEL') then
        call cyclicsymmetrymodels(inpc,textpart,set,
     &       istartset,iendset,
     &       ialset,nset,tieset,tietol,co,nk,ipompc,nodempc,
     &       coefmpc,nmpc,nmpc_,ikmpc,ilmpc,mpcfree,dcs(lprev+1),
     &       dcs(ncs_+lprev+1),ics(lprev+1),ics(ncs_+lprev+1),
     &       ics(2*ncs_+lprev+1),dcs(2*ncs_+lprev+1),
     &       dcs(3*ncs_+lprev+1),ncs_,cs,labmpc,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,ntie,mcs,lprev,ithermal,
     &       dcs(4*ncs_+1),dcs(6*ncs_+1),dcs(8*ncs_+1),dcs(10*ncs_+1),
     &       ics(3*ncs_+1),ics(5*ncs_+1),ics(7*ncs_+1),ics(8*ncs_+1),
     &       dcs(12*ncs_+1),ne,ipkon,kon,lakon,ics(14*ncs_+1),
     &       ics(16*ncs_+1),ics(18*ncs_+1),ipoinpc,
     &       maxsectors,trab,ntrans,ntrans_,jobnamec,vold,nef,mi,
     &       iaxial,ier)
!     
      elseif(textpart(1)(1:8).eq.'*DAMPING') then
        call dampings(inpc,textpart,xmodal,istep,
     &       istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,irstrt,ier,
     &       dacon,nmat)
!     
      elseif(textpart(1)(1:8).eq.'*DASHPOT') then
        call dashpots(inpc,textpart,nelcon,nmat,ntmat_,npmat_,
     &       plicon,nplicon,
     &       ncmat_,elcon,matname,irstrt,istep,istat,n,iline,ipol,
     &       inl,ipoinp,inp,nmat_,set,istartset,iendset,ialset,
     &       nset,ielmat,ielorien,ipoinpc,mi,ier)
!     
      elseif(textpart(1)(1:22).eq.'*DEFORMATIONPLASTICITY') then
        call deformationplasticitys(inpc,textpart,elcon,nelcon,
     &       imat,ntmat_,ncmat_,irstrt,istep,istat,n,iperturb,
     &       iline,ipol,inl,ipoinp,inp,ipoinpc,ier,iplas)
!     
      elseif(textpart(1)(1:8).eq.'*DENSITY') then
        call densitys(inpc,textpart,rhcon,nrhcon,
     &       imat,ntmat_,irstrt,istep,istat,n,iline,ipol,
     &       inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:7).eq.'*DEPVAR') then
        call depvars(inpc,textpart,nelcon,imat,
     &       nstate_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &       ncocon,ipoinpc,ier)
!     
      elseif(textpart(1)(1:16).eq.'*DESIGNVARIABLES') then
        call designvariabless(inpc,textpart,tieset,tietol,istep,
     &       istat,n,iline,ipol,inl,ipoinp,inp,ntie,ntie_,ipoinpc,
     &       set,nset,ier)
!     
!     energy is calculated in case it is an objective function
!     (is not clear at the time of reading *DESIGN VARIABLES)
!     
        nener=1
!     
      elseif(textpart(1)(1:15).eq.'*DESIGNRESPONSE') then
        call designresponses(inpc,textpart,istat,n,iline,ipol,
     &       inl,ipoinp,inp,ipoinpc,nobject,objectset,
     &       set,nset,ntie,tieset,ier,nmethod) 
!     
      elseif(textpart(1)(1:6).eq.'*DFLUX') then
        call dfluxs(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,nelemload,sideload,xload,nload,nload_,
     &       ielmat,ntmat_,iamload,amname,nam,lakon,ne,dflux_flag,
     &       istep,istat,n,iline,ipol,inl,ipoinp,inp,nam_,namtot_,
     &       namta,amta,ipoinpc,mi,idefload,iamplitudedefault,
     &       namtot,ier)
        dflux_flag=.true.
!     
      elseif(textpart(1)(1:21).eq.'*DISTRIBUTINGCOUPLING') then
        call distributingcouplings(inpc,textpart,ipompc,nodempc,
     &       coefmpc,nmpc,nmpc_,mpcfree,nk,ikmpc,ilmpc,
     &       labmpc,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &       lakon,kon,ipkon,set,nset,istartset,iendset,ialset,co,
     &       ier)
!     
!     *DISTRIBUTION TABLE: nothing to do
!     however, check has to precede *DISTRIBUTION
!     
      elseif(textpart(1)(1:18).eq.'*DISTRIBUTIONTABLE') then
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
!     
      elseif(textpart(1)(1:13).eq.'*DISTRIBUTION') then
        call distributions(inpc,textpart,orname,orab,norien,
     &       norien_,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &       ier,set,istartset,iendset,ialset,nset,ne,ielorien,mi)
!     
      elseif((textpart(1)(1:6).eq.'*DLOAD').or.
     &       (textpart(1)(1:7).eq.'*DSLOAD')) then
        call dloads(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,nelemload,sideload,xload,nload,nload_,
     &       ielmat,iamload,
     &       amname,nam,lakon,ne,dload_flag,istep,istat,n,
     &       iline,ipol,inl,ipoinp,inp,cbody,ibody,xbody,nbody,nbody_,
     &       xbodyold,iperturb,physcon,nam_,namtot_,namta,amta,nmethod,
     &       ipoinpc,maxsectors,mi,idefload,idefbody,ipkon,thicke,
     &       iamplitudedefault,namtot,ier)
        dload_flag=.true.
!     
      elseif(textpart(1)(1:8).eq.'*DYNAMIC') then
        call dynamics(inpc,textpart,nmethod,iperturb,tinc,tper,
     &       tmin,tmax,idrct,alpha,iexpl,isolver,istep,
     &       istat,n,iline,ipol,inl,ipoinp,inp,ithermal,ipoinpc,nef,
     &       ctrl,tincf,nener,ier,mortar)
!     
      elseif(textpart(1)(1:8).eq.'*ELASTIC') then
        call elastics(inpc,textpart,elcon,nelcon,
     &       imat,ntmat_,ncmat_,irstrt,istep,istat,n,
     &       iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:23).eq.'*ELECTRICALCONDUCTIVITY') then
        call electricalconductivitys(inpc,textpart,alcon,nalcon,
     &       alzero,imat,ntmat_,irstrt,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:17).eq.'*ELECTROMAGNETICS') then
        call electromagneticss(inpc,textpart,nmethod,iperturb,
     &       isolver,istep,istat,n,tinc,tper,tmin,tmax,idrct,
     &       ithermal,iline,ipol,inl,ipoinp,inp,alpha,mei,fei,
     &       ipoinpc,ctrl,ttime,ier)
!     
      elseif((textpart(1)(1:8).eq.'*ELEMENT').and.
     &       (textpart(1)(1:14).ne.'*ELEMENTOUTPUT')) then
        call elements(inpc,textpart,kon,ipkon,lakon,nkon,
     &       ne,ne_,set,istartset,iendset,ialset,nset,nset_,nalset,
     &       nalset_,mi(1),ixfree,iponor,xnor,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,iaxial,ipoinpc,solid,
     &       network,filab,nlabel,out3d,iuel,nuel_,ier,iparentel)
!     
      elseif((textpart(1)(1:7).eq.'*ELFILE').or.
     &       (textpart(1)(1:14).eq.'*ELEMENTOUTPUT')) then
        if(textpart(1)(1:7).eq.'*ELFILE') then
          output='asc '
        else
          output='bin '
        endif
        ifile_output=2
        call noelfiles(inpc,textpart,jout,filab,nmethod,
     &       nodefile_flag,elfile_flag,ifile_output,nener,ithermal,
     &       istep,istat,n,iline,ipol,inl,ipoinp,inp,out3d,nlabel,
     &       amname,nam,itpamp,idrct,ipoinpc,nef,contactfile_flag,
     &       set,nset,xmodal,ier,physcon,output)
        elfile_flag=.true.
!     
      elseif(textpart(1)(1:8).eq.'*ELPRINT') then
        call elprints(inpc,textpart,set,
     &       nset,nprint,nprint_,jout,
     &       prlab,prset,nmethod,elprint_flag,nener,ithermal,
     &       istep,istat,n,iline,ipol,inl,ipoinp,inp,amname,nam,itpamp,
     &       idrct,ipoinpc,nef,ier)
        elprint_flag=.true.
!     
      elseif(textpart(1)(1:6).eq.'*ELSET') then
        call noelsets(inpc,textpart,set,istartset,iendset,ialset,
     &       nset,nset_,nalset,nalset_,nk,ne,irstrt,istep,istat,n,
     &       iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:8).eq.'*ENDSTEP') then
        exit
!     
      elseif(textpart(1)(1:10).eq.'*EQUATIONF') then
        M_or_SPC=1
        call equationfs(inpc,textpart,ipompc,nodempc,coefmpc,
     &       nmpc,nmpc_,mpcfree,co,trab,ntrans,ikmpc,ilmpc,
     &       labmpc,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &       lakon,ne,nload,sideload,ipkon,kon,nelemload,ier)
!     
      elseif(textpart(1)(1:9).eq.'*EQUATION') then
        M_or_SPC=1
        call equations(inpc,textpart,ipompc,nodempc,coefmpc,
     &       nmpc,nmpc_,mpcfree,nk,co,trab,inotr,ntrans,ikmpc,ilmpc,
     &       labmpc,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &       set,istartset,iendset,ialset,nset,nodempcref,coefmpcref,
     &       ikmpcref,memmpcref_,mpcfreeref,maxlenmpcref,memmpc_,
     &       maxlenmpc,ier)
!     
      elseif(textpart(1)(1:10).eq.'*EXPANSION') then
        call expansions(inpc,textpart,alcon,nalcon,
     &       alzero,imat,ntmat_,irstrt,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:13).eq.'*SECTIONPRINT') then
        call sectionprints(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,nset_,nalset,nprint,nprint_,jout,prlab,
     &       prset,sectionprint_flag,ithermal,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,amname,nam,itpamp,idrct,ipoinpc,nef,
     &       ier)
        sectionprint_flag=.true.
!     
      elseif(textpart(1)(1:18).eq.'*FEASIBLEDIRECTION') then
        call feasibledirections(inpc,textpart,istat,n,key,iline,ipol,
     &       inl,ipoinp,inp,ipoinpc,nmethod,istep,ier,tmax)
!     
      elseif(textpart(1)(1:5).eq.'*FILM') then
        call films(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,nelemload,sideload,xload,nload,nload_,
     &       ielmat,ntmat_,iamload,amname,nam,lakon,ne,film_flag,
     &       istep,istat,n,iline,ipol,inl,ipoinp,inp,nam_,namtot_,
     &       namta,amta,ipoinpc,mi,iamplitudedefault,namtot,ier)
        film_flag=.true.
!     
      elseif(textpart(1)(1:7).eq.'*FILTER') then
        call filters(inpc,textpart,istep,istat,n,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc,objectset,ier,nobject,nmethod)
!     
      elseif(textpart(1)(1:15).eq.'*FLUIDCONSTANTS') then
        call fluidconstantss(inpc,textpart,shcon,nshcon,
     &       imat,ntmat_,irstrt,istep,istat,n,iline,ipol,
     &       inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:13).eq.'*FLUIDSECTION') then
        call fluidsections(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,ielmat,matname,nmat,irstrt,istep,istat,n,
     &       iline,ipol,inl,ipoinp,inp,lakon,ielprop,nprop,
     &       nprop_,prop,ipoinpc,mi,ier)
!     
      elseif(textpart(1)(1:10).eq.'*FREQUENCY') then
        call frequencys(inpc,textpart,nmethod,
     &       mei,fei,iperturb,istep,istat,n,iline,ipol,
     &       inl,ipoinp,inp,ithermal,isolver,xboun,nboun,ipoinpc,
     &       ipompc,labmpc,fmpc,ikmpc,ilmpc,nmpc,ier)
!     
      elseif(textpart(1)(1:9).eq.'*FRICTION') then
        call frictions(inpc,textpart,elcon,nelcon,
     &       imat,ntmat_,ncmat_,irstrt,istep,istat,n,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc,nstate_,ichangefriction,
     &       mortar,ier)
!     
      elseif(textpart(1)(1:5).eq.'*GAP ') then
        call gaps(inpc,textpart,nelcon,nmat,ntmat_,npmat_,
     &       plicon,nplicon,
     &       ncmat_,elcon,matname,irstrt,istep,istat,n,iline,ipol,
     &       inl,ipoinp,inp,nmat_,set,istartset,iendset,ialset,
     &       nset,ielmat,ielorien,ipoinpc,mi,ier)
!     
      elseif(textpart(1)(1:15).eq.'*GAPCONDUCTANCE') then
        call gapconductances(inpc,textpart,nelcon,nmat,ntmat_,
     &       npmat_,plkcon,nplkcon,iperturb,irstrt,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:18).eq.'*GAPHEATGENERATION') then
        call gapheatgenerations(inpc,textpart,elcon,nelcon,
     &       imat,ntmat_,ncmat_,irstrt,istep,istat,n,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc,nstate_,ier)
!     
      elseif(textpart(1)(1:19).eq.'*GEOMETRICTOLERANCE') then
         call geometrictolerances(inpc,textpart,set,istartset,iendset,
     &        ialset,nset,nk,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &        ipoinpc,ier,irobustdesign,irandomtype,randomval)
!     
      elseif(textpart(1)(1:20).eq.'*GEOMETRICCONSTRAINT') then
         call geometricconstraints(inpc,textpart,istat,n,
     &        iline,ipol,inl,ipoinp,inp,ipoinpc,nobject,objectset,
     &        set,nset,ier,nmethod)
!     
      elseif(textpart(1)(1:6).eq.'*GREEN') then
        call greens(inpc,textpart,nmethod,
     &       mei,iperturb,istep,istat,n,iline,ipol,inl,
     &       ipoinp,inp,ithermal,isolver,xboun,nboun,ipoinpc,
     &       ier)
!     
      elseif(textpart(1)(1:4).eq.'*HCF') then
        call hcfs(inpc,textpart,istep,istat,n,iline,ipol,inl,ipoinp,
     &       inp,ipoinpc,ier,jobnamec,mei,tincf)
!     
      elseif(textpart(1)(1:8).eq.'*HEADING') then
        call headings(inpc,textpart,istat,n,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc,heading,istep,irstrt,ier)
!     
      elseif(textpart(1)(1:13).eq.'*HEATTRANSFER') then
        call heattransfers(inpc,textpart,nmethod,iperturb,isolver,
     &       istep,istat,n,tinc,tper,tmin,tmax,idrct,ithermal,iline,
     &       ipol,inl,ipoinp,inp,alpha,mei,fei,ipoinpc,ctrl,ttime,
     &       ier)
!     
      elseif(textpart(1)(1:13).eq.'*HYPERELASTIC') then
        call hyperelastics(inpc,textpart,elcon,nelcon,
     &       imat,ntmat_,ncmat_,irstrt,istep,istat,n,iperturb,
     &       iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:10).eq.'*HYPERFOAM') then
        call hyperfoams(inpc,textpart,elcon,nelcon,
     &       imat,ntmat_,ncmat_,irstrt,istep,istat,n,iperturb,iline,
     &       ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:18).eq.'*INITIALCONDITIONS') then
        call initialconditionss(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,t0,t1,prestr,iprestr,ithermal,veold,inoelfree,
     &       nk_,mi(1),istep,istat,n,iline,ipol,inl,ipoinp,inp,lakon,
     &       kon,co,ne,ipkon,vold,ipoinpc,xstate,nstate_,nk,t0g,
     &       t1g,iaxial,ielprop,prop,ier,nuel_)
!     
      elseif(textpart(1)(1:12).eq.'*INITIALMESH') then
         call initialmeshs(inpc,textpart,coini,nk,nk_,set,istat,n,
     &        iline,ipol,inl,ipoinp,inp,ipoinpc,ier)              
!     
      elseif(textpart(1)(1:22).eq.'*INITIALSTRAININCREASE') then
        call initialstrainincreases(inpc,textpart,prestr,iprestr,
     &       mi,istep,istat,n,iline,ipol,inl,ipoinp,inp,ne,ipoinpc,
     &       ier)
!     
      elseif(textpart(1)(1:21).eq.'*MAGNETICPERMEABILITY') then
        call magneticpermeabilitys(inpc,textpart,elcon,nelcon,
     &       imat,ntmat_,ncmat_,irstrt,istep,istat,n,
     &       iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:9).eq.'*MASSFLOW') then
        call massflows(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,nelemload,sideload,xload,nload,nload_,
     &       iamload,lakon,ne,istep,istat,n,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc,idefload,nam,ier)
!     
      elseif(textpart(1)(1:5).eq.'*MASS') then
        call masss(inpc,textpart,nrhcon,nmat,ntmat_,
     &       rhcon,matname,irstrt,istep,istat,n,iline,ipol,
     &       inl,ipoinp,inp,nmat_,set,istartset,iendset,ialset,
     &       nset,ielmat,ielorien,ipoinpc,mi,iaxial,ier)
!     
      elseif(textpart(1)(1:9).eq.'*MATERIAL') then
        call materials(inpc,textpart,matname,nmat,nmat_,
     &       irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &       ipoinpc,imat,ier)
!     
      elseif(textpart(1)(1:16).eq.'*MEMBRANESECTION') then
        call membranesections(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,ielmat,matname,nmat,ielorien,orname,
     &       norien,thicke,kon,ipkon,offset,irstrt,istep,istat,n,
     &       iline,ipol,inl,ipoinp,inp,lakon,iaxial,ipoinpc,mi,
     &       icomposite,nelcon,ier)
!     
      elseif(textpart(1)(1:13).eq.'*MODALDAMPING') then
        call modaldampings(inpc,textpart,nmethod,xmodal,istep,
     &       istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:13).eq.'*MODALDYNAMIC') then
        call modaldynamics(inpc,textpart,nmethod,tinc,tper,iexpl,
     &       istep,istat,n,iline,ipol,inl,ipoinp,inp,iperturb,
     &       isolver,cs,mcs,ipoinpc,idrct,ctrl,tmin,tmax,
     &       nforc,nload,nbody,iprestr,t0,t1,ithermal,nk,vold,veold,
     &       xmodal,set,nset,mi,cyclicsymmetry,ier)
!     
      elseif(textpart(1)(1:12).eq.'*MODELCHANGE') then
        call modelchanges(inpc,textpart,tieset,istat,n,iline,
     &       ipol,inl,ipoinp,inp,ntie,ipoinpc,istep,ipkon,nset,
     &       istartset,iendset,set,ialset,ne,mi,ielmat,iprestr,
     &       iperturb,ier,tietol)
!     
      elseif(textpart(1)(1:4).eq.'*MPC') then
        call mpcs(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,nset_,nalset,nalset_,ipompc,nodempc,
     &       coefmpc,labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,lakon,
     &       ipkon,kon,nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &       nboun,nboun_,iperturb,ne_,co,xboun,ctrl,typeboun,
     &       istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &       ier)
!     
      elseif(textpart(1)(1:11).eq.'*NETWORKMPC') then
        M_or_SPC=1
        call networkmpcs(inpc,textpart,ipompc,nodempc,coefmpc,
     &       nmpc,nmpc_,mpcfree,nk,ikmpc,ilmpc,
     &       labmpc,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &       ier)
!     
      elseif(textpart(1)(1:11).eq.'*NOANALYSIS') then
        call noanalysiss(inpc,textpart,nmethod,iperturb,istep,
     &       istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,tper,ier)
!     
      elseif(textpart(1)(1:15).eq.'*NODALTHICKNESS') then
        call nodalthicknesss(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,thickn,nk,istep,istat,n,iline,ipol,inl,
     &       ipoinp,inp,iaxial,ipoinpc,ier)
!     
      elseif((textpart(1)(1:5).eq.'*NODE').and.
     &       (textpart(1)(1:10).ne.'*NODEPRINT').and.
     &       (textpart(1)(1:11).ne.'*NODEOUTPUT').and.
     &       (textpart(1)(1:9).ne.'*NODEFILE')) then
        call nodes(inpc,textpart,co,nk,nk_,set,istartset,iendset,
     &       ialset,nset,nset_,nalset,nalset_,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif((textpart(1)(1:9).eq.'*NODEFILE').or.
     &       (textpart(1)(1:11).eq.'*NODEOUTPUT')) then
        if(textpart(1)(1:9).eq.'*NODEFILE') then
          output='asc '
        else
          output='bin '
        endif
        ifile_output=1
        call noelfiles(inpc,textpart,jout,filab,nmethod,
     &       nodefile_flag,elfile_flag,ifile_output,nener,ithermal,
     &       istep,istat,n,iline,ipol,inl,ipoinp,inp,out3d,nlabel,
     &       amname,nam,itpamp,idrct,ipoinpc,nef,contactfile_flag,
     &       set,nset,xmodal,ier,physcon,output)
        nodefile_flag=.true.
!     
      elseif(textpart(1)(1:10).eq.'*NODEPRINT') then
        call nodeprints(inpc,textpart,set,istartset,iendset,ialset,
     &       nset,nset_,nalset,nprint,nprint_,jout,
     &       prlab,prset,nodeprint_flag,ithermal,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,amname,nam,itpamp,idrct,ipoinpc,nef,
     &       ier)
        nodeprint_flag=.true.
!     
      elseif(textpart(1)(1:7).eq.'*NORMAL') then
        call normals(inpc,textpart,iponor,xnor,ixfree,
     &       ipkon,kon,nk,nk_,ne,lakon,istep,istat,n,iline,ipol,
     &       inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:5).eq.'*NSET') then
        call noelsets(inpc,textpart,set,istartset,iendset,ialset,
     &       nset,nset_,nalset,nalset_,nk,ne,irstrt,istep,istat,n,
     &       iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:10).eq.'*OBJECTIVE') then
        call objectives(inpc,textpart,istat,n,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc,nobject,objectset,ier,nmethod,
     &       objective_flag)
        objective_flag=.true.
!     
      elseif(textpart(1)(1:12).eq.'*ORIENTATION') then
        call orientations(inpc,textpart,orname,orab,norien,
     &       norien_,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &       ier)
!     
      elseif(textpart(1)(1:7).eq.'*OUTPUT') then
        call outputs(inpc,textpart,jout,itpamp,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:18).eq.'*PHYSICALCONSTANTS') then
        call physicalconstantss(inpc,textpart,physcon,
     &       istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:8).eq.'*PLASTIC') then
        call plastics(inpc,textpart,nelcon,imat,ntmat_,npmat_,
     &       plicon,nplicon,plkcon,nplkcon,iplas,iperturb,nstate_,
     &       ncmat_,elcon,matname,irstrt,istep,istat,n,iline,ipol,
     &       inl,ipoinp,inp,ipoinpc,ianisoplas,ier)
!     
      elseif(textpart(1)(1:19).eq.'*PRE-TENSIONSECTION') then
        call pretensionsections(inpc,textpart,ipompc,nodempc,
     &       coefmpc,nmpc,nmpc_,mpcfree,nk,ikmpc,ilmpc,
     &       labmpc,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &       lakon,kon,ipkon,set,nset,istartset,iendset,ialset,co,
     &       ics,dcs,t0,ithermal,ne,ier)
        pretension=.true.
!     
      elseif(textpart(1)(1:8).eq.'*RADIATE') then
        call radiates(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,nelemload,sideload,xload,nload,nload_,
     &       ielmat,ntmat_,iamload,amname,nam,lakon,ne,radiate_flag,
     &       istep,istat,n,iline,ipol,inl,ipoinp,inp,physcon,nam_,
     &       namtot_,namta,amta,ipoinpc,mi,iamplitudedefault,
     &       namtot,ier)
        radiate_flag=.true.
!     
      elseif(textpart(1)(1:12).eq.'*RANDOMFIELD') then
        call randomfields(inpc,textpart,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,ipoinpc,nener,physcon,ier)        
!     
      elseif(textpart(1)(1:14).eq.'*RATEDEPENDENT') then
        call ratedependents(inpc,textpart,nelcon,nmat,ntmat_,
     &     iplas,iperturb,nstate_,ncmat_,elcon,matname,irstrt,istep,
     &     istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:11).eq.'*REFINEMESH') then
        call refinemeshs(inpc,textpart,filab,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,ipoinpc,ier,set,nset)
!     
      elseif(textpart(1)(1:8).eq.'*RESTART') then
        call restarts(istep,nset,nload,nforc, nboun,nk,ne,
     &       nmpc,nalset,nmat,ntmat_,npmat_,norien,nam,nprint,
     &       mi(1),ntrans,ncs_,namtot,ncmat_,mpcfree,
     &       maxlenmpc,ne1d,
     &       ne2d,nflow,nlabel,iplas,nkon,ithermal,nmethod,
     &       iperturb,nstate_,nener,set,istartset,iendset,ialset,co,
     &       kon,ipkon,lakon,nodeboun,ndirboun,iamboun,xboun,ikboun,
     &       ilboun,ipompc,nodempc,coefmpc,labmpc,ikmpc,ilmpc,
     &       nodeforc,ndirforc,iamforc,xforc,ikforc,ilforc,
     &       nelemload,iamload,sideload,xload,
     &       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,plicon,
     &       nplicon,plkcon,nplkcon,orname,orab,ielorien,trab,inotr,
     &       amname,amta,namta,t0,t1,iamt1,veold,ielmat,
     &       matname,prlab,prset,filab,vold,nodebounold,
     &       ndirbounold,xbounold,xforcold,xloadold,t1old,eme,
     &       iponor,xnor,knor,thicke,offset,iponoel,
     &       inoel,rig,shcon,nshcon,cocon,
     &       ncocon,ics,sti,ener,xstate,jobnamec,infree,
     &       irstrt,inpc,textpart,istat,n,key,prestr,iprestr,
     &       cbody,ibody,xbody,nbody,xbodyold,ttime,qaold,
     &       cs,mcs,output,physcon,ctrl,typeboun,iline,ipol,inl,
     &       ipoinp,inp,fmpc,tieset,ntie,tietol,ipoinpc,nslavs,
     &       t0g,t1g,nprop,ielprop,prop,mortar,nintpoint,ifacecount,
     &       islavsurf,pslavsurf,clearini,ier,vel,nef,velo,veloo,
     &       ne2boun,heading,network,irestartread,nfc,ndc,coeffc,
     &       ikdc,edc)
!     
      elseif(textpart(1)(1:18).eq.'*RETAINEDNODALDOFS') then
        call retainednodaldofss(inpc,textpart,set,istartset,
     &       iendset,ialset,nset,nodeboun,ndirboun,xboun,nboun,
     &       nboun_,nk,iamboun,nam,ipompc,nodempc,coefmpc,nmpc,
     &       nmpc_,mpcfree,inotr,trab,ikboun,ilboun,ikmpc,ilmpc,
     &       nk_,co,labmpc,typeboun,istat,n,iline,ipol,
     &       inl,ipoinp,inp,nmethod,iperturb,
     &       ipoinpc,vold,mi,istep,ier)
!     
      elseif(textpart(1)(1:10).eq.'*RIGIDBODY') then
        call rigidbodys(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,nset_,nalset,nalset_,ipompc,nodempc,
     &       coefmpc,labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,lakon,
     &       ipkon,kon,nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &       nboun,nboun_,iperturb,ne_,ctrl,typeboun,
     &       istep,istat,n,iline,ipol,inl,ipoinp,inp,co,ipoinpc,
     &       ier)
!     
      elseif(textpart(1)(1:13).eq.'*ROBUSTDESIGN') then
        call robustdesigns(inpc,textpart,nmethod,
     &       istep,istat,n,iline,ipol,inl,
     &       ipoinp,inp,tieset,ipoinpc,ntie,tinc,tper,tmin,tmax,tincf,
     &       isens,ier,physcon,irobustdesign)
!     
      elseif(textpart(1)(1:26).eq.'*SELECTCYCLICSYMMETRYMODES') then
        call selectcyclicsymmetrymodess(inpc,textpart,cs,ics,
     &       tieset,istartset,
     &       iendset,ialset,ipompc,nodempc,coefmpc,nmpc,nmpc_,ikmpc,
     &       ilmpc,mpcfree,mcs,set,nset,labmpc,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,nmethod,key,ipoinpc,ier)
!     
      elseif(textpart(1)(1:12).eq.'*SENSITIVITY') then
        call sensitivitys(inpc,textpart,nmethod,
     &       istep,istat,n,iline,ipol,inl,
     &       ipoinp,inp,tieset,ipoinpc,ntie,tinc,tper,tmin,tmax,tincf,
     &       isens,objectset,ier)
!     
      elseif(textpart(1)(1:13).eq.'*SHELLSECTION') then
        call shellsections(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,ielmat,matname,nmat,ielorien,orname,
     &       norien,thicke,kon,ipkon,offset,irstrt,istep,istat,n,
     &       iline,ipol,inl,ipoinp,inp,lakon,iaxial,ipoinpc,mi,
     &       icomposite,nelcon,ier)
!     
      elseif(textpart(1)(1:13).eq.'*SOLIDSECTION') then
        call solidsections(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,ielmat,matname,nmat,ielorien,orname,
     &       norien,lakon,thicke,kon,ipkon,irstrt,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,cs,mcs,iaxial,ipoinpc,mi,co,
     &       ixfree,xnor,iponor,ier,orab)
!     
      elseif(textpart(1)(1:20).eq.'*SPECIFICGASCONSTANT') then
        call specificgasconstants(inpc,textpart,shcon,nshcon,
     &       imat,ntmat_,irstrt,istep,istat,n,iline,ipol,
     &       inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:13).eq.'*SPECIFICHEAT') then
        call specificheats(inpc,textpart,shcon,nshcon,
     &       imat,ntmat_,irstrt,istep,istat,n,iline,ipol,
     &       inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:7).eq.'*SPRING') then
        call springs(inpc,textpart,nelcon,nmat,ntmat_,npmat_,
     &       plicon,nplicon,
     &       ncmat_,elcon,matname,irstrt,istep,istat,n,iline,ipol,
     &       inl,ipoinp,inp,nmat_,set,istartset,iendset,ialset,
     &       nset,ielmat,ielorien,ipoinpc,mi,norien,orname,ier)
!     
      elseif(textpart(1)(1:7).eq.'*STATIC') then
        call statics(inpc,textpart,nmethod,iperturb,isolver,istep,
     &       istat,n,tinc,tper,tmin,tmax,idrct,iline,ipol,inl,ipoinp,
     &       inp,ithermal,cs,ics,tieset,istartset,
     &       iendset,ialset,ipompc,nodempc,coefmpc,nmpc,nmpc_,ikmpc,
     &       ilmpc,mpcfree,mcs,set,nset,labmpc,ipoinpc,iexpl,nef,ttime,
     &       iaxial,nelcon,nmat,tincf,ier)
!     
      elseif(textpart(1)(1:20).eq.'*STEADYSTATEDYNAMICS') then
        call steadystatedynamicss(inpc,textpart,nmethod,
     &       iexpl,istep,istat,n,iline,ipol,inl,ipoinp,inp,iperturb,
     &       isolver,xmodal,cs,mcs,ipoinpc,nforc,nload,nbody,iprestr,
     &       t0,t1,ithermal,nk,set,nset,cyclicsymmetry,ibody,ier)
!     
      elseif(textpart(1)(1:5).eq.'*STEP') then
        call steps(inpc,textpart,iperturb,iprestr,nbody,nforc,
     &       nload,ithermal,t0,t1,nk,irstrt,istep,istat,n,
     &       jmax,ctrl,iline,ipol,inl,ipoinp,inp,newstep,
     &       ipoinpc,network,iamplitudedefault,amname,nam,
     &       nam_,namta,amta,namtot,nstam,ier,namtot_,
     &       physcon)
        boun_flag=.false.
!     
      elseif(textpart(1)(1:9).eq.'*SUBMODEL') then
        call submodels(inpc,textpart,set,istartset,iendset,ialset,
     &       nset,nset_,nalset,nalset_,nk,istep,istat,n,iline,ipol,
     &       inl,ipoinp,inp,ipoinpc,nsubmodel,tieset,tietol,ntie,
     &       ntie_,jobnamec,amta,namta,nam,nam_,namtot_,ier)
!     
      elseif(textpart(1)(1:21).eq.'*SUBSTRUCTUREGENERATE') then
        call substructuregenerates(inpc,textpart,nmethod,iperturb,
     &       isolver,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &       ithermal,ipoinpc,filab,ier)
!     
      elseif(textpart(1)(1:25).eq.'*SUBSTRUCTUREMATRIXOUTPUT') then
        call substructurematrixoutputs(textpart,istep,inpc,
     &       istat,n,key,iline,ipol,inl,ipoinp,inp,jobnamec,ipoinpc,
     &       ier)
!     
      elseif(textpart(1)(1:9).eq.'*SURFACE ') then
        call surfaces(inpc,textpart,set,istartset,iendset,ialset,
     &       nset,nset_,nalset,nalset_,nk,ne,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,lakon,ipoinpc,ier)
!     
      elseif(textpart(1)(1:16).eq.'*SURFACEBEHAVIOR') then
        call surfacebehaviors(inpc,textpart,elcon,nelcon,
     &       imat,ntmat_,ncmat_,irstrt,istep,istat,n,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc,npmat_,plicon,nplicon,nstate_,
     &       ichangesurfacebehavior,ier)
!     
      elseif(textpart(1)(1:19).eq.'*SURFACEINTERACTION') then
        call surfaceinteractions(inpc,textpart,matname,nmat,nmat_,
     &       irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,nrhcon,
     &       ipoinpc,imat,ier)
!     
      elseif(textpart(1)(1:12).eq.'*TEMPERATURE') then
        call temperatures(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,t0,t1,nk,ithermal,iamt1,amname,nam,
     &       inoelfree,nk_,nmethod,temp_flag,istep,istat,n,iline,
     &       ipol,inl,ipoinp,inp,nam_,namtot_,namta,amta,ipoinpc,t1g,
     &       iamplitudedefault,namtot,ier,itempuser,jobnamec,nuel_)
        temp_flag=.true.
!     
      elseif(textpart(1)(1:4).eq.'*TIE') then
        call ties(inpc,textpart,tieset,tietol,istep,
     &       istat,n,iline,ipol,inl,ipoinp,inp,ntie,ntie_,ipoinpc,
     &       ier)
!     
      elseif(textpart(1)(1:11).eq.'*TIMEPOINTS') then
        call timepointss(inpc,textpart,amname,amta,namta,nam,
     &       nam_,namtot_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,
     &       inp,ipoinpc,namtot,ier)
!     
      elseif(textpart(1)(1:11).eq.'*TRANSFORM ') then
        if(M_or_SPC.eq.1) then
          write(*,*) '*WARNING in calinput: SPCs or MPCs have'
          write(*,*) '         been defined before the definition'
          write(*,*) '         of a transformation'
          write(*,*)
        endif
        call transforms(inpc,textpart,trab,ntrans,ntrans_,
     &       inotr,set,istartset,iendset,ialset,nset,istep,istat,
     &       n,iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:11).eq.'*TRANSFORMF') then
        if(M_or_SPC.eq.1) then
          write(*,*) '*WARNING in calinput: SPCs or MPCs have'
          write(*,*) '         been defined before the definition'
          write(*,*) '         of a transformation'
          write(*,*)
        endif
        call transformfs(inpc,textpart,trab,ntrans,ntrans_,
     &       set,istartset,iendset,ialset,nset,istep,istat,
     &       n,iline,ipol,inl,ipoinp,inp,ipoinpc,xload,sideload,
     &       nelemload,idefload,nload,nload_,ne,nam,iamload,ier)
!     
      elseif(textpart(1)(1:34).eq.
     &       '*UNCOUPLEDTEMPERATURE-DISPLACEMENT') then
        call uncouptempdisps(inpc,textpart,
     &       nmethod,iperturb,isolver,
     &       istep,istat,n,tinc,tper,tmin,tmax,idrct,ithermal,iline,
     &       ipol,inl,ipoinp,inp,ipoinpc,alpha,ctrl,ttime,nener,ier)
!     
      elseif(textpart(1)(1:12).eq.'*USERELEMENT') then
        if(istep.gt.0) then
          write(*,*) '*ERROR reading *USER ELEMENT:'
          write(*,*) '       *USER ELEMENT should be placed'
          write(*,*) '  before all step definitions'
          call exit(201)
        endif
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
!     
      elseif(textpart(1)(1:13).eq.'*USERMATERIAL') then
        call usermaterials(inpc,textpart,elcon,nelcon,
     &       imat,ntmat_,ncmat_,iperturb,irstrt,istep,istat,n,
     &       iline,ipol,inl,ipoinp,inp,cocon,ncocon,ipoinpc,ier)
!     
      elseif(textpart(1)(1:12).eq.'*USERSECTION') then
        call usersections(inpc,textpart,set,istartset,iendset,
     &       ialset,nset,ielmat,matname,nmat,irstrt,istep,istat,n,
     &       iline,ipol,inl,ipoinp,inp,lakon,ielprop,nprop,
     &       nprop_,prop,ipoinpc,mi,ier)
!     
      elseif(textpart(1)(1:17).eq.'*VALUESATINFINITY') then
        call valuesatinfinitys(inpc,textpart,physcon,
     &       istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
      elseif(textpart(1)(1:11).eq.'*VIEWFACTOR') then
        call viewfactors(textpart,iviewfile,istep,inpc,
     &       istat,n,key,iline,ipol,inl,ipoinp,inp,jobnamec,ipoinpc,
     &       ier,irestartread)
!     
      elseif(textpart(1)(1:7).eq.'*VISCO') then
        call viscos(inpc,textpart,nmethod,iperturb,isolver,istep,
     &       istat,n,tinc,tper,tmin,tmax,idrct,iline,ipol,inl,ipoinp,
     &       inp,ipoinpc,nelcon,nmat,ncmat_,ctrl,ier)
!     
!     check for zero-character-lines?
!     
      elseif(ipoinpc(iline-1).eq.ipoinpc(iline)) then
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
!     
      else
        write(*,*) '*WARNING in calinput. Card image cannot be inter
     &preted:'
        call inputwarning(inpc,ipoinpc,iline,
     &       "the input file%")
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      endif
!     
      if(ier.eq.1) then
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if(key.eq.1) exit
        enddo
        ier=2
      endif
!     
      enddo loop
!     
!     check whether the *END STEP card was preceded by a *STEP card    
!     
      if(newstep.eq.0) then
        write(*,*) '*ERROR in calinput: *END STEP card in step ',
     &       istep+1
        write(*,*) '       was not preceded by a *STEP card'
      endif
!     
!     reorganizing the input in field inpc
!     
      j=1
      do
        if(j.eq.1) then
          inp(1,j)=iline+1
        else
          inp(1,j)=inp(1,inl)
        endif
        inp(2,j)=inp(2,inl)
        if(inp(3,inl).eq.0) then
          inp(3,j)=0
          ipoinp(2,nentries)=j
          exit
        else
          inl=inp(3,inl)
          inp(3,j)=j+1
          j=j+1
        endif
      enddo
      do j=1,nentries-1
        ipoinp(1,j)=0
      enddo
      ipoinp(1,nentries)=1
!     
!     expanding the 1-D and 2-D elements to volume elements
!     treating the incompressibility constraint
!     
      do j=1,nforc_
        idefforc(j)=0
      enddo
!     
!     copying field t1 and iamt1 in newly created nodes for
!     pretension
!     
      if(pretension.and.(ithermal(1).eq.1)) then
        do i=1,nmpc
          if(labmpc(i)(1:11).eq.'THERMALPRET') then
            nodedep=nodempc(1,ipompc(i))
            nodeind=nodempc(1,nodempc(3,ipompc(i)))
            t1(nodedep)=t1(nodeind)
            if(nam.gt.0) then
              iamt1(nodedep)=iamt1(nodeind)
            endif
          endif
        enddo
      endif
!     
      call gen3delem(kon,ipkon,lakon,ne,ipompc,nodempc,coefmpc,
     &     nmpc,nmpc_,mpcfree,ikmpc,ilmpc,labmpc,ikboun,ilboun,nboun,
     &     nboun_,nodeboun,ndirboun,xboun,iamboun,nam,
     &     inotr,trab,nk,nk_,iponoel,inoel,iponor,xnor,thicke,thickn,
     &     knor,istep,offset,t0,t1,ikforc,ilforc,rig,nforc,
     &     nforc_,nodeforc,ndirforc,xforc,iamforc,sideload,
     &     nload,ithermal,ntrans,co,ixfree,ikfree,inoelfree,iponoelmax,
     &     iperturb,tinc,tper,tmin,tmax,ctrl,typeboun,nmethod,nset,set,
     &     istartset,iendset,ialset,prop,ielprop,vold,mi,nkon,ielmat,
     &     icomposite,t0g,t1g,idefforc,iamt1,orname,orab,norien,norien_,
     &     ielorien,jobnamec,ne2boun)
!     
!     New multistage Routine Call
!     
      call multistages(nkon,set,istartset,iendset,
     &     ialset,nset,tieset,tietol,co,nk,ipompc,nodempc,
     &     coefmpc,nmpc,nmpc_,ikmpc,ilmpc,mpcfree,dcs(lprev+1),
     &     dcs(ncs_+lprev+1),ics(lprev+1),ics(ncs_+lprev+1),
     &     ics(2*ncs_+lprev+1),dcs(2*ncs_+lprev+1),
     &     dcs(3*ncs_+lprev+1),ncs_,cs,labmpc,ntie,mcs,
     &     dcs(4*ncs_+1),dcs(6*ncs_+1),dcs(8*ncs_+1),dcs(10*ncs_+1),
     &     ics(3*ncs_+1),ics(5*ncs_+1),ics(7*ncs_+1),ics(8*ncs_+1),
     &     dcs(12*ncs_+1),ne,ipkon,kon,lakon,ics(14*ncs_+1),
     &     ics(16*ncs_+1),ics(18*ncs_+1),jobnamec)
!     
      infree(1)=ixfree
      infree(2)=ikfree
      infree(3)=inoelfree
      infree(4)=iponoelmax
!     
!     check of the selected options
!     
      if((filab(11)(1:2).eq.'PU').or.(filab(18)(1:3).eq.'PHS').or.
     &     (filab(19)(1:4).eq.'MAXU').or.(filab(20)(1:4).eq.'MAXS'))
     &     then
        if((nmethod.eq.1).or.(nmethod.eq.3).or.(nmethod.eq.4).or.
     &       ((nmethod.eq.2).and.((mcs.eq.0).or.(cs(2,1).lt.0)))) then
          write(*,*) '*WARNING in calinput: PU, PHS, MAXU or MAXS'
          write(*,*) '         was selected for a static, a non-'
          write(*,*) '         cyclic-symmetric frequency, a'
          write(*,*) '         buckling or a modal dynamic'
          write(*,*) '         calculation; the output option is'
          write(*,*) '         removed.'
          write(*,*)
          filab(11)(1:4)='    '
          filab(18)(1:4)='    '
          filab(19)(1:4)='    '
          filab(20)(1:4)='    '
        endif
      endif
!     
      if(((iplas.eq.0).and.(ianisoplas.eq.0)).or.(nmethod.eq.2)) then
        if(filab(6)(1:4).eq.'PEEQ') then
          write(*,*) '*WARNING in calinput: PEEQ-output requested'
          write(*,*) '         yet no (visco)plastic calculation'
          write(*,*)
          filab(6)='      '
        endif
        ii=0
        do i=1,nprint
          if(prlab(i)(1:4).eq.'PEEQ') then
            write(*,*) '*WARNING in calinput: PEEQ-output requested'
            write(*,*) '         yet no (visco)plastic calculation'
            write(*,*)
            cycle
          endif
          ii=ii+1
          prlab(ii)=prlab(i)
          prset(ii)=prset(i)
        enddo
        nprint=ii
      endif
!     
!     for frequency calculations no kinetic energy is calculated
!     
      if(nmethod.eq.2) then
        ii=0
        do i=1,nprint
          if(prlab(i)(1:4).eq.'ELKE') then
            write(*,*) '*WARNING in calinput: ELKE-output requested'
            write(*,*) '         this is not available in a'
            write(*,*) '         frequency calculation: ELKE'
            write(*,*) '         is deactivated in this and all'
            write(*,*) '         subsequent steps (unless reactived'
            write(*,*) '         explicitly)'
            write(*,*)
            cycle
          endif
          ii=ii+1
          prlab(ii)=prlab(i)
          prset(ii)=prset(i)
        enddo
        nprint=ii
      endif
!     
      if((ithermal(2).eq.0).and.(nmethod.le.7)) then
        if(filab(2)(1:2).eq.'NT') then
          write(*,*) '*WARNING in calinput: temperature output'
          write(*,*) '         requested, yet no thermal loading'
          write(*,*) '         active'
          write(*,*)
          filab(2)='      '
        endif
        ii=0
        do i=1,nprint
          if(prlab(i)(1:4).eq.'NT  ') then
            write(*,*) '*WARNING in calinput: temperature output'
            write(*,*) '         requested, yet no thermal loading'
            write(*,*) '         active'
            write(*,*)
            cycle
          endif
          ii=ii+1
          prlab(ii)=prlab(i)
          prset(ii)=prset(i)
        enddo
        nprint=ii
      endif
!     
      if((ithermal(2).le.1).and.(nmethod.le.7)) then
        if(filab(9)(1:3).eq.'HFL') then
          write(*,*) '*WARNING in calinput: heat flux output'
          write(*,*) '         requested, yet no heat transfer'
          write(*,*) '         calculation'
          write(*,*)
          filab(9)='      '
        endif
        if(filab(10)(1:3).eq.'RFL') then
          write(*,*) '*WARNING in calinput: heat source output'
          write(*,*) '         requested, yet no heat transfer'
          write(*,*) '         calculation'
          write(*,*)
          filab(10)='      '
        endif
        if(filab(14)(1:2).eq.'MF') then
          write(*,*) '*WARNING in calinput: mass flow output'
          write(*,*) '         requested, yet no heat transfer'
          write(*,*) '         calculation'
          write(*,*)
          filab(14)='      '
        endif
        if(filab(15)(1:2).eq.'PT') then
          write(*,*) '*WARNING in calinput: total pressure output'
          write(*,*) '         requested, yet no heat transfer'
          write(*,*) '         calculation'
          write(*,*)
          filab(15)='      '
        endif
        if(filab(16)(1:2).eq.'TT') then
          write(*,*) '*WARNING in calinput: total temperature output'
          write(*,*) '         requested, yet no heat transfer'
          write(*,*) '         calculation'
          write(*,*)
          filab(16)='      '
        endif
        ii=0
        do i=1,nprint
          if(prlab(i)(1:4).eq.'HFL ') then
            write(*,*) '*WARNING in calinput: heat flux output'
            write(*,*) '         requested, yet no heat transfer'
            write(*,*) '         calculation'
            write(*,*)
            cycle
          elseif(prlab(i)(1:4).eq.'RFL ') then
            write(*,*) '*WARNING in calinput: heat source output'
            write(*,*) '         requested, yet no heat transfer'
            write(*,*) '         calculation'
            write(*,*)
            cycle
          elseif(prlab(i)(1:4).eq.'MF  ') then
            write(*,*) '*WARNING in calinput: mass flow output'
            write(*,*) '         requested, yet no heat transfer'
            write(*,*) '         calculation'
            write(*,*)
            cycle
          elseif(prlab(i)(1:4).eq.'PN  ') then
            write(*,*) '*WARNING in calinput: pressure output'
            write(*,*) '         requested, yet no heat transfer'
            write(*,*) '         calculation'
            write(*,*)
            cycle
          elseif(prlab(i)(1:4).eq.'TS  ') then
            write(*,*) 
     &           '*WARNING in calinput: total temperature output'
            write(*,*) '         requested, yet no heat transfer'
            write(*,*) '         calculation'
            write(*,*)
            cycle
          endif
          ii=ii+1
          prlab(ii)=prlab(i)
          prset(ii)=prset(i)
        enddo
        nprint=ii
      endif
!     
!     check whether a material was assigned to each active element
!     
      ierror=0
      do i=1,ne
        if(ipkon(i).lt.0) cycle
!     
!     gaps and DCOUP3D-elements do not need a material assignment
!     
        if(lakon(i)(1:1).eq.'G') cycle
        if(lakon(i)(1:7).eq.'DCOUP3D') cycle
        if(ielmat(1,i).eq.0) then
          ierror=1
          write(*,*) '*ERROR in calinput: no material was assigned'
          write(*,*) '       to element ',i
        endif
      enddo
      if(ierror.eq.1) call exit(201)
!     
!     check whether for mechanical calculations with temperature
!     conditions an initial and final temperature value was assigned
!     to each node belonging to an element
!     
      if((ithermal(1).eq.1).and.(istep.eq.1).and.(itempuser(1).eq.0))
     &     then
        call checktemp(t0,t1,lakon,ne,ipkon,kon)
      endif
!     
!     check whether the density was defined for dynamic calculations
!     and transient thermal calculations
!     
      if((((nbody.gt.0).or.
     &     (nmethod.eq.2).or.(nmethod.eq.4)).and.(nef.eq.0)).or.
     &     ((nef.gt.0).and.(iexpl.ne.1))) then
        ierror=0
        do i=1,nmat
          if((nrhcon(i).ne.0).or.(matname(i)(1:6).eq.'SPRING').or.
     &         (matname(i)(1:7).eq.'DASHPOT')) then
            if(nrhcon(i).gt.0) ierror=ierror+1
          else
            write(*,*)'*WARNING in calinput: no density was assigned'
            write(*,*) '         to material ',
     &           matname(i)(1:index(matname(i),' ')-1),
     &           ' in a dynamic'
            write(*,*) '         calculation or a calculation with'
            write(*,*) '         centrifugal or gravitational loads'
            write(*,*) '         or an incompressible CFD calculation'
            write(*,*)
          endif
        enddo
        if(ierror.eq.0) then
          write(*,*) '*ERROR in calinput: no density was assigned'
          write(*,*) '       to any material ',
     &         ' in a dynamic'
          write(*,*) '       calculation or a calculation with'
          write(*,*) '       centrifugal or gravitational loads'
          call exit(201)
        endif
      endif
!     
!     check whether the specific heat was defined for 
!     transient thermal calculations
!     
      if((nmethod.eq.2).or.(nmethod.eq.4)) then
        if(ithermal(1).ge.2) then
          ierror=0
          do i=1,nmat
            if(nshcon(i).ne.0) then
              ierror=ierror+1
            else
              write(*,*) '*WARNING in calinput: no specific heat '
              write(*,*) '         was assigned to material ',
     &             matname(i)(1:index(matname(i),' ')-1),
     &             ' in a transient'
              write(*,*) '         heat transfer calculation'
              write(*,*)
            endif
          enddo
          if(ierror.eq.0) then
            write(*,*) '*ERROR in calinput: no specific heat was'
            write(*,*) '       assigned to any material ',
     &           ' in a transient'
            write(*,*) '       heat transfer calculation'
            write(*,*)
            call exit(201)
          endif
        endif
      endif
!     
!     check whether a *FLUID CONSTANTS card was used for 
!     3D compressible fluid calculations
!     
      if(((nef.gt.0).or.(network.gt.0)).and.(ithermal(1).ge.2)) then
        ierror=0
        do i=1,nmat
          if(nshcon(i).ne.0) then
            ierror=ierror+1
          else
            write(*,*) '*WARNING in calinput: no specific heat '
            write(*,*) '         was assigned to material ',
     &           matname(i)(1:index(matname(i),' ')-1),
     &           ' in a transient'
            write(*,*) '         CFD or network calculation'
            write(*,*)
          endif
        enddo
        if(ierror.eq.0) then
          write(*,*) '*ERROR in calinput: no specific heat was'
          write(*,*) '       assigned to any material ',
     &         ' in a transient'
          write(*,*) '       CFD or network calculation'
          write(*,*)
          call exit(201)
        endif
      endif
!     
!     check whether the elastic constants were defined for 
!     mechanical calculations
!     
      if((ithermal(1).ne.2).and.solid.and.(nmethod.le.7)) then
        ierror=0
        do i=1,nmat
          if(nelcon(1,i).ne.0) then
            ierror=ierror+1
          else
            write(*,*)'*WARNING in calinput: no elastic constants '
            write(*,*)' were assigned to material ',
     &           matname(i)(1:index(matname(i),' ')-1)
            write(*,*) '         in a (thermo)mechanical calculation'
            write(*,*)
          endif
        enddo
        if(ierror.eq.0) then
          write(*,*) '*ERROR in calinput: no elastic constants'
          write(*,*) '       were assigned to any material in a'
          write(*,*) '       (thermo)mechanical calculation'
          write(*,*)
          call exit(201)
        endif
      endif
!     
!     check whether the magnetic constants were defined for 
!     electromagnetic calculations
!     
      if((ithermal(1).ne.2).and.solid.and.(nmethod.ge.8)) then
        ierror=0
        do i=1,nmat
          if(nelcon(1,i).ne.0) then
            ierror=ierror+1
          else
            write(*,*)'*WARNING in calinput: no magnetic constants '
            write(*,*)' were assigned to material ',
     &           matname(i)(1:index(matname(i),' ')-1)
            write(*,*) '         in an electromagnetic calculation'
            write(*,*)
          endif
        enddo
        if(ierror.eq.0) then
          write(*,*) '*ERROR in calinput: no magnetic constants'
          write(*,*) '       were assigned to any material in a'
          write(*,*) '       electromagnetic calculation'
          write(*,*)
          call exit(201)
        endif
      endif
!     
!     check whether the conductivity was defined for thermal calculations
!     
      if((ithermal(1).ge.2).and.(nef.eq.0)) then
        ierror=0
        do i=1,nmat
          if(ncocon(1,i).ne.0) then
            ierror=ierror+1
          else
            write(*,*) '*WARNING in calinput: no conductivity '
            write(*,*) 
     &           '         constants were assigned to material ',
     &           matname(i)(1:index(matname(i),' ')-1)
            write(*,*) '         in a thermo(mechanical) calculation'
            write(*,*)
          endif
        enddo
      endif
!     
      if(nef.gt.0) then
        if(iperturb(1).eq.0) then
          iperturb(1)=2
        elseif(iperturb(1).eq.1) then
          write(*,*) '*ERROR in calinput: PERTURBATION and fluids'
          write(*,*) '       are mutually exclusive; '
          ier=1
        endif
      endif
!     
      write(*,*)
      write(*,*) 'STEP ',istep
      write(*,*)
      if(nmethod.eq.-1) then
        write(*,*) 'Visco analysis was selected'
      elseif(nmethod.eq.0) then
        write(*,*) 'No analysis was selected'
      elseif(nmethod.eq.1) then
        write(*,*) 'Static analysis was selected'
      elseif(nmethod.eq.2) then
        write(*,*) 'Frequency analysis was selected'
      elseif(nmethod.eq.3) then
        write(*,*) 'Buckling analysis was selected'
      elseif(nmethod.eq.4) then
        write(*,*) 'Dynamic analysis was selected'
      endif
      write(*,*)
      if(iperturb(1).eq.1) then
        write(*,*) 'Perturbation parameter is active'
        write(*,*)
      elseif(iperturb(1).eq.3) then
        write(*,*) 'Nonlinear material laws are taken into account'
        write(*,*)
      endif
!     
      if(iperturb(1).ge.2) then
        write(*,*) 'Newton-Raphson iterative procedure is active'
        write(*,*)
      endif
!     
      if(iperturb(2).eq.1) then
        write(*,*) 'Nonlinear geometric effects are taken into account'
        write(*,*)
      endif
!     
      timepar(1)=tinc
      timepar(2)=tper
      timepar(3)=tmin
      timepar(4)=tmax
      timepar(5)=tincf
!     
      if(istep.eq.1) ncs_=lprev
!     
      if(ier.ge.1) then
        write(*,*) '*ERROR in calinput: at least one fatal'
        write(*,*) '       error message while reading the'
        write(*,*) '       input deck: CalculiX stops.'
        write(*,*)
        call exit(201)
      endif
!
!     the step number is written in the .dat-file if data output
!     was requested or if the procedure always leads to some
!     output in the .dat file; this is the case for the following
!     procedures:
!     *FREQUENCY
!     *BUCKLE
!     *STEADY STATE DYNAMICS
!     *COMPLEX FREQUENCY,CORIOLIS
!     *COMPLEX FREQUENCY,FLUTTER
!     *SENSITIVITY
!     *ROBUST DESIGN
!     *FEASIBLE DIRECTION
!
      if(((nmethod.eq.2).or.(nmethod.eq.3).or.(nmethod.eq.5).or.
     &     (nmethod.eq.6).or.(nmethod.eq.7).or.(nmethod.eq.12).or.
     &     (nmethod.eq.14).or.(nmethod.eq.16)).or.(nprint.gt.0)) then
        write(5,*)
        write(5,100) istep
 100    format('                        S T E P ',i7)
        write(5,*)
      endif
!     
      return
      end
