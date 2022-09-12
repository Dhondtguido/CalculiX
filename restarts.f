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
      subroutine restarts(istep,nset,nload,nforc, nboun,nk,ne,
     &  nmpc,nalset,nmat,ntmat_,npmat_,norien,nam,nprint,mi,
     &  ntrans,ncs_,namtot,ncmat_,mpcfree,maxlenmpc,
     &  ne1d,ne2d,nflow,nlabel,iplas,
     &  nkon,ithermal,nmethod,iperturb,nstate_,nener,set,istartset,
     &  iendset,ialset,co,kon,ipkon,lakon,nodeboun,ndirboun,iamboun,
     &  xboun,ikboun,ilboun,ipompc,nodempc,coefmpc,labmpc,ikmpc,ilmpc,
     &  nodeforc,ndirforc,iamforc,xforc,ikforc,ilforc,nelemload,iamload,
     &  sideload,xload,elcon,nelcon,rhcon,nrhcon,
     &  alcon,nalcon,alzero,plicon,nplicon,plkcon,nplkcon,orname,orab,
     &  ielorien,trab,inotr,amname,amta,namta,t0,t1,iamt1,veold,
     &  ielmat,matname,prlab,prset,filab,vold,nodebounold,
     &  ndirbounold,xbounold,xforcold,xloadold,t1old,eme,
     &  iponor,xnor,knor,thicke,offset,iponoel,inoel,rig,
     &  shcon,nshcon,cocon,ncocon,ics,sti,
     &  ener,xstate,jobnamec,infree,irstrt,inpc,textpart,istat,n,
     &  key,prestr,iprestr,cbody,ibody,xbody,nbody,xbodyold,
     &  ttime,qaold,cs,mcs,output,physcon,ctrl,typeboun,iline,ipol,inl,
     &  ipoinp,inp,fmpc,tieset,ntie,tietol,ipoinpc,nslavs,t0g,t1g,nprop,
     &  ielprop,prop,mortar,nintpoint,ifacecount,islavsurf,pslavsurf,
     &  clearini,ier,vel,nef,velo,veloo,ne2boun,heading,network,
     &  irestartread,nfc,ndc,coeffc,ikdc,edc)
!
      implicit none
!
      character*1 typeboun(*),inpc(*)
      character*4 output
      character*6 prlab(*)
      character*8 lakon(*)
      character*20 labmpc(*),sideload(*)
      character*66 heading(*)
      character*80 orname(*),amname(*),matname(*)
      character*81 set(*),prset(*),tieset(3,*),cbody(*)
      character*87 filab(*)
      character*132 jobnamec(*),textpart(16)
!
      integer istep,nset,nload,nforc,nboun,nk,ne,nmpc,nalset,nmat,
     &  ntmat_,npmat_,norien,nam,nprint,mi(*),ntrans,ncs_,
     &  namtot,ncmat_,mpcfree,ne1d,ne2d,nflow,nlabel,iplas,nkon,
     &  ithermal(*),nmethod,iperturb(*),nstate_,istartset(*),iendset(*),
     &  ialset(*),kon(*),ipkon(*),nodeboun(*),ndirboun(*),iamboun(*),
     &  ikboun(*),ilboun(*),ipompc(*),nodempc(*),ikmpc(*),ilmpc(*),
     &  nodeforc(*),ndirforc(*),iamforc(*),ikforc(*),ilforc(*),
     &  nelemload(*),iamload(*),nelcon(*),ipoinpc(0:*),
     &  nrhcon(*),nalcon(*),nplicon(*),nplkcon(*),ielorien(*),
     &  inotr(*),nprop,ielprop(*),mortar,nintpoint,ifacecount,
     &  namta(*),iamt1(*),ielmat(*),nodebounold(*),ndirbounold(*),
     &  iponor(*),knor(*),iponoel(*),inoel(*),rig(*),islavsurf(*),
     &  nshcon(*),ncocon(*),ics(*),infree(*),ier,network,
     &  nener,irestartstep,irestartread,irstrt(*),istat,n,i,key,
     &  iprestr,mcs,maxlenmpc,iline,ipol,inl,nfc,ndc,ikdc(*),
     &  ipoinp(2,*),inp(3,*),ntie,ibody(*),nbody,nslavs,nef,
     &  ne2boun(2,*)
!
      real*8 co(*),xboun(*),coefmpc(*),xforc(*),xload(*),elcon(*),
     &  rhcon(*),alcon(*),alzero(*),plicon(*),plkcon(*),orab(*),
     &  trab(*),amta(*),t0(*),t1(*),prestr(*),veold(*),tietol(*),
     &  vold(*),xbounold(*),xforcold(*),xloadold(*),t1old(*),eme(*),
     &  xnor(*),thicke(*),offset(*),t0g(*),t1g(*),clearini(*),
     &  shcon(*),cocon(*),sti(*),ener(*),xstate(*),prop(*),edc(*),
     &  ttime,qaold(2),cs(17,*),physcon(*),pslavsurf(*),coeffc(*),
     &  ctrl(*),fmpc(*),xbody(*),xbodyold(*),vel(*),velo(*),veloo(*)
!
      irestartstep=0
!
      do i=2,n
         if(textpart(i)(1:4).eq.'READ') then
            irestartread=2
         elseif(textpart(i)(1:5).eq.'STEP=') then
            read(textpart(i)(6:15),'(i10)',iostat=istat) irestartstep
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*RESTART%",ier)
               return
            endif
         elseif(textpart(i)(1:5).eq.'WRITE') then
            irstrt(1)=1
         elseif(textpart(i)(1:10).eq.'FREQUENCY=') then
            read(textpart(i)(11:20),'(i10)',iostat=istat) irstrt(1)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*RESTART%",ier)
               return
            endif
         elseif(textpart(i)(1:7).eq.'OVERLAY') then
            irstrt(2)=1
         else
            write(*,*) 
     &        '*WARNING in restarts: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*RESTART%")
         endif
      enddo
!
      if(irestartread.eq.2) then
        irestartread=1
        call restartread(istep,nset,nload,nforc, nboun,nk,ne,
     &  nmpc,nalset,nmat,ntmat_,npmat_,norien,nam,nprint,mi,
     &  ntrans,ncs_,namtot,ncmat_,mpcfree,maxlenmpc,
     &  ne1d,ne2d,nflow,nlabel,iplas,
     &  nkon,ithermal,nmethod,iperturb,nstate_,nener,set,istartset,
     &  iendset,ialset,co,kon,ipkon,lakon,nodeboun,ndirboun,iamboun,
     &  xboun,ikboun,ilboun,ipompc,nodempc,coefmpc,labmpc,ikmpc,ilmpc,
     &  nodeforc,ndirforc,iamforc,xforc,ikforc,ilforc,nelemload,iamload,
     &  sideload,xload,elcon,nelcon,rhcon,nrhcon,
     &  alcon,nalcon,alzero,plicon,nplicon,plkcon,nplkcon,orname,orab,
     &  ielorien,trab,inotr,amname,amta,namta,t0,t1,iamt1,veold,
     &  ielmat,matname,prlab,prset,filab,vold,nodebounold,
     &  ndirbounold,xbounold,xforcold,xloadold,t1old,eme,
     &  iponor,xnor,knor,thicke,offset,iponoel,inoel,rig,
     &  shcon,nshcon,cocon,ncocon,ics,sti,
     &  ener,xstate,jobnamec,infree,irestartstep,prestr,iprestr,
     &  cbody,ibody,xbody,nbody,xbodyold,ttime,qaold,cs,mcs,
     &  output,physcon,ctrl,typeboun,fmpc,tieset,ntie,tietol,nslavs,
     &  t0g,t1g,nprop,ielprop,prop,mortar,nintpoint,ifacecount,
     &  islavsurf,pslavsurf,clearini,irstrt,vel,nef,velo,veloo,
     &  ne2boun,heading,network,nfc,ndc,coeffc,ikdc,edc)
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end


