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
      subroutine noelfiles(inpc,textpart,jout,filab,nmethod,
     &  nodefile_flag,elfile_flag,ifile_output,nener,ithermal,
     &  istep,istat,n,iline,ipol,inl,ipoinp,inp,out3d,nlabel,
     &  amname,nam,itpamp,idrct,ipoinpc,nef,contactfile_flag,
     &  set,nset,xmodal,ier,physcon,output)
!
!     reading the *NODE FILE, *EL FILE and *CONTACT FILE cards in the 
!     input deck
!
      implicit none
!
      logical nodefile_flag,elfile_flag,out3d,sectionforces,
     &  contactfile_flag
!
      character*1 nodesys,elemsys,inpc(*)
      character*4 output
      character*80 amname(*),timepointsname
      character*81 noset,set(*)
      character*87 filab(*)
      character*132 textpart(16)
!
      integer istep,istat,n,key,ii,jout(2),joutl,nmethod,nener,
     &  ithermal(*),ier,id,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),j,nlabel,nam,itpamp,i,
     &  idrct,ipoinpc(0:*),nef,ifile_output,ipos,nset
!
      real*8 xmodal(*),physcon(*)
!
      save sectionforces
!
!     default ist the global system
!
      nodesys='G'
      elemsys='G'
!
      ipos=0
            noset='
     &                           '
!
      if(istep.lt.1) then
         write(*,*) 
     &'*ERROR reading *NODE/EL/CONTACT FILE: *NODE FILE, *EL FILE'
         write(*,*) '       *CONTACT FILE'
         write(*,*) '       should only be used within a *STEP' 
         write(*,*) '       definition'
         ier=1
         return
      endif
!
      if(ifile_output.eq.1) then
!
!        reset the nodal output requests
!
         if(.not.nodefile_flag) then
            filab(1)(1:2)='  '
            if((.not.elfile_flag).and.(.not.contactfile_flag)) then
               filab(1)(3:4)='  '
            endif
            filab(2)(1:4)='    '
            filab(5)(1:4)='    '
            do j=10,12
               filab(j)(1:4)='    '
            enddo
            do j=14,17
               filab(j)(1:4)='    '
            enddo
            filab(19)(1:4)='    '
            do j=21,25
               filab(j)(1:4)='    '
            enddo
            filab(28)(1:4)='    '
            filab(29)(1:4)='    '
            filab(31)(1:4)='    '
            do j=34,38
               filab(j)(1:4)='    '
            enddo
            filab(43)(1:4)='    '
            filab(47)(1:4)='    '
            filab(49)(1:4)='    '
            filab(50)(1:4)='    '
            filab(54)(1:4)='    '
!
            filab(1)(6:87)=' '
            filab(2)(6:87)=' '
            filab(5)(6:87)=' '
            do j=10,12
               filab(j)(6:87)=' '
            enddo
            do j=14,17
               filab(j)(6:87)=' '
            enddo
            filab(19)(6:87)=' '
            do j=21,25
               filab(j)(6:87)=' '
            enddo
            filab(28)(6:87)=' '
            filab(29)(6:87)=' '
            filab(31)(6:87)=' '
            do j=34,38
               filab(j)(6:87)=' '
            enddo
            filab(43)(6:87)=' '
            filab(47)(6:87)=' '
            filab(49)(6:87)=' '
            filab(50)(6:87)=' '
            filab(54)(6:87)=' '
         endif
      elseif(ifile_output.eq.2) then
!
!        reset the element output requests
!
         if(.not.elfile_flag) then
!
!           reset "last iterations" and "contact elements"
!
            if((.not.nodefile_flag).and.(.not.contactfile_flag)) then
               filab(1)(3:4)='  '
            endif
            filab(3)(1:4)='    '
            filab(4)(1:4)='    '
            do j=6,9
               filab(j)(1:4)='    '
            enddo
            filab(13)(1:4)='    '
            filab(18)(1:4)='    '
            filab(20)(1:4)='    '
            filab(30)(1:4)='    '
            filab(32)(1:4)='    '
            do j=39,42
               filab(j)(1:4)='    '
            enddo
            filab(44)(1:4)='    '
            filab(45)(1:4)='    '
            do j=51,53
               filab(j)(1:4)='    '
            enddo
            filab(55)(1:4)='    '
!
            filab(3)(6:87)=' '
            filab(4)(6:87)=' '
            do j=6,9
               filab(j)(6:87)=' '
            enddo
            filab(13)(6:87)=' '
            filab(18)(6:87)=' '
            filab(20)(6:87)=' '
            filab(30)(6:87)=' '
            filab(32)(6:87)=' '
            do j=39,42
               filab(j)(6:87)=' '
            enddo
            filab(44)(6:87)=' '
            filab(45)(6:87)=' '
            do j=51,53
               filab(j)(6:87)=' '
            enddo
            filab(55)(6:87)=' '
!
            sectionforces=.false.
         endif
      elseif(ifile_output.eq.3)then
!
!     reset the contact print requests
!         
         if(.not.contactfile_flag) then
!
!           reset "last iterations" and "contact elements"
!
            if((.not.nodefile_flag).and.(.not.elfile_flag)) then
               filab(1)(3:4)='  '
            endif
            filab(26)(1:4)='    '
            filab(26)(6:87)='    '
            filab(27)(1:4)='    '
            filab(27)(6:87)='    '
            filab(46)(1:4)='    '
            filab(46)(6:87)='    '
         endif
      endif
!
      do ii=2,n
        if(textpart(ii)(1:10).eq.'FREQUENCY=') then
           read(textpart(ii)(11:20),'(i10)',iostat=istat) joutl
           if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &"*NODE FILE/OUTPUT or *EL FILE/OUTPUT or *CONTACT FILE/OUTPUT %",
     &             ier)
              return
           endif
           if(joutl.eq.0) then
              do
                 call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                inl,ipoinp,inp,ipoinpc)
                 if((key.eq.1).or.(istat.lt.0)) return
              enddo
           endif
           if(joutl.gt.0) then
              jout(1)=joutl
              itpamp=0
           endif
        elseif(textpart(ii)(1:11).eq.'FREQUENCYF=') then
           read(textpart(ii)(12:21),'(i10)',iostat=istat) joutl
           if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &"*NODE FILE/OUTPUT or *EL FILE/OUTPUT or *CONTACT FILE/OUTPUT %",
     &             ier)
              return
           endif
           if(joutl.eq.0) then
              do
                 call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                inl,ipoinp,inp,ipoinpc)
                 if((key.eq.1).or.(istat.lt.0)) return
              enddo
           endif
           if(joutl.gt.0) then
              jout(2)=joutl
              itpamp=0
           endif
        elseif(textpart(ii)(1:10).eq.'GLOBAL=YES') then
           nodesys='G'
           elemsys='G'
        elseif(textpart(ii)(1:9).eq.'GLOBAL=NO') then
           nodesys='L'
           elemsys='L'
        elseif((textpart(ii)(1:9).eq.'OUTPUT=2D').or.
     &         (textpart(ii)(1:9).eq.'OUTPUT=2d'))then
           if(istep.eq.1) then
              out3d=.false.
              do j=1,nlabel
                 if(filab(j)(5:5).eq.'E') filab(j)(5:5)='I'
              enddo
           elseif(out3d) then
              write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: OUTPUT=2D has no'
              write(*,*) '         effect in all but the first step'
           endif
        elseif((textpart(ii)(1:9).eq.'OUTPUT=3D').or.
     &         (textpart(ii)(1:9).eq.'OUTPUT=3d'))then
           if(istep.eq.1) then
              out3d=.true.
              do j=1,nlabel
                 filab(j)(5:5)='E'
              enddo
           elseif(.not.out3d) then
              write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: OUTPUT=3D has no'
              write(*,*) '         effect in all but the first step'
           endif
        elseif(textpart(ii)(1:13).eq.'SECTIONFORCES') then
              filab(3)(5:5)='M'
           elseif(textpart(ii)(1:9).eq.'OUTPUTALL') then
              output(4:4)='a'
        elseif(textpart(ii)(1:11).eq.'TIMEPOINTS=') then
           timepointsname=textpart(ii)(12:91)
           do i=1,nam
              if(amname(i).eq.timepointsname) then
                 itpamp=i
                 exit
              endif
           enddo
           if(i.gt.nam) then
              write(*,*) '*ERROR reading *NODE/EL/CONTACT FILE: time'
              write(*,*) '       points definition',
     &               timepointsname,' is unknown'
              ier=1
              return
           endif
           if(idrct.eq.1) then
              write(*,*) 
     &'*ERROR reading *NODE/EL/CONTACT FILE: the DIRECT option'
              write(*,*) '       collides with a TIME POINTS '
              write(*,*) '       specification'
              ier=1
              return
           endif
           jout(1)=1
           jout(2)=1
        elseif(textpart(ii)(1:5).eq.'NSET=') then
           noset=textpart(ii)(6:85)
           noset(81:81)=' '
           ipos=index(noset,' ')
           noset(ipos:ipos)='N'
        elseif(textpart(ii)(1:14).eq.'LASTITERATIONS') then
           filab(1)(4:4)='I'
        elseif(textpart(ii)(1:15).eq.'CONTACTELEMENTS') then
           filab(1)(3:3)='C'
        elseif(textpart(ii)(1:6).eq.'DOUBLE') then
           output(1:3)='dbi '
        else
            write(*,*) 
     &             '*WARNING reading *NODE/EL/CONTACT FILE:' 
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(ii)(1:index(textpart(ii),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*NODE FILE/OUTPUT or *EL FILE/OUTPUT or *CONTACT FILE/OUTPUT %")
        endif
      enddo
!
!     check whether SECTION FORCES and OUTPUT=3D are both active
!
      if((filab(3)(5:5).eq.'M').and.(out3d)) then
         write(*,*) 
     &'*ERROR reading *NODE/EL/CONTACT FILE: SECTION FORCES and'
         write(*,*) '       OUTPUT=3D are mutually exclusive'
         call inputerror(inpc,ipoinpc,iline,
     &"*NODE FILE/OUTPUT or *EL FILE/OUTPUT or *CONTACT FILE/OUTPUT %",
     &        ier)
         return
      endif
!
!     check the existence of the node set (if any was specified)
!
      if(ipos.ne.0) then
c         do i=1,nset
c            if(set(i).eq.noset) exit
c         enddo
        call cident81(set,noset,nset,id)
        i=nset+1
        if(id.gt.0) then
          if(noset.eq.set(id)) then
            i=id
          endif
        endif
         if(i.gt.nset) then
            noset(ipos:ipos)=' '
            write(*,*) 
     &        '*ERROR reading *NODE/EL/CONTACT FILE: node set ',noset
            write(*,*) '  has not yet been defined.'
            ier=1
            return
         endif
      endif
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((key.eq.1).or.(istat.lt.0)) return
         do ii=1,n
            if(textpart(ii)(1:4).eq.'U   ') then
               filab(1)(1:2)='U '
               filab(1)(6:6)=nodesys
               filab(1)(7:87)=noset
            elseif(textpart(ii)(1:4).eq.'NT  ') then
               filab(2)(1:4)='NT  '
               filab(2)(6:6)=nodesys
               filab(2)(7:87)=noset
            elseif(textpart(ii)(1:4).eq.'S   ') then
               filab(3)(1:4)='S   '
               filab(3)(6:6)=elemsys
               filab(3)(7:87)=noset
!
!              next 3 lines introduced on 01.12.2017: stress
!              request automatically induces stress error request
!
               filab(13)(1:4)='ERR '
               filab(13)(6:6)=elemsys
               filab(13)(7:87)=noset
            elseif(textpart(ii)(1:4).eq.'E   ') then
               filab(4)(1:4)='E   '
               filab(4)(6:6)=elemsys
               filab(4)(7:87)=noset
            elseif(textpart(ii)(1:4).eq.'RF  ') then
               filab(5)(1:4)='RF  '
               filab(5)(6:6)=nodesys
               filab(5)(7:87)=noset
            elseif(textpart(ii)(1:4).eq.'PEEQ') then
               if((nmethod.eq.2).or.(nmethod.eq.3)) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: selection of PEEQ'
                  write(*,*) '         does not make sense for a'
                  write(*,*) '         frequency or bucking calculation'
               else
                  filab(6)(1:4)='PEEQ'
                  filab(6)(6:6)=elemsys
                  filab(6)(7:87)=noset
               endif
            elseif(((textpart(ii)(1:4).eq.'CEEQ').or.
     &              (textpart(ii)(1:2).eq.'CE').or.
     &              (textpart(ii)(1:2).eq.'PE')).and.
     &             (textpart(ii)(1:4).ne.'CELS')) then
               textpart(ii)(1:4)='PEEQ'
               if((nmethod.eq.2).or.(nmethod.eq.3)) then
                  write(*,*) 
     &                '*WARNING reading *NODE/EL/CONTACT FILE:' 
                  write(*,*) '         selection of CEEQ or CE or PE'
                  write(*,*) '         does not make sense for a'
                  write(*,*) '         frequency or bucking calculation'
               else
                  write(*,*) 
     &            '*WARNING in elprints: selection of CEEQ or CE or PE'
                  write(*,*) '         is converted into PEEQ; no distin
     &ction'
                  write(*,*) 
     &             '        is made between PEEQ, CEEQ, CE and PE'
                  filab(6)(1:4)='PEEQ'
                  filab(6)(6:6)=elemsys
                  filab(6)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'ENER') then
               filab(7)(1:4)='ENER'
               filab(7)(6:6)=elemsys
               filab(7)(7:87)=noset
              nener=1
            elseif(textpart(ii)(1:4).eq.'SDV ') then
               if((nmethod.eq.2).or.(nmethod.eq.3)) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: selection of SDV'
                  write(*,*) '         does not make sense for a'
                  write(*,*) '         frequency or bucking calculation'
               else
                  filab(8)(1:4)='SDV '
                  filab(8)(6:6)=elemsys
                  filab(8)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'HFL ') then
               if(ithermal(1).le.1) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: HFL only makes '
                  write(*,*) '         sense for heat transfer '
                  write(*,*) '          calculations'
               else
                  filab(9)(1:4)='HFL '
                  filab(9)(6:6)=elemsys
                  filab(9)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'RFL ') then
               if(ithermal(1).le.1) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: RFL only makes '
                  write(*,*) '         sense for heat transfer '
                  write(*,*) '          calculations'
               else
                  filab(10)(1:4)='RFL '
                  filab(10)(6:6)=nodesys
                  filab(10)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'PU  ') then
               if((nmethod.ne.2).and.(nmethod.ne.5).and.
     &            (nmethod.ne.6).and.(nmethod.ne.7)) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: PU only makes'
                  write(*,*) '         sense for frequency and steady'
                  write(*,*) '         state dynamics calculations'
               elseif((nmethod.eq.5).and.(xmodal(7).gt.0.d0)) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: PU does not make'
                  write(*,*) '         sense for nonharmonic periodic'
                  write(*,*) '         excitations; use U instead'
               else
                  filab(11)(1:4)='PU  '
                  filab(11)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'PNT ') then
               filab(12)(1:4)='PNT '
               filab(12)(7:87)=noset
            elseif(textpart(ii)(1:3).eq.'ZZS') then
               filab(13)(1:4)='ZZS '
               filab(13)(6:6)=elemsys
               filab(13)(7:87)=noset
            elseif(textpart(ii)(1:3).eq.'ERR') then
               filab(13)(1:4)='ERR '
               filab(13)(6:6)=elemsys
               filab(13)(7:87)=noset
!
!              next three lines introduced on 01.12.2017
!
            elseif(textpart(ii)(1:3).eq.'NOE') then
               filab(13)(1:4)='    '
               filab(13)(6:87)=' '
            elseif(textpart(ii)(1:4).eq.'TT  ') then
               filab(14)(1:4)='TT  '
               filab(14)(6:6)=nodesys
               filab(14)(7:87)=noset
            elseif(textpart(ii)(1:4).eq.'MF  ') then
               filab(15)(1:4)='MF  '
               filab(15)(6:6)=nodesys
               filab(15)(7:87)=noset
            elseif(textpart(ii)(1:4).eq.'PT  ') then
               filab(16)(1:4)='PT  '
               filab(16)(6:6)=nodesys
               filab(16)(7:87)=noset
            elseif(textpart(ii)(1:4).eq.'TS  ') then
               filab(17)(1:4)='TS  '
               filab(17)(6:6)=nodesys
               filab(17)(7:87)=noset
            elseif(textpart(ii)(1:4).eq.'PHS ') then
               if((nmethod.ne.2).and.(nmethod.ne.5)) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: PHS only makes'
                  write(*,*) '         sense for frequency and steady'
                  write(*,*) '         state dynamics calculations'
               else
                  filab(18)(1:4)='PHS '
                  filab(18)(6:6)=elemsys
                  filab(18)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'MAXU') then
               if(nmethod.ne.2) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: MAXU only makes'
                  write(*,*) '         sense for frequency calculations'
               else
                  filab(19)(1:4)='MAXU'
                  filab(19)(6:6)=nodesys
                  filab(19)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'MAXS') then
               if(nmethod.ne.2) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: MAXS only makes'
                  write(*,*) '         sense for frequency calculations'
               else
                  filab(20)(1:4)='MAXS'
                  filab(20)(6:6)=elemsys
                  filab(20)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'V   ') then
                if((nmethod.eq.4).or.(nef.gt.0)) then
                   filab(21)(1:4)='V   '
                   filab(21)(6:6)=nodesys
                   filab(21)(7:87)=noset
                else
                   write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: V only available'
                   write(*,*) '         for dynamic calculations and'
                   write(*,*) '         3-D fluid calculations'
                endif
            elseif(textpart(ii)(1:4).eq.'PS  ') then
               filab(22)(1:4)='PS  '
               filab(22)(6:6)=nodesys
               filab(22)(7:87)=noset
            elseif(textpart(ii)(1:4).eq.'MACH') then
               if(nef.eq.0) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: MACH only makes'
                  write(*,*) '         sense for 3D fluid calculations'
               else
                  filab(23)(1:4)='MACH'
                  filab(23)(6:6)=nodesys
                  filab(23)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'CP  ') then
               if(nef.eq.0) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: CP only makes'
                  write(*,*) '         sense for 3D fluid calculations'
               else
                  filab(24)(1:4)='CP  '
                  filab(24)(6:6)=nodesys
                  filab(24)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'TURB') then
               if(nef.eq.0) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: TURB only makes'
                  write(*,*) '         sense for 3D fluid calculations'
                elseif((physcon(9).lt.1.d0).or.
     &                 ((physcon(9).gt.10.d0).and.(physcon(9).lt.11.d0))
     &                 )then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: TURB only makes'
                  write(*,*) '         sense for 3D fluid calculations'
                  write(*,*) '         with an active turbulence model'
               else
                  filab(25)(1:4)='TURB'
                  filab(25)(6:6)=nodesys
                  filab(25)(7:87)=noset
               endif
            elseif((textpart(ii)(1:4).eq.'CSTR').or.
     &             (textpart(ii)(1:4).eq.'CDIS')) then
               filab(26)(1:4)='CONT'
               filab(26)(6:6)=nodesys
               filab(26)(7:87)=noset
            elseif(textpart(ii)(1:4).eq.'CELS') then
               filab(27)(1:4)='CELS'
               filab(27)(6:6)=nodesys
               filab(27)(7:87)=noset
               nener=1
            elseif(textpart(ii)(1:4).eq.'DEPT') then
                  filab(28)(1:4)='DEPT'
                  filab(28)(6:6)=nodesys
                  filab(28)(7:87)=noset
            elseif(textpart(ii)(1:4).eq.'HCRI') then
                  filab(29)(1:4)='HCRI'
                  filab(29)(6:6)=nodesys
                  filab(29)(7:87)=noset
            elseif(textpart(ii)(1:4).eq.'MAXE') then
               if(nmethod.ne.2) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: MAXE only makes'
                  write(*,*) '         sense for frequency calculations'
               else
                  filab(30)(1:4)='MAXE'
                  filab(30)(6:6)=elemsys
                  filab(30)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'PRF ') then
               if((nmethod.ne.2).and.(nmethod.ne.5).and.
     &            (nmethod.ne.6).and.(nmethod.ne.7)) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: PRF only makes'
                  write(*,*) '         sense for frequency and steady'
                  write(*,*) '         state dynamics calculations'
               elseif((nmethod.eq.5).and.(xmodal(7).gt.0.d0)) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: PRF does not make'
                  write(*,*) '         sense for nonharmonic periodic'
                  write(*,*) '         excitations; use RF instead'
               else
                  filab(31)(1:4)='PRF '
                  filab(31)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'ME  ') then
               filab(32)(1:4)='ME  '
               filab(32)(6:6)=elemsys
               filab(32)(7:87)=noset
            elseif(textpart(ii)(1:3).eq.'HER') then
               if(ithermal(1).le.1) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: HER only makes '
                  write(*,*) '         sense for heat transfer '
                  write(*,*) '          calculations'
               else
                 filab(33)(1:4)='HER '
                 filab(33)(6:6)=elemsys
                 filab(33)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'VF  ') then
               if(nef.eq.0) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: VF only makes'
                  write(*,*) '         sense for 3D fluid calculations'
               else
                  filab(34)(1:4)='VF  '
                  filab(34)(6:6)=nodesys
                  filab(34)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'PSF ') then
               if(nef.eq.0) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: PSF only makes'
                  write(*,*) '         sense for 3D fluid calculations'
               else
                  filab(35)(1:4)='PSF '
                  filab(35)(6:6)=nodesys
                  filab(35)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'TSF ') then
               if(nef.eq.0) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: TSF only makes'
                  write(*,*) '         sense for 3D fluid calculations'
               else
                  filab(36)(1:4)='TSF '
                  filab(36)(6:6)=nodesys
                  filab(36)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'PTF ') then
               if(nef.eq.0) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: PTF only makes'
                  write(*,*) '         sense for 3D fluid calculations'
               else
                  filab(37)(1:4)='PTF '
                  filab(37)(6:6)=nodesys
                  filab(37)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'TTF ') then
               if(nef.eq.0) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: TTF only makes'
                  write(*,*) '         sense for 3D fluid calculations'
               else
                  filab(38)(1:4)='TTF '
                  filab(38)(6:6)=nodesys
                  filab(38)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'SF  ') then
               if(nef.eq.0) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: SF only makes'
                  write(*,*) '         sense for 3D fluid calculations'
               else
                  filab(39)(1:4)='SF  '
                  filab(39)(6:6)=elemsys
                  filab(39)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'HFLF') then
               if(nef.eq.0) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: HFLF only makes'
                  write(*,*) '         sense for 3D fluid calculations'
               else
                  filab(40)(1:4)='HFLF'
                  filab(40)(6:6)=elemsys
                  filab(40)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'SVF ') then
               if(nef.eq.0) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: SVF only makes'
                  write(*,*) '         sense for 3D fluid calculations'
               else
                  filab(41)(1:4)='SVF '
                  filab(41)(6:6)=elemsys
                  filab(41)(7:87)=noset
               endif
            elseif(textpart(ii)(1:3).eq.'ECD') then
               if(nmethod.lt.8) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: ECD only makes'
                  write(*,*) 
     &         '         sense for electromagnetic calculations'
               else
                  filab(42)(1:4)='ECD '
                  filab(42)(6:6)=elemsys
                  filab(42)(7:87)=noset
               endif
            elseif(textpart(ii)(1:3).eq.'POT') then
               if(nmethod.lt.8) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: POT only makes'
                  write(*,*) 
     &         '         sense for electromagnetic calculations'
               else
                  filab(43)(1:4)='POT '
                  filab(43)(6:6)=nodesys
                  filab(43)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'EMFE') then
               if(nmethod.lt.8) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: EMFE only makes'
                  write(*,*) 
     &         '         sense for electromagnetic calculations'
               else
                  filab(44)(1:4)='EMFE'
                  filab(44)(6:6)=elemsys
                  filab(44)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'EMFB') then
               if(nmethod.lt.8) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: EMFB only makes'
                  write(*,*) 
     &         '         sense for electromagnetic calculations'
               else
                  filab(45)(1:4)='EMFB'
                  filab(45)(6:6)=elemsys
                  filab(45)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'PCON') then
               if((nmethod.ne.2).and.(nmethod.ne.5)) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: PCON only makes'
                  write(*,*) '         sense for frequency and steady'
                  write(*,*) '         state dynamics calculations'
               else
                  filab(46)(1:4)='PCON'
                  filab(46)(6:6)=elemsys
                  filab(46)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'DEPF') then
               if(nef.eq.0) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: DEPF only makes'
                  write(*,*) '         sense for 3D fluid calculations'
               else
                  filab(49)(1:4)='DEPF'
                  filab(49)(6:6)=elemsys
                  filab(49)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'DTF ') then
               if(nef.eq.0) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: DTF only makes'
                  write(*,*) '         sense for 3D fluid calculations'
               else
                  filab(50)(1:4)='DTF '
                  filab(50)(6:6)=elemsys
                  filab(50)(7:87)=noset
               endif
            elseif(textpart(ii)(1:4).eq.'SNEG') then
               filab(51)(1:4)='SNEG'
               filab(51)(6:6)=elemsys
               filab(51)(7:87)=noset
            elseif(textpart(ii)(1:4).eq.'SMID') then
               filab(52)(1:4)='SMID'
               filab(52)(6:6)=elemsys
               filab(52)(7:87)=noset
            elseif(textpart(ii)(1:4).eq.'SPOS') then
               filab(53)(1:4)='SPOS'
               filab(53)(6:6)=elemsys
               filab(53)(7:87)=noset
            elseif(textpart(ii)(1:4).eq.'KEQ ') then
               if(nmethod.ne.15) then
                  write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: KEQ only makes'
                  write(*,*) '         sense for crack propagation'
                  write(*,*) '         calculations'
               else
                  filab(54)(1:4)='KEQ '
                  filab(54)(6:6)=elemsys
                  filab(54)(7:87)=noset
               endif
             elseif(textpart(ii)(1:4).eq.'THE ') then
!
               write(*,*) '*INFO reading *NODE/EL/CONTACT FILE:'
               write(*,*)
     &'      the selection of the thermal strains triggers the'
               write(*,*)
     &'      extrapolation to the nodes and storage of the total'
               write(*,*)
     &'      and the mechanical strains as well; this is needed,'
               write(*,*)
     &'      since the thermal strains at the nodes are computed as'
               write(*,*)
     &'      the difference of the extrapolated total strains to the'
               write(*,*)
     &'      nodes and the extrapolated mechanical strains to the'
               write(*,*)
     &'      nodes'            
               write(*,*)
!
               filab(4)(1:4)='E   '
               filab(4)(6:6)=elemsys
               filab(4)(7:87)=noset
               filab(32)(1:4)='ME  '
               filab(32)(6:6)=elemsys
               filab(32)(7:87)=noset
               filab(55)(1:4)='THE '
               filab(55)(6:6)=elemsys
               filab(55)(7:87)=noset
             else
               write(*,*) 
     &'*WARNING reading *NODE/EL/CONTACT FILE: label not applicable'
               write(*,*) '         or unknown; '
               call inputwarning(inpc,ipoinpc,iline,
     &"*NODE FILE/OUTPUT or *EL FILE/OUTPUT or *CONTACT FILE/OUTPUT %")
            endif
         enddo
      enddo
!
      return
      end






