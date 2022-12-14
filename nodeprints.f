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
      subroutine nodeprints(inpc,textpart,set,istartset,iendset,ialset,
     &  nset,nset_,nalset,nprint,nprint_,jout,prlab,prset,
     &  nodeprint_flag,ithermal,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,amname,nam,itpamp,idrct,ipoinpc,nef,ier)
!
!     reading the *NODE PRINT cards in the input deck
!
      implicit none
!
      logical nodeprint_flag
!
      character*1 total,nodesys,inpc(*)
      character*6 prlab(*)
      character*80 amname(*),timepointsname
      character*81 set(*),prset(*),noset
      character*132 textpart(16)
!
      integer istartset(*),iendset(*),ialset(*),ii,i,nam,itpamp,id,
     &  jout(2),joutl,ithermal(*),nset,nset_,nalset,nprint,nprint_,
     &  istat,n,key,ipos,iline,ipol,inl,ipoinp(2,*),inp(3,*),idrct,
     &  ipoinpc(0:*),nef,ier,istep
!
      if(istep.lt.1) then
         write(*,*) 
     &        '*ERROR reading *NODE PRINT: *NODE PRINT should only be'
         write(*,*) '  used within a *STEP definition'
         ier=1
         return
      endif
!
      nodesys='L'
!
!     reset the nodal print requests (element print requests, if any,
!     are kept)
!
      if(.not.nodeprint_flag) then
         ii=0
         do i=1,nprint
            if((prlab(i)(1:4).eq.'U   ').or.
     &         (prlab(i)(1:4).eq.'NT  ').or.
     &         (prlab(i)(1:4).eq.'TS  ').or.
     &         (prlab(i)(1:4).eq.'RF  ').or.
     &         (prlab(i)(1:4).eq.'RFL ').or.
     &         (prlab(i)(1:4).eq.'PS  ').or.
     &         (prlab(i)(1:4).eq.'PN  ').or.
     &         (prlab(i)(1:4).eq.'MF  ').or.
     &         (prlab(i)(1:4).eq.'VF  ').or.
     &         (prlab(i)(1:4).eq.'PSF ').or.
     &         (prlab(i)(1:4).eq.'TSF ').or.
     &         (prlab(i)(1:4).eq.'MACH').or.
     &         (prlab(i)(1:4).eq.'DEPF').or.
     &         (prlab(i)(1:4).eq.'TTF ').or.
     &         (prlab(i)(1:4).eq.'PTF ').or.
     &         (prlab(i)(1:4).eq.'CP  ').or.
     &         (prlab(i)(1:4).eq.'TURB').or.
     &         (prlab(i)(1:4).eq.'V   ')) cycle
            ii=ii+1
            prlab(ii)=prlab(i)
            prset(ii)=prset(i)
         enddo
         nprint=ii
      endif
!
      do ii=1,81
         noset(ii:ii)=' '
      enddo
      total=' '
!
      do ii=2,n
        if(textpart(ii)(1:5).eq.'NSET=') then
          noset(1:80)=textpart(ii)(6:85)
          ipos=index(noset,' ')
          noset(ipos:ipos)='N'
          call cident81(set,noset,nset,id)
          i=nset+1
          if(id.gt.0) then
            if(noset.eq.set(id)) then
              i=id
            endif
          endif
          if(i.gt.nset) then
             write(*,*) '*WARNING reading *NODE PRINT: node set ',
     &            noset(1:ipos-1),' does not exist'
             call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &            ipoinp,inp,ipoinpc)
             return
          endif
        elseif(textpart(ii)(1:10).eq.'FREQUENCY=') then
           read(textpart(ii)(11:20),'(i10)',iostat=istat) joutl
           if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*NODE PRINT%",ier)
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
     &             "*NODE PRINT%",ier)
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
        elseif(textpart(ii)(1:10).eq.'TOTALS=YES') then
           total='T'
        elseif(textpart(ii)(1:11).eq.'TOTALS=ONLY') then
           total='O'
        elseif(textpart(ii)(1:10).eq.'GLOBAL=YES') then
           nodesys='G'
        elseif(textpart(ii)(1:9).eq.'GLOBAL=NO') then
           nodesys='L'
        elseif(textpart(ii)(1:11).eq.'TIMEPOINTS=') then
           timepointsname=textpart(ii)(12:91)
           do i=1,nam
              if(amname(i).eq.timepointsname) then
                 itpamp=i
                 exit
              endif
           enddo
           if(i.gt.nam) then
              ipos=index(timepointsname,' ')
              write(*,*) 
     &          '*ERROR reading *NODE PRINT: time points definition '
     &               ,timepointsname(1:ipos-1),' is unknown or empty'
              ier=1
              return
           endif
           if(idrct.eq.1) then
              write(*,*) '*ERROR reading *NODE PRINT: the DIRECT option'
              write(*,*) '       collides with a TIME POINTS '
              write(*,*) '       specification'
              ier=1
              return
           endif
           jout(1)=1
           jout(2)=1
         else
            write(*,*) 
     &        '*WARNING in modaldynamics: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(ii)(1:index(textpart(ii),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*NODE PRINT%")
        endif
      enddo
!
!     check whether a set was defined
!
      if(noset(1:1).eq.' ') then
         write(*,*) '*WARNING reading *NODE PRINT: no set was defined'
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         return
      endif
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if(key.eq.1) exit
         do ii=1,n
            if((textpart(ii)(1:4).ne.'U   ').and.
     &         (textpart(ii)(1:4).ne.'NT  ').and.
     &         (textpart(ii)(1:4).ne.'TS  ').and.
     &         (textpart(ii)(1:4).ne.'RF  ').and.
     &         (textpart(ii)(1:4).ne.'RFL ').and.
     &         (textpart(ii)(1:4).ne.'PS  ').and.
     &         (textpart(ii)(1:4).ne.'PN  ').and.
     &         (textpart(ii)(1:4).ne.'MF  ').and.
     &         (textpart(ii)(1:4).ne.'V   ').and.
     &         (textpart(ii)(1:4).ne.'VF  ').and.
     &         (textpart(ii)(1:4).ne.'PSF ').and.
     &         (textpart(ii)(1:4).ne.'TSF ').and.
     &         (textpart(ii)(1:4).ne.'MACH').and.
     &         (textpart(ii)(1:4).ne.'DEPF').and.
     &         (textpart(ii)(1:4).ne.'TTF ').and.
     &         (textpart(ii)(1:4).ne.'PTF ').and.
     &         (textpart(ii)(1:4).ne.'CP  ').and.
     &         (textpart(ii)(1:4).ne.'TURB')) then
               write(*,*) 
     &            '*WARNING reading *NODE PRINT: label not applicable'
               write(*,*) '         or unknown; '
               call inputwarning(inpc,ipoinpc,iline,
     &"*NODE PRINT%")
               cycle
            endif
            if(textpart(ii)(1:4).eq.'RFL ') then
               if(ithermal(1).lt.2) then
                  write(*,*) 
     &              '*WARNING reading *NODE PRINT: RFL only makes '
                  write(*,*) '         sense for heat transfer '
                  write(*,*) '          calculations'
                  cycle
               endif
            elseif((textpart(ii)(1:4).eq.'VF  ').or.
     &         (textpart(ii)(1:4).eq.'PSF ').or.
     &         (textpart(ii)(1:4).eq.'TSF ').or.
     &         (textpart(ii)(1:4).eq.'MACH').or.
     &         (textpart(ii)(1:4).eq.'DEPF').or.
     &         (textpart(ii)(1:4).eq.'TTF ').or.
     &         (textpart(ii)(1:4).eq.'PTF ').or.
     &         (textpart(ii)(1:4).eq.'CP  ').or.
     &         (textpart(ii)(1:4).eq.'TURB')) then
               if(nef.eq.0) then
                  write(*,*) 
     &               '*WARNING reading *NODE PRINT: VF, PSF, TSF,'
                  write(*,*) '         MACH, DEPF, TTF, PTF, CP or '
                  write(*,*) '         TURB only make sense for '
                  write(*,*) '         3D-fluid calculations'
                  cycle
               endif
            endif
            nprint=nprint+1
            if(nprint.gt.nprint_) then
               write(*,*) '*ERROR reading *NODE PRINT: increase nprint_'
               ier=1
               return
            endif
            prset(nprint)=noset
            prlab(nprint)(1:4)=textpart(ii)(1:4)
            prlab(nprint)(5:5)=total
            prlab(nprint)(6:6)=nodesys
         enddo
      enddo
!
      return
      end

