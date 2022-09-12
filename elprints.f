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
      subroutine elprints(inpc,textpart,set,
     &  nset,nprint,nprint_,jout,prlab,prset,
     &  nmethod,elprint_flag,nener,ithermal,istep,istat,n,iline,ipol,
     &  inl,ipoinp,inp,amname,nam,itpamp,idrct,ipoinpc,nef,ier)
!
!     reading the *ELEMENT PRINT cards in the input deck
!
      implicit none
!
      logical elprint_flag
!
      character*1 total,elemsys,inpc(*)
      character*6 prlab(*)
      character*80 amname(*),timepointsname
      character*81 set(*),elset,prset(*)
      character*132 textpart(16)
!
      integer nset,nprint,nprint_,istep,istat,n,i,ii,key,ier,id,
     &  jout(2),joutl,ipos,nmethod,nener,ithermal(*),iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),nam,itpamp,idrct,ipoinpc(0:*),nef
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *EL PRINT: *EL PRINT should only be'
         write(*,*) '       used within a *STEP definition'
         ier=1
         return
      endif
!
      elemsys='L'
!
!     reset the element print requests
!
      if(.not.elprint_flag) then
         ii=0
         do i=1,nprint
            if((prlab(i)(1:4).eq.'S   ').or.
     &         (prlab(i)(1:4).eq.'E   ').or.
     &         (prlab(i)(1:4).eq.'ME  ').or.
     &         (prlab(i)(1:4).eq.'PEEQ').or.
     &         (prlab(i)(1:4).eq.'ENER').or.
     &         (prlab(i)(1:4).eq.'SDV ').or.
     &         (prlab(i)(1:4).eq.'ELSE').or.
     &         (prlab(i)(1:4).eq.'ELKE').or.
     &         (prlab(i)(1:4).eq.'EVOL').or.
     &         (prlab(i)(1:4).eq.'EMAS').or.
     &         (prlab(i)(1:4).eq.'EBHE').or.
     &         (prlab(i)(1:4).eq.'SVF ').or.
     &         (prlab(i)(1:4).eq.'HFLF').or.
     &         (prlab(i)(1:4).eq.'COOR').or.
     &         (prlab(i)(1:4).eq.'CENT').or.
     &         (prlab(i)(1:4).eq.'HFL ')) cycle
            ii=ii+1
            prlab(ii)=prlab(i)
            prset(ii)=prset(i)
         enddo
         nprint=ii
      endif
!
      do ii=1,81
         elset=' '
      enddo
      total=' '
!
      do ii=2,n
        if(textpart(ii)(1:6).eq.'ELSET=') then
          elset(1:80)=textpart(ii)(7:86)
          ipos=index(elset,' ')
          elset(ipos:ipos)='E'
c          do i=1,nset
c            if(set(i).eq.elset) exit
c          enddo
          call cident81(set,elset,nset,id)
          i=nset+1
          if(id.gt.0) then
            if(elset.eq.set(id)) then
              i=id
            endif
          endif
          if(i.gt.nset) then
             write(*,*) '*WARNING reading *EL PRINT: elementset ',
     &           elset(1:ipos-1),' does not exist'
             write(*,*)
             call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &            ipoinp,inp,ipoinpc)
             return
          endif
        elseif(textpart(ii)(1:10).eq.'FREQUENCY=') then
           read(textpart(ii)(11:20),'(i10)',iostat=istat) joutl
           if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*EL PRINT%",ier)
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
     &             "*EL PRINT%",ier)
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
           elemsys='G'
        elseif(textpart(ii)(1:9).eq.'GLOBAL=NO') then
           elemsys='L'
        elseif(textpart(ii)(1:11).eq.'TIMEPOINTS=') then
           timepointsname=textpart(ii)(12:91)
           do i=1,nam
              if(amname(i).eq.timepointsname) then
                 itpamp=i
                 exit
              endif
           enddo
           if(i.gt.nam) then
              write(*,*) '*ERROR elprints: time'
              write(*,*) '       points definition',
     &               timepointsname,' is unknown'
              ier=1
              return
           endif
           if(idrct.eq.1) then
              write(*,*) '*ERROR reading *EL PRINT: the DIRECT option'
              write(*,*) '       collides with a TIME POINTS '
              write(*,*) '       specification'
              ier=1
              return
           endif
           jout(1)=1
           jout(2)=1
         else
            write(*,*) 
     &        '*WARNING reading *EL PRINT: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(ii)(1:index(textpart(ii),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*EL PRINT%")
        endif
      enddo
!
!     check whether a set was defined
!
      if(elset.eq.'                     ') then
         write(*,*) '*WARNING reading *EL PRINT: no set was defined'
         write(*,*)
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
            if(textpart(ii)(1:4).eq.'PEEQ') then
               if((nmethod.eq.2).or.(nmethod.eq.3)) then
                  write(*,*) 
     &               '*WARNING reading *EL PRINT: selection of PEEQ'
                  write(*,*) '         does not make sense for a'
                  write(*,*) '         frequency or bucking calculation'
                  write(*,*)
                  cycle
               endif
            elseif((textpart(ii)(1:4).eq.'CEEQ').or.
     &              ((textpart(ii)(1:2).eq.'CE').and.
     &               (textpart(ii)(1:4).ne.'CENT')).or.
     &             (textpart(ii)(1:2).eq.'PE')) then
               if((nmethod.eq.2).or.(nmethod.eq.3)) then
                  write(*,*) 
     &       '*WARNING reading *EL PRINT: selection of CEEQ or CE or PE'
                  write(*,*) '         does not make sense for a'
                  write(*,*) '         frequency or bucking calculation'
                  write(*,*)
                  cycle
               endif
               textpart(ii)(1:4)='PEEQ'
               write(*,*) 
     &       '*WARNING reading *EL PRINT: selection of CEEQ or CE or PE'
               write(*,*)
     &            '         is converted into PEEQ; no distinction'
               write(*,*) 
     &          '        is made between PEEQ, CEEQ, CE and PE'
               write(*,*)
            elseif(textpart(ii)(1:3).eq.'SDV') then
               if((nmethod.eq.2).or.(nmethod.eq.3)) then
                  write(*,*) 
     &              '*WARNING reading *EL PRINT: selection of SDV'
                  write(*,*) '         does not make sense for a'
                  write(*,*) '         frequency or bucking calculation'
                  write(*,*)
                  cycle
               endif
            elseif((textpart(ii)(1:4).eq.'ENER').or.
     &             (textpart(ii)(1:4).eq.'ELSE').or.
     &             (textpart(ii)(1:4).eq.'ELKE')) then
               nener=1
            elseif(textpart(ii)(1:4).eq.'HFL ') then
               if(ithermal(1).lt.2) then
                  write(*,*) 
     &               '*WARNING reading *EL PRINT: HFL only makes '
                  write(*,*) '         sense for heat transfer '
                  write(*,*) '         calculations'
                  write(*,*)
                  cycle
               endif
            elseif((textpart(ii)(1:4).eq.'SVF ').or.
     &             (textpart(ii)(1:4).eq.'SF  ').or.
     &             (textpart(ii)(1:4).eq.'HFLF')) then
               if(nef.eq.0) then
                  write(*,*) 
     &               '*ERROR reading *EL PRINT: SVF or HFLF only'
                  write(*,*) '       make sense for 3D fluid'
                  write(*,*) '       calculations; '
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*EL PRINT%",ier)
                  return
               endif
            elseif((textpart(ii)(1:4).ne.'S   ').and.
     &             (textpart(ii)(1:4).ne.'E   ').and.
     &             (textpart(ii)(1:4).ne.'ME  ').and.
     &             (textpart(ii)(1:4).ne.'EVOL').and.
     &             (textpart(ii)(1:4).ne.'EMAS').and.
     &             (textpart(ii)(1:4).ne.'COOR').and.
     &             (textpart(ii)(1:4).ne.'CENT').and.
     &             (textpart(ii)(1:4).ne.'EBHE')) then
               write(*,*) 
     &             '*ERROR reading *EL PRINT: label not applicable'
               write(*,*) '       or unknown; '
               call inputerror(inpc,ipoinpc,iline,
     &              "*EL PRINT%",ier)
               return
            endif
            nprint=nprint+1
            if(nprint.gt.nprint_) then
               write(*,*) '*ERROR reading *EL PRINT: increase nprint_'
               ier=1
               return
            endif
            prset(nprint)=elset
            prlab(nprint)(1:4)=textpart(ii)(1:4)
            prlab(nprint)(5:5)=total
            prlab(nprint)(6:6)=elemsys
         enddo
      enddo
!
      return
      end

