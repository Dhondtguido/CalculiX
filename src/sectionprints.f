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
      subroutine sectionprints(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,nset_,nalset,nprint,nprint_,jout,prlab,prset,
     &  sectionprint_flag,ithermal,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,amname,nam,itpamp,idrct,ipoinpc,nef,ier)
!
!     reading the *NODE PRINT cards in the input deck
!
      implicit none
!
      logical sectionprint_flag
!
      character*1 total,nodesys,inpc(*)
      character*6 prlab(*)
      character*80 amname(*),timepointsname
      character*81 set(*),prset(*),noset
      character*132 textpart(16),name
!
      integer istartset(*),iendset(*),ialset(*),ii,i,nam,itpamp,id,
     &  jout(2),joutl,ithermal(*),nset,nset_,nalset,nprint,nprint_,
     &  istat,n,key,ipos,iline,ipol,inl,ipoinp(2,*),inp(3,*),idrct,
     &  ipoinpc(0:*),nef,ier,istep
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *SECTION PRINT: *SECTION PRINT'
         write(*,*) '       should only be used within a *STEP' 
         write(*,*) '       definition'
         ier=1
         return
      endif
!
      nodesys='G'
!
!     reset the facial print requests (nodal and element print requests, 
!     if any,are kept)
!
      if(.not.sectionprint_flag) then
         ii=0
         do i=1,nprint
            if((prlab(i)(1:4).eq.'DRAG').or.(prlab(i)(1:4).eq.'FLUX')
     &     .or.(prlab(i)(1:3).eq.'SOF').or.(prlab(i)(1:3).eq.'SOM')
     &     .or.(prlab(i)(1:6).eq.'SOAREA'))
     &           cycle
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
      name(1:1)=' '
      do ii=2,n
        if(textpart(ii)(1:8).eq.'SURFACE=') then
          noset(1:80)=textpart(ii)(9:88)
          ipos=index(noset,' ')
          noset(ipos:ipos)='T'
c          do i=1,nset
c            if(set(i).eq.noset) exit
c          enddo
          call cident81(set,noset,nset,id)
          i=nset+1
          if(id.gt.0) then
            if(noset.eq.set(id)) then
              i=id
            endif
          endif
          if(i.gt.nset) then
             write(*,*) 
     &           '*WARNING reading *SECTION PRINT: element surface ',
     &            noset(1:ipos-1),' does not exist'
             call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &            ipoinp,inp,ipoinpc)
             return
          endif
        elseif(textpart(ii)(1:11).eq.'FREQUENCYF=') then
           read(textpart(ii)(12:21),'(i10)',iostat=istat) joutl
           if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*SECTION PRINT%",ier)
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
     &          '*ERROR reading *SECTION PRINT: time points definition '
     &               ,timepointsname(1:ipos-1),' is unknown or empty'
              ier=1
              return
           endif
           if(idrct.eq.1) then
              write(*,*) 
     &            '*ERROR reading *SECTION PRINT: the DIRECT option'
              write(*,*) '       collides with a TIME POINTS '
              write(*,*) '       specification'
              ier=1
              return
           endif
           jout(1)=1
           jout(2)=1
        elseif(textpart(ii)(1:5).eq.'NAME=') then
            name(1:127)=textpart(ii)(6:132)
        else
            write(*,*) 
     &      '*WARNING reading *SECTION PRINT: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(ii)(1:index(textpart(ii),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*SECTION PRINT%")
        endif
      enddo
!
      if(name(1:1).eq.' ') then
         write(*,*)
     &        '*ERROR reading *SECTION PRINT: no NAME given'
         write(*,*) '  '
         call inputerror(inpc,ipoinpc,iline,
     &        "*SECTION PRINT%",ier)
         return
      endif
!
!     check whether a set was defined
!
      if(noset(1:1).eq.' ') then
         write(*,*) 
     &       '*WARNING reading *SECTION PRINT: no set was defined'
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         return
      endif
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if(key.eq.1) exit
         loop: do ii=1,n
            if((textpart(ii)(1:4).ne.'DRAG').and.
     &         (textpart(ii)(1:4).ne.'FLUX').and.
     &         (textpart(ii)(1:3).ne.'SOF').and.
     &         (textpart(ii)(1:3).ne.'SOM').and.
     &         (textpart(ii)(1:6).ne.'SOAREA')) then
               write(*,*) 
     &           '*WARNING reading *SECTION PRINT: label not applicable'
               write(*,*) '         or unknown; '
               call inputwarning(inpc,ipoinpc,iline,
     &"*SECTION PRINT%")
               cycle
            endif
            if((nef.eq.0).and.(textpart(ii)(1:4).eq.'DRAG')) then
               write(*,*) 
     &           '*WARNING reading *SECTION PRINT: DRAG only makes '
               write(*,*) '         sense for 3D fluid '
               write(*,*) '         calculations'
               cycle
            endif
!
!           SOF, SOM and SOAREA generate the same output
!
            if(textpart(ii)(1:3).eq.'SOM') textpart(ii)(1:3)='SOF'
            if(textpart(ii)(1:6).eq.'SOAREA') textpart(ii)(1:6)='SOF   '
            do i=1,nprint
               if(prlab(i)(1:3).eq.'SOF') then
                  if(prset(i).eq.noset) cycle loop
               endif
            enddo
!
            nprint=nprint+1
            if(nprint.gt.nprint_) then
               write(*,*) 
     &            '*ERROR reading *SECTION PRINT: increase nprint_'
               ier=1
               return
            endif
            prset(nprint)=noset
            prlab(nprint)(1:4)=textpart(ii)(1:4)
            prlab(nprint)(5:5)=total
            prlab(nprint)(6:6)=nodesys
         enddo loop
      enddo
!
      return
      end

