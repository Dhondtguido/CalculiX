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
      subroutine contactprints(inpc,textpart,
     &     nprint,nprint_,jout,prlab,prset,
     &     contactprint_flag,ithermal,istep,istat,n,iline,ipol,inl,
     &     ipoinp,inp,amname,nam,itpamp,idrct,ipoinpc,nener,ier,
     &     ntie,tieset)
!
!     reading the *CONTACT PRINT cards in the input deck
!
      implicit none
!
      logical contactprint_flag
!
      character*1 total,nodesys,inpc(*)
      character*6 prlab(*)
      character*80 amname(*),timepointsname
      character*81 prset(*),noset,mastersurface,slavesurface,tieset(3,*)
      character*132 textpart(16)
!
      integer ii,i,nam,itpamp,ier,ntie,iposslave,iposmaster,
     &  jout(2),joutl,ithermal(*),nprint,nprint_,istep,itie,
     &  istat,n,key,ipos,iline,ipol,inl,ipoinp(2,*),inp(3,*),idrct,
     &  ipoinpc(0:*),nener
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *CONTACT PRINT: *CONTACT PRINT 
     &        should only be'
         write(*,*) '  used within a *STEP definition'
         ier=1
         return
      endif
!
      nodesys='L'
      do i=1,81
         slavesurface(i:i)=' '
         mastersurface(i:i)=' '
      enddo
!
!     reset the contact print requests (node and element print requests,
!     if any, are kept)
!
      if(.not.contactprint_flag) then
         ii=0
         do i=1,nprint
            if((prlab(i)(1:4).eq.'CSTR').or.
     &         (prlab(i)(1:4).eq.'CDIS').or.
     &         (prlab(i)(1:4).eq.'CNUM').or.
     &         (prlab(i)(1:4).eq.'CELS').or.
     &         (prlab(i)(1:4).eq.'CF  ').or.
     &         (prlab(i)(1:4).eq.'CFN ').or.
     &         (prlab(i)(1:4).eq.'CFS ')) cycle
            ii=ii+1
            prlab(ii)=prlab(i)
            prset(ii)=prset(i)
         enddo
         nprint=ii
      endif
!
c      jout=max(jout,1)
      do ii=1,81
         noset(ii:ii)=' '
      enddo
      total=' '
!
      do ii=2,n
         if(textpart(ii)(1:10).eq.'FREQUENCY=') then
            read(textpart(ii)(11:20),'(i10)',iostat=istat) joutl
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*CONTACT PRINT%",ier)
               return
            endif
            if(joutl.eq.0) then
               do
                  call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                 inl,ipoinp,inp,ipoinpc)
                  if((key.eq.1).or.(istat.lt.0)) return
               enddo
            endif
           if(joutl.gt.0) then
              jout(1)=joutl
              itpamp=0
           endif
        elseif(textpart(ii)(1:10).eq.'TOTALS=YES') then
           total='T'
        elseif(textpart(ii)(1:11).eq.'TOTALS=ONLY') then
           total='O'
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
              write(*,*) '*ERROR reading *CONTACT PRINT: time points de
     &finition '
     &               ,timepointsname(1:ipos-1),' is unknown or empty'
              ier=1
              return
           endif
           if(idrct.eq.1) then
              write(*,*) 
     &            '*ERROR reading *CONTACT PRINT: the DIRECT option'
              write(*,*) '       collides with a TIME POINTS '
              write(*,*) '       specification'
              ier=1
              return
           endif
           jout(1)=1
           jout(2)=1
        elseif(textpart(ii)(1:6).eq.'SLAVE=') then
           read(textpart(ii)(7:86),'(a80)',iostat=istat) slavesurface
           if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*CONTACT PRINT%",ier)
              return
           endif
        elseif(textpart(ii)(1:7).eq.'MASTER=') then
           read(textpart(ii)(8:87),'(a80)',iostat=istat) mastersurface
           if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*CONTACT PRINT%",ier)
              return
           endif
        else
            write(*,*) 
     &      '*WARNING reading *CONTACT PRINT: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(ii)(1:index(textpart(ii),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*CONTACT PRINT%")
        endif
      enddo
!
!     determining the appropriate contact tie
!
      itie=0
      iposslave=index(slavesurface(1:80),' ')
      iposmaster=index(mastersurface(1:80),' ')
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.'C') cycle
         ipos=index(tieset(2,i),' ')-1
         if(ipos.ne.iposslave) cycle
         if(tieset(2,i)(1:ipos-1).ne.slavesurface(1:ipos-1)) cycle
         ipos=index(tieset(3,i),' ')-1
         if(ipos.ne.iposmaster) cycle
         if(tieset(3,i)(1:ipos-1).ne.mastersurface(1:ipos-1)) cycle
         itie=i
         exit
      enddo
!     
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if(key.eq.1) exit
         do ii=1,n
            if((textpart(ii)(1:4).ne.'CSTR').and.
     &         (textpart(ii)(1:4).ne.'CELS').and.
     &         (textpart(ii)(1:4).ne.'CNUM').and.
     &         (textpart(ii)(1:4).ne.'CDIS').and.
     &         (textpart(ii)(1:4).ne.'CF  ').and.
     &         (textpart(ii)(1:4).ne.'CFN ').and.
     &         (textpart(ii)(1:4).ne.'CFS')) then
               write(*,*) '*WARNING reading *CONTACT PRINT: label not ap
     &plicable'
               write(*,*) '         or unknown; '
               call inputwarning(inpc,ipoinpc,iline,
     &"*CONTACT PRINT%")
               cycle
            endif
!
            if(textpart(ii)(1:4).eq.'CELS') nener=1
!
            nprint=nprint+1
            if(nprint.gt.nprint_) then
               write(*,*) '*ERROR in contatcprints: increase nprint_'
               ier=1
               return
            endif
!
!           for contact tie sets the tie number is stored,
!           else the surface name
!
            if(textpart(ii)(1:2).eq.'CF') then
               if(itie.eq.0) then
                  write(*,*) '*ERROR reading *CONTACT PRINT: no'
                  write(*,*) '       existing contact pair defined'
                  ier=1
                  return
               endif
               write(prset(nprint)(1:10),'(i10)') itie
               do i=11,81
                  prset(nprint)(i:i)=' '
               enddo
            else
               prset(nprint)=noset
            endif
!            
            prlab(nprint)(1:4)=textpart(ii)(1:4)
            prlab(nprint)(5:5)=total
            prlab(nprint)(6:6)=nodesys
         enddo
      enddo
!
      return
      end

