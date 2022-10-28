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
      subroutine temperatures(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,t0,t1,nk,ithermal,iamt1,amname,nam,inoelfree,nk_,
     &  nmethod,temp_flag,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &  nam_,namtot_,namta,amta,ipoinpc,t1g,iamplitudedefault,
     &  namtot,ier,itempuser,jobnamec,nuel_)
!
!     reading the input deck: *TEMPERATURE
!
!     itempuser(1): flag: 0: temperatures in the input deck
!                         1: user subroutine utemp is used
!                         2: temperatures read from frd-file
!                            (user subroutine may be used depending
!                             on the values in the file)
!
      implicit none
!
      logical temp_flag,user,submodel,utempusesteps
!
      character*1 inpc(*)
      character*80 amname(*),amplitude
      character*81 set(*),noset
      character*132 textpart(16),jobnamec(*)
!
      integer istartset(*),iendset(*),ialset(*),iamt1(*),nmethod,id,
     &  nset,nk,ithermal(*),istep,istat,n,key,i,j,k,l,nam,ipoinpc(0:*),
     &  iamplitude,ipos,inoelfree,nk_,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),nam_,namtot,namtot_,namta(3,*),idelay,iglobstep,
     &  iamplitudedefault,ier,itempuser(*),nuel_
!
      real*8 t0(*),t1(*),temperature,tempgrad1,tempgrad2,amta(2,*),
     &  t1g(2,*)
!
      iamplitude=iamplitudedefault
      idelay=0
      user=.false.
      utempusesteps=.false.
      iglobstep=0
      submodel=.false.
!
!     defaults: 1) temperatures read from the input deck
!               2) if read from file, then from the first step      
!
      itempuser(1)=0
      itempuser(3)=1
!
      if(nmethod.eq.3) then
         write(*,*) '*ERROR reading *TEMPERATURE: temperature'
         write(*,*) '       loading is not allowed in a linear'
         write(*,*) '       buckling step; perform a static'
         write(*,*) '       nonlinear calculation instead'
         ier=1
         return
      endif
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *TEMPERATURE: *TEMPERATURE'
         write(*,*) '  should only be used within a STEP'
         ier=1
         return
      endif
!
      if(ithermal(1).ne.1) then
         write(*,*) '*ERROR reading *TEMPERATURE: a *TEMPERATURE'
         write(*,*) '  card is detected but no thermal'
         write(*,*) '  *INITIAL CONDITIONS are given'
         ier=1
         return
      endif
!
      do i=2,n
         if((textpart(i).eq.'OP=NEW').and.(.not.temp_flag)) then
            do j=1,nk
               t1(j)=t0(j)
            enddo
         elseif(textpart(i)(1:10).eq.'AMPLITUDE=') then
            read(textpart(i)(11:90),'(a80)') amplitude
            do j=nam,1,-1
               if(amname(j).eq.amplitude) then
                  iamplitude=j
                  exit
               endif
            enddo
            if(j.eq.0) then
               write(*,*)
     &           '*ERROR reading *TEMPERATURE: nonexistent amplitude'
               write(*,*) '  '
               call inputerror(inpc,ipoinpc,iline,
     &              "*TEMPERATURE%",ier)
               return
            endif
            iamplitude=j
         elseif(textpart(i)(1:10).eq.'TIMEDELAY=') THEN
            if(idelay.ne.0) then
               write(*,*) 
     &           '*ERROR reading *TEMPERATURE: the parameter TIME'
               write(*,*) '       DELAY is used twice in the same'
               write(*,*) '       keyword; '
               call inputerror(inpc,ipoinpc,iline,
     &              "*TEMPERATURE%",ier)
               return
            else
               idelay=1
            endif
            nam=nam+1
            if(nam.gt.nam_) then
               write(*,*) '*ERROR reading *TEMPERATURE: increase nam_'
               ier=1
               return
            endif
            amname(nam)='
     &                                 '
            if(iamplitude.eq.0) then
               write(*,*) 
     &           '*ERROR reading *TEMPERATURE: time delay must be'
               write(*,*) '       preceded by the amplitude parameter'
               ier=1
               return
            endif
            namta(3,nam)=sign(iamplitude,namta(3,iamplitude))
            iamplitude=nam
            namtot=namtot+1
            if(namtot.gt.namtot_) then
               write(*,*) '*ERROR temperatures: increase namtot_'
               ier=1
               return
            endif
            namta(1,nam)=namtot
            namta(2,nam)=namtot
            read(textpart(i)(11:30),'(f20.0)',iostat=istat) 
     &           amta(1,namtot)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*TEMPERATURE%",ier)
               return
            endif
         elseif(textpart(i)(1:4).eq.'USER') then
            user=.true.
         elseif(textpart(i)(1:4).eq.'FILE') then
            itempuser(1)=2
            jobnamec(6)(1:127)=textpart(i)(6:132)
            jobnamec(6)(128:132)='     '
            loop1: do j=1,127
               if(jobnamec(6)(j:j).eq.'"') then
                  do k=j+1,127
                     if(jobnamec(6)(k:k).eq.'"') then
                        do l=k-1,127
                           jobnamec(6)(l:l)=' '
                           exit loop1
                        enddo
                     endif
                     jobnamec(6)(k-1:k-1)=jobnamec(6)(k:k)
                  enddo
                  jobnamec(6)(127:127)=' '
               endif
            enddo loop1
         user=.true.
         elseif(textpart(i)(1:6).eq.'BSTEP=') then
            read(textpart(i)(7:16),'(i10)',iostat=istat) itempuser(3)
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*TEMPERATURE%",ier)
            utempusesteps=.true.
         elseif(textpart(i)(1:8).eq.'SUBMODEL') then
            submodel=.true.
         elseif(textpart(i)(1:5).eq.'STEP=') then
            read(textpart(i)(6:15),'(i10)',iostat=istat) iglobstep
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*TEMPERATURE%",ier)
               return
            endif
         elseif(textpart(i)(1:8).eq.'DATASET=') then
            read(textpart(i)(9:18),'(i10)',iostat=istat) iglobstep
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*TEMPERATURE%",ier)
               return
            endif
!
!           the mode number for submodels
!           is stored as a negative global step
!
            iglobstep=-iglobstep
         else
            write(*,*) 
     &        '*WARNING reading *TEMPERATURE: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*TEMPERATURE%")
         endif
      enddo
!
      if(utempusesteps.and.(itempuser(1).ne.2)) then
         write(*,*) '*ERROR reading *TEMPERATURE: the BSTEP parameter '
         write(*,*) '       can only be used in combination '
         write(*,*) '       with parameter FILE.'
         ier=1
         return
      endif
!
!     check whether global step was specified for submodel
!
      if((submodel).and.(iglobstep.eq.0)) then
         write(*,*) '*ERROR reading *TEMPERATURE: no global step'
         write(*,*) '       step specified for the submodel'
         call inputerror(inpc,ipoinpc,iline,
     &        "*TEMPERATURE%",ier)
         return
      endif
!
!     storing the step for submodels in iamboun
!
      if(submodel) then
         if(iamplitude.ne.0) then
            write(*,*) '*WARNING reading *TEMPERATURE:'
            write(*,*) '         no amplitude definition is allowed'
            write(*,*) '         in combination with a submodel'
         endif
         iamplitude=iglobstep
      endif
!
      if(user.and.(iamplitude.ne.0)) then
         write(*,*) 
     &     '*WARNING reading *TEMPERATURE: no amplitude definition is'
         write(*,*) '          allowed for temperatures defined by a'
         write(*,*) '          user routine'
         iamplitude=0
      endif
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
         read(textpart(2)(1:20),'(f20.0)',iostat=istat) temperature
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*TEMPERATURE%",ier)
            return
         endif
!
!        dummy temperature consisting of the first primes
!
         if(user) temperature=1.2357111317d0
         if(submodel) temperature=1.9232931374d0
!
         if((inoelfree.ne.0).or.(nuel_.gt.0)) then
            tempgrad1=0.d0
            tempgrad2=0.d0
            if(n.gt.2) then
               read(textpart(3)(1:20),'(f20.0)',iostat=istat) tempgrad1
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*TEMPERATURE%",ier)
                  return
               endif
            endif
            if(n.gt.3) then
               read(textpart(4)(1:20),'(f20.0)',iostat=istat) tempgrad2
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*TEMPERATURE%",ier)
                  return
               endif
            endif
         endif
!            
         read(textpart(1)(1:10),'(i10)',iostat=istat) l
         if(istat.eq.0) then
            if(l.gt.nk) then
               write(*,*) '*WARNING reading *TEMPERATURE: node ',l
               write(*,*) '         exceeds the largest defined ',
     &            'node number'
               cycle
            endif
            t1(l)=temperature
            if(nam.gt.0) iamt1(l)=iamplitude
            if((inoelfree.ne.0).or.(nuel_.gt.0)) then
               t1g(1,l)=tempgrad1
               t1g(2,l)=tempgrad2
            endif
         else
            read(textpart(1)(1:80),'(a80)',iostat=istat) noset
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
c            do i=1,nset
c               if(set(i).eq.noset) exit
c            enddo
            call cident81(set,noset,nset,id)
            i=nset+1
            if(id.gt.0) then
              if(noset.eq.set(id)) then
                i=id
              endif
            endif
            if(i.gt.nset) then
               noset(ipos:ipos)=' '
               write(*,*) '*ERROR reading *TEMPERATURE: node set ',noset
               write(*,*) '       has not yet been defined. '
               call inputerror(inpc,ipoinpc,iline,
     &              "*TEMPERATURE%",ier)
               return
            endif
            do j=istartset(i),iendset(i)
               if(ialset(j).gt.0) then
                  t1(ialset(j))=temperature
                  if(nam.gt.0) iamt1(ialset(j))=iamplitude
                  if((inoelfree.ne.0).or.(nuel_.gt.0)) then
                     t1g(1,ialset(j))=tempgrad1
                     t1g(2,ialset(j))=tempgrad2
                  endif
               else
                  k=ialset(j-2)
                  do
                     k=k-ialset(j)
                     if(k.ge.ialset(j-1)) exit
                     t1(k)=temperature
                     if(nam.gt.0) iamt1(k)=iamplitude
                     if((inoelfree.ne.0).or.(nuel_.gt.0)) then
                        t1g(1,k)=tempgrad1
                        t1g(2,k)=tempgrad2
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
!
      return
      end

