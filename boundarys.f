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
      subroutine boundarys(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,nodeboun,ndirboun,xboun,nboun,nboun_,nk,
     &  iamboun,amname,nam,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &  mpcfree,inotr,trab,ntrans,ikboun,ilboun,ikmpc,ilmpc,nk_,
     &  co,labmpc,boun_flag,typeboun,istat,n,iline,ipol,
     &  inl,ipoinp,inp,nam_,namtot_,namta,amta,nmethod,iperturb,
     &  iaxial,ipoinpc,vold,mi,iamplitudedefault,namtot,ier)
!
!     reading the input deck: *BOUNDARY
!
      implicit none
!
      logical boun_flag,user,massflowrate,fixed,submodel
!
      character*1 typeboun(*),type,inpc(*)
      character*20 labmpc(*),label
      character*80 amname(*),amplitude
      character*81 set(*),noset
      character*132 textpart(16)
!
      integer istartset(*),iendset(*),ialset(*),nodeboun(*),
     &  ndirboun(*),nbounnew,kflag,
     &  nset,nboun,nboun_,istat,n,i,j,k,l,ibounstart,ibounend,
     &  key,nk,iamboun(*),nam,iamplitude,ipompc(*),nodempc(3,*),
     &  nmpc,nmpc_,mpcfree,inotr(2,*),ikboun(*),ilboun(*),ikmpc(*),
     &  ilmpc(*),nmpcold,id,idof,index1,ntrans,nk_,ipos,m,node,is,ie,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),nam_,namtot,namtot_,
     &  namta(3,*),idelay,nmethod,iperturb(*),lc,iaxial,ipoinpc(0:*),
     &  ktrue,mi(*),iglobstep,iamplitudedefault,ier
!
      real*8 xboun(*),bounval,coefmpc(*),trab(7,*),co(3,*),amta(2,*),
     &  vold(0:mi(2),*)
!
      type='B'
      label='                    '
      iamplitude=iamplitudedefault
      idelay=0
      user=.false.
      massflowrate=.false.
      fixed=.false.
      lc=1
      iglobstep=0
      submodel=.false.
!
      do i=2,n
         if((textpart(i)(1:6).eq.'OP=NEW').and.(.not.boun_flag)) then
!
!           spc's in nonglobal coordinates result in mpc's
!           removing these mpc's
!           necessary and sufficient condition for a MPC to be removed:
!           - on the dependent side a node "a" corresponding to SPC "b"
!             (no matter which DOF); SPC "b" is applied in direction "c" of
!             node "a" and corresponds to node "d" to account for the
!             inhomogeneous term
!           - on the independent side a term for node "d" in direction "c".
!            
            if(ntrans.gt.0) then
               nmpcold=nmpc
               do j=1,nk
                  if(inotr(2,j).gt.0) then
                     do k=1,3
                        idof=8*(inotr(2,j)-1)+k
                        call nident(ikboun,idof,nboun,id)
                        if(id.gt.0) then
                           if(ikboun(id).eq.idof) then
!
!           if a SPC is defined in direction k for a node j for which a 
!           local coordinate system applies, then the coordinate system
!           number is stored in inotr(1,j) and the additional node
!           for the inhomogeneous term is stored in inotr(2,j). The
!           SPC DOF is (inotr(2,j)-1)*3+k, however, the independent
!           MPC DOF is (j-1)*3+l, where l can be different from k,
!           since (j-1)*3+k might already be taken by another MPC, or
!           the coefficient for this direction might be zero.
!
                              loop: do l=1,3
                                 idof=8*(j-1)+l
                                 call nident(ikmpc,idof,nmpc,id)
                                 if(id.gt.0) then
                                    if(ikmpc(id).eq.idof) then
                                       index1=ipompc(ilmpc(id))
                                       if(index1.eq.0) cycle
                                       do
                                          if((nodempc(1,index1).eq.
     &                                         inotr(2,j)).and.
     &                                       (nodempc(2,index1).eq.k))
     &                                         then
                                             nodempc(3,index1)=mpcfree
                                             mpcfree=ipompc(ilmpc(id))
                                             ipompc(ilmpc(id))=0
                                             do m=id,nmpc-1
                                                ikmpc(m)=ikmpc(m+1)
                                                ilmpc(m)=ilmpc(m+1)
                                             enddo
                                             ikmpc(nmpc)=0
                                             ilmpc(nmpc)=0
                                             nmpc=nmpc-1
                                             exit
                                          endif
                                          index1=nodempc(3,index1)
                                          if(index1.eq.0) exit
                                       enddo
                                    endif
                                 endif
                              enddo loop
!
                           endif
                        endif
                     enddo
                  endif
               enddo
!
!              getting rid of the superfluous lines in ipompc and labmpc
!
               k=0
               do j=1,nmpcold
                  if(ipompc(j).ne.0) then
                     k=k+1
                     ipompc(k)=ipompc(j)
                     labmpc(k)=labmpc(j)
                     index1=ipompc(j)
                     idof=8*(nodempc(1,index1)-1)+nodempc(2,index1)
                     call nident(ikmpc,idof,nmpc,id)
                     if(id.eq.0) then
                        write(*,*) '*ERROR reading *BOUNDARY'
                        ier=1
                        return
                     elseif(ikmpc(id).ne.idof) then
                        write(*,*) '*ERROR reading *BOUNDARY'
                        ier=1
                        return
                     endif
                     ilmpc(id)=k
                  endif
               enddo
            endif
!
            nbounnew=0
            do j=1,nboun
               if(typeboun(j).ne.'B') then
                  nbounnew=nbounnew+1
                  nodeboun(nbounnew)=nodeboun(j)
                  ndirboun(nbounnew)=ndirboun(j)
                  xboun(nbounnew)=xboun(j)
                  typeboun(nbounnew)=typeboun(j)
                  if(nam.gt.0) iamboun(nbounnew)=iamboun(j)
                  ikboun(nbounnew)=8*(nodeboun(j)-1)+ndirboun(j)
                  ilboun(nbounnew)=nbounnew
               endif
            enddo
            nboun=nbounnew
            kflag=2
            call isortii(ikboun,ilboun,nboun,kflag)
         elseif(textpart(i)(1:6).eq.'OP=NEW') then
               write(*,*)
     &           '*INFO reading *BOUNDARY: OP=NEW only makes sense'
               write(*,*)
     &              '      on the first *BOUNDARY card within a step'
               write(*,*)
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
     &           '*ERROR reading *BOUNDARY: nonexistent amplitude'
               write(*,*) '  '
               call inputerror(inpc,ipoinpc,iline,
     &              "*BOUNDARY%",ier)
               return
            endif
            iamplitude=j
         elseif(textpart(i)(1:10).eq.'TIMEDELAY=') THEN
            if(idelay.ne.0) then
               write(*,*) '*ERROR reading *BOUNDARY: the parameter TIME'
               write(*,*) '       DELAY is used twice in the same'
               write(*,*) '       keyword; '
               call inputerror(inpc,ipoinpc,iline,
     &              "*BOUNDARY%",ier)
               return
            else
               idelay=1
            endif
            nam=nam+1
            if(nam.gt.nam_) then
               write(*,*) '*ERROR reading *BOUNDARY: increase nam_'
               ier=1
               return
            endif
            amname(nam)='
     &                                 '
            if(iamplitude.eq.0) then
               write(*,*) '*ERROR reading *BOUNDARY: time delay must be'
               write(*,*) '       preceded by the amplitude parameter'
               ier=1
               return
            endif
            namta(3,nam)=sign(iamplitude,namta(3,iamplitude))
            iamplitude=nam
            namtot=namtot+1
            if(namtot.gt.namtot_) then
               write(*,*) '*ERROR boundaries: increase namtot_'
               ier=1
               return
            endif
            namta(1,nam)=namtot
            namta(2,nam)=namtot
            read(textpart(i)(11:30),'(f20.0)',iostat=istat) 
     &           amta(1,namtot)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*BOUNDARY%",ier)
               return
            endif
         elseif(textpart(i)(1:9).eq.'LOADCASE=') then
            read(textpart(i)(10:19),'(i10)',iostat=istat) lc
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*BOUNDARY%",ier)
               return
            endif
            if(nmethod.ne.5) then
               write(*,*) '*ERROR reading *BOUNDARY: the parameter LOAD'
               write(*,*) '       CASE is only allowed in STEADY STATE'
               write(*,*) '       DYNAMICS calculations'
               ier=1
               return
            endif
         elseif(textpart(i)(1:4).eq.'USER') then
            user=.true.
         elseif(textpart(i)(1:8).eq.'MASSFLOW') then
            massflowrate=.true.
         elseif(textpart(i)(1:5).eq.'FIXED') then
            fixed=.true.
         elseif(textpart(i)(1:8).eq.'SUBMODEL') then
            submodel=.true.
         elseif(textpart(i)(1:5).eq.'STEP=') then
            read(textpart(i)(6:15),'(i10)',iostat=istat) iglobstep
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*BOUNDARY%",ier)
               return
            endif
         elseif(textpart(i)(1:8).eq.'DATASET=') then
            read(textpart(i)(9:18),'(i10)',iostat=istat) iglobstep
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*BOUNDARY%",ier)
               return
            endif
!
!           the mode number for submodels
!           is stored as a negative global step
!
            iglobstep=-iglobstep
         else
            write(*,*) 
     &        '*WARNING reading *BOUNDARY: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*BOUNDARY%")
         endif
      enddo
!
!     check whether global step was specified for submodel
!
      if((submodel).and.(iglobstep.eq.0)) then
         write(*,*) '*ERROR reading *BOUNDARY: no global step'
         write(*,*) '       step specified for the submodel'
         call inputerror(inpc,ipoinpc,iline,
     &        "*BOUNDARY%",ier)
         return
      endif
!
!     storing the step for submodels in iamboun
!
      if(submodel) then
         if(iamplitude.ne.0) then
            write(*,*) '*WARNING reading *BOUNDARY:'
            write(*,*) '         no amplitude definition is allowed'
            write(*,*) '         in combination with a submodel'
         endif
         iamplitude=iglobstep
      endif
!
      if(user.and.(iamplitude.ne.0)) then
         write(*,*) '*WARNING reading *BOUNDARY: '
         write(*,*) '         no amplitude definition is allowed'
         write(*,*) '         for boundary conditions defined by a'
         write(*,*) '         user routine'
         iamplitude=0
      endif
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
!
         read(textpart(2)(1:10),'(i10)',iostat=istat) ibounstart
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*BOUNDARY%",ier)
            return
         endif
!     
         if(textpart(3)(1:1).eq.' ') then
            ibounend=ibounstart
         else
            read(textpart(3)(1:10),'(i10)',iostat=istat) ibounend
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*BOUNDARY%",ier)
               return
            endif
         endif
!     
         if(textpart(4)(1:1).eq.' ') then
            bounval=0.d0
         else
            read(textpart(4)(1:20),'(f20.0)',iostat=istat) bounval
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*BOUNDARY%",ier)
            endif
         endif
         if((massflowrate).and.(iaxial.eq.180)) bounval=bounval/iaxial
!     
!        dummy boundary condition consisting of the first primes
!
         if(user) bounval=1.2357111317d0
         if(submodel) bounval=1.9232931374d0
!
         read(textpart(1)(1:10),'(i10)',iostat=istat) l
         if(istat.eq.0) then
            if((l.gt.nk).or.(l.le.0)) then
               write(*,*) '*ERROR reading *BOUNDARY:'
               write(*,*) '       node ',l,' is not defined'
               ier=1
               return
            endif
            if(submodel) then
               if(ntrans.gt.0) then
                  if(inotr(1,l).gt.0) then
                     write(*,*) '*ERROR reading *BOUNDARY: in submodel'
                     write(*,*) '       node',l,' a local coordinate'
                     write(*,*) '       system was defined. This is not'
                     write(*,*) '       allowed'
                     ier=1
                     return
                  endif
               endif
            endif
            ktrue=l
            if(lc.ne.1) l=l+nk
            call bounadd(l,ibounstart,ibounend,bounval,
     &        nodeboun,ndirboun,xboun,nboun,nboun_,
     &        iamboun,iamplitude,nam,ipompc,nodempc,
     &        coefmpc,nmpc,nmpc_,mpcfree,inotr,trab,
     &        ntrans,ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &        type,typeboun,nmethod,iperturb,fixed,vold,ktrue,
     &        mi,label)
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
               write(*,*) '*ERROR reading *BOUNDARY: node set ',noset
               write(*,*) '       has not yet been defined. '
               call inputerror(inpc,ipoinpc,iline,
     &              "*BOUNDARY%",ier)
               return
            endif
            do j=istartset(i),iendset(i)
               if(ialset(j).gt.0) then
                  k=ialset(j)
                  if(submodel) then
                     if(ntrans.gt.0) then
                        if(inotr(1,k).gt.0) then
                           write(*,*) 
     &                       '*ERROR reading *BOUNDARY: in submodel'
                           write(*,*) '       node',k,
     &                       ' a local coordinate'
                           write(*,*) 
     &                       '       system was defined. This is not'
                           write(*,*) '       allowed'
                           ier=1
                           return
                        endif
                     endif
                  endif
                  ktrue=k
                  if(lc.ne.1) k=k+nk
                  call bounadd(k,ibounstart,ibounend,bounval,
     &               nodeboun,ndirboun,xboun,nboun,nboun_,
     &               iamboun,iamplitude,nam,ipompc,nodempc,
     &               coefmpc,nmpc,nmpc_,mpcfree,inotr,trab,
     &               ntrans,ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &               type,typeboun,nmethod,iperturb,fixed,vold,ktrue,
     &               mi,label)
               else
                  k=ialset(j-2)
                  do
                     k=k-ialset(j)
                     if(k.ge.ialset(j-1)) exit
                     if(submodel) then
                        if(ntrans.gt.0) then
                           if(inotr(1,k).gt.0) then
                              write(*,*) 
     &                          '*ERROR reading *BOUNDARY: in submodel'
                              write(*,*) '       node',k,
     &                          ' a local coordinate'
                              write(*,*) 
     &                          '       system was defined. This is not'
                              write(*,*) '       allowed'
                              ier=1
                              return
                           endif
                        endif
                     endif
                     ktrue=k
                     if(lc.ne.1) k=k+nk
                     call bounadd(k,ibounstart,ibounend,bounval,
     &                 nodeboun,ndirboun,xboun,nboun,nboun_,
     &                 iamboun,iamplitude,nam,ipompc,nodempc,
     &                 coefmpc,nmpc,nmpc_,mpcfree,inotr,trab,
     &                 ntrans,ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,
     &                 labmpc,type,typeboun,nmethod,iperturb,fixed,
     &                 vold,ktrue,mi,label)
                  enddo
               endif
            enddo
         endif
      enddo
!
      return
      end

