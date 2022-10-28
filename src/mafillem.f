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
      subroutine mafillem(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
     &  xboun,nboun,
     &  ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
     &  nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
     &  ad,au,fext,nactdof,icol,jq,irow,neq,nzl,nmethod,
     &  ikmpc,ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,
     &  nrhcon,alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
     &  t0,t1,ithermal,prestr,
     &  iprestr,vold,iperturb,sti,nzs,stx,adb,aub,iexpl,plicon,
     &  nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &  matname,mi,ncmat_,mass,stiffness,buckling,rhsi,intscheme,
     &  physcon,shcon,nshcon,cocon,ncocon,ttime,time,istep,iinc,
     &  coriolis,ibody,xloadold,reltime,veold,springarea,nstate_,
     &  xstateini,xstate,thicke,integerglob,doubleglob,
     &  tieset,istartset,iendset,ialset,ntie,nasym,iactive,h0,
     &  pslavsurf,pmastsurf,mortar,clearini,ielprop,prop,
     &  iponoel,inoel,network)
!
!     filling the stiffness matrix in spare matrix format (sm)
!
!     domain 1: phi-domain (air)
!     domain 2: A,V-domain (body)
!     domain 3: A-domain (air, the union of domain 2 and 3 should
!               be simple connected)
!
      implicit none
!
      logical mass(2),stiffness,buckling,rhsi,stiffonly(2),coriolis
!
      character*8 lakon(*)
      character*20 sideload(*)
      character*80 matname(*)
      character*81 tieset(3,*)
!
      integer kon(*),nodeboun(*),ndirboun(*),ipompc(*),nodempc(3,*),
     &  nodeforc(2,*),ndirforc(*),nelemload(2,*),icol(*),jq(*),ikmpc(*),
     &  ilmpc(*),ikboun(*),ilboun(*),mi(*),nstate_,ne0,nasym,
     &  nactdof(0:mi(2),*),konl(26),irow(*),icolumn,ialset(*),
     &  nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),ntie,
     &  ielorien(mi(3),*),integerglob(*),istartset(*),iendset(*),
     &  ipkon(*),intscheme,ncocon(2,*),nshcon(*),ipobody(2,*),nbody,
     &  ibody(3,*),nk,ne,nboun,nmpc,nforc,nload,neq(2),nzl,nmethod,
     &  ithermal(*),iprestr,iperturb(*),nzs(3),i,j,k,l,m,idist,jj,
     &  ll,id,id1,id2,ist,ist1,ist2,index,jdof1,jdof2,idof1,idof2,
     &  mpc1,mpc2,index1,index2,jdof,node1,node2,kflag,icalccg,
     &  ntmat_,indexe,nope,norien,iexpl,i0,ncmat_,istep,iinc,mortar,
     &  nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,iactive(3),
     &  ielprop(*),iponoel(*),inoel(2,*),network
!
      real*8 co(3,*),xboun(*),coefmpc(*),xforc(*),xload(2,*),p1(3),
     &  p2(3),ad(*),au(*),bodyf(3),fext(*),xloadold(2,*),reltime,
     &  t0(*),t1(*),prestr(6,mi(1),*),vold(0:mi(2),*),s(100,100),
     &  sti(6,mi(1),*),sm(100,100),stx(6,mi(1),*),adb(*),aub(*),
     &  elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),springarea(2,*),
     &  alcon(0:6,ntmat_,*),physcon(*),cocon(0:6,ntmat_,*),ff(100),
     &  xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &  shcon(0:3,ntmat_,*),alzero(*),orab(7,*),xbody(7,*),cgr(4,*),
     &  plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &  xstiff(27,mi(1),*),veold(0:mi(2),*),om,valu2,value,dtime,ttime,
     &  time,thicke(mi(3),*),doubleglob(*),h0(3,*),st(60,60),smt(60,60),
     &  pslavsurf(3,*),pmastsurf(6,*),clearini(3,9,*),prop(*)
!
      kflag=2
      i0=0
      icalccg=0
!
      if(stiffness.and.(.not.mass(1))) then
         stiffonly(1)=.true.
      else
         stiffonly(1)=.false.
      endif
      if(stiffness.and.(.not.mass(2))) then
         stiffonly(2)=.true.
      else
         stiffonly(2)=.false.
      endif
!
!     determining nzl
!
      nzl=0
      do i=neq(2),1,-1
         if(icol(i).gt.0) then
            nzl=i
            exit
         endif
      enddo
!
!     initializing the matrices
!
      do i=1,neq(2)
         ad(i)=0.d0
      enddo
      do i=1,nzs(3)
         au(i)=0.d0
      enddo
!     
      if(rhsi) then
         do i=1,neq(2)
            fext(i)=0.d0
         enddo
      endif
!
      if(mass(1)) then
         do i=1,neq(1)
            adb(i)=0.d0
         enddo
         do i=1,nzs(1)
            aub(i)=0.d0
         enddo
      endif
      if(mass(2)) then
         do i=neq(1)+1,neq(2)
            adb(i)=0.d0
         enddo
         do i=nzs(1)+1,nzs(2)
            aub(i)=0.d0
         enddo
      endif
!
!     electromagnetic force should always be taken into account
!
      if(rhsi) idist=1
!
      if((ithermal(1).le.1).or.(ithermal(1).eq.3)) then
!
!     electromagnetic analysis: loop over all elements
!
      ne0=0
      do i=1,ne
!
        if(ipkon(i).lt.0) cycle
        indexe=ipkon(i)
        if(lakon(i)(4:5).eq.'20') then
           nope=20
        elseif(lakon(i)(4:4).eq.'8') then
           nope=8
        elseif(lakon(i)(4:5).eq.'10') then
           nope=10
        elseif(lakon(i)(4:4).eq.'4') then
           nope=4
        elseif(lakon(i)(4:5).eq.'15') then
           nope=15
        elseif(lakon(i)(4:4).eq.'6') then
           nope=6
        else
           cycle
        endif
!
        do j=1,nope
          konl(j)=kon(indexe+j) 
        enddo
!
        call e_c3d_em(co,konl,lakon(i),s,sm,ff,i,nmethod,
     &          ielmat,ntmat_,t1,ithermal,vold,
     &          idist,matname,mi,mass(1),rhsi,
     &          ncmat_,elcon,nelcon,h0,iactive,
     &          alcon,nalcon,istartset,iendset,ialset)
!
        do jj=1,5*nope
!
          j=(jj-1)/5+1
          k=jj-5*(j-1)
!
          node1=kon(indexe+j)
          jdof1=nactdof(k,node1)
!
          do ll=jj,5*nope
!
            l=(ll-1)/5+1
            m=ll-5*(l-1)
!
            node2=kon(indexe+l)
            jdof2=nactdof(m,node2)
!
!           check whether one of the DOF belongs to a SPC or MPC
!
            if((jdof1.gt.0).and.(jdof2.gt.0)) then
               if(stiffonly(1)) then
                  call add_sm_st(au,ad,jq,irow,jdof1,jdof2,
     &                 s(jj,ll),jj,ll)
               else
                  call add_sm_ei(au,ad,aub,adb,jq,irow,jdof1,jdof2,
     &                 s(jj,ll),sm(jj,ll),jj,ll)
               endif
            elseif((jdof1.gt.0).or.(jdof2.gt.0)) then
!
!              idof1: genuine DOF
!              idof2: nominal DOF of the SPC/MPC
!
               if(jdof1.le.0) then
                  idof1=jdof2
c                  idof2=(node1-1)*8+k
                  idof2=jdof1
               else
                  idof1=jdof1
c                  idof2=(node2-1)*8+m
                  idof2=jdof2
               endif
               if(nmpc.gt.0) then
c                  call nident(ikmpc,idof2,nmpc,id)
c                  if((id.gt.0).and.(ikmpc(id).eq.idof2)) then
                  if(idof2.ne.2*(idof2/2)) then
!
!                    regular DOF / MPC
!
c                     id=ilmpc(id)
                     id=(-idof2+1)/2
                     ist=ipompc(id)
                     index=nodempc(3,ist)
                     if(index.eq.0) cycle
                     do
                        idof2=nactdof(nodempc(2,index),nodempc(1,index))
                        value=-coefmpc(index)*s(jj,ll)/coefmpc(ist)
                        if(idof1.eq.idof2) value=2.d0*value
                        if(idof2.gt.0) then
                           if(stiffonly(1)) then
                              call add_sm_st(au,ad,jq,irow,idof1,
     &                             idof2,value,i0,i0)
                           else
                              valu2=-coefmpc(index)*sm(jj,ll)/
     &                               coefmpc(ist)
c
                              if(idof1.eq.idof2) valu2=2.d0*valu2
c
                              call add_sm_ei(au,ad,aub,adb,jq,irow,
     &                             idof1,idof2,value,valu2,i0,i0)
                           endif
                        endif
                        index=nodempc(3,index)
                        if(index.eq.0) exit
                     enddo
                     cycle
                  endif
               endif
!
!              regular DOF / SPC
!
               if(rhsi) then
               elseif(nmethod.eq.2) then
                  value=s(jj,ll)
c                  call nident(ikboun,idof2,nboun,id)
c                  icolumn=neq(2)+ilboun(id)
                  icolumn=neq(2)-idof2/2
                  call add_bo_st(au,jq,irow,idof1,icolumn,value)
               endif
            else
c               idof1=(node1-1)*8+k
c               idof2=(node2-1)*8+m
               idof1=jdof1
               idof2=jdof2
               mpc1=0
               mpc2=0
               if(nmpc.gt.0) then
c                  call nident(ikmpc,idof1,nmpc,id1)
c                  if((id1.gt.0).and.(ikmpc(id1).eq.idof1)) mpc1=1
c                  call nident(ikmpc,idof2,nmpc,id2)
c                  if((id2.gt.0).and.(ikmpc(id2).eq.idof2)) mpc2=1
                  if(idof1.ne.2*(idof1/2)) mpc1=1
                  if(idof2.ne.2*(idof2/2)) mpc2=1
               endif
               if((mpc1.eq.1).and.(mpc2.eq.1)) then
c                  id1=ilmpc(id1)
c                  id2=ilmpc(id2)
                  id1=(-idof1+1)/2
                  id2=(-idof2+1)/2
                  if(id1.eq.id2) then
!
!                    MPC id1 / MPC id1
!
                     ist=ipompc(id1)
                     index1=nodempc(3,ist)
                     if(index1.eq.0) cycle
                     do
                        idof1=nactdof(nodempc(2,index1),
     &                                nodempc(1,index1))
                        index2=index1
                        do
                           idof2=nactdof(nodempc(2,index2),
     &                                   nodempc(1,index2))
                           value=coefmpc(index1)*coefmpc(index2)*
     &                          s(jj,ll)/coefmpc(ist)/coefmpc(ist)
                           if((idof1.gt.0).and.(idof2.gt.0)) then
                              if(stiffonly(1)) then
                                 call add_sm_st(au,ad,jq,irow,
     &                             idof1,idof2,value,i0,i0)
                              else
                                 valu2=coefmpc(index1)*coefmpc(index2)*
     &                             sm(jj,ll)/coefmpc(ist)/coefmpc(ist)
                                 call add_sm_ei(au,ad,aub,adb,jq,
     &                             irow,idof1,idof2,value,valu2,i0,i0)
                              endif
                           endif
!
                           index2=nodempc(3,index2)
                           if(index2.eq.0) exit
                        enddo
                        index1=nodempc(3,index1)
                        if(index1.eq.0) exit
                     enddo
                   else
!
!                    MPC id1 / MPC id2
!
                     ist1=ipompc(id1)
                     index1=nodempc(3,ist1)
                     if(index1.eq.0) cycle
                     do
                        idof1=nactdof(nodempc(2,index1),
     &                                nodempc(1,index1))
                        ist2=ipompc(id2)
                        index2=nodempc(3,ist2)
                        if(index2.eq.0) then
                           index1=nodempc(3,index1)
                           if(index1.eq.0) then
                              exit
                           else
                              cycle
                           endif
                        endif
                        do
                           idof2=nactdof(nodempc(2,index2),
     &                                   nodempc(1,index2))
                           value=coefmpc(index1)*coefmpc(index2)*
     &                          s(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
                           if(idof1.eq.idof2) value=2.d0*value
                           if((idof1.gt.0).and.(idof2.gt.0)) then
                              if(stiffonly(1)) then
                                 call add_sm_st(au,ad,jq,irow,
     &                             idof1,idof2,value,i0,i0)
                              else
                                 valu2=coefmpc(index1)*coefmpc(index2)*
     &                             sm(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
c
                                 if(idof1.eq.idof2) valu2=2.d0*valu2
c
                                 call add_sm_ei(au,ad,aub,adb,jq,
     &                             irow,idof1,idof2,value,valu2,i0,i0)
                              endif
                           endif
!
                           index2=nodempc(3,index2)
                           if(index2.eq.0) exit
                        enddo
                        index1=nodempc(3,index1)
                        if(index1.eq.0) exit
                     enddo
                  endif
               endif
            endif
          enddo
!
          if(rhsi) then
!
!            distributed forces
!
             if(idist.ne.0) then
                if(jdof1.le.0) then
                   if(nmpc.ne.0) then
c                      idof1=(node1-1)*8+k
                      idof1=jdof1
c                      call nident(ikmpc,idof1,nmpc,id)
c                      if((id.gt.0).and.(ikmpc(id).eq.idof1)) then
                      if(idof1.ne.2*(idof1/2)) then
c                         id=ilmpc(id)
                         id=(-idof1+1)/2
                         ist=ipompc(id)
                         index=nodempc(3,ist)
                         if(index.eq.0) cycle
                         do
                            jdof1=nactdof(nodempc(2,index),
     &                           nodempc(1,index))
                            if(jdof1.gt.0) then
                               fext(jdof1)=fext(jdof1)
     &                              -coefmpc(index)*ff(jj)
     &                              /coefmpc(ist)
                            endif
                            index=nodempc(3,index)
                            if(index.eq.0) exit
                         enddo
                      endif
                   endif
                   cycle
                endif
                fext(jdof1)=fext(jdof1)+ff(jj)
             endif
          endif
!
        enddo
      enddo
!
      endif
      if(ithermal(1).gt.1) then
!
!     thermal analysis: loop over all elements
!
      do i=1,ne
!
        if(ipkon(i).lt.0) cycle
!
!       only elements belonging to the A-V-domain should be
!       included in the thermal analysis
!
        if(int(elcon(2,1,ielmat(1,i))).ne.2) cycle
        indexe=ipkon(i)
        if(lakon(i)(4:5).eq.'20') then
           nope=20
        elseif(lakon(i)(4:4).eq.'8') then
           nope=8
        elseif(lakon(i)(4:5).eq.'10') then
           nope=10
        elseif(lakon(i)(4:4).eq.'4') then
           nope=4
        elseif(lakon(i)(4:5).eq.'15') then
           nope=15
        elseif(lakon(i)(4:4).eq.'6') then
           nope=6
         elseif((lakon(i)(1:1).eq.'E').and.(lakon(i)(7:7).ne.'A')) then
!
!          advection elements
!
           read(lakon(i)(8:8),'(i1)') nope
           nope=nope+1
        elseif((lakon(i)(1:2).eq.'D ').or.
     &         ((lakon(i)(1:1).eq.'D').and.(network.eq.1))) then
!
!          asymmetrical contribution -> mafillsmas.f
!
           cycle
        else
           cycle
        endif
!
        call e_c3d_th(co,nk,kon,lakon(i),st,smt,
     &  ff,i,nmethod,rhcon,nrhcon,ielmat,ielorien,norien,orab,
     &  ntmat_,t0,t1,ithermal,vold,iperturb,nelemload,
     &  sideload,xload,nload,idist,iexpl,dtime,
     &  matname,mi(1),mass(2),stiffness,buckling,rhsi,intscheme,
     &  physcon,shcon,nshcon,cocon,ncocon,ttime,time,istep,iinc,
     &  xstiff,xloadold,reltime,ipompc,nodempc,coefmpc,nmpc,ikmpc,
     &  ilmpc,springarea,plkcon,nplkcon,npmat_,ncmat_,elcon,nelcon,
     &  lakon,pslavsurf,pmastsurf,mortar,clearini,plicon,nplicon,
     &  ipkon,ielprop,prop,iponoel,inoel,sti,xstateini,xstate,
     &  nstate_,network,ipobody,xbody,ibody)
!
        do jj=1,nope
!
          j=jj
!
          node1=kon(indexe+j)
          jdof1=nactdof(0,node1)
!
          do ll=jj,nope
!
            l=ll
!
            node2=kon(indexe+l)
            jdof2=nactdof(0,node2)
!
!           check whether one of the DOF belongs to a SPC or MPC
!
            if((jdof1.gt.0).and.(jdof2.gt.0)) then
               if(stiffonly(2)) then
                  call add_sm_st(au,ad,jq,irow,jdof1,jdof2,
     &                 st(jj,ll),jj,ll)
               else
                  call add_sm_ei(au,ad,aub,adb,jq,irow,jdof1,jdof2,
     &                 st(jj,ll),smt(jj,ll),jj,ll)
               endif
            elseif((jdof1.gt.0).or.(jdof2.gt.0)) then
!
!              idof1: genuine DOF
!              idof2: nominal DOF of the SPC/MPC
!
               if(jdof1.le.0) then
                  idof1=jdof2
c                  idof2=(node1-1)*8
                  idof2=jdof1
               else
                  idof1=jdof1
c                  idof2=(node2-1)*8
                  idof2=jdof2
               endif
               if(nmpc.gt.0) then
c                  call nident(ikmpc,idof2,nmpc,id)
c                  if((id.gt.0).and.(ikmpc(id).eq.idof2)) then
                  if(idof2.ne.2*(idof2/2)) then
!
!                    regular DOF / MPC
!
c                     id=ilmpc(id)
                     id=(-idof2+1)/2
                     ist=ipompc(id)
                     index=nodempc(3,ist)
                     if(index.eq.0) cycle
                     do
                        idof2=nactdof(nodempc(2,index),nodempc(1,index))
                        value=-coefmpc(index)*st(jj,ll)/coefmpc(ist)
                        if(idof1.eq.idof2) value=2.d0*value
                        if(idof2.gt.0) then
                           if(stiffonly(2)) then
                              call add_sm_st(au,ad,jq,irow,idof1,
     &                             idof2,value,i0,i0)
                           else
                              valu2=-coefmpc(index)*smt(jj,ll)/
     &                               coefmpc(ist)
!
                              if(idof1.eq.idof2) valu2=2.d0*valu2
!
                              call add_sm_ei(au,ad,aub,adb,jq,irow,
     &                             idof1,idof2,value,valu2,i0,i0)
                           endif
                        endif
                        index=nodempc(3,index)
                        if(index.eq.0) exit
                     enddo
                     cycle
                  endif
               endif
!
!              regular DOF / SPC
!
               if(rhsi) then
               elseif(nmethod.eq.2) then
                  value=st(jj,ll)
c                  call nident(ikboun,idof2,nboun,id)
c                  icolumn=neq(2)+ilboun(id)
                  icolumn=neq(2)-idof2/2
                  call add_bo_st(au,jq,irow,idof1,icolumn,value)
               endif
            else
c               idof1=(node1-1)*8
c               idof2=(node2-1)*8
               idof1=jdof1
               idof2=jdof2
               mpc1=0
               mpc2=0
               if(nmpc.gt.0) then
c                  call nident(ikmpc,idof1,nmpc,id1)
c                  if((id1.gt.0).and.(ikmpc(id1).eq.idof1)) mpc1=1
c                  call nident(ikmpc,idof2,nmpc,id2)
c                  if((id2.gt.0).and.(ikmpc(id2).eq.idof2)) mpc2=1
                  if(idof1.ne.2*(idof1/2)) mpc1=1
                  if(idof2.ne.2*(idof2/2)) mpc2=1
               endif
               if((mpc1.eq.1).and.(mpc2.eq.1)) then
c                  id1=ilmpc(id1)
c                  id2=ilmpc(id2)
                  id1=(-idof1+1)/2
                  id2=(-idof2+1)/2
                  if(id1.eq.id2) then
!
!                    MPC id1 / MPC id1
!
                     ist=ipompc(id1)
                     index1=nodempc(3,ist)
                     if(index1.eq.0) cycle
                     do
                        idof1=nactdof(nodempc(2,index1),
     &                                nodempc(1,index1))
                        index2=index1
                        do
                           idof2=nactdof(nodempc(2,index2),
     &                                   nodempc(1,index2))
                           value=coefmpc(index1)*coefmpc(index2)*
     &                          st(jj,ll)/coefmpc(ist)/coefmpc(ist)
                           if((idof1.gt.0).and.(idof2.gt.0)) then
                              if(stiffonly(2)) then
                                 call add_sm_st(au,ad,jq,irow,
     &                             idof1,idof2,value,i0,i0)
                              else
                                 valu2=coefmpc(index1)*coefmpc(index2)*
     &                             smt(jj,ll)/coefmpc(ist)/coefmpc(ist)
                                 call add_sm_ei(au,ad,aub,adb,jq,
     &                             irow,idof1,idof2,value,valu2,i0,i0)
                              endif
                           endif
!
                           index2=nodempc(3,index2)
                           if(index2.eq.0) exit
                        enddo
                        index1=nodempc(3,index1)
                        if(index1.eq.0) exit
                     enddo
                   else
!
!                    MPC id1 / MPC id2
!
                     ist1=ipompc(id1)
                     index1=nodempc(3,ist1)
                     if(index1.eq.0) cycle
                     do
                        idof1=nactdof(nodempc(2,index1),
     &                                nodempc(1,index1))
                        ist2=ipompc(id2)
                        index2=nodempc(3,ist2)
                        if(index2.eq.0) then
                           index1=nodempc(3,index1)
                           if(index1.eq.0) then
                              exit
                           else
                              cycle
                           endif
                        endif
                        do
                           idof2=nactdof(nodempc(2,index2),
     &                                   nodempc(1,index2))
                           value=coefmpc(index1)*coefmpc(index2)*
     &                          st(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
                           if(idof1.eq.idof2) value=2.d0*value
                           if((idof1.gt.0).and.(idof2.gt.0)) then
                              if(stiffonly(2)) then
                                 call add_sm_st(au,ad,jq,irow,
     &                             idof1,idof2,value,i0,i0)
                              else
                                 valu2=coefmpc(index1)*coefmpc(index2)*
     &                            smt(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
!
                                 if(idof1.eq.idof2) valu2=2.d0*valu2
!
                                 call add_sm_ei(au,ad,aub,adb,jq,
     &                             irow,idof1,idof2,value,valu2,i0,i0)
                              endif
                           endif
!
                           index2=nodempc(3,index2)
                           if(index2.eq.0) exit
                        enddo
                        index1=nodempc(3,index1)
                        if(index1.eq.0) exit
                     enddo
                  endif
               endif
            endif
          enddo
!
          if(rhsi) then
!
!            distributed forces
!
             if(idist.ne.0) then
                if(jdof1.le.0) then
                   if(nmpc.ne.0) then
c                      idof1=(node1-1)*8
                      idof1=jdof1
c                      call nident(ikmpc,idof1,nmpc,id)
c                      if((id.gt.0).and.(ikmpc(id).eq.idof1)) then
                      if(idof1.ne.2*(idof1/2)) then
c                         id=ilmpc(id)
                         id=(-idof1+1)/2
                         ist=ipompc(id)
                         index=nodempc(3,ist)
                         if(index.eq.0) cycle
                         do
                            jdof1=nactdof(nodempc(2,index),
     &                           nodempc(1,index))
                            if(jdof1.gt.0) then
                               fext(jdof1)=fext(jdof1)
     &                              -coefmpc(index)*ff(jj)
     &                              /coefmpc(ist)
                            endif
                            index=nodempc(3,index)
                            if(index.eq.0) exit
                         enddo
                      endif
                   endif
                   cycle
                endif
                fext(jdof1)=fext(jdof1)+ff(jj)
             endif
          endif
!
        enddo
      enddo
!
      endif
!
      return
      end
