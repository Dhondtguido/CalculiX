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
      subroutine mafillsm(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
     &     xboun,nboun,
     &     ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
     &     nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
     &     ad,au,fext,nactdof,icol,jq,irow,neq,nzl,nmethod,
     &     ikmpc,ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,
     &     nrhcon,alcon,nalcon,alzero,ielmat,ielorien,norien,orab,
     &     ntmat_,t0,t1,ithermal,prestr,
     &     iprestr,vold,iperturb,sti,nzs,stx,adb,aub,iexpl,plicon,
     &     nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &     matname,mi,ncmat_,mass,stiffness,buckling,rhsi,intscheme,
     &     physcon,shcon,nshcon,cocon,ncocon,ttime,time,istep,iinc,
     &     coriolis,ibody,xloadold,reltime,veold,springarea,nstate_,
     &     xstateini,xstate,thicke,integerglob,doubleglob,tieset,
     &     istartset,iendset,ialset,ntie,nasym,pslavsurf,pmastsurf,
     &     mortar,clearini,ielprop,prop,ne0,fnext,nea,neb,kscale,
     &     iponoel,inoel,network,smscale,mscalmethod,set,nset,
     &     islavelinv,autloc,irowtloc,jqtloc,mortartrafoflag)
!     
!     filling the stiffness matrix in spare matrix format (sm)
!     
      implicit none
!     
      integer mass(2),stiffness,buckling,rhsi,stiffonly(2),coriolis
!     
      character*8 lakon(*)
      character*20 sideload(*)
      character*80 matname(*)
      character*81 tieset(3,*),set(*)
!     
      integer kon(*),nodeboun(*),ndirboun(*),ipompc(*),nodempc(3,*),
     &     nodeforc(2,*),ndirforc(*),nelemload(2,*),icol(*),jq(*),
     &     ilmpc(*),ikboun(*),ilboun(*),mi(*),nstate_,ne0,nasym,
     &     nactdof(0:mi(2),*),irow(*),icolumn,ialset(*),ielprop(*),
     &     nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),ntie,
     &     ielorien(mi(3),*),integerglob(*),istartset(*),iendset(*),
     &     ipkon(*),intscheme,ncocon(2,*),nshcon(*),ipobody(2,*),nbody,
     &     ibody(3,*),nk,ne,nboun,nmpc,nforc,nload,neq(2),nzl,nmethod,
     &     ithermal(*),iprestr,iperturb(*),nzs(3),i,j,k,l,m,idist,jj,
     &     ll,id,id1,id2,ist,ist1,ist2,index,jdof1,jdof2,idof1,idof2,
     &     mpc1,mpc2,index1,index2,jdof,node1,node2,kflag,icalccg,
     &     ntmat_,indexe,nope,norien,iexpl,i0,ncmat_,istep,iinc,
     &     nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,mortar,
     &     nea,neb,kscale,iponoel(*),inoel(2,*),network,ndof,
     &     nset,islavelinv(*),jqtloc(*),irowtloc(*),ii,jqtloc1(21),
     &     irowtloc1(96),i1,j1,j2,konl(26),mortartrafoflag,ikmpc(*),
     &     mscalmethod
!     
      real*8 co(3,*),xboun(*),coefmpc(*),xforc(*),xload(2,*),p1(3),
     &     p2(3),ad(*),au(*),bodyf(3),fext(*),xloadold(2,*),reltime,
     &     t0(*),t1(*),prestr(6,mi(1),*),vold(0:mi(2),*),s(60,60),
     &     ff(60),fnext(0:mi(2),*),
     &     sti(6,mi(1),*),sm(60,60),stx(6,mi(1),*),adb(*),aub(*),
     &     elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),springarea(2,*),
     &     alcon(0:6,ntmat_,*),physcon(*),cocon(0:6,ntmat_,*),prop(*),
     &     xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &     shcon(0:3,ntmat_,*),alzero(*),orab(7,*),xbody(7,*),cgr(4,*),
     &     plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &     xstiff(27,mi(1),*),veold(0:mi(2),*),om,valu2,value,dtime,
     &     time,thicke(mi(3),*),doubleglob(*),clearini(3,9,*),ttime,
     &     pslavsurf(3,*),pmastsurf(6,*),smscale(*),autloc(*),
     &     autloc1(96)
!     
      kflag=2
      i0=0
      icalccg=0
c     write(*,*) loc(kflag)
c     write(*,*) loc(s)
c     write(*,*) loc(sm)
c     write(*,*) loc(ff)
c     write(*,*) loc(index1)
!     
      if((stiffness.eq.1).and.(mass(1).eq.0).and.(buckling.eq.0)) then
        stiffonly(1)=1
      else
        stiffonly(1)=0
      endif
      if((stiffness.eq.1).and.(mass(2).eq.0).and.(buckling.eq.0)) then
        stiffonly(2)=1
      else
        stiffonly(2)=0
      endif
!     
      if(rhsi.eq.1) then
!     
!     distributed forces (body forces or thermal loads or
!     residual stresses or distributed face loads)
!     
        if((nbody.ne.0).or.(ithermal(1).ne.0).or.
     &       (iprestr.ne.0).or.(nload.ne.0)) then
          idist=1
        else
          idist=0
        endif
!     
      endif
!     
      if((ithermal(1).le.1).or.(ithermal(1).eq.3)) then
!     
!     mechanical analysis: loop over all elements
!     
        do i=nea,neb
!     
          if((ipkon(i).lt.0).or.(lakon(i)(1:1).eq.'F')) cycle
          indexe=ipkon(i)
c     Bernhardi start
          if(lakon(i)(1:5).eq.'C3D8I') then
            nope=11
            ndof=3
          elseif(lakon(i)(4:5).eq.'20') then
c     Bernhardi end
            nope=20
            ndof=3
          elseif(lakon(i)(4:4).eq.'8') then
            nope=8
            ndof=3
          elseif(lakon(i)(4:5).eq.'10') then
            nope=10
            ndof=3
          elseif(lakon(i)(4:4).eq.'4') then
            nope=4
            ndof=3
          elseif(lakon(i)(4:5).eq.'15') then
            nope=15
            ndof=3
          elseif(lakon(i)(4:4).eq.'6') then
            nope=6
            ndof=3
          elseif((lakon(i)(1:2).eq.'ES').and.(lakon(i)(7:7).ne.'F'))
     &           then
!     
!     spring and contact spring elements (NO dashpot elements
!     = ED... elements)
!     
            nope=ichar(lakon(i)(8:8))-47
            ndof=3
!     
!     local contact spring number
!     if friction is involved, the contact spring element
!     matrices are determined in mafillsmas.f
!     
            if(lakon(i)(7:7).eq.'C') then
              if(nasym.eq.1) cycle
              if(mortar.eq.1) nope=kon(indexe)
            endif
          elseif(lakon(i)(1:4).eq.'MASS') then
            nope=1
            ndof=3
          elseif(lakon(i)(1:1).eq.'U') then
            ndof=ichar(lakon(i)(7:7))
            nope=ichar(lakon(i)(8:8))
          else
            cycle
          endif
!     
c     mortar start
          if(mortartrafoflag.eq.1) then
            do j=1,nope
              konl(j)=kon(indexe+j)
            enddo
          endif
c     mortar end
!     
          om=0.d0
!     
          if((nbody.gt.0).and.(lakon(i)(1:1).ne.'E')) then
!     
!     assigning centrifugal forces
!     
            bodyf(1)=0.d0
            bodyf(2)=0.d0
            bodyf(3)=0.d0
!     
            index=i
            do
              j=ipobody(1,index)
              if(j.eq.0) exit
              if(ibody(1,j).eq.1) then
                om=xbody(1,j)
                p1(1)=xbody(2,j)
                p1(2)=xbody(3,j)
                p1(3)=xbody(4,j)
                p2(1)=xbody(5,j)
                p2(2)=xbody(6,j)
                p2(3)=xbody(7,j)
!     
!     assigning gravity forces
!     
              elseif(ibody(1,j).eq.2) then
                bodyf(1)=bodyf(1)+xbody(1,j)*xbody(2,j)
                bodyf(2)=bodyf(2)+xbody(1,j)*xbody(3,j)
                bodyf(3)=bodyf(3)+xbody(1,j)*xbody(4,j)
!     
!     assigning newton gravity forces
!     
              elseif(ibody(1,j).eq.3) then
                call newton(icalccg,ne,ipkon,lakon,kon,t0,co,rhcon,
     &               nrhcon,ntmat_,physcon,i,cgr,bodyf,ielmat,ithermal,
     &               vold,mi)
              endif
              index=ipobody(2,index)
              if(index.eq.0) exit
            enddo
          endif
!     
!     mortar start
!     
!     generate local transformation matrix for current element
!     
          if(mortartrafoflag.eq.1) then
            if(islavelinv(i).gt.0) then
              if((nope.eq.20).or.(nope.eq.10).or.(nope.eq.15)) then
                jqtloc1(1)=1
                ii=1
                do i1=1,nope
                  node1=konl(i1)
                  do j1=jqtloc(node1),jqtloc(node1+1)-1
                    node2=irowtloc(j1)
                    do j2=1,nope
                      if(konl(j2).eq.node2) then
                        autloc1(ii)=autloc(j1)
                        irowtloc1(ii)=j2
                        ii=ii+1
                      endif
                    enddo
                  enddo
                  jqtloc1(i1+1)=ii
                enddo
              else
                jqtloc1(1)=1
                ii=1
                do i1=1,nope
                  jqtloc1(i1+1)=ii
                enddo
              endif
            endif
          endif
!     
!     mortar end
!     
          if(lakon(i)(1:1).ne.'U') then
            call e_c3d(co,kon,lakon(i),p1,p2,om,bodyf,nbody,s,sm,ff,i,
     &           nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &           alzero,ielmat,ielorien,norien,orab,ntmat_,
     &           t0,t1,ithermal,vold,iperturb,nelemload,sideload,xload,
     &           nload,idist,sti,stx,iexpl,plicon,
     &           nplicon,plkcon,nplkcon,xstiff,npmat_,
     &           dtime,matname,mi(1),ncmat_,mass(1),stiffness,buckling,
     &           rhsi,intscheme,ttime,time,istep,iinc,coriolis,xloadold,
     &           reltime,ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,veold,
     &           springarea,nstate_,xstateini,xstate,ne0,ipkon,thicke,
     &           integerglob,doubleglob,tieset,istartset,
     &           iendset,ialset,ntie,nasym,pslavsurf,pmastsurf,mortar,
     &           clearini,ielprop,prop,kscale,smscale(i),mscalmethod,
     &           set,nset,islavelinv,autloc1,irowtloc1,jqtloc1,
     &           mortartrafoflag)
          else
            call e_c3d_u(co,kon,lakon(i),p1,p2,om,bodyf,nbody,s,sm,ff,i,
     &           nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &           alzero,ielmat,ielorien,norien,orab,ntmat_,
     &           t0,t1,ithermal,vold,iperturb,nelemload,sideload,xload,
     &           nload,idist,sti,stx,iexpl,plicon,
     &           nplicon,plkcon,nplkcon,xstiff,npmat_,
     &           dtime,matname,mi(1),ncmat_,mass(1),stiffness,buckling,
     &           rhsi,intscheme,ttime,time,istep,iinc,coriolis,xloadold,
     &           reltime,ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,veold,
     &           ne0,ipkon,thicke,
     &           integerglob,doubleglob,tieset,istartset,
     &           iendset,ialset,ntie,nasym,
     &           ielprop,prop)
          endif
!     
          do jj=1,ndof*nope
!     
            j=(jj-1)/ndof+1
            k=jj-ndof*(j-1)
!     
            node1=kon(indexe+j)
            jdof1=nactdof(k,node1)
!     
            do ll=jj,ndof*nope
!     
              l=(ll-1)/ndof+1
              m=ll-ndof*(l-1)
!     
              node2=kon(indexe+l)
              jdof2=nactdof(m,node2)
!     
!     check whether one of the DOF belongs to a SPC or MPC
!     
              if((jdof1.gt.0).and.(jdof2.gt.0)) then
                if(stiffonly(1).eq.1) then
                  call add_sm_st(au,ad,jq,irow,jdof1,jdof2,
     &                 s(jj,ll),jj,ll)
                else
                  call add_sm_ei(au,ad,aub,adb,jq,irow,jdof1,jdof2,
     &                 s(jj,ll),sm(jj,ll),jj,ll)
                endif
              elseif((jdof1.gt.0).or.(jdof2.gt.0)) then
!     
!     idof1: genuine DOF
!     idof2: nominal DOF of the SPC/MPC
!     
                if(jdof1.le.0) then
                  idof1=jdof2
                  idof2=jdof1
                else
                  idof1=jdof1
                  idof2=jdof2
                endif
                if(nmpc.gt.0) then
                  if(idof2.ne.2*(idof2/2)) then
!     
!     regular DOF / MPC
!     
                    id=(-idof2+1)/2
                    ist=ipompc(id)
                    index=nodempc(3,ist)
                    if(index.eq.0) cycle
                    do
                      idof2=nactdof(nodempc(2,index),nodempc(1,index))
                      value=-coefmpc(index)*s(jj,ll)/coefmpc(ist)
                      if(idof1.eq.idof2) value=2.d0*value
                      if(idof2.gt.0) then
                        if(stiffonly(1).eq.1) then
                          call add_sm_st(au,ad,jq,irow,idof1,
     &                         idof2,value,i0,i0)
                        else
                          valu2=-coefmpc(index)*sm(jj,ll)/
     &                         coefmpc(ist)
!     
                          if(idof1.eq.idof2) valu2=2.d0*valu2
!     
                          call add_sm_ei(au,ad,aub,adb,jq,irow,
     &                         idof1,idof2,value,valu2,i0,i0)
                        endif
                      endif
                      index=nodempc(3,index)
                      if(index.eq.0) exit
                    enddo
                    cycle
                  endif
                endif
!     
!     regular DOF / SPC
!     
                if(rhsi.eq.1) then
                elseif(nmethod.eq.2) then
                  value=s(jj,ll)
                  icolumn=neq(2)-idof2/2
                  call add_bo_st(au,jq,irow,idof1,icolumn,value)
                endif
              else
                idof1=jdof1
                idof2=jdof2
                mpc1=0
                mpc2=0
                if(nmpc.gt.0) then
                  if(idof1.ne.2*(idof1/2)) mpc1=1
                  if(idof2.ne.2*(idof2/2)) mpc2=1
                endif
                if((mpc1.eq.1).and.(mpc2.eq.1)) then
                  id1=(-idof1+1)/2
                  id2=(-idof2+1)/2
                  if(id1.eq.id2) then
!     
!     MPC id1 / MPC id1
!     
                    ist=ipompc(id1)
                    index1=nodempc(3,ist)
                    if(index1.eq.0) cycle
                    do
                      idof1=nactdof(nodempc(2,index1),
     &                     nodempc(1,index1))
                      index2=index1
                      do
                        idof2=nactdof(nodempc(2,index2),
     &                       nodempc(1,index2))
                        value=coefmpc(index1)*coefmpc(index2)*
     &                       s(jj,ll)/coefmpc(ist)/coefmpc(ist)
                        if((idof1.gt.0).and.(idof2.gt.0)) then
                          if(stiffonly(1).eq.1) then
                            call add_sm_st(au,ad,jq,irow,
     &                           idof1,idof2,value,i0,i0)
                          else
                            valu2=coefmpc(index1)*coefmpc(index2)*
     &                           sm(jj,ll)/coefmpc(ist)/coefmpc(ist)
                            call add_sm_ei(au,ad,aub,adb,jq,
     &                           irow,idof1,idof2,value,valu2,i0,i0)
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
!     MPC id1 / MPC id2
!     
                    ist1=ipompc(id1)
                    index1=nodempc(3,ist1)
                    if(index1.eq.0) cycle
                    do
                      idof1=nactdof(nodempc(2,index1),
     &                     nodempc(1,index1))
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
     &                       nodempc(1,index2))
                        value=coefmpc(index1)*coefmpc(index2)*
     &                       s(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
                        if(idof1.eq.idof2) value=2.d0*value
                        if((idof1.gt.0).and.(idof2.gt.0)) then
                          if(stiffonly(1).eq.1) then
                            call add_sm_st(au,ad,jq,irow,
     &                           idof1,idof2,value,i0,i0)
                          else
                            valu2=coefmpc(index1)*coefmpc(index2)*
     &                           sm(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
!     
                            if(idof1.eq.idof2) valu2=2.d0*valu2
!     
                            call add_sm_ei(au,ad,aub,adb,jq,
     &                           irow,idof1,idof2,value,valu2,i0,i0)
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
            if(rhsi.eq.1) then
!     
!     distributed forces
!     
              if(idist.ne.0) then
!     
!     updating the external force vector for dynamic
!     calculations
!     
                if(nmethod.eq.4) fnext(k,node1)=fnext(k,node1)+ff(jj)
!     
                if(jdof1.le.0) then
                  if(nmpc.ne.0) then
                    if(jdof1.ne.2*(jdof1/2)) then
                      id=(-jdof1+1)/2
                      ist=ipompc(id)
                      index=nodempc(3,ist)
                      if(index.eq.0) cycle
                      do
                        jdof1=nactdof(nodempc(2,index),
     &                       nodempc(1,index))
                        if(jdof1.gt.0) then
                          fext(jdof1)=fext(jdof1)
     &                         -coefmpc(index)*ff(jj)
     &                         /coefmpc(ist)
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
        do i=nea,neb
!     
          if((ipkon(i).lt.0).or.(lakon(i)(1:1).eq.'F')) cycle
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
!     contact spring and advection elements
!     
            nope=ichar(lakon(i)(8:8))-47
!     
!     local contact spring number
!     
            if(lakon(i)(7:7).eq.'C') then
              if(mortar.eq.1) nope=kon(indexe)
            endif
          elseif((lakon(i)(1:2).eq.'D ').or.
     &           ((lakon(i)(1:1).eq.'D').and.(network.eq.1))) then
!     
!     asymmetrical contribution -> mafillsmas.f
!     
            cycle
          else
            cycle
          endif
!     
          call e_c3d_th(co,nk,kon,lakon(i),s,sm,
     &         ff,i,nmethod,rhcon,nrhcon,ielmat,ielorien,norien,orab,
     &         ntmat_,t0,t1,ithermal,vold,iperturb,nelemload,
     &         sideload,xload,nload,idist,iexpl,dtime,matname,mi(1),
     &         mass(2),stiffness,buckling,rhsi,intscheme,physcon,shcon,
     &         nshcon,cocon,ncocon,ttime,time,istep,iinc,xstiff,
     &         xloadold,reltime,ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,
     &         springarea,plkcon,nplkcon,npmat_,ncmat_,elcon,nelcon,
     &         lakon,pslavsurf,pmastsurf,mortar,clearini,plicon,nplicon,
     &         ipkon,ielprop,prop,iponoel,inoel,sti,xstateini,xstate,
     &         nstate_,network,ipobody,xbody,ibody)
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
!     check whether one of the DOF belongs to a SPC or MPC
!     
              if((jdof1.gt.0).and.(jdof2.gt.0)) then
                if(stiffonly(2).eq.1) then
                  call add_sm_st(au,ad,jq,irow,jdof1,jdof2,
     &                 s(jj,ll),jj,ll)
                else
                  call add_sm_ei(au,ad,aub,adb,jq,irow,jdof1,jdof2,
     &                 s(jj,ll),sm(jj,ll),jj,ll)
                endif
              elseif((jdof1.gt.0).or.(jdof2.gt.0)) then
!     
!     idof1: genuine DOF
!     idof2: nominal DOF of the SPC/MPC
!     
                if(jdof1.le.0) then
                  idof1=jdof2
                  idof2=jdof1
                else
                  idof1=jdof1
                  idof2=jdof2
                endif
                if(nmpc.gt.0) then
                  if(idof2.ne.2*(idof2/2)) then
!     
!     regular DOF / MPC
!     
                    id=(-idof2+1)/2
                    ist=ipompc(id)
                    index=nodempc(3,ist)
                    if(index.eq.0) cycle
                    do
                      idof2=nactdof(nodempc(2,index),nodempc(1,index))
                      value=-coefmpc(index)*s(jj,ll)/coefmpc(ist)
                      if(idof1.eq.idof2) value=2.d0*value
                      if(idof2.gt.0) then
                        if(stiffonly(2).eq.1) then
                          call add_sm_st(au,ad,jq,irow,idof1,
     &                         idof2,value,i0,i0)
                        else
                          valu2=-coefmpc(index)*sm(jj,ll)/
     &                         coefmpc(ist)
!     
                          if(idof1.eq.idof2) valu2=2.d0*valu2
!     
                          call add_sm_ei(au,ad,aub,adb,jq,irow,
     &                         idof1,idof2,value,valu2,i0,i0)
                        endif
                      endif
                      index=nodempc(3,index)
                      if(index.eq.0) exit
                    enddo
                    cycle
                  endif
                endif
!     
!     regular DOF / SPC
!     
                if(rhsi.eq.1) then
                elseif(nmethod.eq.2) then
                  value=s(jj,ll)
                  icolumn=neq(2)-idof2/2
                  call add_bo_st(au,jq,irow,idof1,icolumn,value)
                endif
              else
                idof1=jdof1
                idof2=jdof2
                mpc1=0
                mpc2=0
                if(nmpc.gt.0) then
                  if(idof1.ne.2*(idof1/2)) mpc1=1
                  if(idof2.ne.2*(idof2/2)) mpc2=1
                endif
                if((mpc1.eq.1).and.(mpc2.eq.1)) then
                  id1=(-idof1+1)/2
                  id2=(-idof2+1)/2
                  if(id1.eq.id2) then
!     
!     MPC id1 / MPC id1
!     
                    ist=ipompc(id1)
                    index1=nodempc(3,ist)
                    if(index1.eq.0) cycle
                    do
                      idof1=nactdof(nodempc(2,index1),
     &                     nodempc(1,index1))
                      index2=index1
                      do
                        idof2=nactdof(nodempc(2,index2),
     &                       nodempc(1,index2))
                        value=coefmpc(index1)*coefmpc(index2)*
     &                       s(jj,ll)/coefmpc(ist)/coefmpc(ist)
                        if((idof1.gt.0).and.(idof2.gt.0)) then
                          if(stiffonly(2).eq.1) then
                            call add_sm_st(au,ad,jq,irow,
     &                           idof1,idof2,value,i0,i0)
                          else
                            valu2=coefmpc(index1)*coefmpc(index2)*
     &                           sm(jj,ll)/coefmpc(ist)/coefmpc(ist)
                            call add_sm_ei(au,ad,aub,adb,jq,
     &                           irow,idof1,idof2,value,valu2,i0,i0)
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
!     MPC id1 / MPC id2
!     
                    ist1=ipompc(id1)
                    index1=nodempc(3,ist1)
                    if(index1.eq.0) cycle
                    do
                      idof1=nactdof(nodempc(2,index1),
     &                     nodempc(1,index1))
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
     &                       nodempc(1,index2))
                        value=coefmpc(index1)*coefmpc(index2)*
     &                       s(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
                        if(idof1.eq.idof2) value=2.d0*value
                        if((idof1.gt.0).and.(idof2.gt.0)) then
                          if(stiffonly(2).eq.1) then
                            call add_sm_st(au,ad,jq,irow,
     &                           idof1,idof2,value,i0,i0)
                          else
                            valu2=coefmpc(index1)*coefmpc(index2)*
     &                           sm(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
!     
                            if(idof1.eq.idof2) valu2=2.d0*valu2
!     
                            call add_sm_ei(au,ad,aub,adb,jq,
     &                           irow,idof1,idof2,value,valu2,i0,i0)
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
            if(rhsi.eq.1) then
!     
!     distributed forces
!     
              if(idist.ne.0) then
                if(jdof1.le.0) then
                  if(nmpc.ne.0) then
                    if(jdof1.ne.2*(jdof1/2)) then
                      id=(-jdof1+1)/2
                      ist=ipompc(id)
                      index=nodempc(3,ist)
                      if(index.eq.0) cycle
                      do
                        jdof1=nactdof(nodempc(2,index),
     &                       nodempc(1,index))
                        if(jdof1.gt.0) then
                          fext(jdof1)=fext(jdof1)
     &                         -coefmpc(index)*ff(jj)
     &                         /coefmpc(ist)
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
