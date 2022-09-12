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
      subroutine mafillsmas(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
     &     xboun,nboun,
     &     ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
     &     nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
     &     ad,au,bb,nactdof,icol,jq,irow,neq,nzl,nmethod,
     &     ikmpc,ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,
     &     nrhcon,alcon,nalcon,alzero,ielmat,ielorien,norien,orab,
     &     ntmat_,t0,t1,ithermal,prestr,
     &     iprestr,vold,iperturb,sti,nzs,stx,adb,aub,iexpl,plicon,
     &     nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &     matname,mi,ncmat_,mass,stiffness,buckling,rhsi,intscheme,
     &     physcon,shcon,nshcon,cocon,ncocon,ttime,time,istep,iinc,
     &     coriolis,ibody,xloadold,reltime,veold,springarea,nstate_,
     &     xstateini,xstate,thicke,
     &     integerglob,doubleglob,tieset,istartset,iendset,ialset,
     &     ntie,nasym,pslavsurf,pmastsurf,mortar,clearini,ielprop,
     &     prop,ne0,kscale,iponoel,inoel,network,neam,nebm,neat,nebt,
     &     set,nset)
!     
!     filling the stiffness matrix in spare matrix format (sm)
!     asymmetric contributions
!     
      implicit none
!     
      integer mass(2),stiffness,buckling,rhsi,coriolis
!     
      character*8 lakon(*)
      character*20 sideload(*)
      character*80 matname(*)
      character*81 tieset(3,*),set(*)
!     
      integer kon(*),nodeboun(*),ndirboun(*),ipompc(*),nodempc(3,*),
     &     nodeforc(2,*),ndirforc(*),nelemload(2,*),icol(*),jq(*),
     &     ilmpc(*),ikboun(*),ilboun(*),mi(*),integerglob(*),ist,mpc1,
     &     nactdof(0:mi(2),*),irow(*),istartset(*),iendset(*),
     &     nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),index,
     &     ielorien(mi(3),*),ialset(*),ntie,ne0,nstate_,index1,index2,
     &     ipkon(*),intscheme,ncocon(2,*),nshcon(*),ipobody(2,*),nbody,
     &     ibody(3,*),nk,ne,nboun,nmpc,nforc,nload,neq,nzl,nmethod,
     &     ithermal(*),iprestr,iperturb(*),nzs(3),i,j,k,l,m,idist,jj,
     &     ll,jdof1,jdof2,node1,node2,id,i0,id1,id2,idof1,idof2,nasym,
     &     ntmat_,indexe,nope,norien,iexpl,ncmat_,istep,iinc,mpc2,
     &     nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,ist1,ist2,
     &     mortar,ielprop(*),kscale,iponoel(*),inoel(2,*),network,
     &     neam,nebm,neat,nebt,ikmpc(*),nset,islavelinv(1),irowtloc1(1),
     &     jqtloc1(1),mortartrafoflag,mscalmethod
!     
      real*8 co(3,*),xboun(*),coefmpc(*),xforc(*),xload(2,*),p1(3),
     &     p2(3),ad(*),au(*),bodyf(3),bb(*),xloadold(2,*),value,
     &     t0(*),t1(*),prestr(6,mi(1),*),vold(0:mi(2),*),s(60,60),
     &     ff(60),smscale(1),
     &     sti(6,mi(1),*),sm(60,60),stx(6,mi(1),*),adb(*),aub(*),
     &     elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),reltime,
     &     alcon(0:6,ntmat_,*),physcon(*),cocon(0:6,ntmat_,*),
     &     shcon(0:3,ntmat_,*),alzero(*),orab(7,*),xbody(7,*),cgr(4,*),
     &     xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &     springarea(2,*),thicke(mi(3),*),clearini(3,9,*),prop(*),
     &     plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &     xstiff(27,mi(1),*),veold(0:mi(2),*),doubleglob(*),
     &     om,dtime,ttime,time,pslavsurf(3,*),pmastsurf(6,*),
     &     autloc1(1)
!     
      i0=0
      mortartrafoflag=0
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
!     mechanical analysis: asymmetric contributions
!     only contact friction
!     
        do i=neam,nebm
!     
          if(ipkon(i).lt.0) cycle
          if(lakon(i)(1:2).ne.'ES') cycle
          indexe=ipkon(i)
          read(lakon(i)(8:8),'(i1)') nope
          nope=nope+1
!     
!     local contact spring number
!     
          if(lakon(i)(7:7).eq.'C') then
            if(mortar.eq.1) nope=kon(indexe)
          else
            cycle
          endif
!     
          call e_c3d(co,kon,lakon(i),p1,p2,om,bodyf,nbody,s,sm,ff,i,
     &         nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &         alzero,ielmat,ielorien,norien,orab,ntmat_,
     &         t0,t1,ithermal,vold,iperturb,nelemload,sideload,xload,
     &         nload,idist,sti,stx,iexpl,plicon,
     &         nplicon,plkcon,nplkcon,xstiff,npmat_,
     &         dtime,matname,mi(1),ncmat_,mass,stiffness,buckling,rhsi,
     &         intscheme,ttime,time,istep,iinc,coriolis,xloadold,
     &         reltime,ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,veold,
     &         springarea,nstate_,xstateini,xstate,ne0,ipkon,thicke,
     &         integerglob,
     &         doubleglob,tieset,istartset,iendset,ialset,ntie,nasym,
     &         pslavsurf,pmastsurf,mortar,clearini,ielprop,prop,kscale,
     &         smscale(1),mscalmethod,
     &         set,nset,islavelinv,autloc1,irowtloc1,jqtloc1,
     &         mortartrafoflag)
!     
          do jj=1,3*nope
!     
            j=(jj-1)/3+1
            k=jj-3*(j-1)
!     
            node1=kon(indexe+j)
            jdof1=nactdof(k,node1)
!     
            do ll=1,3*nope
!     
              l=(ll-1)/3+1
              m=ll-3*(l-1)
!     
              node2=kon(indexe+l)
              jdof2=nactdof(m,node2)
!     
!     check whether one of the DOF belongs to a SPC or MPC
!     
              if((jdof1.gt.0).and.(jdof2.gt.0)) then
                call add_sm_st_as(au,ad,jq,irow,jdof1,jdof2,
     &               s(jj,ll),jj,ll,nzs)
              elseif((jdof1.gt.0).or.(jdof2.gt.0)) then
!     
!     idof1: genuine DOF
!     idof2: nominal DOF of the SPC/MPC
!     
                if(jdof1.le.0) then
                  idof1=jdof1
                  idof2=jdof2
                  if(nmpc.gt.0) then
                    if(idof1.ne.2*(idof1/2)) then
!     
!     regular DOF / MPC
!     
                      id=(-idof1+1)/2
                      ist=ipompc(id)
                      index=nodempc(3,ist)
                      if(index.eq.0) cycle
                      do
                        idof1=nactdof(nodempc(2,index),
     &                       nodempc(1,index))
                        value=-coefmpc(index)*s(jj,ll)/coefmpc(ist)
                        if(idof1.gt.0) then
                          call add_sm_st_as(au,ad,jq,irow,idof1,
     &                         idof2,value,i0,i0,nzs)
                        endif
                        index=nodempc(3,index)
                        if(index.eq.0) exit
                      enddo
                      cycle
                    endif
                  endif
                else
                  idof1=jdof1
                  idof2=jdof2
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
                        idof2=nactdof(nodempc(2,index),
     &                       nodempc(1,index))
                        value=-coefmpc(index)*s(jj,ll)/coefmpc(ist)
                        if(idof2.gt.0) then
                          call add_sm_st_as(au,ad,jq,irow,idof1,
     &                         idof2,value,i0,i0,nzs)
                        endif
                        index=nodempc(3,index)
                        if(index.eq.0) exit
                      enddo
                      cycle
                    endif
                  endif
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
                        if((idof1.gt.0).and.(idof2.gt.0)) then
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1,idof2,value,i0,i0,nzs)
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
                        if((idof1.gt.0).and.(idof2.gt.0)) then
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1,idof2,value,i0,i0,nzs)
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
          enddo
        enddo
!     
      endif
!     
      if(ithermal(1).gt.1) then
!     
!     thermal analysis: asymmetric contributions
!     
        do i=neat,nebt
!     
          if(ipkon(i).lt.0) cycle
          if((lakon(i)(1:2).ne.'D ').and.
     &         ((lakon(i)(1:1).ne.'D').or.(network.ne.1))) cycle
          indexe=ipkon(i)
!     
!     no entry or exit elements
!     
          if((kon(indexe+1).eq.0).or.(kon(indexe+3).eq.0)) cycle
          nope=3
!     
          call e_c3d_th(co,nk,kon,lakon(i),s,sm,
     &         ff,i,nmethod,rhcon,nrhcon,ielmat,ielorien,norien,orab,
     &         ntmat_,t0,t1,ithermal,vold,iperturb,nelemload,
     &         sideload,xload,nload,idist,iexpl,dtime,
     &         matname,mi(1),mass(2),stiffness,buckling,rhsi,intscheme,
     &         physcon,shcon,nshcon,cocon,ncocon,ttime,time,istep,iinc,
     &         xstiff,xloadold,reltime,ipompc,nodempc,coefmpc,nmpc,
     &         ikmpc,ilmpc,springarea,plkcon,nplkcon,npmat_,ncmat_,
     &         elcon,nelcon,lakon,pslavsurf,pmastsurf,mortar,clearini,
     &         plicon,nplicon,ipkon,ielprop,prop,iponoel,inoel,sti,
     &         xstateini,xstate,nstate_,network,ipobody,xbody,ibody)
!     
          do jj=1,nope
!     
            j=jj
!     
            node1=kon(indexe+j)
            jdof1=nactdof(0,node1)
!     
            do ll=1,nope
!     
              l=ll
!     
              node2=kon(indexe+l)
              jdof2=nactdof(0,node2)
!     
!     check whether one of the DOF belongs to a SPC or MPC
!     
              if((jdof1.gt.0).and.(jdof2.gt.0)) then
                call add_sm_st_as(au,ad,jq,irow,jdof1,jdof2,
     &               s(jj,ll),jj,ll,nzs)
              elseif((jdof1.gt.0).or.(jdof2.gt.0)) then
!     
!     idof1: genuine DOF
!     idof2: nominal DOF of the SPC/MPC
!     
                if(jdof1.le.0) then
                  idof1=jdof1
                  idof2=jdof2
                  if(nmpc.gt.0) then
                    if(idof1.ne.2*(idof1/2)) then
!     
!     regular DOF / MPC
!     
                      id=(-idof1+1)/2
                      ist=ipompc(id)
                      index=nodempc(3,ist)
                      if(index.eq.0) cycle
                      do
                        idof1=nactdof(nodempc(2,index),
     &                       nodempc(1,index))
                        value=-coefmpc(index)*s(jj,ll)/coefmpc(ist)
                        if(idof1.gt.0) then
                          call add_sm_st_as(au,ad,jq,irow,idof1,
     &                         idof2,value,i0,i0,nzs)
                        endif
                        index=nodempc(3,index)
                        if(index.eq.0) exit
                      enddo
                      cycle
                    endif
                  endif
                else
                  idof1=jdof1
                  idof2=jdof2
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
                        idof2=nactdof(nodempc(2,index),
     &                       nodempc(1,index))
                        value=-coefmpc(index)*s(jj,ll)/coefmpc(ist)
                        if(idof2.gt.0) then
                          call add_sm_st_as(au,ad,jq,irow,idof1,
     &                         idof2,value,i0,i0,nzs)
                        endif
                        index=nodempc(3,index)
                        if(index.eq.0) exit
                      enddo
                      cycle
                    endif
                  endif
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
                        if((idof1.gt.0).and.(idof2.gt.0)) then
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1,idof2,value,i0,i0,nzs)
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
                        if((idof1.gt.0).and.(idof2.gt.0)) then
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1,idof2,value,i0,i0,nzs)
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
          enddo
        enddo
!     
      endif
!     
      return
      end
