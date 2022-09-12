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
      subroutine mafillsmcsas(co,nk,kon,ipkon,lakon,ne,nodeboun,
     &     ndirboun,xboun,nboun,
     &     ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
     &     nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
     &     ad,au,fext,nactdof,icol,jq,irow,neq,nzl,nmethod,
     &     ikmpc,ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,
     &     nrhcon,alcon,nalcon,alzero,ielmat,ielorien,norien,orab,
     &     ntmat_,t0,t1,ithermal,prestr,
     &     iprestr,vold,iperturb,sti,nzs,stx,adb,aub,iexpl,plicon,
     &     nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &     matname,mi,ics,cs,nm,ncmat_,labmpc,mass,stiffness,buckling,
     &     rhsi,intscheme,mcs,coriolis,ibody,xloadold,reltime,ielcs,
     &     veold,springarea,thicke,integerglob,doubleglob,
     &     tieset,istartset,iendset,ialset,ntie,nasym,nstate_,xstateini,
     &     xstate,pslavsurf,pmastsurf,mortar,clearini,ielprop,prop,ne0,
     &     kscale,set,nset)
!     
!     filling the stiffness matrix in spare matrix format (sm)
!     for cyclic symmetry calculations
!     
      implicit none
!     
      logical mass,stiffness,buckling,rhsi,coriolis
!     
      character*8 lakon(*)
      character*20 labmpc(*),sideload(*)
      character*80 matname(*)
      character*81 tieset(3,*),set(*)
!     
      integer kon(*),nodeboun(*),ndirboun(*),ipompc(*),nodempc(3,*),
     &     nodeforc(2,*),ndirforc(*),nelemload(2,*),icol(*),jq(*),
     &     ilmpc(*),ikboun(*),ilboun(*),mi(*),nstate_,ne0,ielprop(*),
     &     nactdof(0:mi(2),*),irow(*),istartset(*),iendset(*),kscale,
     &     nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),ikmpc(*),
     &     ielorien(mi(3),*),integerglob(*),ialset(*),ntie,
     &     ipkon(*),ics(*),ij,ilength,lprev,ipobody(2,*),nbody,
     &     ibody(3,*),nk,ne,nboun,nmpc,nforc,nload,neq,nzl,nmethod,
     &     ithermal(*),iprestr,iperturb(*),nzs(3),i,j,k,l,m,idist,jj,
     &     ll,id,id1,id2,ist,ist1,ist2,index,jdof1,jdof2,idof1,idof2,
     &     mpc1,mpc2,index1,index2,node1,node2,kflag,nasym,mortar,
     &     ntmat_,indexe,nope,norien,iexpl,i0,nm,inode,icomplex,
     &     inode1,icomplex1,inode2,icomplex2,ner,ncmat_,intscheme,istep,
     &     iinc,mcs,ielcs(*),nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),
     &     npmat_,nset,islavelinv(1),irowtloc(1),jqtloc(1),
     &     mortartrafoflag,mscalmethod
!     
      real*8 co(3,*),xboun(*),coefmpc(*),xforc(*),xload(2,*),p1(3),
     &     p2(3),ad(*),au(*),bodyf(3),fext(*),xbody(7,*),cgr(4,*),
     &     t0(*),t1(*),prestr(6,mi(1),*),vold(0:mi(2),*),s(60,60),
     &     ff(60),smscale(1),
     &     sti(6,mi(1),*),sm(60,60),stx(6,mi(1),*),adb(*),aub(*),
     &     elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),xloadold(2,*),
     &     alcon(0:6,ntmat_,*),cs(17,*),alzero(*),orab(7,*),reltime,
     &     springarea(2,*),plicon(0:2*npmat_,ntmat_,*),prop(*),
     &     plkcon(0:2*npmat_,ntmat_,*),thicke(mi(3),*),doubleglob(*),
     &     xstiff(27,mi(1),*),pi,theta,ti,tr,veold(0:mi(2),*),om,
     &     xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &     value,dtime,walue,time,ttime,clearini(3,9,*),
     &     pslavsurf(3,*),pmastsurf(6,*),autloc(1)
!
      mortartrafoflag=0
!     
      pi=4.d0*datan(1.d0)
!     
      do i=1,mcs
        theta=nm*2.d0*pi/cs(1,i)
        cs(15,i)=dcos(theta)
        cs(16,i)=dsin(theta)
      enddo
!     
      kflag=2
      i0=0
!     
      ner=neq/2
!     
!     storing the symmetric matrix in asymmetric format
!     
      do i=1,nzs(2)
        au(nzs(3)+i)=au(i)
        aub(nzs(3)+i)=aub(i)
      enddo
!     
!     loop over all elements
!     
      do i=1,ne
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
     &       nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &       alzero,ielmat,ielorien,norien,orab,ntmat_,
     &       t0,t1,ithermal,vold,iperturb,nelemload,sideload,xload,
     &       nload,idist,sti,stx,iexpl,plicon,
     &       nplicon,plkcon,nplkcon,xstiff,npmat_,
     &       dtime,matname,mi(1),ncmat_,mass,stiffness,buckling,rhsi,
     &       intscheme,ttime,time,istep,iinc,coriolis,xloadold,
     &       reltime,ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,veold,
     &       springarea,nstate_,xstateini,xstate,ne0,ipkon,thicke,
     &       integerglob,doubleglob,tieset,istartset,
     &       iendset,ialset,ntie,nasym,pslavsurf,pmastsurf,mortar,
     &       clearini,ielprop,prop,kscale,smscale(1),mscalmethod,
     &       set,nset,islavelinv,autloc,
     &       irowtloc,jqtloc,mortartrafoflag)
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
            if (mcs.gt.1)then
              if(ielcs(i).gt.0) then
                s(jj,ll)=(cs(1,(ielcs(i)+1))/cs(1,1))*s(jj,ll)
              endif
            endif
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
     &             s(jj,ll),jj,ll,nzs)
              call add_sm_st_as(au,ad,jq,irow,jdof1+ner,jdof2+ner,
     &             s(jj,ll),jj,ll,nzs)
            elseif((jdof1.gt.0).or.(jdof2.gt.0)) then
!     
!     idof1: genuine DOF
!     idof2: nominal DOF of the SPC/MPC
!     
              if(jdof1.le.0) then
                idof1=jdof1
                idof2=jdof2
!     
                if(nmpc.gt.0) then
                  if(idof1.ne.2*(idof1/2)) then
!     
!     regular DOF / MPC
!     
                    id1=(-idof1+1)/2
                    ist=ipompc(id1)
                    index=nodempc(3,ist)
                    if(index.eq.0) cycle
                    do
                      inode=nodempc(1,index)
                      icomplex=0
                      if(labmpc(id1)(1:6).eq.'CYCLIC') then
                        read(labmpc(id1)(7:20),'(i14)') icomplex
                      elseif(labmpc(id1)(1:9).eq.'SUBCYCLIC') then
                        do ij=1,mcs
                          ilength=int(cs(4,ij))
                          lprev=int(cs(14,ij))
                          call nident(ics(lprev+1),inode,ilength,
     &                         id)
                          if(id.gt.0) then
                            if(ics(lprev+id).eq.inode) then
                              icomplex=ij
                              exit
                            endif
                          endif
                        enddo
                      endif
                      idof1=nactdof(nodempc(2,index),inode)
                      if(idof1.gt.0) then
                        value=-coefmpc(index)*s(jj,ll)/
     &                       coefmpc(ist)
                        if(icomplex.eq.0) then
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1,idof2,value,i0,i0,nzs)
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1+ner,idof2+ner,value,i0,i0,nzs)
                        else
                          walue=value*cs(15,icomplex)
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1,idof2,walue,i0,i0,nzs)
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1+ner,idof2+ner,walue,i0,i0,nzs)
                          if(idof1.ne.idof2) then
                            walue=value*cs(16,icomplex)
                            call add_sm_st_as(au,ad,jq,irow,
     &                           idof1,idof2+ner,walue,i0,i0,nzs)
                            walue=-walue
                            call add_sm_st_as(au,ad,jq,irow,
     &                           idof1+ner,idof2,walue,i0,i0,nzs)
                          endif
                        endif
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
!     
                if(nmpc.gt.0) then
                  if(idof2.ne.2*(idof2/2)) then
!     
!     regular DOF / MPC
!     
                    id1=(-idof2+1)/2
                    ist=ipompc(id1)
                    index=nodempc(3,ist)
                    if(index.eq.0) cycle
                    do
                      inode=nodempc(1,index)
                      icomplex=0
                      if(labmpc(id1)(1:6).eq.'CYCLIC') then
                        read(labmpc(id1)(7:20),'(i14)') icomplex
                      elseif(labmpc(id1)(1:9).eq.'SUBCYCLIC') then
                        do ij=1,mcs
                          ilength=int(cs(4,ij))
                          lprev=int(cs(14,ij))
                          call nident(ics(lprev+1),inode,ilength,
     &                         id)
                          if(id.gt.0) then
                            if(ics(lprev+id).eq.inode) then
                              icomplex=ij
                              exit
                            endif
                          endif
                        enddo
                      endif
                      idof2=nactdof(nodempc(2,index),inode)
                      if(idof2.gt.0) then
                        value=-coefmpc(index)*s(jj,ll)/
     &                       coefmpc(ist)
                        if(icomplex.eq.0) then
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1,idof2,value,i0,i0,nzs)
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1+ner,idof2+ner,value,i0,i0,nzs)
                        else
                          walue=value*cs(15,icomplex)
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1,idof2,walue,i0,i0,nzs)
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1+ner,idof2+ner,walue,i0,i0,nzs)
                          if(idof1.ne.idof2) then
                            walue=value*cs(16,icomplex)
                            call add_sm_st_as(au,ad,jq,irow,
     &                           idof1,idof2+ner,walue,i0,i0,nzs)
                            walue=-walue
                            call add_sm_st_as(au,ad,jq,irow,
     &                           idof1+ner,idof2,walue,i0,i0,nzs)
                          endif
                        endif
                      endif
                      index=nodempc(3,index)
                      if(index.eq.0) exit
                    enddo
                    cycle
                  endif
                endif
              endif
!     
            else
              idof1=jdof1
              idof2=jdof2
!     
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
                    inode1=nodempc(1,index1)
!     
                    icomplex1=0
                    if(labmpc(id1)(1:6).eq.'CYCLIC') then
                      read(labmpc(id1)(7:20),'(i14)') icomplex1
                    elseif(labmpc(id1)(1:9).eq.'SUBCYCLIC') then
                      do ij=1,mcs
                        ilength=int(cs(4,ij))
                        lprev=int(cs(14,ij))
                        call nident(ics(lprev+1),inode1,
     &                       ilength,id)
                        if(id.gt.0) then
                          if(ics(lprev+id).eq.inode1) then
                            icomplex1=ij
                            exit
                          endif
                        endif
                      enddo
                    endif
                    idof1=nactdof(nodempc(2,index1),inode1)
!     
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
!     
                    do
                      inode2=nodempc(1,index2)
                      icomplex2=0
                      if(labmpc(id1)(1:6).eq.'CYCLIC') then
                        read(labmpc(id1)(7:20),'(i14)') icomplex2
                      elseif(labmpc(id1)(1:9).eq.'SUBCYCLIC') then
                        do ij=1,mcs
                          ilength=int(cs(4,ij))
                          lprev=int(cs(14,ij))
                          call nident(ics(lprev+1),inode2,
     &                         ilength,id)
                          if(id.gt.0) then
                            if(ics(lprev+id).eq.inode2) then
                              icomplex2=ij
                              exit
                            endif
                          endif
                        enddo
                      endif
                      idof2=nactdof(nodempc(2,index2),inode2)
                      if((idof1.gt.0).and.(idof2.gt.0)) then
                        value=coefmpc(index1)*coefmpc(index2)*
     &                       s(jj,ll)/coefmpc(ist)/coefmpc(ist)
                        if((icomplex1.eq.0).and.
     &                       (icomplex2.eq.0)) then
                          call add_sm_st_as(au,ad,jq,
     &                         irow,idof1,idof2,value,i0,i0,nzs)
                          call add_sm_st_as(au,ad,jq,
     &                         irow,idof1+ner,idof2+ner,value,
     &                         i0,i0,nzs)
                        elseif((icomplex1.ne.0).and.
     &                         (icomplex2.ne.0)) then
                          if(icomplex1.eq.icomplex2) then
                            call add_sm_st_as(au,ad,jq,
     &                           irow,idof1,idof2,value,i0,i0,nzs)
                            call add_sm_st_as(au,ad,jq,
     &                           irow,idof1+ner,idof2+ner,value,
     &                           i0,i0,nzs)
                          else
                            tr=cs(15,icomplex1)*cs(15,icomplex2)
     &                           +cs(16,icomplex1)*cs(16,icomplex2)
                            walue=value*tr
                            call add_sm_st_as(au,ad,jq,
     &                           irow,idof1,idof2,walue,i0,i0,nzs)
                            call add_sm_st_as(au,ad,jq,
     &                           irow,idof1+ner,idof2+ner,walue,
     &                           i0,i0,nzs)
                            ti=cs(15,icomplex1)*cs(16,icomplex2)
     &                           -cs(15,icomplex2)*cs(16,icomplex1)
                            walue=value*ti
                            call add_sm_st_as(au,ad,jq,irow
     &                           ,idof1,idof2+ner,walue,i0,i0,nzs)
                            walue=-walue
                            call add_sm_st_as(au,ad,jq,irow
     &                           ,idof1+ner,idof2,walue,i0,i0,nzs)
                          endif
                        elseif((icomplex1.eq.0).or.
     &                         (icomplex2.eq.0)) then
                          if(icomplex2.ne.0) then
                            walue=value*cs(15,icomplex2)
                          else
                            walue=value*cs(15,icomplex1)
                          endif
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1,idof2,walue,i0,i0,nzs)
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1+ner,idof2+ner,walue,i0,i0,nzs)
                          if(icomplex2.ne.0) then
                            walue=value*cs(16,icomplex2)
                          else
                            walue=-value*cs(16,icomplex1)
                          endif
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1,idof2+ner,walue,i0,i0,nzs)
                          walue=-walue
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1+ner,idof2,walue,i0,i0,nzs)
                        endif
                      endif
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
                    inode1=nodempc(1,index1)
                    icomplex1=0
                    if(labmpc(id1)(1:6).eq.'CYCLIC') then
                      read(labmpc(id1)(7:20),'(i14)') icomplex1
                    elseif(labmpc(id1)(1:9).eq.'SUBCYCLIC') then
                      do ij=1,mcs
                        ilength=int(cs(4,ij))
                        lprev=int(cs(14,ij))
                        call nident(ics(lprev+1),inode1,
     &                       ilength,id)
                        if(id.gt.0) then
                          if(ics(lprev+id).eq.inode1) then
                            icomplex1=ij
                            exit
                          endif
                        endif
                      enddo
                    endif
                    idof1=nactdof(nodempc(2,index1),inode1)
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
                      inode2=nodempc(1,index2)
                      icomplex2=0
                      if(labmpc(id2)(1:6).eq.'CYCLIC') then
                        read(labmpc(id2)(7:20),'(i14)') icomplex2
                      elseif(labmpc(id2)(1:9).eq.'SUBCYCLIC') then
                        do ij=1,mcs
                          ilength=int(cs(4,ij))
                          lprev=int(cs(14,ij))
                          call nident(ics(lprev+1),inode2,
     &                         ilength,id)
                          if(id.gt.0) then
                            if(ics(lprev+id).eq.inode2) then
                              icomplex2=ij
                              exit
                            endif
                          endif
                        enddo
                      endif
                      idof2=nactdof(nodempc(2,index2),inode2)
                      if((idof1.gt.0).and.(idof2.gt.0)) then
                        value=coefmpc(index1)*coefmpc(index2)*
     &                       s(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
                        if((icomplex1.eq.0).and.
     &                       (icomplex2.eq.0)) then
                          call add_sm_st_as(au,ad,jq,
     &                         irow,idof1,idof2,value,i0,i0,nzs)
                          call add_sm_st_as(au,ad,jq,
     &                         irow,idof1+ner,idof2+ner,value,
     &                         i0,i0,nzs)
                        elseif((icomplex1.ne.0).and.
     &                         (icomplex2.ne.0)) then
                          if(icomplex1.eq.icomplex2) then
                            call add_sm_st_as(au,ad,jq,
     &                           irow,idof1,idof2,value,i0,i0,nzs)
                            call add_sm_st_as(au,ad,jq,
     &                           irow,idof1+ner,idof2+ner,value,
     &                           i0,i0,nzs)
                          else
                            tr=cs(15,icomplex1)*cs(15,icomplex2)
     &                           +cs(16,icomplex1)*cs(16,icomplex2)
                            walue=value*tr
                            call add_sm_st_as(au,ad,jq,
     &                           irow,idof1,idof2,walue,i0,i0,nzs)
                            call add_sm_st_as(au,ad,jq,
     &                           irow,idof1+ner,idof2+ner,walue,
     &                           i0,i0,nzs)
                            ti=cs(15,icomplex1)*cs(16,icomplex2)
     &                           -cs(15,icomplex2)*cs(16,icomplex1)
                            walue=value*ti
                            call add_sm_st_as(au,ad,jq,irow
     &                           ,idof1,idof2+ner,walue,i0,i0,nzs)
                            walue=-walue
                            call add_sm_st_as(au,ad,jq,irow
     &                           ,idof1+ner,idof2,walue,i0,i0,nzs)
                          endif
                        elseif((icomplex1.eq.0).or.
     &                         (icomplex2.eq.0)) then
                          if(icomplex2.ne.0) then
                            walue=value*cs(15,icomplex2)
                          else
                            walue=value*cs(15,icomplex1)
                          endif
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1,idof2,walue,i0,i0,nzs)
                          call add_sm_st_as(au,ad,jq,irow,
     &                         idof1+ner,idof2+ner,walue,i0,i0,nzs)
                          if(idof1.ne.idof2) then
                            if(icomplex2.ne.0) then
                              walue=value*cs(16,icomplex2)
                            else
                              walue=-value*cs(16,icomplex1)
                            endif
                            call add_sm_st_as(au,ad,jq,
     &                           irow,idof1,idof2+ner,walue,
     &                           i0,i0,nzs)
                            walue=-walue
                            call add_sm_st_as(au,ad,jq,
     &                           irow,idof1+ner,idof2,walue,
     &                           i0,i0,nzs)
                          endif
                        endif
                      endif
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
        enddo
      enddo
!     
      return
      end
