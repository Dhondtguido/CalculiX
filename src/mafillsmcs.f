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
      subroutine mafillsmcs(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
     &     xboun,nboun,
     &     ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
     &     nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
     &     ad,au,bb,nactdof,icol,jq,irow,neq,nzl,nmethod,
     &     ikmpc,ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,
     &     nrhcon,alcon,nalcon,alzero,ielmat,ielorien,norien,orab,
     &     ntmat_,t0,t1,ithermal,prestr,
     &     iprestr,vold,iperturb,sti,nzs,stx,adb,aub,iexpl,plicon,
     &     nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &     matname,mi,ics,cs,nm,ncmat_,labmpc,mass,stiffness,buckling,
     &     rhsi,intscheme,mcs,coriolis,ibody,xloadold,reltime,ielcs,
     &     veold,springarea,thicke,integerglob,doubleglob,tieset,
     &     istartset,iendset,ialset,ntie,nasym,pslavsurf,pmastsurf,
     &     mortar,clearini,ielprop,prop,ne0,kscale,xstateini,xstate,
     &     nstate_,set,nset)
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
     &     nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),nset,
     &     ielorien(mi(3),*),integerglob(*),ialset(*),ntie,ikmpc(*),
     &     ipkon(*),ics(*),ij,ilength,lprev,ipobody(2,*),nbody,
     &     ibody(3,*),nk,ne,nboun,nmpc,nforc,nload,neq,nzl,nmethod,
     &     ithermal(*),iprestr,iperturb(*),nzs,i,j,k,l,m,idist,jj,
     &     ll,id,id1,id2,ist,ist1,ist2,index,jdof1,jdof2,idof1,idof2,
     &     mpc1,mpc2,index1,index2,node1,node2,kflag,nasym,mortar,
     &     ntmat_,indexe,nope,norien,iexpl,i0,nm,inode,icomplex,
     &     inode1,icomplex1,inode2,icomplex2,ner,ncmat_,intscheme,istep,
     &     iinc,mcs,ielcs(*),nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),
     &     npmat_,islavelinv(1),irowtloc(1),jqtloc(1),mortartrafoflag,
     &     mscalmethod
!     
      real*8 co(3,*),xboun(*),coefmpc(*),xforc(*),xload(2,*),p1(3),
     &     p2(3),ad(*),au(*),bodyf(3),bb(*),xbody(7,*),cgr(4,*),prop(*),
     &     t0(*),t1(*),prestr(6,mi(1),*),vold(0:mi(2),*),s(60,60),
     &     xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),ff(60),
     &     sti(6,mi(1),*),sm(60,60),stx(6,mi(1),*),adb(*),aub(*),
     &     elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),xloadold(2,*),
     &     alcon(0:6,ntmat_,*),cs(17,*),alzero(*),orab(7,*),reltime,
     &     springarea(2,*),plicon(0:2*npmat_,ntmat_,*),smscale(1),
     &     plkcon(0:2*npmat_,ntmat_,*),thicke(mi(3),*),doubleglob(*),
     &     xstiff(27,mi(1),*),pi,theta,ti,tr,veold(0:mi(2),*),om,valu2,
     &     value,dtime,walue,walu2,time,ttime,clearini(3,9,*),
     &     pslavsurf(3,*),pmastsurf(6,*),autloc(1)
!
      mortartrafoflag=0
!     
!     calculating the scaling factors for the cyclic symmetry calculation
!     
      pi=4.d0*datan(1.d0)
!     
      do i=1,mcs
        write(*,*) '*INFO in mafillsmcs: calculating nodal diameter',nm
        write(*,*) '      for',cs(1,i),'sectors'
        write(*,*) '      (cyclic symmetry definition number',i,')'
        write(*,*)
        theta=nm*2.d0*pi/cs(1,i)
        cs(15,i)=dcos(theta)
        cs(16,i)=dsin(theta)
      enddo
!     
      kflag=2
      i0=0
!     
!     determining nzl
!     
      nzl=0
      do i=neq,1,-1
        if(icol(i).gt.0) then
          nzl=i
          exit
        endif
      enddo
!     
!     initializing the matrices
!     
      do i=1,neq
        ad(i)=0.d0
      enddo
      do i=1,nzs
        au(i)=0.d0
      enddo
!     
      do i=1,neq
        adb(i)=0.d0
      enddo
      do i=1,nzs
        aub(i)=0.d0
      enddo
!     
      ner=neq/2
!     
!     loop over all elements
!     
!     initialisation of the error parameter
!     
c     ne0=0
      do i=1,ne
!     
        if(ipkon(i).lt.0) cycle
        indexe=ipkon(i)
c     Bernhardi start
        if(lakon(i)(4:5).eq.'8I') then
          nope=11
c     Bernhardi end
        elseif(lakon(i)(4:5).eq.'20') then
          nope=20
        elseif(lakon(i)(4:4).eq.'2') then
          nope=26
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
        elseif(lakon(i)(1:2).eq.'ES') then
c     begin 16.07.2020             
          nope=ichar(lakon(i)(8:8))-47
c     read(lakon(i)(8:8),'(i1)') nope
c     nope=nope+1
c     end 16.07.2020             
!     
!     local contact spring number
!     
c     write(*,*) 'nope before= ',nope
          if(lakon(i)(7:7).eq.'C') then
            if(nasym.eq.1) cycle
c     begin 16.07.2020             
            if(mortar.eq.1) nope=kon(indexe)
c     end 16.07.2020             
          endif
c     write(*,*) 'nope after= ',nope
        elseif(lakon(i)(1:4).eq.'MASS') then
          nope=1
        else
          cycle
        endif
!     
        om=0.d0
!     
        if((nbody.gt.0).and.(lakon(i)(1:1).ne.'E')) then
!     
!     assigning centrifugal forces
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
            endif
            index=ipobody(2,index)
            if(index.eq.0) exit
          enddo
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
          do ll=jj,3*nope
            if (mcs.gt.1)then
              if(ielcs(i).gt.0) then
                s(jj,ll)=(cs(1,(ielcs(i)+1))/cs(1,1))*s(jj,ll)
                sm(jj,ll)=(cs(1,(ielcs(i)+1))/cs(1,1))*sm(jj,ll)
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
              call add_sm_ei(au,ad,aub,adb,jq,irow,jdof1,jdof2,
     &             s(jj,ll),sm(jj,ll),jj,ll)
              call add_sm_ei(au,ad,aub,adb,jq,irow,jdof1+ner,jdof2+ner,
     &             s(jj,ll),sm(jj,ll),jj,ll)
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
                        call nident(ics(lprev+1),inode,ilength,id)
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
                      value=-coefmpc(index)*s(jj,ll)/coefmpc(ist)
                      valu2=-coefmpc(index)*sm(jj,ll)/
     &                     coefmpc(ist)
                      if(idof1.eq.idof2) then
                        value=2.d0*value
                        valu2=2.d0*valu2
                      endif
                      if(icomplex.eq.0) then
                        call add_sm_ei(au,ad,aub,adb,jq,irow,
     &                       idof1,idof2,value,valu2,i0,i0)
                        call add_sm_ei(au,ad,aub,adb,jq,irow,
     &                       idof1+ner,idof2+ner,value,valu2,i0,i0)
                      else
                        walue=value*cs(15,icomplex)
                        walu2=valu2*cs(15,icomplex)
                        call add_sm_ei(au,ad,aub,adb,jq,irow,
     &                       idof1,idof2,walue,walu2,i0,i0)
                        call add_sm_ei(au,ad,aub,adb,jq,irow,
     &                       idof1+ner,idof2+ner,walue,walu2,i0,i0)
                        if(idof1.ne.idof2) then
                          walue=value*cs(16,icomplex)
                          walu2=valu2*cs(16,icomplex)
                          call add_sm_ei(au,ad,aub,adb,jq,irow,
     &                         idof1,idof2+ner,walue,walu2,i0,i0)
                          walue=-walue
                          walu2=-walu2
                          call add_sm_ei(au,ad,aub,adb,jq,irow,
     &                         idof1+ner,idof2,walue,walu2,i0,i0)
                        endif
                      endif
                    endif
                    index=nodempc(3,index)
                    if(index.eq.0) exit
                  enddo
                  cycle
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
                    index2=index1
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
                        valu2=coefmpc(index1)*coefmpc(index2)*
     &                       sm(jj,ll)/coefmpc(ist)/coefmpc(ist)
                        if((icomplex1.eq.0).and.
     &                       (icomplex2.eq.0)) then
                          call add_sm_ei(au,ad,aub,adb,jq,
     &                         irow,idof1,idof2,value,valu2,i0,i0)
                          call add_sm_ei(au,ad,aub,adb,jq,
     &                         irow,idof1+ner,idof2+ner,value,
     &                         valu2,i0,i0)
                        elseif((icomplex1.ne.0).and.
     &                         (icomplex2.ne.0)) then
                          if(icomplex1.eq.icomplex2) then
                            call add_sm_ei(au,ad,aub,adb,jq,
     &                           irow,idof1,idof2,value,valu2,i0,i0)
                            call add_sm_ei(au,ad,aub,adb,jq,
     &                           irow,idof1+ner,idof2+ner,value,
     &                           valu2,i0,i0)
                          else
                            tr=cs(15,icomplex1)*cs(15,icomplex2)
     &                           +cs(16,icomplex1)*cs(16,icomplex2)
                            walue=value*tr
                            walu2=valu2*tr
                            call add_sm_ei(au,ad,aub,adb,jq,
     &                           irow,idof1,idof2,walue,walu2,i0,i0)
                            call add_sm_ei(au,ad,aub,adb,jq,
     &                           irow,idof1+ner,idof2+ner,walue,
     &                           walu2,i0,i0)
                            ti=cs(15,icomplex1)*cs(16,icomplex2)
     &                           -cs(15,icomplex2)*cs(16,icomplex1)
                            walue=value*ti
                            walu2=valu2*ti
                            call add_sm_ei(au,ad,aub,adb,jq,irow
     &                           ,idof1,idof2+ner,walue,walu2,i0,i0)
                            walue=-walue
                            walu2=-walu2
                            call add_sm_ei(au,ad,aub,adb,jq,irow
     &                           ,idof1+ner,idof2,walue,walu2,i0,i0)
                          endif
                        elseif((icomplex1.eq.0).or.
     &                         (icomplex2.eq.0)) then
                          if(icomplex2.ne.0) then
                            walue=value*cs(15,icomplex2)
                            walu2=valu2*cs(15,icomplex2)
                          else
                            walue=value*cs(15,icomplex1)
                            walu2=valu2*cs(15,icomplex1)
                          endif
                          call add_sm_ei(au,ad,aub,adb,jq,irow,
     &                         idof1,idof2,walue,walu2,i0,i0)
                          call add_sm_ei(au,ad,aub,adb,jq,irow,
     &                         idof1+ner,idof2+ner,walue,walu2,i0,i0)
                          if(icomplex2.ne.0) then
                            walue=value*cs(16,icomplex2)
                            walu2=valu2*cs(16,icomplex2)
                          else
                            walue=-value*cs(16,icomplex1)
                            walu2=-valu2*cs(16,icomplex1)
                          endif
                          call add_sm_ei(au,ad,aub,adb,jq,irow,
     &                         idof1,idof2+ner,walue,walu2,i0,i0)
                          walue=-walue
                          walu2=-walu2
                          call add_sm_ei(au,ad,aub,adb,jq,irow,
     &                         idof1+ner,idof2,walue,walu2,i0,i0)
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
                        valu2=coefmpc(index1)*coefmpc(index2)*
     &                       sm(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
                        if(idof1.eq.idof2) then
                          value=2.d0*value
                          valu2=2.d0*valu2
                        endif
                        if((icomplex1.eq.0).and.
     &                       (icomplex2.eq.0)) then
                          call add_sm_ei(au,ad,aub,adb,jq,
     &                         irow,idof1,idof2,value,valu2,i0,i0)
                          call add_sm_ei(au,ad,aub,adb,jq,
     &                         irow,idof1+ner,idof2+ner,value,
     &                         valu2,i0,i0)
                        elseif((icomplex1.ne.0).and.
     &                         (icomplex2.ne.0)) then
                          if(icomplex1.eq.icomplex2) then
                            call add_sm_ei(au,ad,aub,adb,jq,
     &                           irow,idof1,idof2,value,valu2,i0,i0)
                            call add_sm_ei(au,ad,aub,adb,jq,
     &                           irow,idof1+ner,idof2+ner,value,
     &                           valu2,i0,i0)
                          else
                            tr=cs(15,icomplex1)*cs(15,icomplex2)
     &                           +cs(16,icomplex1)*cs(16,icomplex2)
c     write(*,*) 'tr= ',tr
                            walue=value*tr
                            walu2=valu2*tr
                            call add_sm_ei(au,ad,aub,adb,jq,
     &                           irow,idof1,idof2,walue,walu2,i0,i0)
                            call add_sm_ei(au,ad,aub,adb,jq,
     &                           irow,idof1+ner,idof2+ner,walue,
     &                           walu2,i0,i0)
                            ti=cs(15,icomplex1)*cs(16,icomplex2)
     &                           -cs(15,icomplex2)*cs(16,icomplex1)
                            walue=value*ti
                            walu2=valu2*ti
                            call add_sm_ei(au,ad,aub,adb,jq,irow
     &                           ,idof1,idof2+ner,walue,walu2,i0,i0)
                            walue=-walue
                            walu2=-walu2
                            call add_sm_ei(au,ad,aub,adb,jq,irow
     &                           ,idof1+ner,idof2,walue,walu2,i0,i0)
                          endif
                        elseif((icomplex1.eq.0).or.
     &                         (icomplex2.eq.0)) then
                          if(icomplex2.ne.0) then
                            walue=value*cs(15,icomplex2)
                            walu2=valu2*cs(15,icomplex2)
                          else
                            walue=value*cs(15,icomplex1)
                            walu2=valu2*cs(15,icomplex1)
                          endif
                          call add_sm_ei(au,ad,aub,adb,jq,irow,
     &                         idof1,idof2,walue,walu2,i0,i0)
                          call add_sm_ei(au,ad,aub,adb,jq,irow,
     &                         idof1+ner,idof2+ner,walue,walu2,i0,i0)
                          if(idof1.ne.idof2) then
                            if(icomplex2.ne.0) then
                              walue=value*cs(16,icomplex2)
                              walu2=valu2*cs(16,icomplex2)
                            else
                              walue=-value*cs(16,icomplex1)
                              walu2=-valu2*cs(16,icomplex1)
                            endif
                            call add_sm_ei(au,ad,aub,adb,jq,
     &                           irow,idof1,idof2+ner,walue,
     &                           walu2,i0,i0)
                            walue=-walue
                            walu2=-walu2
                            call add_sm_ei(au,ad,aub,adb,jq,
     &                           irow,idof1+ner,idof2,walue,
     &                           walu2,i0,i0)
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
