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
      subroutine mafilldmss(co,nk,kon,ipkon,lakon,ne,
     &     ipompc,nodempc,coefmpc,nmpc,
     &     nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
     &     ad,au,nactdof,jq,irow,neq,nmethod,
     &     ikmpc,ilmpc,elcon,nelcon,rhcon,
     &     nrhcon,alcon,nalcon,alzero,ielmat,ielorien,norien,orab,
     &     ntmat_,t0,t1,ithermal,vold,iperturb,sti,stx,iexpl,plicon,
     &     nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &     matname,mi,ncmat_,
     &     physcon,ttime,time,istep,iinc,
     &     ibody,xloadold,reltime,veold,springarea,nstate_,
     &     xstateini,xstate,thicke,integerglob,doubleglob,
     &     tieset,istartset,iendset,ialset,ntie,nasym,pslavsurf,
     &     pmastsurf,mortar,clearini,ielprop,prop,ne0,nea,neb,
     &     freq,ndamp,dacon,set,nset)
!     
!     filling the stiffness matrix in spare matrix format (sm)
!     
      implicit none
!     
      integer mass(2),stiffness(2),buckling,rhsi,coriolis
!     
      character*8 lakon(*)
      character*20 sideload(*)
      character*80 matname(*)
      character*81 tieset(3,*),set(*)
!     
      integer kon(*),ipompc(*),nodempc(3,*),nelemload(2,*),jq(*),
     &     ikmpc(*),ilmpc(*),mi(*),nstate_,ne0,nasym,konl(20),
     &     nactdof(0:mi(2),*),irow(*),icolumn,ialset(*),ielprop(*),one,
     &     nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),ntie,
     &     jfaces,ielorien(mi(3),*),integerglob(*),istartset(*),
     &     iendset(*),ipkon(*),intscheme,ipobody(2,*),nbody,nset,
     &     ibody(3,*),nk,ne,nmpc,nload,neq(2),nmethod,mscalmethod,
     &     ithermal(*),iperturb(*),i,j,k,l,m,idist,jj,
     &     ll,id,id1,id2,ist,ist1,ist2,index,jdof1,jdof2,idof1,idof2,
     &     mpc1,mpc2,index1,index2,node1,node2,kflag,icalccg,ndamp,
     &     ntmat_,indexe,nope,norien,iexpl,i0,ncmat_,istep,iinc,imat,
     &     nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,mortar,kode,
     &     nea,neb,kscale,ndof,ii,igauss,islavelinv(1),irowtloc(1),
     &     jqtloc(1),mortartrafoflag
!     
      real*8 co(3,*),coefmpc(*),xload(2,*),p1(3),smscale(1),
     &     p2(3),ad(*),au(*),bodyf(3),xloadold(2,*),reltime,
     &     t0(*),t1(*),vold(0:mi(2),*),s(60,60),
     &     ff(60),xl(3,26),voldl(0:mi(2),26),
     &     sti(6,mi(1),*),sm(60,60),stx(6,mi(1),*),
     &     elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),springarea(2,*),
     &     alcon(0:6,ntmat_,*),physcon(*),prop(*),
     &     xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &     alzero(*),orab(7,*),xbody(7,*),cgr(4,*),
     &     plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &     xstiff(27,mi(1),*),veold(0:mi(2),*),om,value,dtime,ttime,
     &     time,thicke(mi(3),*),doubleglob(*),clearini(3,9,*),damping,
     &     pslavsurf(3,*),pmastsurf(6,*),freq,elas(21),dacon(*),
     &     autloc(1)
!     
      mortartrafoflag=0
      one=1
!     
      kflag=2
      i0=0
      icalccg=0
!     
      mass(1)=0
      stiffness(1)=1
      buckling=0
      rhsi=0
      intscheme=0
      coriolis=0
      kscale=1
!     
      do i=nea,neb
!     
        if(ipkon(i).lt.0) cycle
!     
        if(ndamp.gt.0) then
          damping=dacon(ielmat(1,i))
          if(mi(3).gt.1) then
            do j=2,mi(3)
              if(ielmat(j,i).gt.0) then
                if(dacon(ielmat(j,i)).ne.damping) then
                  write(*,*) '*ERROR in mafilldmss: element',i
                  write(*,*) '       contains several layers'
                  write(*,*) '       with different damping '
                  write(*,*) '       coefficients:'
                  write(*,*) '       for layer ',one,':',damping
                  write(*,*) '       for layer ',j,':',
     &                 dacon(ielmat(j,i))
                  call exit(201)
                endif
              endif
            enddo
          endif
        else
          damping=0.d0
        endif
        
        if(lakon(i)(1:2).eq.'ED') then
          indexe=ipkon(i)
          read(lakon(i)(8:8),'(i1)') nope
          nope=nope+1
!     
          do j=1,nope
            konl(j)=kon(indexe+j) 
          enddo
!     
          call e_damp(co,nk,konl,lakon(i),p1,p2,om,bodyf,nbody,s,
     &         sm,ff,i,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &         alzero,ielmat,ielorien,norien,orab,ntmat_,
     &         t0,t1,ithermal,vold,iperturb,
     &         nelemload,sideload,xload,nload,idist,sti,stx,
     &         iexpl,plicon,nplicon,plkcon,nplkcon,xstiff,npmat_,
     &         dtime,matname,mi(1),ncmat_,ttime,time,istep,iinc,
     &         nmethod)
          ndof=3
          do jj=1,ndof*nope
            do ii=1,jj
              s(ii,jj)=s(ii,jj)*freq
            enddo
          enddo
        elseif((lakon(i)(1:2).eq.'ES').and.
     &         (lakon(i)(7:7).eq.'C')) then
          indexe=ipkon(i)
          if(mortar.eq.0) then
            read(lakon(i)(8:8),'(i1)') nope
            nope=nope+1
            konl(nope+1)=kon(indexe+nope+1)
          elseif(mortar.eq.1) then
            nope=kon(indexe)
          endif
          imat=ielmat(1,i)
!     
!     computation of the coordinates of the local nodes
!     
          do k=1,nope
            konl(k)=kon(indexe+k)
            do j=1,3
              xl(j,k)=co(j,konl(k))
              voldl(j,k)=vold(j,konl(k))
            enddo
          enddo
!     
!     initialisation of s
!     
          do k=1,3*nope
            do j=1,3*nope
              s(k,j)=0.d0
            enddo
          enddo
!     
          kode=nelcon(1,imat)
!     
!     as soon as the first contact element is discovered ne0 is
!     determined and saved
!     
          if(mortar.eq.0) then
            call springdamp_n2f(xl,elas,voldl,s,imat,elcon,
     &           ncmat_,ntmat_,nope,iperturb,
     &           springarea(1,konl(nope+1)),nmethod,
     &           mi,reltime,nasym)
          elseif(mortar.eq.1) then
            jfaces=kon(indexe+nope+2)
            igauss=kon(indexe+nope+1) 
            call springdamp_f2f(xl,elas,voldl,s,imat,elcon,
     &           ncmat_,ntmat_,nope,lakon(i),iperturb,
     &           springarea(1,igauss),
     &           nmethod,mi,reltime,nasym,jfaces,igauss,pslavsurf,
     &           pmastsurf,clearini)
          endif
          ndof=3
          do jj=1,ndof*nope
            do ii=1,jj
              s(ii,jj)=s(ii,jj)*freq
            enddo
          enddo
        elseif(damping.ne.0.d0) then
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
          elseif((lakon(i)(1:2).eq.'ES').and.
     &           (lakon(i)(7:7).ne.'F')) then
!     
!     spring and contact spring elements (NO dashpot elements
!     = ED... elements)
!     
            nope=ichar(lakon(i)(8:8))-47
            ndof=3
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
     &               nrhcon,ntmat_,physcon,i,cgr,bodyf,ielmat,
     &               ithermal,vold,mi)
              endif
              index=ipobody(2,index)
              if(index.eq.0) exit
            enddo
          endif
!     
          if(lakon(i)(1:1).ne.'U') then
            call e_c3d(co,kon,lakon(i),p1,p2,om,bodyf,nbody,s,sm,ff,
     &           i,nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
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
     &           clearini,ielprop,prop,kscale,smscale(1),mscalmethod,
     &           set,nset,islavelinv,
     &           autloc,irowtloc,jqtloc,mortartrafoflag)
          else
            call e_c3d_u(co,kon,lakon(i),p1,p2,om,bodyf,nbody,s,sm,
     &           ff,i,nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
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
          do jj=1,ndof*nope
            do ii=1,jj
              s(ii,jj)=s(ii,jj)*damping
            enddo
          enddo
        else
          cycle
        endif
!     
        if(nasym.eq.0) then
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
                call add_sm_st(au,ad,jq,irow,jdof1,jdof2,
     &               s(jj,ll),jj,ll)
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
                      idof2=nactdof(nodempc(2,index),
     &                     nodempc(1,index))
                      value=-coefmpc(index)*s(jj,ll)/
     &                     coefmpc(ist)
                      if(idof1.eq.idof2) value=2.d0*value
                      if(idof2.gt.0) then
                        call add_sm_st(au,ad,jq,irow,idof1,
     &                       idof2,value,i0,i0)
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
                          call add_sm_st(au,ad,jq,irow,
     &                         idof1,idof2,value,i0,i0)
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
                          call add_sm_st(au,ad,jq,irow,
     &                         idof1,idof2,value,i0,i0)
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
!     
          enddo
        endif
      enddo
!     
      return
      end
      
