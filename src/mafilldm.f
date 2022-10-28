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
      subroutine mafilldm(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
     &  xboun,nboun,
     &  ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
     &  nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
     &  ad,au,nactdof,icol,jq,irow,neq,nzl,nmethod,
     &  ikmpc,ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,
     &  nrhcon,alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
     &  t0,t1,ithermal,prestr,
     &  iprestr,vold,iperturb,sti,nzs,stx,adb,aub,iexpl,plicon,
     &  nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &  matname,mi,ncmat_,ttime,time,istep,iinc,ibody,
     &  clearini,mortar,springarea,pslavsurf,pmastsurf,reltime,nasym)
!
!     filling the damping matrix in spare matrix format (sm)
!
      implicit none
!
      character*8 lakon(*)
      character*20 sideload(*)
      character*80 matname(*)
!
      integer kon(*),nodeboun(*),ndirboun(*),ipompc(*),nodempc(3,*),
     &  nodeforc(2,*),ndirforc(*),nelemload(2,*),icol(*),jq(*),ikmpc(*),
     &  ilmpc(*),ikboun(*),ilboun(*),mi(*),nactdof(0:mi(2),*),konl(20),
     &  irow(*),nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),
     &  ielorien(mi(3),*),ipkon(*),ipobody(2,*),nbody,ibody(3,*),
     &  nk,ne,nboun,nmpc,nforc,nload,neq(2),nzl,nmethod,icolumn,
     &  ithermal(*),iprestr,iperturb(*),nzs(3),i,j,k,l,m,idist,jj,
     &  ll,id,id1,id2,ist,ist1,ist2,index,jdof1,jdof2,idof1,idof2,
     &  mpc1,mpc2,index1,index2,node1,node2,ne0,igauss,imat,
     &  ntmat_,indexe,nope,norien,iexpl,i0,ncmat_,istep,iinc,
     &  nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,jfaces,
     &  kode,mortar,nasym
!
      real*8 co(3,*),xboun(*),coefmpc(*),xforc(*),xload(2,*),p1(3),
     &  p2(3),ad(*),au(*),bodyf(3),voldl(0:mi(2),26),clearini(3,9,*),
     &  t0(*),t1(*),prestr(6,mi(1),*),vold(0:mi(2),*),s(60,60),
     &  ff(60),
     &  sti(6,mi(1),*),sm(60,60),stx(6,mi(1),*),adb(*),aub(*),
     &  elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),elas(21),
     &  alcon(0:6,ntmat_,*),alzero(*),orab(7,*),xbody(7,*),cgr(4,*),
     &  plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &  xstiff(27,mi(1),*),om,value,dtime,ttime,time,elconloc(ncmat_),
     &  xl(3,26),springarea(2,*),pslavsurf(3,*),reltime,pmastsurf(6,*)
!
      i0=0
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
      if((ithermal(1).le.1).or.(ithermal(1).eq.3)) then
!     
!     mechanical analysis: loop over all elements
!     
         ne0=0
         do i=1,ne
!     
            if(ipkon(i).lt.0) cycle
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
     &           sm,ff,i,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &           alzero,ielmat,ielorien,norien,orab,ntmat_,
     &           t0,t1,ithermal,vold,iperturb,
     &           nelemload,sideload,xload,nload,idist,sti,stx,
     &           iexpl,plicon,nplicon,plkcon,nplkcon,xstiff,npmat_,
     &           dtime,matname,mi(1),ncmat_,ttime,time,istep,iinc,
     &           nmethod)
            elseif((lakon(i)(1:2).eq.'ES').and.
     &             (lakon(i)(7:7).eq.'C')) then
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
!              computation of the coordinates of the local nodes
!
               do k=1,nope
                  konl(k)=kon(indexe+k)
                  do j=1,3
                     xl(j,k)=co(j,konl(k))
                     voldl(j,k)=vold(j,konl(k))
                  enddo
               enddo
!
!              initialisation of s
!
               do k=1,3*nope
                  do j=1,3*nope
                     s(k,j)=0.d0
                  enddo
               enddo
!
               kode=nelcon(1,imat)
!
!              as soon as the first contact element is discovered ne0 is
!              determined and saved
!
               if(ne0.eq.0) ne0=i-1
               if(mortar.eq.0) then
                  call springdamp_n2f(xl,elas,voldl,s,imat,elcon,
     &                 ncmat_,ntmat_,nope,iperturb,
     &                 springarea(1,konl(nope+1)),nmethod,
     &                 mi,reltime,nasym)
               elseif(mortar.eq.1) then
                  jfaces=kon(indexe+nope+2)
                  igauss=kon(indexe+nope+1) 
                  call springdamp_f2f(xl,elas,voldl,s,imat,elcon,
     &                 ncmat_,ntmat_,nope,lakon(i),iperturb,
     &                 springarea(1,igauss),
     &                 nmethod,mi,reltime,nasym,jfaces,igauss,pslavsurf,
     &                 pmastsurf,clearini)
               endif
            else
               cycle
            endif
!
            if(nasym.eq.0) then
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
                     call add_sm_st(au,ad,jq,irow,jdof1,jdof2,
     &                    s(jj,ll),jj,ll)
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
     &                                        nodempc(1,index))
                              value=-coefmpc(index)*s(jj,ll)/
     &                                        coefmpc(ist)
                              if(idof1.eq.idof2) value=2.d0*value
                              if(idof2.gt.0) then
                                 call add_sm_st(au,ad,jq,irow,idof1,
     &                                idof2,value,i0,i0)
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
     &                             nodempc(1,index1))
                              index2=index1
                              do
                                 idof2=nactdof(nodempc(2,index2),
     &                                nodempc(1,index2))
                                 value=coefmpc(index1)*coefmpc(index2)*
     &                                s(jj,ll)/coefmpc(ist)/coefmpc(ist)
                                 if((idof1.gt.0).and.(idof2.gt.0)) then
                                    call add_sm_st(au,ad,jq,irow,
     &                                   idof1,idof2,value,i0,i0)
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
     &                             nodempc(1,index1))
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
     &                                nodempc(1,index2))
                                 value=coefmpc(index1)*coefmpc(index2)*
     &                                s(jj,ll)/coefmpc(ist1)/
     &                                         coefmpc(ist2)
                                 if(idof1.eq.idof2) value=2.d0*value
                                 if((idof1.gt.0).and.(idof2.gt.0)) then
                                    call add_sm_st(au,ad,jq,irow,
     &                                   idof1,idof2,value,i0,i0)
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
            enddo
!
!           asymmetrical calculations (due to friction or tangential damping...)
!
            else
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
     &                s(jj,ll),jj,ll,nzs)
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
     &                              nodempc(1,index))
                             value=-coefmpc(index)*s(jj,ll)/coefmpc(ist)
                             if(idof1.gt.0) then
                                call add_sm_st_as(au,ad,jq,irow,idof1,
     &                               idof2,value,i0,i0,nzs)
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
     &                             nodempc(1,index))
                             value=-coefmpc(index)*s(jj,ll)/coefmpc(ist)
                             if(idof2.gt.0) then
                                call add_sm_st_as(au,ad,jq,irow,idof1,
     &                               idof2,value,i0,i0,nzs)
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
     &                         nodempc(1,index1))
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
     &                            nodempc(1,index2))
                             value=coefmpc(index1)*coefmpc(index2)*
     &                            s(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
                             if((idof1.gt.0).and.(idof2.gt.0)) then
                                call add_sm_st_as(au,ad,jq,irow,
     &                               idof1,idof2,value,i0,i0,nzs)
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
     &                         nodempc(1,index1))
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
     &                            nodempc(1,index2))
                             value=coefmpc(index1)*coefmpc(index2)*
     &                            s(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
                             if((idof1.gt.0).and.(idof2.gt.0)) then
                                call add_sm_st_as(au,ad,jq,irow,
     &                               idof1,idof2,value,i0,i0,nzs)
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
            endif
         enddo
!     
      endif
!     
      return
      end
      
