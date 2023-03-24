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
      subroutine mafillsmse(co,kon,ipkon,lakon,ne,ipompc,nodempc,
     &  coefmpc,nmpc,nelemload,sideload,xload,nload,xbody,ipobody,
     &  nbody,cgr,nactdof,neq,nmethod,ikmpc,ilmpc,elcon,nelcon,rhcon,
     &  nrhcon,alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
     &  t0,t1,ithermal,iprestr,vold,iperturb,sti,stx,iexpl,plicon,
     &  nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,matname,mi,
     &  ncmat_,mass,stiffness,buckling,rhsi,intscheme,physcon,ttime,
     &  time,istep,iinc,coriolis,ibody,xloadold,reltime,veold,
     &  springarea,nstate_,xstateini,xstate,thicke,integerglob,
     &  doubleglob,tieset,istartset,iendset,ialset,ntie,nasym,
     &  pslavsurf,pmastsurf,mortar,clearini,ielprop,prop,ne0,nea,
     &  neb,distmin,ndesi,nodedesi,df,jqs,irows,dfl,
     &  icoordinate,dxstiff,xdesi,istartelem,ialelem,v,sigma,
     &  ieigenfrequency,set,nset,sigmak)
!
      implicit none
!
      character*8 lakon(*)
      character*20 sideload(*)
      character*80 matname(*)
      character*81 tieset(3,*),set(*)
!
      integer kon(*),ipompc(*),nodempc(3,*),nelemload(2,*),ikmpc(*),
     &  ilmpc(*),mi(*),nstate_,ne0,nasym,nactdof(0:mi(2),*),ialset(*),
     &  ielprop(*),nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),
     &  ntie,ielorien(mi(3),*),integerglob(*),istartset(*),iendset(*),
     &  ipkon(*),intscheme,ipobody(2,*),nbody,ibody(3,*),ne,nmpc,
     &  nload,neq(2),nmethod,ithermal(*),iprestr,iperturb(*),i,j,k,
     &  idist,jj,id,ist,index,jdof1,idof1,node1,kflag,icalccg,ntmat_,
     &  indexe,nope,norien,iexpl,i0,ncmat_,istep,iinc,jqs(*),irows(*),
     &  nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,mortar,nea,
     &  neb,ndesi,nodedesi(*),idesvar,istartelem(*),ialelem(*),
     &  icoordinate,ii,ieigenfrequency,mass(2),stiffness,buckling,rhsi,
     &  stiffonly(2),coriolis,idesloc,nset
!
      real*8 co(3,*),coefmpc(*),xload(2,*),p1(3),p2(3),bodyf(3),
     &  xloadold(2,*),reltime,t0(*),t1(*),vold(0:mi(2),*),
     &  s(60,60),ff(60),sti(6,mi(1),*),sm(60,60),xdesi(3,*),
     &  stx(6,mi(1),*),elcon(0:ncmat_,ntmat_,*),val,sigma,
     &  rhcon(0:1,ntmat_,*),springarea(2,*),alcon(0:6,ntmat_,*),
     &  physcon(*),prop(*),xstate(nstate_,mi(1),*),sigmak,
     &  xstateini(nstate_,mi(1),*),alzero(*),orab(7,*),
     &  xbody(7,*),cgr(4,*),plicon(0:2*npmat_,ntmat_,*),
     &  plkcon(0:2*npmat_,ntmat_,*),xstiff(27,mi(1),*),
     &  veold(0:mi(2),*),om,dtime,ttime,time,thicke(mi(3),*),
     &  doubleglob(*),clearini(3,9,*),pslavsurf(3,*),
     &  pmastsurf(6,*),distmin,dfl(20,60),
     &  df(*),dxstiff(27,mi(1),ne,*),v(0:mi(2),*)
!
      kflag=2
      i0=0
      icalccg=0
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
!        distributed forces (body forces or thermal loads or
!        residual stresses or distributed face loads)
!
         if((nbody.ne.0).or.(ithermal(1).ne.0).or.
     &      (iprestr.ne.0).or.(nload.ne.0)) then
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
            if(ipkon(i).lt.0) cycle
            indexe=ipkon(i)
c     Bernhardi start
            if(lakon(i)(1:5).eq.'C3D8I') then
               nope=11
            elseif(lakon(i)(4:5).eq.'20') then
c     Bernhardi end
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
            elseif((lakon(i)(1:2).eq.'ES').and.(lakon(i)(7:7).ne.'F'))
     &              then
!     
!     spring and contact spring elements (NO dashpot elements
!     = ED... elements)
!     
               nope=ichar(lakon(i)(8:8))-47
!     
!     local contact spring number
!     if friction is involved, the contact spring element
!     matrices are determined in mafillsmas.f
!     
               if(lakon(i)(7:7).eq.'C') then
                  if(nasym.eq.1) cycle
                  if(mortar.eq.1) nope=kon(indexe)
               endif
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
     &                    nrhcon,ntmat_,physcon,i,cgr,bodyf,ielmat,
     &                    ithermal,vold,mi)
                  endif
                  index=ipobody(2,index)
                  if(index.eq.0) exit
               enddo
            endif
!     
            call e_c3d_se(co,kon,lakon(i),p1,p2,om,bodyf,nbody,s,sm,ff,
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
     &           clearini,ielprop,prop,distmin,ndesi,nodedesi,
     &           dfl,icoordinate,dxstiff,ne,xdesi,istartelem,
     &           ialelem,v,sigma,ieigenfrequency,set,nset,sigmak)
!     
            do ii=istartelem(i),istartelem(i+1)-1
               idesvar=ialelem(ii)
               if(idesvar.eq.0) cycle
               idesloc=ii-istartelem(i)+1
!
               do jj=1,3*nope
!     
                  j=(jj-1)/3+1
                  k=jj-3*(j-1)
!     
                  node1=kon(indexe+j)
                  jdof1=nactdof(k,node1)
!     
                  if(jdof1.le.0) then
                     if(nmpc.ne.0) then
                        idof1=jdof1
                        if(idof1.ne.2*(idof1/2)) then
                           id=(-idof1+1)/2
                           ist=ipompc(id)
                           index=nodempc(3,ist)
                           if(index.eq.0) cycle
                           do
                              jdof1=nactdof(nodempc(2,index),
     &                             nodempc(1,index))
                              if(jdof1.gt.0) then
                                 val=-coefmpc(index)*dfl(idesloc,jj)
     &                                /coefmpc(ist)
                                 call add_bo_st(df,jqs,irows,jdof1,
     &                                          idesvar,val)
                              endif
                              index=nodempc(3,index)
                              if(index.eq.0) exit
                           enddo
                        endif
                     endif
                     cycle
                  endif  
                  call add_bo_st(df,jqs,irows,jdof1,idesvar,
     &                 dfl(idesloc,jj))
               enddo
            enddo
         enddo
      endif
!     
      return
      end
      
