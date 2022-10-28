!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine gencontelem_f2f(tieset,ntie,itietri,ne,ipkon,kon,lakon,
     &     cg,straight,ifree,koncont,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,
     &     ielmat,elcon,istep,iinc,iit,ncmat_,ntmat_,mi,imastop,
     &     islavsurf,itiefac,springarea,tietol,reltime,filab,nasym,
     &     pslavsurf,pmastsurf,clearini,theta,xstateini,xstate,nstate_,
     &     ne0,icutb,ialeatoric,nmethod,jobnamef,alea)
!     
!     generate contact elements for the slave contact nodes
!     
      implicit none
!     
      character*8 lakon(*)
      character*33 cfile
      character*81 tieset(3,*),setname
      character*87 filab(*)
      character*132 jobnamef(*)
!     
      integer ntie,ifree,nasym,nstate_,ne0,
     &     itietri(2,ntie),ipkon(*),kon(*),koncont(4,*),ne,nel,
     &     neigh(1),iflag,kneigh,i,j,k,l,jj,nn,isol,icutb,
     &     itri,ll,kflag,n,nx(*),ny(*),istep,iinc,mi(*),
     &     nz(*),nstart,ielmat(mi(3),*),imat,ifaceq(8,6),ifacet(6,4),
     &     ifacew1(4,5),ifacew2(8,5),nelemm,jfacem,indexe,iit,
     &     nface,nope,nodefm(9),ncmat_,ntmat_,number(4),lenset,
     &     iteller,ifaces,jfaces,ifacem,indexel,iloop,iprev,iact,
     &     imastop(3,*), itriangle(100),ntriangle,ntriangle_,itriold,
     &     itrinew,id,islavsurf(2,*),itiefac(2,*),nelems,m,mint2d,nopes,
     &     igauss,nopem,nodefs(9),indexf,ialeatoric,nmethod
!     
      real*8 cg(3,*),straight(16,*),co(3,*),vold(0:mi(2),*),p(3),
     &     dist,xo(*),yo(*),zo(*),x(*),y(*),z(*),clearini(3,9,*),
     &     elcon(0:ncmat_,ntmat_,*),weight,theta,harvest,alea,
     &     springarea(2,*),xl2(3,9),area,xi,et,shp2(7,9),
     &     xs2(3,2),xsj2(3),tietol(4,*),reltime,xstate(nstate_,mi(1),*),
     &     clear,ratio(9),pl(3,9),xstateini(nstate_,mi(1),*),
     &     pproj(3),al(3),xn(3),xm(3),dm,pslavsurf(3,*),pmastsurf(6,*)
!     
!     nodes per face for hex elements
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
!     
!     nodes per face for tet elements
!     
      data ifacet /1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/
!     
!     nodes per face for linear wedge elements
!     
      data ifacew1 /1,3,2,0,
     &     4,5,6,0,
     &     1,2,5,4,
     &     2,3,6,5,
     &     3,1,4,6/
!     
!     nodes per face for quadratic wedge elements
!     
      data ifacew2 /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     3,1,4,6,9,13,12,15/
!     
!     flag for shape functions
!     
      data iflag /2/
      data indexel /0/
!     
      save indexel
!     
      include "gauss.f"
!     
      if((iinc.eq.1).and.(iit.le.0)) call random_seed()
!     
!     opening a file to store the contact spring elements
!     
      if(filab(1)(3:3).eq.'C') then
        iteller=iteller+1
!     
        do i=1,132
          if(jobnamef(1)(i:i).eq.' ') exit
        enddo
        i=i-1
        cfile=jobnamef(1)(1:i)//'.cel'
        open(27,file=cfile,status='unknown',position='append')
!     
        setname(1:15)='contactelements'
        lenset=15
        number(1)=istep
        number(2)=iinc
        number(3)=icutb+1
        number(4)=iit
        if(number(4).le.0) number(4)=1
        do i=1,4
          setname(lenset+1:lenset+1)='_'
          if(i.eq.1) then
            setname(lenset+2:lenset+3)='st'
          elseif(i.eq.2) then
            setname(lenset+2:lenset+3)='in'
          elseif(i.eq.3) then
            setname(lenset+2:lenset+3)='at'
          else
            setname(lenset+2:lenset+3)='it'
          endif
          if(number(i).lt.10) then
            write(setname(lenset+4:lenset+4),'(i1)') number(i)
            lenset=lenset+4
          elseif(number(i).lt.100) then
            write(setname(lenset+4:lenset+5),'(i2)') number(i)
            lenset=lenset+5
          elseif(number(i).lt.1000) then
            write(setname(lenset+4:lenset+6),'(i3)') number(i)
            lenset=lenset+6
          else
            write(*,*) '*ERROR in gencontelem_f2f: no more than 1000'
            write(*,*) '       steps/increments/cutbacks/iterations'
            write(*,*) '       allowed (for output in '
            write(*,*) '       contactelements.inp)'
            stop
          endif
        enddo
      endif
!     
      igauss=0
!     
!     loop over all active contact ties
!     
      do i=1,ntie
        if(tieset(1,i)(81:81).ne.'C') cycle
        imat=int(tietol(2,i))
!     
!     sorting the centers of gravity of the triangulation
!     of the master surface in order to find the master
!     location corresponding to each slave integration point
!     (only at the start of an increment)
!     
        if(iit.le.0) then
          kneigh=1
          nstart=itietri(1,i)-1
          n=itietri(2,i)-nstart
          if(n.lt.kneigh) kneigh=n
          do j=1,n
            xo(j)=cg(1,nstart+j)
            x(j)=xo(j)
            nx(j)=j
            yo(j)=cg(2,nstart+j)
            y(j)=yo(j)
            ny(j)=j
            zo(j)=cg(3,nstart+j)
            z(j)=zo(j)
            nz(j)=j
          enddo
          kflag=2
          call dsort(x,nx,n,kflag)
          call dsort(y,ny,n,kflag)
          call dsort(z,nz,n,kflag)
        endif
!     
!     loop over all slave faces
!     
!     iprev is the number of contact elements at the end
!     of the previous increment at the start of a
!     new increment, else zero
!     iact is the number of actual contact elements
!     
        do iloop=1,2
          iprev=0
          iact=0
          do jj=itiefac(1,i), itiefac(2,i)
            ifaces=islavsurf(1,jj)
            nelems=int(ifaces/10)
            jfaces=ifaces-nelems*10            
!     
            if(lakon(nelems)(4:5).eq.'8R') then
              mint2d=1
              nopes=4
              nope=8
            elseif(lakon(nelems)(4:4).eq.'8') then
              mint2d=4
              nopes=4
              nope=8
            elseif(lakon(nelems)(4:6).eq.'20R') then
              mint2d=4
              nopes=8
              nope=20
            elseif(lakon(nelems)(4:5).eq.'20') then
              mint2d=9
              nopes=8
              nope=20
            elseif(lakon(nelems)(4:5).eq.'10') then
              mint2d=3
              nopes=6
              nope=10
            elseif(lakon(nelems)(4:4).eq.'4') then
              mint2d=1
              nopes=3
              nope=4
!     
!     treatment of wedge faces
!     
            elseif(lakon(nelems)(4:4).eq.'6') then
              mint2d=1
              nope=6
              if(jfaces.le.2) then
                nopes=3
              else
                nopes=4
              endif
            elseif(lakon(nelems)(4:5).eq.'15') then
              nope=15
              if(jfaces.le.2) then
                mint2d=3
                nopes=6
              else
                mint2d=4
                nopes=8
              endif
            endif
!     
!     actual position of the nodes belonging to the
!     slave surface
!     
            if((nope.eq.20).or.(nope.eq.8)) then
              do m=1,nopes
                nodefs(m)=kon(ipkon(nelems)+ifaceq(m,jfaces))
                do j=1,3
                  xl2(j,m)=co(j,nodefs(m))+clearini(j,m,jj)
     &                 *reltime+vold(j,nodefs(m))
                enddo
              enddo
            elseif((nope.eq.10).or.(nope.eq.4)) then
              do m=1,nopes
                nodefs(m)=kon(ipkon(nelems)+ifacet(m,jfaces))
                do j=1,3
                  xl2(j,m)=co(j,nodefs(m))+clearini(j,m,jj)
     &                 *reltime+vold(j,nodefs(m))
                enddo
              enddo
            elseif(nope.eq.15) then
              do m=1,nopes
                nodefs(m)=kon(ipkon(nelems)+ifacew2(m,jfaces))
                do j=1,3
                  xl2(j,m)=co(j,nodefs(m))+clearini(j,m,jj)
     &                 *reltime+vold(j,nodefs(m))
                enddo
              enddo
            else
              do m=1,nopes
                nodefs(m)=kon(ipkon(nelems)+ifacew1(m,jfaces))
                do j=1,3
                  xl2(j,m)=co(j,nodefs(m))+clearini(j,m,jj)
     &                 *reltime+vold(j,nodefs(m))
                enddo
              enddo
            endif
!     
!     loop over all integration points in the slave surface
!     
!     actions which are only done in the first iteration of an
!     increment (for each slave face integration point)
!     
!     - determine an opposite master face
!     - determine the local coordinates of the opposite master
!     location
!     - determine the local normal on the opposite master
!     location
!     - determine the representative slave area for the slave
!     integration point
!     
            area=0.d0
            mint2d=islavsurf(2,jj+1)-islavsurf(2,jj)
            if(mint2d.eq.0) cycle
            indexf=islavsurf(2,jj)
!     
            do m=1,mint2d
              igauss=indexf+m
              xi=pslavsurf(1,indexf+m)
              et=pslavsurf(2,indexf+m)
              weight=pslavsurf(3,indexf+m)
!     
c     if(nopes.eq.9) then
c     call shape9q(xi,et,xl2,xsj2,xs2,shp2,iflag)
              if(nopes.eq.8) then
                call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
              elseif(nopes.eq.4) then
                call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
              elseif(nopes.eq.6) then
                call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
c     elseif(nopes.eq.7) then
c     call shape7tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
              else
                call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
              endif
              if(iit.le.0) then
                area=dsqrt(xsj2(1)**2+xsj2(2)**2+xsj2(3)**2)*weight
              endif
!     
!     search a master face for each gauss point and generate a contact
!     spring element if successful
!     
              do k=1,3
                p(k)=0.d0
                do j=1,nopes
                  p(k)=p(k)+shp2(4,j)*xl2(k,j)
                enddo
              enddo
              if(iit.gt.0) then
                if(int(pmastsurf(3,igauss)).eq.0) then
                  isol=0
                else
                  isol=1
                endif
              else
!     
!     determining the kneigh neighboring master contact
!     triangle centers of gravity
!     
                call near3d(xo,yo,zo,x,y,z,nx,ny,nz,p(1),p(2),p(3),
     &               n,neigh,kneigh)
!     
                isol=0
!     
                itriold=0
                itri=neigh(1)+itietri(1,i)-1
                ntriangle=0
                ntriangle_=100
!     
                loop1: do
                  do l=1,3
                    ll=4*l-3
                    dist=straight(ll,itri)*p(1)+
     &                   straight(ll+1,itri)*p(2)+
     &                   straight(ll+2,itri)*p(3)+
     &                   straight(ll+3,itri)
!     
!     1.d-6 was increased to 1.d-3 on 19/04/2012
!     this is important for 2d-calculations or
!     calculations for which structures fit exactly
!     at their boundaries
!     
                    if(dist.gt.1.d-3*dsqrt(area)) then
                      itrinew=imastop(l,itri)
                      if(itrinew.eq.0) then
c     write(*,*) '**border reached'
                        exit loop1
                      elseif((itrinew.lt.itietri(1,i)).or.
     &                       (itrinew.gt.itietri(2,i))) then
c     write(*,*) '**border reached'
                        exit loop1
                      elseif(itrinew.eq.itriold) then
c     write(*,*) '**solution in between triangles'
                        isol=itri
                        exit loop1
                      else
                        call nident(itriangle,itrinew,
     &                       ntriangle,id)
                        if(id.gt.0) then
                          if(itriangle(id).eq.itrinew) then
c     write(*,*) '**circular path;no solution'
                            exit loop1
                          endif
                        endif
                        ntriangle=ntriangle+1
                        if(ntriangle.gt.ntriangle_) then
c     write(*,*) '**too many iterations'
                          exit loop1
                        endif
                        do k=ntriangle,id+2,-1
                          itriangle(k)=itriangle(k-1)
                        enddo
                        itriangle(id+1)=itrinew
                        itriold=itri
                        itri=itrinew
                        cycle loop1
                      endif
                    elseif(l.eq.3) then
c     write(*,*) '**regular solution'
                      isol=itri
                      exit loop1
                    endif
                  enddo
                enddo loop1
              endif
!     
!     integration point is catalogued if opposite master
!     face was detected
!     
              if(isol.eq.0) then
                if(iit.le.0) then
                  pmastsurf(3,igauss)=0.5d0
                endif
              else
!     
!     determining the clearance
!     
!     identifying the element face to which the
!     triangle belongs
!     
                if(iit.le.0) then
                  springarea(1,igauss)=area
                  springarea(2,igauss)=0.d0
                  ifacem=koncont(4,itri)
                  pmastsurf(3,igauss)=ifacem+0.5d0
                else
                  ifacem=int(pmastsurf(3,igauss))
                endif
                nelemm=int(ifacem/10.d0)
                jfacem=ifacem-10*nelemm
!     
                indexe=ipkon(nelemm)
                if(lakon(nelemm)(4:5).eq.'20') then
                  nopem=8
                  nface=6
                elseif(lakon(nelemm)(4:4).eq.'8') then
                  nopem=4
                  nface=6
                elseif(lakon(nelemm)(4:5).eq.'10') then
                  nopem=6
                  nface=4
                elseif(lakon(nelemm)(4:4).eq.'4') then
                  nopem=3
                  nface=4
                elseif(lakon(nelemm)(4:5).eq.'15') then
                  if(jfacem.le.2) then
                    nopem=6
                  else
                    nopem=8
                  endif
                  nface=5
                  nope=15
                elseif(lakon(nelemm)(4:4).eq.'6') then
                  if(jfacem.le.2) then
                    nopem=3
                  else
                    nopem=4
                  endif
                  nface=5
                  nope=6
                else
                  cycle
                endif
!     
!     determining the nodes of the master face
!     
                if(nface.eq.4) then
                  do k=1,nopem
                    nodefm(k)=kon(indexe+ifacet(k,jfacem))
                  enddo
                elseif(nface.eq.5) then
                  if(nope.eq.6) then
                    do k=1,nopem
                      nodefm(k)=kon(indexe+ifacew1(k,jfacem))
                    enddo
                  elseif(nope.eq.15) then
                    do k=1,nopem
                      nodefm(k)=kon(indexe+ifacew2(k,jfacem))
                    enddo
                  endif
                elseif(nface.eq.6) then
                  do k=1,nopem
                    nodefm(k)=kon(indexe+ifaceq(k,jfacem))
                  enddo
                endif
!     
!     total number of nodes belonging to the contact 
!     element
!     
                nope=nopem+nopes
!     
!     orthogonal projection on the master element face
!     
                do k=1,nopem
                  do nn=1,3
                    pl(nn,k)=co(nn,nodefm(k))+vold(nn,nodefm(k))
                  enddo
                enddo
!     
                if(iit.le.0) then 
                  do nn=1,3
                    pproj(nn)=p(nn)
                  enddo
                  call attach_2d(pl,pproj,nopem,ratio,dist,xi,et)
                  pmastsurf(1,indexf+m)=xi
                  pmastsurf(2,indexf+m)=et
                else
                  xi=pmastsurf(1,indexf+m)
                  et=pmastsurf(2,indexf+m)
                endif
!     
c     if(nopem.eq.9) then
c     call shape9q(xi,et,pl,xm,xs2,shp2,iflag)
                if(nopem.eq.8) then
                  call shape8q(xi,et,pl,xm,xs2,shp2,iflag)
                elseif(nopem.eq.4) then
                  call shape4q(xi,et,pl,xm,xs2,shp2,iflag)
                elseif(nopem.eq.6) then
                  call shape6tri(xi,et,pl,xm,xs2,shp2,iflag)
c     elseif(nopem.eq.7) then
c     call shape7tri(xi,et,pl,xm,xs2,shp2,iflag)
                else
                  call shape3tri(xi,et,pl,xm,xs2,shp2,iflag)
                endif
!     
                if(iit.gt.0) then
                  do nn=1,3
                    pproj(nn)=0.d0
                    do k=1,nopem
                      pproj(nn)=pproj(nn)+shp2(4,k)*pl(nn,k)
                    enddo
                  enddo
                endif
!     
                do nn=1,3
                  al(nn)=p(nn)-pproj(nn)
                enddo
!     
!     normal on the surface
!     
                if(iit.le.0) then
                  dm=dsqrt(xm(1)*xm(1)+xm(2)*xm(2)+xm(3)*xm(3))
                  do nn=1,3
                    xn(nn)=xm(nn)/dm
                  enddo
                  pmastsurf(4,igauss)=xn(1)
                  pmastsurf(5,igauss)=xn(2)
                  pmastsurf(6,igauss)=xn(3)
                else
                  xn(1)=pmastsurf(4,igauss)
                  xn(2)=pmastsurf(5,igauss)
                  xn(3)=pmastsurf(6,igauss)
                endif
!     
!     distance from surface along normal (= clearance)
!     
                clear=al(1)*xn(1)+al(2)*xn(2)+al(3)*xn(3)
!     
                if((nmethod.eq.4).or.(nmethod.eq.2)) then
!     
!     dynamic calculation
!     
                  if((clear.gt.0.d0).and.
     &                 (int(elcon(3,1,imat)).ne.4)) then
                    isol=0
                  endif
                else
!     
!     static calculation: more sophisticated
!     algorithm to detect contact
!     
ccccc replace next line on 4 Dec 2021
c     if((istep.eq.1).and.(iit.le.0)) then
ccccc by
                  if((iit.le.0).and.
     &                 ((istep.eq.1).or.(tietol(4,i).gt.0.d0))) then
ccccc end of replacement on 4 Dec 2021                 
                    if(clear.lt.0.d0) then
cccc  start change 12 nov 2020
c     if(iinc.eq.1) then
cccc  end change 12 nov 2020
                      springarea(2,igauss)=clear/(1.d0-theta)
cccc  start change 12 nov 2020
c     endif
cccc  end change 12 nov 2020
                    elseif(clear.lt.1.d0/elcon(2,1,imat)) then
                      clear=0.d0
                    endif
                  endif
                  clear=clear-springarea(2,igauss)*(1.d0-reltime)
!     
!     if iloop=1 AND a new increment was started the
!     number of contact elements at the end of the
!     previous increment is counted (using xstateini)
!     
!     furthermore, the regular penetration criterion is
!     used to decide on the creation of a contact element
!     
                  if(iloop.eq.1) then
                    if((isol.ne.0).and.
     &                   (int(elcon(3,1,imat)).ne.4))then
                      if(((istep.gt.1).or.(iinc.gt.1)).and.
     &                     (iit.le.0).and.(ncmat_.ge.7).and.
     &                     (elcon(6,1,imat).gt.0.d0)) then
                        if(dsqrt(xstateini(4,1,ne0+igauss)**2+
     &                       xstateini(5,1,ne0+igauss)**2+
     &                       xstateini(6,1,ne0+igauss)**2)
     &                       .ge.1.d-30) then
                          iprev=iprev+1
                        endif
                      endif
!     
!     if a local minimum was detected in checkconvergence:
!     remove in an aleatoric way 10 % of the contact elements
!     
                      if(ialeatoric.eq.1) then
                        call random_number(harvest)
c     if(harvest.gt.0.9d0) isol=0
                        if(harvest.gt.(1.d0-alea)) isol=0
                      endif
!     
                      if(isol.ne.0) then
                        if((icutb.eq.0).or.(ncmat_.lt.7).or.
     &                       (elcon(6,1,imat).le.0.d0)) then
!     
!     this is no cut-back: only use a negative
!     clearance as contact criterion
!     (this also applies if no friction is
!     defined)
!     
                          if(clear.gt.0.d0) then
                            isol=0
                          else
                            iact=iact+1
                          endif
                        else
!     
!     this is the first iteration in a cut-back:
!     all contact elements from the end of the
!     previous increment are included as well
!     
                          if((dsqrt(
     &                         xstateini(4,1,ne0+igauss)**2+
     &                         xstateini(5,1,ne0+igauss)**2+
     &                         xstateini(6,1,ne0+igauss)**2)
     &                         .lt.1.d-30).and.
     &                         (clear.gt.0.d0)) then
                            isol=0
                          else
                            iact=iact+1
                          endif
                        endif
                      endif
                    endif
                  else
!     
!     if iloop=2 it means that at the start of a new
!     increment no contact elements were generated
!     based on penetration, although there were contact
!     elements at the end of the last increment. In that
!     case the contact elements from the end of the last
!     increment are taken.
!     
                    if((isol.ne.0).and.
     &                   (int(elcon(3,1,imat)).ne.4)) then
                      if(dsqrt(xstateini(4,1,ne0+igauss)**2+
     &                     xstateini(5,1,ne0+igauss)**2+
     &                     xstateini(6,1,ne0+igauss)**2)
     &                     .lt.1.d-30) isol=0
                    endif
                  endif
!     
                endif
!     
              endif
!     
              if(isol.ne.0) then
!     
!     generation of a contact spring element
!     
                ne=ne+1
                nel=ne
!     
                id=ifree+1
                ifree=ifree+nopem+nopes+4
!     
                ipkon(nel)=id
                lakon(nel)(1:7)='ESPRNGC'
                lakon(nel)(8:8)=char(nopem+48)
                ielmat(1,nel)=imat
!     
!     nasym indicates whether at least one contact
!     spring elements exhibits friction in the present
!     step. If so, nasym=1, else nasym=0; nasym=1
!     triggers the asymmetric equation solver
!     
                if(ncmat_.ge.7) then
                  if((elcon(6,1,imat).gt.0).and.
     &                 (int(elcon(3,1,imat)).ne.4)) then
                    nasym=1
                  endif
                endif
                if(ncmat_.ge.8) then
                  if(elcon(8,1,imat).gt.0) then
                    nasym=1
                  endif
                endif
!     
                kon(id)=nopes+nopem
!     
                do k=1,nopem
                  kon(id+k)=nodefm(k) 
                enddo
                id=id+nopem
                do k=1,nopes
                  kon(id+k)=nodefs(k)
                enddo
                id=id+nopes
                kon(id+1)=igauss
                kon(id+2)=jj
                kon(id+3)=indexf+m
!     
c     write(lakon(ne)(8:8),'(i1)') nopem
c     lakon(ne)(8:8)=char(nopem+48)
!     
c     indexel=indexel+1
!     
                if(filab(1)(3:3).eq.'C') then
                  indexel=indexel+1
                  if((nopem.eq.4).or.(nopem.eq.8)) then
                    if((nopes.eq.4).or.(nopes.eq.8)) then
c     if(filab(1)(3:3).eq.'C') then
                      write(27,100) setname(1:lenset)
 100                  format('*ELEMENT,TYPE=C3D8,ELSET=',A)
                      write(27,*) ne0+indexel,',',nodefm(1),',',
     &                     nodefm(2),',',nodefm(3),',',nodefm(4),
     &                     ',',nodefs(2),',',nodefs(1),',',
     &                     nodefs(4),',',nodefs(3)
c     endif
                    endif
                    if((nopes.eq.3).or.(nopes.eq.6)) then
c     if(filab(1)(3:3).eq.'C') then
                      write(27,101) setname(1:lenset)
 101                  format('*ELEMENT,TYPE=C3D8,ELSET=',A)
                      write(27,*) ne0+indexel,',',nodefm(1),',',
     &                     nodefm(2),',',nodefm(3),',',nodefm(4),
     &                     ',',nodefs(2),',',nodefs(1),',',
     &                     nodefs(3),',',nodefs(3)
c     endif
                    endif
                  endif
                  if((nopem.eq.3).or.(nopem.eq.6)) then
                    if((nopes.eq.4).or.(nopes.eq.8)) then
c     if(filab(1)(3:3).eq.'C') then
                      write(27,102) setname(1:lenset)
 102                  format('*ELEMENT,TYPE=C3D8,ELSET=',A)
                      write(27,*) ne0+indexel,',',nodefm(1),',',
     &                     nodefm(2),',',nodefm(3),',',nodefm(3),
     &                     ',',nodefs(2),',',nodefs(1),',',
     &                     nodefs(4),',',nodefs(3)
c     endif
                    endif
                    if((nopes.eq.3).or.(nopes.eq.6)) then
c     if(filab(1)(3:3).eq.'C') then
                      write(27,103) setname(1:lenset)
 103                  format('*ELEMENT,TYPE=C3D6,ELSET=',A)
                      write(27,*) ne0+indexel,',',nodefm(1),',',
     &                     nodefm(2),',',nodefm(3),',',nodefs(2),
     &                     ',',nodefs(1),',',nodefs(3)
c     endif
                    endif
                  endif
                endif
              else
!     
!     no contact: set internal variables at end of increment
!     to zero 
!     
!     next line is inserted because xstate is not yet
!     allocated for iit <= 0
!     
                if(iit.gt.0) then
                  do k=1,nstate_
                    xstate(k,1,ne0+igauss)=0.d0
                  enddo
                endif
!     
              endif
!     
            enddo               ! m
          enddo                 ! jj
          if((iact.ne.0).or.(iprev.eq.0).or.(nmethod.eq.4)) exit
          write(*,*)'*INFO in gencontelem_f2f: contact lost at the'
          write(*,*)'      start of a new increment; contact'
          write(*,*)'      elements from end of previous increment'
          write(*,*)'      are kept.'
          write(*,*)'      number of previous contact elements:',iprev
          write(*,*)'      number of actual contact elements:',iact
        enddo                   ! iloop
        tietol(4,i)=0.d0
      enddo                     ! ntie
!     
!     closing the file containing the contact elements
!     
      if(filab(1)(3:3).eq.'C') then
        close(27)
      endif
!     
      return
      end
