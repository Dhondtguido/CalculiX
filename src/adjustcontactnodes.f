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
      subroutine adjustcontactnodes(tieset,ntie,itietri,cg,straight,
     &     co,vold,xo,yo,zo,x,y,z,nx,ny,nz,istep,iinc,iit,
     &     mi,imastop,nslavnode,islavnode,set,nset,istartset,
     &     iendset,ialset,tietol,clearini,clearslavnode,itiefac,
     &     ipkon,kon,lakon,islavsurf)
!     
!     Adjusting contact nodes if requested by the user
!     
      implicit none
!     
      character*8 lakon(*)
      character*81 tieset(3,*),slavset,set(*),noset
!     
      integer ntie,
     &     itietri(2,ntie),node,neigh(1),kneigh,i,j,k,l,m,isol,iset,
     &     idummy,itri,ll,kflag,n,nx(*),ny(*),istep,iinc,mi(*),
     &     nz(*),nstart,iit,imastop(3,*),itriangle(100),ntriangle,
     &     ntriangle_,itriold,itrinew,id,nslavnode(*),islavnode(*),
     &     ipos,nset,istartset(*),iendset(*),ialset(*),iclear,
     &     istart,ilength,nope,nopes,nelems,jj,ifaces,itiefac(2,*),
     &     jfaces,islavsurf(2,*),ifaceq(8,6),ifacet(6,4),ifacew1(4,5),
     &     ifacew2(8,5),ipkon(*),kon(*)
!     
      real*8 cg(3,*),straight(16,*),co(3,*),vold(0:mi(2),*),p(3),
     &     dist,xo(*),yo(*),zo(*),x(*),y(*),z(*),tietol(4,*),adjust,
     &     clearini(3,9,*),clearslavnode(3,*),clear
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
      do i=1,ntie
cccc  replaced on 06.12.2021 next line      
c     if(tieset(1,i)(81:81).ne.'C') cycle
cccc  by
c     if(tieset(1,i)(81:81).ne.'C') cycle
        if((tieset(1,i)(81:81).ne.'C').and.(tieset(1,i)(81:81).ne.'-'))
     &       cycle
cccc  end replace        
        kneigh=1
        slavset=tieset(2,i)
!     
!     check whether a clearance has been specified by the user
!     
c     tietol(3,i)=1.2357111317d0
        if((tietol(3,i).gt.1.2357111316d0).and.
     &       (tietol(3,i).lt.1.2357111318d0)) then
          iclear=0
        else
          iclear=1
          clear=tietol(3,i)
        endif
!     
!     check whether an adjust node set has been defined
!     only checked in the first increment of the first step
!     
        if(istep.eq.1) then
          iset=0
          if(tieset(1,i)(1:1).ne.' ') then
            noset(1:80)=tieset(1,i)(1:80)
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
c     do iset=1,nset
c     if(set(iset).eq.noset) exit
c     enddo
            call cident81(set,noset,nset,id)
            iset=nset+1
            if(id.gt.0) then
              if(noset.eq.set(id)) then
                iset=id
              endif
            endif
            kflag=1
            call isortii(ialset(istartset(iset)),idummy,
     &           iendset(iset)-istartset(iset)+1,kflag)
          endif
        endif
!     
!     search a master face for each slave node; determine the
!     distance and adjust it if requested by the user
!     
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
!     
        do j=nslavnode(i)+1,nslavnode(i+1)
          node=islavnode(j)
!     
          do k=1,3
            p(k)=co(k,node)
c     p(k)=co(k,node)+vold(k,node)
          enddo
!     
!     determining the kneigh neighboring master contact
!     triangle centers of gravity
!     
          call near3d(xo,yo,zo,x,y,z,nx,ny,nz,p(1),p(2),p(3),
     &         n,neigh,kneigh)
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
     &             straight(ll+1,itri)*p(2)+
     &             straight(ll+2,itri)*p(3)+
     &             straight(ll+3,itri)
              if(dist.gt.1.d-6) then
                itrinew=imastop(l,itri)
                if(itrinew.eq.0) then
c     write(*,*) '**border reached'
                  exit loop1
                elseif(itrinew.eq.itriold) then
c     write(*,*) '**solution in between triangles'
                  isol=itri
                  exit loop1
                else
                  call nident(itriangle,itrinew,ntriangle,id)
                  if(id.gt.0) then
                    if(itriangle(id).eq.itrinew) then
c     write(*,*) '**circular path; no solution'
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
!     
!     calculating the distance
!     
          if(isol.ne.0) then
            dist=straight(13,itri)*p(1)+
     &           straight(14,itri)*p(2)+
     &           straight(15,itri)*p(3)+
     &           straight(16,itri)
!     
!     check for an adjust parameter (only in the first
!     increment of the first step)
!     
            if(istep.eq.1) then
              if(iset.ne.0) then
!     
!     check whether node belongs to the adjust node
!     set
!     
                call nident(ialset(istartset(iset)),node,
     &               iendset(iset)-istartset(iset)+1,id)
                if(id.gt.0) then
                  if(ialset(istartset(iset)+id-1).eq.node) then
                    do k=1,3
                      co(k,node)=co(k,node)-
     &                     dist*straight(12+k,itri)
                    enddo
                    dist=0.d0
                  endif
                endif
              elseif(dabs(tietol(1,i)).ge.2.d0) then
!     
!     adjust parameter
!     
                adjust=dabs(tietol(1,i))-2.d0
                
                if(dist.le.adjust) then
                  do k=1,3
                    co(k,node)=co(k,node)-
     &                   dist*straight(12+k,itri)
                  enddo
                  dist=0.d0
                endif
              endif
            endif
!     
!     clearance in each slave node
!     
            if(iclear.eq.1) then
              do k=1,3
                clearslavnode(k,j)=(clear-dist)*straight(12+k,itri)
              enddo
            endif
!     
          endif
        enddo
!     
!     loop over all slave faces
!     
        if(iclear.eq.1) then
          do jj=itiefac(1,i), itiefac(2,i)
            ifaces=islavsurf(1,jj)
            nelems=int(ifaces/10)
            jfaces=ifaces - nelems*10            
!     
            if(lakon(nelems)(4:5).eq.'8R') then
              nopes=4
              nope=8
            elseif(lakon(nelems)(4:4).eq.'8') then
              nopes=4
              nope=8
            elseif(lakon(nelems)(4:6).eq.'20R') then
              nopes=8
              nope=20
            elseif(lakon(nelems)(4:5).eq.'20') then
              nopes=8
              nope=20
            elseif(lakon(nelems)(4:5).eq.'10') then
              nopes=6
              nope=10
            elseif(lakon(nelems)(4:4).eq.'4') then
              nopes=3
              nope=4
!     
!     treatment of wedge faces
!     
            elseif(lakon(nelems)(4:4).eq.'6') then
              nope=6
              if(jfaces.le.2) then
                nopes=3
              else
                nopes=4
              endif
            elseif(lakon(nelems)(4:5).eq.'15') then
              nope=15
              if(jfaces.le.2) then
                nopes=6
              else
                nopes=8
              endif
            endif
!     
!     actual position of the nodes belonging to the
!     slave surface
!     
            istart=nslavnode(i)+1
            ilength=nslavnode(i+1)-nslavnode(i)
!     
            do m=1,nopes
              if((nope.eq.20).or.(nope.eq.8))then
                node=kon(ipkon(nelems)+ifaceq(m,jfaces))
              elseif((nope.eq.10).or.(nope.eq.4).or.(nope.eq.14)) 
     &               then
                node=kon(ipkon(nelems)+ifacet(m,jfaces))
              elseif(nope.eq.15) then
                node=kon(ipkon(nelems)+ifacew2(m,jfaces))
              else
                node=kon(ipkon(nelems)+ifacew1(m,jfaces))
              endif
              call nident(islavnode(istart),node,ilength,id)
              do k=1,3
                clearini(k,m,jj)=
     &               clearslavnode(k,nslavnode(i)+id)
              enddo
            enddo
          enddo
        endif
!     
!     transfer the clearance from the slave nodes in islavnode
!     to the nodes per slave face in islavsurf
!     
      enddo
!     
      return
      end
      
