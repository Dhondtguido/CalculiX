!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
!     
!     subroutine generating the first actice set 
!     and calculating the areas for all faces on slaveside
!     
!     [in]  imastop        	(l,i) for edge l in triagle i neightbouring triangle
!     [in]  islavsurf      	islavsurf(1,i) slaveface i islavsurf(2,i) pointer into imastsurf and pmastsurf
!     [in]  itiefac     	pointer into field islavsurf: (1,i) beginning slave_i (2,i) end of slave_i
!     [out]    areaslav    	(i)area of face i, ONLY HELP FIELD
!     [in,out] ifree       	in=0 out=# active nodes
!     
      subroutine genfirstactif(tieset,ntie,itietri,ipkon,kon,lakon,cg,
     &     straight,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,istep,iinc,iit,
     &     mi,imastop,nslavnode,islavnode,islavsurf,itiefac,areaslav,
     &     set,nset,istartset,iendset,ialset,islavact,ifree,tietol)
!     
!     Initialization of the Active slave nodes set
!     
!     Author: Samoela Rakotonanahary, Saskia Sitzmann      
!
!     guido: areaslav not needed??
!      
      implicit none
!     
      character*8 lakon(*)
      character*81 tieset(3,*),slavset,set(*),noset
!     
      integer ntie,
     &     itietri(2,ntie),ipkon(*),kon(*),node,neigh(1),iflag,kneigh,
     &     i,j,k,l,isol,iset,idummy,itri,ll,kflag,n,nx(*),ny(*),istep,
     &     iinc,mi(*),nz(*),nstart,iit,nope,iteller,ifaces,jfaces,
     &     imastop(3,*), itriangle(100),ntriangle,ntriangle_,itriold,
     &     itrinew,id,nslavnode(*),islavnode(*),islavsurf(2,*),
     &     itiefac(2,*),konl(20),nelems,m,mint2d,nopes,
     &     ipos,nset,istartset(*),iendset(*),
     &     ialset(*),islavact(*),ifree,ifac,getlocno
!     
      real*8 cg(3,*),straight(16,*),co(3,*),vold(0:mi(2),*),p(3),
     &     dist,xo(*),yo(*),zo(*),x(*),y(*),z(*),c0,weight,
     &     areaslav(*),xl2(3,8),area,xi,et,shp2(7,8),
     &     xs2(3,2),xsj2(3),tietol(4,*),adjust
!     
!     flag for shape functions
!     
      data iflag /2/
!     
      data iteller /0/
      save iteller
!     
      include "gauss.f"
!     
      do i=1,ntie
        if(tieset(1,i)(81:81).ne.'C') cycle
        kneigh=1
        slavset=tieset(2,i)
!     
!     check whether an adjust node set has been defined
!     only checked in the first increment of the first step
!     
        if((istep.eq.1).and.(iinc.eq.1).and.(iit.le.1)) then
          iset=0
          if(tieset(1,i)(1:1).ne.' ') then
            noset(1:80)=tieset(1,i)(1:80)
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
c            do iset=1,nset
c              if(set(iset).eq.noset) exit
c            enddo
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
!     determine the area of the slave surfaces
!     
        do l=itiefac(1,i),itiefac(2,i)
          ifaces=islavsurf(1,l)
          nelems=int(ifaces/10)
          jfaces=ifaces-nelems*10
!     
!     Decide on the max integration points number, just consider 2D situation 
!     
          call getnumberofnodes(nelems,jfaces,lakon,nope,nopes,mint2d)
!     
!     actual position of the nodes belonging to the
!     slave surface
!     
          do j=1,nope
            konl(j)=kon(ipkon(nelems)+j)
          enddo
          do m=1,nopes
            ifac=getlocno(m,jfaces,nope)
            do j=1,3                  
              xl2(j,m)=co(j,konl(ifac))+
     &             vold(j,konl(ifac))                
            enddo
          enddo
!     
!     calculating the area of the slave face
!     
          area=0.d0
          do m=1,mint2d
            if((lakon(nelems)(4:5).eq.'8R').or.
     &           ((lakon(nelems)(4:4).eq.'6').and.(nopes.eq.4))) then
              xi=gauss2d1(1,m)
              et=gauss2d1(2,m)
              weight=weight2d1(m)
            elseif((lakon(nelems)(4:4).eq.'8').or.
     &             (lakon(nelems)(4:6).eq.'20R').or.
     &             ((lakon(nelems)(4:5).eq.'15').and.
     &             (nopes.eq.8))) then
              xi=gauss2d2(1,m)
              et=gauss2d2(2,m)
              weight=weight2d2(m)
            elseif(lakon(nelems)(4:4).eq.'2') then
              xi=gauss2d3(1,m)
              et=gauss2d3(2,m)
              weight=weight2d3(m)
            elseif((lakon(nelems)(4:5).eq.'10').or.
     &             ((lakon(nelems)(4:5).eq.'15').and.
     &             (nopes.eq.6))) then
              xi=gauss2d5(1,m)
              et=gauss2d5(2,m)
              weight=weight2d5(m)
            elseif((lakon(nelems)(4:4).eq.'4').or.
     &             ((lakon(nelems)(4:4).eq.'6').and.
     &             (nopes.eq.3))) then
              xi=gauss2d4(1,m)
              et=gauss2d4(2,m)
              weight=weight2d4(m)
            endif
!     
            if(nopes.eq.8) then
              call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.4) then
              call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.6) then
              call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            else
              call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            endif
            area=area+weight*dsqrt(xsj2(1)**2+xsj2(2)**2+
     &           xsj2(3)**2)
          enddo
          areaslav(l)=area
        enddo
!     
!     search a master face for each slave node and generate a contact
!     spring element if successful
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
            p(k)=co(k,node)+vold(k,node)
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
                elseif((itrinew.lt.itietri(1,i)).or.
     &                 (itrinew.gt.itietri(2,i))) then
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
!     check whether distance is larger than c0:
!     no element is generated
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
            if((istep.eq.1).and.(iinc.eq.1).and.(iit.le.1)) then
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
            c0=1.d-10
            if(dabs(tietol(1,i)).ge.2.d0) then
              c0=dabs(tietol(1,i))-2.d0
            endif
            if(dist.gt.c0) then
              isol=0
!     
!     adjusting the bodies at the start of the
!     calculation such that they touch
!     
            endif
          endif
!     
          if(isol.ne.0) then
!     
!     Active node
            islavact(j)=2
            ifree=ifree+1
          else
            islavact(j)=-1
!     
          endif
!     
        enddo
      enddo
!     
      return
      end
