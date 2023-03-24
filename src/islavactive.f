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
!   Sets the value 1 to field islavact for nodes which have 
!   opposite master face.
!
      subroutine islavactive(tieset,ntie,itietri,
     &  cg,straight,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,
     &  mi,imastop,nslavnode,islavnode,islavact)
!
      implicit none
!
      character*81 tieset(3,*)
!
      integer ntie,itietri(2,*),node,neigh(1),iflag,kneigh,i,j,k,l,
     &  isol,itri,ll,kflag,n,nx(*),ny(*),mi(*),nz(*),nstart,
     &  ifacew1(4,5),ifacew2(8,5),imastop(3,*),
     &  itriangle(100),ntriangle,ntriangle_,itriold,itrinew,id,
     &  nslavnode(*),islavnode(*),islavact(*),ifaceq(8,6),
     &  ifacet(6,4)
!
      real*8 cg(3,*),straight(16,*),co(3,*),vold(0:mi(2),*),p(3),
     &  dist,xo(*),yo(*),zo(*),x(*),y(*),z(*)
!
!     nodes per face for hex elements
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
!
!     nodes per face for tet elements
!
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
!
!     nodes per face for linear wedge elements
!
      data ifacew1 /1,3,2,0,
     &             4,5,6,0,
     &             1,2,5,4,
     &             2,3,6,5,
     &             3,1,4,6/
!
!     nodes per face for quadratic wedge elements
!     
      data ifacew2 /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             3,1,4,6,9,13,12,15/
!
      data iflag /2/
!
!
!
!     ***ISLAVACT***
!
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.'C') cycle
         kneigh=1
!
!        search a master face for each slave node
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
c         write(*,*) 'islavactive ntie= ',i
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
c            write(*,*) 'islavactive j= ',j
            call near3d(xo,yo,zo,x,y,z,nx,ny,nz,p(1),p(2),p(3),
     &           n,neigh,kneigh)
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
     &                 straight(ll+1,itri)*p(2)+
     &                 straight(ll+2,itri)*p(3)+
     &                 straight(ll+3,itri)
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
            if(isol.ne.0) then
!     Active node
               islavact(j)=1
            else
               islavact(j)=0
            endif
         enddo
      enddo
!
      return
      end
      
