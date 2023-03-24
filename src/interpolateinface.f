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
!     
!     Subroutine x_interpolate.f 
!     
!     Triangulates the face and interpolates xstate variables according to 
!     the plane equations of these planes
!     
!     by: Jaro Hokkanen
!     
!     
      subroutine interpolateinface(kk,xstate,xstateini,numpts,nstate_,
     &     mi,islavsurf,pslavsurf,
     &     ne0,islavsurfold,pslavsurfold)
!     
      implicit none
!     
      integer numpts,i_int,n_int,i,j,k,kk,l,ll,nn,
     &     koncont(3,2*numpts+1),itri,kflag,neigh(1),kneigh,
     &     imastop(3,2*numpts+1),indexcj,nopespringj,list(numpts),
     &     igauss,mi(*),nstate_,itriangle(100),itriold,
     &     ifaceq(8,6),ip(numpts),ne0,itrinew,ntriangle,
     &     ifacet(6,4),ifacew1(4,5),ifacew2(8,5),n,islavsurf(2,*),
     &     ibin(numpts),ivert1,ntriangle_,nterms,m,islavsurfold(2,*),
     &     nx(2*numpts+1),ny(2*numpts+1),isol,id
!     
      real*8 xstate(nstate_,mi(1),*),p(3),pslavsurfold(3,*),
     &     xstateini(nstate_,mi(1),*),coi(2,numpts+3),pneigh(3,3),
     &     cg(2,2*numpts+1),pslavsurf(3,*),xil,etl,ratio(3),
     &     z(9,3),x(2*numpts+1),xo(2*numpts+1),dis(3),
     &     y(2*numpts+1),yo(2*numpts+1),straight(9,2*numpts+1),dist
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
      kneigh=1
!     
      do i=1,numpts
        list(i)=i
      enddo
!     
!     Loop over the old integration points within the face
!     
      ll=0
      do l=islavsurfold(2,kk)+1,islavsurfold(2,kk+1)
        ll=ll+1
        coi(1,ll)=pslavsurfold(1,l)
        coi(2,ll)=pslavsurfold(2,l)
        ip(ll)=ne0+l
      enddo
!     
!     Calling deltri for triangulation
!     
      kflag=2
!     
      call deltri(numpts,numpts,coi(1,1:numpts+3),coi(2,1:numpts+3),
     &     list,ibin,koncont,imastop,n)
!     
!     rearranging imastop field: imastop(1,i) is the triangle opposite
!     of local node 1 in triangle i etc..      
!     
      nn=0
      do i=1,n
        ivert1=imastop(1,i)
        imastop(1,i)=imastop(2,i)
        imastop(2,i)=imastop(3,i)
        imastop(3,i)=ivert1
      enddo
!     
!     determining the center of gravity and the bounding planes
!     of the triangles
!     
      call updatecont2d(koncont,n,coi,cg,straight)
!     
!     sorting the centers of gravity
!     
      do l=1,n
        xo(l)=cg(1,l)
        x(l)=xo(l)
        nx(l)=l
        yo(l)=cg(2,l)
        y(l)=yo(l)
        ny(l)=l
      enddo
      kflag=2
      call dsort(x,nx,n,kflag)
      call dsort(y,ny,n,kflag)
!     
!     Loop over the new integration points
!
      do igauss=islavsurf(2,kk)+1,islavsurf(2,kk+1)
!     
!     coordinates of the new integration point
!     
        p(1)=pslavsurf(1,igauss)
        p(2)=pslavsurf(2,igauss)
!     
!     closest triangle center of gravity
!     
        call near2d(xo,yo,x,y,nx,ny,p(1),p(2),
     &       n,neigh,kneigh)
!     
        isol=0
!     
        itriold=0
        itri=neigh(1)
        ntriangle=0
        ntriangle_=100
!     
        loop1: do
        do l=1,3
          ll=3*l-2
          dist=straight(ll,itri)*p(1)+
     &         straight(ll+1,itri)*p(2)+
     &         straight(ll+2,itri)
!     
          if(dist.gt.0.d0) then
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
!     
      if(isol.eq.0) then
!     
!     mapping the new integration point onto the nearest 
!     triangle
!     
!     calculating the distance from each edge of the triangle
!     dis(i) is the distance of the new integration point fromt
!     the edge opposite to local node i of the triangle itri
!     
        do l=1,3
          ll=3*l-2
          dis(l)=straight(ll,itri)*p(1)+
     &         straight(ll+1,itri)*p(2)+
     &         straight(ll+2,itri)
        enddo
!
        if(dis(1).gt.0.d0) then
          if(dis(2).gt.0.d0) then
!
!     projection is node 3 of the triangle
!
            p(1)=coi(1,koncont(3,itri))
            p(2)=coi(2,koncont(3,itri))
          elseif(dis(3).gt.0.d0) then
!
!     projection is node 2 of the triangle
!
            p(1)=coi(1,koncont(2,itri))
            p(2)=coi(2,koncont(2,itri))
          else
!
!     projection onto local edge 1 (opposite to local node 1)
!
            p(1)=p(1)-dis(1)*straight(1,itri)
            p(2)=p(2)-dis(1)*straight(2,itri)
          endif
        elseif(dis(2).gt.0.d0) then
          if(dis(3).gt.0.d0) then
!
!     projection is node 1 of the triangle
!
            p(1)=coi(1,koncont(1,itri))
            p(2)=coi(2,koncont(1,itri))
          else
!
!     projection onto local edge 2 (opposite to local node 2)
!
            p(1)=p(1)-dis(2)*straight(4,itri)
            p(2)=p(2)-dis(2)*straight(5,itri)
          endif
        else
!
!     projection onto local edge 3 (opposite to local node 3)
!
          p(1)=p(1)-dis(3)*straight(7,itri)
          p(2)=p(2)-dis(3)*straight(8,itri)
        endif
c     do k=1,3
c     do m=1,2
c     pneigh(m,k)=coi(m,koncont(k,itri))
c     enddo
c     pneigh(3,k)=0.d0
c     enddo
c     p(3)=0.d0
c     nterms=3
c     !
c     call attach_2d(pneigh,p,nterms,ratio,dist,xil,etl)
c     !
c     do m=1,2
c     p(m)=0.d0
c     do k=1,3
c     p(m)=p(m)+ratio(k)*pneigh(m,k)
c     enddo
c     enddo
c        
      endif
!     
!     Assigning xstate values from the vertices of the triangle
!     to the field z (integration point values from xstateini)
!     
      do k=1,3
        do i=1,9
          z(i,k)=xstateini(i,1,ip(koncont(k,itri)))
        enddo
      enddo
!     
!     Calling plane_eq to interpolate values using plane equation
!     
      do i=1,9   
        call plane_eq(
     &       coi(1,koncont(1,itri)),coi(2,koncont(1,itri)),z(i,1),
     &       coi(1,koncont(2,itri)),coi(2,koncont(2,itri)),z(i,2),
     &       coi(1,koncont(3,itri)),coi(2,koncont(3,itri)),z(i,3),
     &       p(1),p(2),xstate(i,1,ne0+igauss))
      enddo
!     
      enddo
!     
      return
      end
      
