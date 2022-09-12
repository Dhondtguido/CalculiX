!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998 Guido Dhondt
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
      subroutine interpoltet(x,y,z,xo,yo,zo,nx,ny,nz,planfa,ifatet,
     &     netet,kontet,cotet,iparent,co,nkfa,nkfb,konl,ratio,ikf)
!     
      implicit none
!     
!     10 nearest nodes: neighbor(10)
!     
      integer ifatet(4,*),neighbor(10),near,nx(*),ny(*),nz(*),
     &     ifs,iface,i,konl(4,*),ielmax,netet,m,ikf(*),j,kontet(*),
     &     iparent(*),nterms,indexe,nelem,ii,nkfa,nkfb,loopa,
     &     jj,inside,konl_opt(4),k
!
      real*8 cotet(3,*),planfa(4,*),dface,dfacemax,tolerance,dist,
     &     ratio(4,*),pneigh(3,4),pnode(3),xi,et,ze,dist_opt,
     &     xp,yp,zp,x(*),y(*),z(*),xo(*),yo(*),zo(*),co(3,*),
     &     ratio_opt(4)
!     
      tolerance=1.d-6
      nterms=4
      loopa=8
c      write(*,*) 'interpoltet'
!     
      do m=nkfa,nkfb
!        
        jj=ikf(m)
        xp=co(1,jj)
        yp=co(2,jj)
        zp=co(3,jj)
!     
        do ii=1,2
!     
          if(ii.eq.1) then
            near=1
          else
            near=10
          endif
          call near3d(xo,yo,zo,x,y,z,nx,ny,nz,xp,yp,zp,netet,neighbor,
     &near)
!     
          inside=0
!
!         a node is inside the tetrahedron if the substitution
!         of its coordinates in the equation of the tet faces
!         (corrected by ifs/abs(ifs)) is for all faces positive
!          
          do i=1,near
            nelem=iparent(neighbor(i))
            dface=0.d0
            do j=1,4
              ifs=ifatet(j,nelem)
              iface=abs(ifs)
              dist=planfa(1,iface)*xp+planfa(2,iface)*yp
     &             +planfa(3,iface)*zp+planfa(4,iface)
              if(dist*ifs.lt.-1.d-10*iface) then
                dface=dface+dist*ifs/iface
              endif
            enddo
            if(dface.gt.-1.d-10) then
              inside=1
c              write(*,*) 'interpoltet inside'
              exit
            endif
            if(i.eq.1) then
              dfacemax=dface
              ielmax=nelem
            else
              if(dface.gt.dfacemax) then
                ielmax=nelem
                dfacemax=dface
              endif
            endif
          enddo
!     
!     if no element was found, the element with the smallest
!     dfacemax (summed distance) is taken 
!     
          if(inside.eq.0) then
            nelem=ielmax
          endif
!     
          indexe=4*(nelem-1)
!     
          do i=1,4
            konl(i,m)=kontet(indexe+i)
          enddo
!     
!     nodes of master element
!     
          do i=1,4
            do j=1,3
              pneigh(j,i)=cotet(j,konl(i,m))
            enddo
          enddo
!     
!     slave node
!     
          pnode(1)=xp
          pnode(2)=yp
          pnode(3)=zp
!     
!     attaching slave node to master element
!     
          call attach_3d(pneigh,pnode,nterms,ratio(1,m),dist,xi,et,ze,
     &         loopa)
c          write(*,*) nelem,dist,tolerance
!     
!     checking the parent elements of the "best" tetrahedra
!     in case the distance between slave node and location of
!     interpolation exceeds "tolerance"
!     
          if(dist.gt.tolerance) then
            if(ii.eq.1) cycle
            do i=1,4
              konl_opt(i)=konl(i,m)
              ratio_opt(i)=ratio(i,m)
            enddo
            dist_opt=dist
!     
            do k=1,near
!     
!     slave node
!     
              pnode(1)=xp
              pnode(2)=yp
              pnode(3)=zp
!     
              nelem=iparent(neighbor(k))
              indexe=4*(nelem-1)
!     
              do i=1,4
                konl(i,m)=kontet(indexe+i)
              enddo
!     
!     nodes of master element
!     
              do i=1,4
                do j=1,3
                  pneigh(j,i)=cotet(j,konl(i,m))
                enddo
              enddo
!     
!     attaching slave node to master element
!     
              call attach_3d(pneigh,pnode,nterms,ratio(1,m),dist,
     &             xi,et,ze,loopa)
!     
!     check whether the present element yields better results
!     
              if(dist.lt.dist_opt) then
                do i=1,4
                  konl_opt(i)=konl(i,m)
                  ratio_opt(i)=ratio(i,m)
                enddo
                dist_opt=dist
              endif
              if(dist.lt.tolerance) exit
            enddo
!     
!     storing the optimal configuration
!     
            do i=1,4
              konl(i,m)=konl_opt(i)
              ratio(i,m)=ratio_opt(i)
            enddo
          endif
!     
          if((ii.eq.2).or.(dist.lt.tolerance)) exit
!
c          if(dist.lt.tolerance) exit
!          
        enddo
c        write(*,*) jj,dist,co(1,jj),ratio(1,m)*pneigh(1,1)+
c     &       ratio(2,m)*pneigh(1,2)+ratio(3,m)*pneigh(1,3)+
c     &       ratio(4,m)*pneigh(1,4)
!     
      enddo
!     
      return
      end
