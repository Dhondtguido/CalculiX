      
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
      subroutine extendmesh(nnfront,istartfront,iendfront,ifront,
     &     ne,nkon,lakon,ipkon,kon,isubsurffront,co,ifronteq,
     &     istartfronteq,iendfronteq,nfront,nfronteq)
!     
!     calculate the extension of the mesh
!     
      implicit none
!     
      character*8 lakon(*)
!     
      integer i,j,nnfront,istartfront(*),iendfront(*),ifront(*),
     &     node1,node2,node1e,node2e,ne,nkon,ipkon(*),kon(*),
     &     isubsurffront(*),k,nx(nfront),ny(nfront),
     &     nz(nfront),kflag,ifronteq(nfronteq),istartfronteq(*),
     &     iendfronteq(*),neigh(1),kneigh,neighbor(nfronteq),
     &     nfronteq,nedges,l,m,nfront
!     
      real*8 co(3,*),x0(nfront),y0(nfront),z0(nfront),
     &     x(nfront),y(nfront),z(nfront),xp,yp,zp
!     
c      write(*,*)'extendmesh',nnfront
!     
!     loop over all fronts
!     
      do i=1,nnfront
         l=istartfronteq(i)
         m=iendfronteq(i)
         k=0
         do j=istartfront(i),iendfront(i)
            k=k+1
            x0(k)=co(1,ifront(j))
            y0(k)=co(2,ifront(j))
            z0(k)=co(3,ifront(j))
            x(k)=x0(k)
            y(k)=y0(k)
            z(k)=z0(k)
            nx(k)=k
            ny(k)=k
            nz(k)=k
         enddo
!     
         kflag=2
         call dsort(x,nx,k,kflag)
         call dsort(y,ny,k,kflag)
         call dsort(z,nz,k,kflag)
!     
!     
         kneigh=1
         neighbor(l)=istartfront(i)
!     
         do j=l+1,m-1           
            xp=co(1,ifronteq(j))
            yp=co(2,ifronteq(j))
            zp=co(3,ifronteq(j))
            call near3d(x0,y0,z0,x,y,z,nx,ny,nz,xp,yp,zp,k,
     &           neigh,kneigh)
            neighbor(j)=istartfront(i)-1+neigh(1)
         enddo    
         neighbor(m)=iendfront(i)
!     
         do j=l,m-1
!     
!     the topology is: first node on the "starting" front and 
!     the other two on the propagated front
!     
            node1=ifront(neighbor(j))
            node1e=ifronteq(j)
            node2e=ifronteq(j+1)
!     
            ne=ne+1
            lakon(ne)='C3D6  L '
            ipkon(ne)=nkon
            kon(nkon+1)=node1
            kon(nkon+2)=node1e
            kon(nkon+3)=node2e
            kon(nkon+4)=node1
            kon(nkon+5)=node1e
            kon(nkon+6)=node2e
            kon(nkon+7)=node1
            kon(nkon+8)=node1e
            kon(nkon+9)=node2e
            nkon=nkon+9 
!     
!     check if two adjacent propagated nodes have the same neighbor:
!     if this is different new triangles are defined 
!     
            if(neighbor(j).ne.neighbor(j+1)) then
!     
!     nedges is the number of edges,on the front, between 
!     neighbor(j) and neighbor(j+1)
!     
               nedges=neighbor(j+1)-neighbor(j)
               do k=1,nedges
                  node1=ifront(neighbor(j)+k-1)
                  node2e=ifronteq(j+1)
                  node2=ifront(neighbor(j)+k)
!     
                  ne=ne+1
                  lakon(ne)='C3D6  L '
                  ipkon(ne)=nkon
                  kon(nkon+1)=node1
                  kon(nkon+2)=node2e
                  kon(nkon+3)=node2
                  kon(nkon+4)=node1
                  kon(nkon+5)=node2e
                  kon(nkon+6)=node2
                  kon(nkon+7)=node1
                  kon(nkon+8)=node2e
                  kon(nkon+9)=node2
                  nkon=nkon+9
               enddo
            endif
         enddo
c         write(*,*)'FRONT NUMBER:',i
c         do j=l,m            
c            write(*,*)'neighbor,ifronteq',ifront(neighbor(j)),
c     &           ifronteq(j)
c         enddo
!     
!     
!     check for a subsurface crack: create two additional
!     shell elements
!     
         if(isubsurffront(i).eq.1) then
c            write(*,*)'isubsurffront=1:',i
            node1=ifront(neighbor(l))
            node1e=ifronteq(m)
            node2e=ifronteq(l)
!     
            ne=ne+1
            lakon(ne)='C3D6  L '
            ipkon(ne)=nkon
            kon(nkon+1)=node1
            kon(nkon+2)=node1e
            kon(nkon+3)=node2e
            kon(nkon+4)=node1
            kon(nkon+5)=node1e
            kon(nkon+6)=node2e
            kon(nkon+7)=node1
            kon(nkon+8)=node1e
            kon(nkon+9)=node2e
            nkon=nkon+9 
!     
            node1=ifront(neighbor(m))
            node2e=ifronteq(l)
            node2=ifront(neighbor(l))
!     
            ne=ne+1
            lakon(ne)='C3D6  L '
            ipkon(ne)=nkon
            kon(nkon+1)=node1
            kon(nkon+2)=node2e
            kon(nkon+3)=node2
            kon(nkon+4)=node1
            kon(nkon+5)=node2e
            kon(nkon+6)=node2
            kon(nkon+7)=node1
            kon(nkon+8)=node2e
            kon(nkon+9)=node2
            nkon=nkon+9
         endif
      enddo
!     
      return
      end

