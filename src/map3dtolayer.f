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
      subroutine map3dtolayer(yn,ipkon,kon,lakon,nfield,
     &  ne,co,ielmat,mi)
!
!     interpolates 3d field nodal values to nodal values in the
!     layers of composite materials
!
      implicit none
!
      character*8 lakon(*)
!
      integer ipkon(*),kon(*),ne,nfield,i,j,k,l,mi(*),nelem,
     &  ielmat(mi(3),*),nlayer,iflag,nbot20(8),nmid20(8),ntop20(8),
     &  nbot15(6),nmid15(6),ntop15(6),ibot,itop,nope,nopes,nopeexp,
     &  indexe,konl(20),imid,node,nodebot,nodetop
!
      real*8 yn(nfield,*),shp(4,20),xsj(3),co(3,*),xl(3,20),
     &  xi,et,ze,dd,dt,xi20(8),et20(8),xi15(6),et15(6)
!
!
!
      data iflag /1/
!
      data nbot20 /1,2,3,4,9,10,11,12/
      data nmid20 /17,18,19,20,0,0,0,0/
      data ntop20 /5,6,7,8,13,14,15,16/
!
      data nbot15 /1,2,3,7,8,9/
      data nmid15 /13,14,15,0,0,0/
      data ntop15 /4,5,6,10,11,12/
!
      data xi20 /-1.d0,1.d0,1.d0,-1.d0,0.d0,1.d0,0.d0,-1.d0/
      data et20 /-1.d0,-1.d0,1.d0,1.d0,-1.d0,0.d0,1.d0,0.d0/
!
      data xi15 /0.d0,1.d0,0.d0,0.5d0,0.5d0,0.d0/
      data et15 /0.d0,0.d0,1.d0,0.d0,0.5d0,0.5d0/
!
      include "gauss.f"
!
      do nelem=1,ne
         if(lakon(nelem)(7:8).eq.'LC') then
!
!        composite materials
!
!        determining the number of layers
!
            nlayer=0
            do k=1,mi(3)
               if(ielmat(k,nelem).ne.0) then
                  nlayer=nlayer+1
               endif
            enddo
!
            indexe=ipkon(nelem)
!
            if(lakon(nelem)(4:5).eq.'20') then
               nope=20
               nopes=8
               nopeexp=28
            elseif(lakon(nelem)(4:5).eq.'15') then
               nope=15
               nopes=6
               nopeexp=21
            endif
!
            do i=1,nope
               konl(i)=kon(indexe+i)
               do j=1,3
                  xl(j,i)=co(j,konl(i))
               enddo
            enddo
!     
            do i=1,nopes
               if(lakon(nelem)(4:5).eq.'20') then
                  ibot=nbot20(i)
                  imid=nmid20(i)
                  itop=ntop20(i)
                  xi=xi20(i)
                  et=et20(i)
               else
                  ibot=nbot15(i)
                  imid=nmid15(i)
                  itop=ntop15(i)
                  xi=xi15(i)
                  et=et15(i)
               endif
!
               nodebot=konl(ibot)
               nodetop=konl(itop)
               dd=sqrt((co(1,nodebot)-co(1,nodetop))**2+
     &                 (co(2,nodebot)-co(2,nodetop))**2+
     &                 (co(3,nodebot)-co(3,nodetop))**2)
               do j=0,nlayer-1
!
!                 bottom node
!
                  node=kon(indexe+nopeexp+j*nope+ibot)
                  dt=sqrt((co(1,nodebot)-co(1,node))**2+
     &                 (co(2,nodebot)-co(2,node))**2+
     &                 (co(3,nodebot)-co(3,node))**2)
                  ze=2.d0*dt/dd-1.d0
c                  write(*,*) 'map3dtolayer',node,xi,et,ze
!
!                 determining the value of the shape functions
!
                  if(lakon(nelem)(4:5).eq.'20') then
                     call shape20h(xi,et,ze,xl,xsj,shp,iflag)
                  elseif(lakon(nelem)(4:5).eq.'15') then
                     call shape15w(xi,et,ze,xl,xsj,shp,iflag)
                  endif
!
                  do k=1,nfield
                     yn(k,node)=0.d0
                     do l=1,nope
                        yn(k,node)=yn(k,node)+
     &                        shp(4,l)*yn(k,konl(l))
c                        write(*,*) 'xi',xi,et,ze,node
c                        write(*,*) l,shp(4,l),yn(k,konl(l))
                     enddo
                  enddo
!
!                 top node
!
                  node=kon(indexe+nopeexp+j*nope+itop)
                  dt=sqrt((co(1,nodebot)-co(1,node))**2+
     &                 (co(2,nodebot)-co(2,node))**2+
     &                 (co(3,nodebot)-co(3,node))**2)
                  ze=2.d0*dt/dd-1.d0
c                  write(*,*) 'map3dtolayer',node,xi,et,ze
!
!                 determining the value of the shape functions
!
                  if(lakon(nelem)(4:5).eq.'20') then
                     call shape20h(xi,et,ze,xl,xsj,shp,iflag)
                  elseif(lakon(nelem)(4:5).eq.'15') then
                     call shape15w(xi,et,ze,xl,xsj,shp,iflag)
                  endif
!
                  do k=1,nfield
                     yn(k,node)=0.d0
                     do l=1,nope
                        yn(k,node)=yn(k,node)+
     &                        shp(4,l)*yn(k,konl(l))
                     enddo
                  enddo
!
!                 middle node, if any
!
                  if(i.gt.nopes/2) cycle
!
                  node=kon(indexe+nopeexp+j*nope+imid)
                  dt=sqrt((co(1,nodebot)-co(1,node))**2+
     &                 (co(2,nodebot)-co(2,node))**2+
     &                 (co(3,nodebot)-co(3,node))**2)
                  ze=2.d0*dt/dd-1.d0
c                  write(*,*) 'map3dtolayer',node,xi,et,ze
!
!                 determining the value of the shape functions
!
                  if(lakon(nelem)(4:5).eq.'20') then
                     call shape20h(xi,et,ze,xl,xsj,shp,iflag)
                  elseif(lakon(nelem)(4:5).eq.'15') then
                     call shape15w(xi,et,ze,xl,xsj,shp,iflag)
                  endif
!
                  do k=1,nfield
                     yn(k,node)=0.d0
                     do l=1,nope
                        yn(k,node)=yn(k,node)+
     &                        shp(4,l)*yn(k,konl(l))
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo
!
      return
      end
