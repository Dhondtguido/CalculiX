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
      subroutine flowbc(ntg,itg,cam,vold,v,nload,sideload,
     &     nelemload,xloadact,nactdog,network,mi,ne,ipkon,lakon,kon)
!     
      implicit none
!    
      character*8 lakon(*),lakonl
      character*20 sideload(*) 
!     
      integer i,j,nload,node,ntg,itg(*),nelemload(2,*),kon(*),
     &     nactdog(0:3,*),network,mi(*),ne,indexe,ipkon(*),
     &     node1,node2,node3
!     
      real*8 cam(5),vold(0:mi(2),*),v(0:mi(2),*),xloadact(2,*)
!     
!     calculating the change of gas temperature: is taken
!     into account in the global convergence for purely
!     thermal networks (reason: for purely thermal networks
!     the network solution is not iterated to speed up
!     the calculation)
!
      if(network.le.1) then
         do i=1,ntg
            node=itg(i)
            if(nactdog(0,node).eq.0) cycle
            if(dabs(vold(0,node)-v(0,node)).gt.cam(2)) then
               cam(2)=dabs(vold(0,node)-v(0,node))
               cam(5)=node+0.5d0
            endif
         enddo
      endif
!     
!     replacing vold by v (including the static temperature for
!     gases and the critical depth for liquid channels)
!
      do i=1,ntg
         node=itg(i)
         do j=0,3
            vold(j,node)=v(j,node)
         enddo
      enddo
!
!     determining all inflowing mass flow in the end nodes and
!     assigning it to the end nodes
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         lakonl=lakon(i)
         if(lakonl(1:1).ne.'D') cycle
         if(lakonl(1:7).eq.'DCOUP3D') cycle
!
         indexe=ipkon(i)
         if(kon(indexe+1).ne.0)  then
            node1=kon(indexe+1)
            vold(1,node1)=0.d0
         endif
         if(kon(indexe+3).ne.0) then
            node3=kon(indexe+3)
            vold(1,node3)=0.d0
         endif
      enddo
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         lakonl=lakon(i)
         if(lakonl(1:1).ne.'D') cycle
         if(lakonl(1:7).eq.'DCOUP3D') cycle
!
         indexe=ipkon(i)
         node2=kon(indexe+2)
         if(kon(indexe+1).ne.0)  then
            node1=kon(indexe+1)
            if(vold(1,node2).lt.0.d0) 
     &         vold(1,node1)=vold(1,node1)-vold(1,node2)
         endif
         if(kon(indexe+3).ne.0) then
            node3=kon(indexe+3)
            if(vold(1,node2).gt.0.d0) 
     &         vold(1,node3)=vold(1,node3)+vold(1,node2)
         endif
      enddo
!     
!     updating the film boundary conditions
!     
      do i=1,nload
         if(sideload(i)(3:4).eq.'FC') then
            node=nelemload(2,i)
            xloadact(2,i)=vold(0,node)
         endif
      enddo
!     
!     updating the pressure boundary conditions
!     
      do i=1,nload
         if(sideload(i)(3:4).eq.'NP') then
            node=nelemload(2,i)
            xloadact(1,i)=vold(2,node)
         endif
      enddo
!      
      return
      end








