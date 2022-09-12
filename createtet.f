!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998 Guido Dhondt
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
      subroutine createtet(kontet,ifatet,ielement,inodfa,
     &  ifreefa,planfa,ipofa,nodes,cotet,iparentelement)
!
      implicit none
!
      integer nodes(4),nodef(3),kontet(4,*),ifatet(4,*),inodfa(4,*),
     &  ipofa(*),ifreetet,ifreefa,ifree,ig(3,4),j,
     &  n1,n2,n3,n4,nxa,nxb,nxc,nxd,nya,nyb,nyc,nyd,nza,nzb,nzc,
     &  nzd,nx1,nx2,ny1,ny2,nz1,nz2,index,j1,j2,j3,i,n,
     &  node,indexold,ielement,ndx,ndy,ndz,iparentelement
!
      real*8 planfa(4,*),cotet(3,*),dd
!
      data ig /2,3,4,3,4,1,4,1,2,1,2,3/
!
!     updating the node per element relationship
!
      do i=1,4
         kontet(i,ielement)=nodes(i)
      enddo
!
!     creating faces
!
      do i=1,4
         nodef(1)=nodes(ig(1,i))
         nodef(2)=nodes(ig(2,i))
         nodef(3)=nodes(ig(3,i))
!
         n=3
         call insertsorti(nodef,n)
c         call isortii(nodef,idum,n,kflag)
!
!        check whether face already exists
!
         node=nodef(1)
         index=ipofa(node)
!
         do
            if(index.eq.0) exit
            if((inodfa(2,index).eq.nodef(2)).and.
     &         (inodfa(3,index).eq.nodef(3))) exit
            indexold=index
            index=inodfa(4,index)
         enddo
!
         if(index.eq.0) then
            index=ifreefa
            ifreefa=inodfa(4,ifreefa)
            if(ifreefa.eq.0) then
               write(*,*) '*ERROR in generatet: increase the dimension'
               write(*,*) '       of inodfa'
            endif
            inodfa(1,index)=nodef(1)
            inodfa(2,index)=nodef(2)
            inodfa(3,index)=nodef(3)
            inodfa(4,index)=0
            if(ipofa(node).eq.0) then
               ipofa(node)=index
            else
               inodfa(4,indexold)=index
            endif
!
            call planeeq(cotet,nodef,planfa(1,index))
!
         endif
!
!        the face number in ifatet is negative, if the equation
!        of the face plane is such, that its value in the 
!        remaining node of the tetrahedron is negative
!
         dd=planfa(1,index)*cotet(1,nodes(i))+
     &      planfa(2,index)*cotet(2,nodes(i))+
     &      planfa(3,index)*cotet(3,nodes(i))+
     &      planfa(4,index)
         if(dabs(dd).lt.1.d-10) then
            write(*,*) '*WARNING in createtet: element ',
     &                     iparentelement
            write(*,*) '         is extremely flat'
            write(*,*) '         the element is deleted'
            ielement=ielement-1
            return
         endif
         if(dd.ge.0.d0) then
            ifatet(i,ielement)=index
         else
            ifatet(i,ielement)=-index
         endif
      enddo
!
      return
      end
