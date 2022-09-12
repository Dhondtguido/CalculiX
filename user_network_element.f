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
      subroutine user_network_element(node1,node2,nodem,nelem,lakon,kon,
     &     ipkon,nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,dvi,numf,set,co,vold,mi,ttime,
     &     time,iaxial,iplausi)
!     
!     user network elements
!
      implicit none
!     
      logical identity
      character*8 lakon(*)
      character*81 set(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(*),idirf(*),kflag,ipkon(*),kon(*),
     &     iaxial,mi(*),iplausi
!     
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(*),R,cp,physcon(*),dvi,
     &     co(3,*),vold(0:mi(2),*),ttime,time
!
!     list of different user network elements
!
!     notice that the input deck is converted into upper case when
!     being read by CalculiX. So even if the user has specified "p1"
!     in his input deck, at the present stage "P1" is stored.
!
      if((lakon(nelem)(3:4).eq.'P0').or.
     &   (lakon(nelem)(3:4).eq.'0 ')) then
!
!        this just contains a skeleton file
!
         call user_network_element_p0(node1,node2,nodem,nelem,lakon,kon,
     &     ipkon,nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,dvi,numf,set,co,vold,mi,ttime,
     &     time,iaxial,iplausi)
      elseif((lakon(nelem)(3:4).eq.'P1').or.
     &   (lakon(nelem)(3:4).eq.'1 ')) then
!
!        this is a working example
!
         call user_network_element_p1(node1,node2,nodem,nelem,lakon,kon,
     &     ipkon,nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,dvi,numf,set,co,vold,mi,ttime,
     &     time,iaxial,iplausi)
      endif
!     
      return
      end
