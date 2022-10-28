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
!
!   Subroutine pre_extrapolate.f
!
!      Interpolates xstate values for the new integration points 
!      at the beginning of the new increment. 
!
!   by: Jaro Hokkanen
!
!
      subroutine interpolatestate(ne,ipkon,kon,lakon,ne0,mi,xstate,
     &  pslavsurf,nstate_,xstateini,islavsurf,islavsurfold,pslavsurfold,
     &  tieset,ntie,itiefac)
!
      implicit none
!
      character*8 lakon(*),lakonl
      character*81 tieset(3,*)
!
      integer ipkon(*),kon(*),ne,i,n,mi(*),indexc,ne0,indexcj,
     &  nstate_,kk,nopespring,iface,ifacej,ielemslave,ll,
     &  numpts,islavsurf(2,*),islavsurfold(2,*),ntie,itiefac(2,*)
!
      real*8 xstate(nstate_,mi(1),*),pslavsurf(3,*),pslavsurfold(3,*),
     &  xstateini(nstate_,mi(1),*)
!
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.'C') cycle
         do kk=itiefac(1,i),itiefac(2,i)
            numpts=islavsurfold(2,kk+1)-islavsurfold(2,kk)
            if(numpts.gt.2) then
               call interpolateinface(kk,xstate,xstateini,numpts,
     &              nstate_,mi,islavsurf,pslavsurf,
     &              ne0,islavsurfold,pslavsurfold)
            endif
         enddo
      enddo
!     
      return
      end
      
