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
      subroutine presgradient(iponoel,inoel,sa,shockcoef,
     &     dtimef,ipkon,kon,lakon,vold,mi,
     &     nactdoh,nka,nkb)
!     
!     determining measure for the pressure gradient
!     
!     Ref: The Finite Element Method for Fluid Dynamics,
!     O.C. Zienkiewicz, R.L. Taylor & P. Nithiarasu
!     6th edition (2006) ISBN 0 7506 6322 7
!     p. 61
!     
      implicit none
!     
      character*8 lakon(*)
!     
      integer iponoel(*),inoel(2,*),i,j,k,index,indexe,nope,
     &     ipkon(*),kon(*),node,ielem,mi(*),nka,nkb,
     &     nactdoh(*)
!     
      real*8 sa(*),shockcoef,dtimef,ca,sum,pa,
     &     vold(0:mi(2),*),sumabs,contribution
!     
      do i=nka,nkb
        if(nactdoh(i).le.0) cycle
        if(iponoel(i).le.0) cycle
        j=nactdoh(i)
!        
        sum=0.d0
        sumabs=0.d0
        pa=vold(4,i)
        index=iponoel(i)
!     
        do
          ielem=inoel(1,index)
          if(ipkon(ielem).lt.0) cycle
          if(lakon(ielem)(1:1).ne.'F') cycle
          if(lakon(ielem)(4:4).eq.'8') then
            nope=8
          elseif(lakon(ielem)(4:4).eq.'4') then
            nope=4
          elseif(lakon(ielem)(4:4).eq.'6') then
            nope=6
          endif
          indexe=ipkon(ielem)
          do k=1,nope
            node=kon(indexe+k)
            if(node.eq.i) cycle
!
!           connecting vector between node i and "node"
!
            contribution=(pa-vold(4,node))
            sum=sum+contribution
            sumabs=sumabs+dabs(contribution)
!            
          enddo
          index=inoel(2,index)
          if(index.eq.0) exit
        enddo
        if(sumabs.lt.1.d-10) then
          sum=0.d0
          sumabs=1.d0
        endif
        sa(j)=dabs(sum)/(sumabs*dtimef)
      enddo
!
      ca=shockcoef*dtimef
      do i=nka,nkb
        sa(i)=ca*sa(i)
      enddo
!     
      return
      end
