!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine networkstiff(voldl,s,imat,konl,mi,ntmat_,shcon,
     &  nshcon,rhcon,nrhcon)
!
!     calculates the stiffness of a generic networkelement
!     element label: D + blank
!
      implicit none
!
      integer konl(20),mi(*),imat,nshcon(*),nrhcon(*),ntmat_
!
      real*8 voldl(0:mi(2),9),s(60,60),gastemp,shcon(0:3,ntmat_,*),
     &  cp,r,dvi,rhcon(0:1,ntmat_,*),rho
!
!
!
      gastemp=(voldl(0,1)+voldl(0,3))/2.d0
!
      call materialdata_tg(imat,ntmat_,gastemp,shcon,nshcon,cp,r,
     &  dvi,rhcon,nrhcon,rho)
!
      if(voldl(1,2).gt.0.d0) then
         s(3,1)=-cp*voldl(1,2)
         s(3,3)=-s(3,1)
      else
         s(1,1)=-cp*voldl(1,2)
         s(1,3)=-s(1,1)
      endif
!
      return
      end

