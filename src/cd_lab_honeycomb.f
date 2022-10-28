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
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!     This subroutine enables to calculate the correction factor for a labyrinth seal
!     wit a honeycomb stator
!     s= gap, hl= width of a honeycomb cell
!     the correction factors are interpolated from a table 
!     
!     H.Zimmermann and K.h. Wolff
!     "Air system correlations part 1 Labyrinth seals"
!     asme 98-GT-206
!
!     author: Yannick Muller
!
      subroutine cd_lab_honeycomb(s,lc,cd_honeycomb)
!
      implicit none
!
      integer id,n11
!
      real*8 s,lc,cd_honeycomb,szlc
!
!     lc=1/8 inch
!     
      real*8 szl(11)
      data szl
     &     /0.05d0,0.06d0,0.075d0,0.081d0,0.1d0,0.13d0,0.15d0,0.16d0,
     &      0.20d0,0.30d0,0.40d0/
!
      real*8 deltamp(11)
      data deltamp
     &     /97.1d0,40d0,32d0,23d0,20d0,0d0,-3.3d0,-5.7d0,-8.5d0,
     &      -11.43d0,-12d0/
!
      data n11 /11/
!
! extrapolation
      szlc=s/lc
!      if (szlc.gt.0.40d0) then
!         cd_honeycomb=deltamp(11)
!      endif
!
!     intrapolation
!
          call ident(szl,szlc,n11,id)
!     call ident(yz,q,11,idy)
            if(id.eq.1) then
              cd_honeycomb=deltamp(1)
            elseif(id.eq.11) then
               cd_honeycomb=deltamp(11)
            else
               cd_honeycomb=deltamp(id)+(deltamp(id+1)-deltamp(id))
     &              *(szlc-szl(id))/(szl(id+1)-szl(id))
            endif
!
            return
!
            end
