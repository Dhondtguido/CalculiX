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
!     This subroutines enables to caclulate a correction term linked to with the radius 
!     of the spike as a function of r/s (radius/gap)
!     the parameter Hst ( height of the step ) enable to select either the table for a
!     straight labyrinth (Hst=0) or for a stepped labyrinth
!
!     H.Zimmermann and K.h. Wolff
!     "Air system correlations part 1 Labyrinth seals"
!     asme 98-GT-206
!
!     author: Yannick Muller
!
      subroutine cd_lab_radius(rad,s,hst,cd_radius)
!
      implicit none
!
      integer id,i,number,n9
!
      real*8 rad,s,cd_radius,rzs_tab(9),cd_sharp(9),rzs,hst
!
      real*8 rzs_tab1(9)
      data rzs_tab1
     &     /0d0,0.05d0,0.100d0,0.150d0,0.200d0,0.250d0,0.300d0,0.350d0,
     &      0.400d0/
!
      real*8 cd_sharp1(9)
      data cd_sharp1
     &     /1d0,1.025d0,1.10d0,1.11d0,1.12d0,1.125d0,1.126d0,1.127d0,
     &      1.127d0/
!
      real*8 rzs_tab2(9)
      data rzs_tab2
     &     /0d0,0.05d0,075d0,0.100d0,0.15d0,0.20d0,0.25d0,0.30d0,0.40d0/
!
      real*8 cd_sharp2(9)
      data cd_sharp2
     &     /1d0,1.10d0,1.15d0,1.20d0,1.26d0,1.31d0,1.34d0,1.36d0,1.37d0/
!
      data n9 /9/
!
      rzs=rad/s
!
!     straight labyrinth
!
      if(hst.eq.0d0) then
         call ident(rzs_tab1,rzs,n9,id)
         number=9
         do i=1,9
            rzs_tab(i)=rzs_tab1(i)
            cd_sharp(i)=cd_sharp1(i)
         enddo
!
!     stepped labyrinth
!
      else
         call ident(rzs_tab2,rzs,n9,id)
         number=9
         do i=1,9
            rzs_tab(i)=rzs_tab2(i)
            cd_sharp(i)=cd_sharp2(i)
         enddo
      endif
!
!     linear interpolation
!
!    
            if(id.eq.1) then
              cd_radius=cd_sharp(1)
            elseif(id.eq.number) then
               cd_radius=cd_sharp(number)
            else
               cd_radius=cd_sharp(id)+(cd_sharp(id+1)-cd_sharp(id))
     &              *(rzs-rzs_tab(id))/(rzs_tab(id+1)-rzs_tab(id))
            endif
!
            return
!
            end
