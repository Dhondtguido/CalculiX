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
!     this subroutines enables to calculate the critical pressure ratio of a straight 
!     labyrinth seal as a function of the number of spikes (n).
!
!     The following table is obtained by solving iteratively the equation :
!     Ps_inf/Pt0=ppkrit=1/dsqrt(1+2.n-ln(ppkrit))
!
!     this equation can be found by using the formula for the ideal mass flow in a straight labyrinth
!     see "Air system Correlations Part 1 : Labyrith Seals" H.Zimmermann and K.H. Wollf ASME98-GT-206
!     and determining the maximum flow for a given number of fin.
!
!     author: Yannick Muller
!
      subroutine lab_straight_ppkrit (n,ppkrit)
!
      implicit none
!
      integer n
!
      real*8 fppkrit(9),ppkrit
!
      data fppkrit
     &     /0.47113022d0,0.37968106d0,0.32930492d0,0.29569704d0,
     &      0.27105479d0,0.25191791d0,0.23646609d0,0.22363192d0,
     &      0.21274011d0/
!
      ppkrit=fppkrit(n)
!     
      return
      end
