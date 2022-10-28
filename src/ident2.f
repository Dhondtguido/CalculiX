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
!     identifies the position id of px in an ordered array
!     x of real numbers; The numbers in x are at positions
!     1, 1+ninc, 1+2*ninc, 1+3*ninc... up to 1+(n-1)*ninc
!  
!     id is such that x(id).le.px and x(id+1).gt.px
!                                                                             
      subroutine ident2(x,px,n,ninc,id)
      implicit none   
! 
      integer n,id,n2,m,ninc
!
      real*8 x(n*ninc),px
!
!
!
      id=0
      if(n.eq.0) return
      n2=n+1
      do                                                            
         m=(n2+id)/2
         if(px.ge.x(1+ninc*(m-1))) then
            id=m      
         else
            n2=m            
         endif
         if((n2-id).eq.1) return
      enddo
      end
                                                                               
