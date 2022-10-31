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
      subroutine anisomaxwavspd(elas,rho,iorth,wavspd)
!
!     Calculates the propagation wave speed in a material, up to its 21 
!     constants. Subroutine for calcmatwavsps.f
!
!     Based on the procedure in:
!     C. Lane. The Development of a 2D Ultrasonic Array Inspection 
!     for Single Crystal Turbine Blades.
!     Switzerland: Springer International Publishing, 2014.
!
!      CARLO MONJARAZ TEC (CMT)
!
!       INPUT:
!       
!       elas: double(21) - The elasticity vector, containing 21 entries. 
!             Non used are zero. If material is orthotropic, values are
!             rearranged to match indexes from anisotropic
!             material card.
!             
!        rho: double - Density of the material
!        
!       iorth: INTEGER - if the value is 1 : material is iorthtropic
!                        for other vaules:   material is anisotropic
!                        
!       OUTPUT:
!       
!       wavspd
!
      implicit none
!     
      integer i,j,k,im,imin,jm,jmin,iorth
!
      real*8 elas(21),c(3,3,3,3),rho,xi(-1:1,-1:1),et(-1:1,-1:1),
     &       wavspd,d1,distmin,a
!     
!     if orthotropic
!      
      if(iorth.eq.1)then
!     
         elas(10)=elas(9)
         elas(15)=elas(8)
         elas(21)=elas(7)
         elas(9)=0.d0
         elas(8)=0.d0
         elas(7)=0.d0
!     
      endif
!     
!     filling matrix c
!       
      call anisotropic(elas,c)
!
      d1=1.d0
!
      xi(0,0)=0.d0
      et(0,0)=0.d0
      call inversewavspd(xi(0,0),et(0,0),c,rho,a)
      distmin=a
      imin=0
      jmin=0
!
      do k=1,8
!
!     initialisation
!
         d1=d1/10.d0
!     
         do i=-1,1
            do j=-1,1
               if((i.eq.0).and.(j.eq.0)) cycle
!
               xi(i,j)=xi(0,0)+i*d1
               et(i,j)=et(0,0)+j*d1
!
!              check whether inside the (-1,1)x(-1,1) domain
!
               if((xi(i,j).le.1.d0).and.
     &              (xi(i,j).ge.-1.d0).and.
     &              (et(i,j).le.1.d0).and.
     &              (et(i,j).ge.-1.d0)) then
                  call inversewavspd(xi(i,j),et(i,j),c,rho,a)
!     
!                 checking for smallest initial distance
!     
                  if(a.lt.distmin) then
                     distmin=a
                     imin=i
                     jmin=j
                  endif
               endif
!
            enddo
         enddo
!     
!     minimizing the distance from the face to the node
!     
         do
!     
!     exit if minimum found
!     
            if((imin.eq.0).and.(jmin.eq.0)) exit
!
!           new center of 3x3 matrix
!
            xi(0,0)=xi(imin,jmin)
            et(0,0)=et(imin,jmin)
!
            im=imin
            jm=jmin
!
            imin=0
            jmin=0
!     
            do i=-1,1
               do j=-1,1
                  if((i+im.lt.-1).or.(i+im.gt.1).or.
     &                 (j+jm.lt.-1).or.(j+jm.gt.1)) then
!
                     xi(i,j)=xi(0,0)+i*d1
                     et(i,j)=et(0,0)+j*d1
!
!              check whether inside the (-1,1)x(-1,1) domain
!
                     if((xi(i,j).le.1.d0).and.
     &                  (xi(i,j).ge.-1.d0).and.
     &                  (et(i,j).le.1.d0).and.
     &                  (et(i,j).ge.-1.d0)) then
                        call inversewavspd(xi(i,j),et(i,j),c,rho,a)
!
!                       check for new minimum
!
                        if(a.lt.distmin) then
                           distmin=a
                           imin=i
                           jmin=j
                        endif
                     endif
!
                  endif
               enddo
            enddo
         enddo
      enddo
!
      call inversewavspd(xi(0,0),et(0,0),c,rho,a)
!
      wavspd=1.d0/a
!     
      return
      end
      
