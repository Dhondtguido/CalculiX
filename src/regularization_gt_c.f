!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
!     regularization function for tangential contact mortar
!     (old, only for stressmortar)
!
!  [in]	lambdatt	lambdatilde_tau=lambda_tau-bar{lambda}_tau
!  [in]	divmode 	indicates whether funtion or derivate 
!                             	should be called
!                    		=0 function called
!                    		=1 derivative called    
!  [in]	regmode        	selects regularization funtion
!                    		=1 perturbed Lagrange
!  [out] gtc	        result regularization function
!  [in]  atauinvloc      stiffness constant for perturbed Lagrange
!
      subroutine regularization_gt_c(lambdatt,divmode,regmode,
     &     gtc,atauinvloc)
!     
!     regularization function for tangential contact
!     Author: Saskia Sitzmann
!     
      implicit none
!     
      integer divmode,regmode,kode
!
      real*8 lambdatt(*),gtc(*),atauinvloc,t1l
!
!
!
      kode=-51
      t1l=0.0
!
!     perturbed Lagrange
!     
      if(regmode.eq.1)then
         if(divmode.eq.0)then
            gtc(1)=atauinvloc*lambdatt(1)
            gtc(2)=atauinvloc*lambdatt(2)
         elseif(divmode.eq.1)then
            gtc(1)=atauinvloc
            gtc(2)=atauinvloc
         else
            write(*,*)'error in regularzation_gt_c.f!'
            call exit(201)
         endif
      else
         gtc(1)=0.0
         gtc(2)=0.0
      endif
!      
      return
      end
