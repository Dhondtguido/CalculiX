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
!cd  incompressible for ASME nozzles eq 4a 4b   
!
!     author: Yannick Muller
!
      subroutine pk_cdi_noz(reynolds,cdi_noz)
!
      implicit none
!      
      real*8 reynolds,cdi_noz,ln_reynolds,cdi_noz_lr,
     &     cdi_noz_hr,e,reynolds_cor
!
      if (reynolds.lt.40000d0) then
!
! formerly pk_cdi_noz_lr : for low Reynolds nsumber
!         
         if (reynolds.eq.0d0) then
            reynolds_cor=1.d0
         else
            reynolds_cor=reynolds
         endif
         e=2.718281828459045d0
         ln_reynolds=log(reynolds_cor)/log(e)
!     
         cdi_noz_lr=0.19436d0+0.152884d0*ln_reynolds
     &        -0.0097785d0*ln_reynolds**2d0+0.00020903d0
     &        *ln_reynolds**3d0
!
         cdi_noz=cdi_noz_lr
!
      elseif (reynolds.lt.50000d0) then
!     
         if (reynolds.eq.0) then
            reynolds_cor=1
         else
            reynolds_cor=reynolds
         endif
!
         e=2.718281828459045d0
         ln_reynolds=log(reynolds_cor)/log(e)
!     
         cdi_noz_lr=0.19436d0+0.152884d0*ln_reynolds
     &        -0.0097785d0*ln_reynolds**2+0.00020903d0
     &        *ln_reynolds**3d0
!     
         cdi_noz_hr=0.9975d0-0.00653d0*dsqrt(1000000d0/50000d0)
         
!     linear interpolation in order to achieve continuity
!     
         cdi_noz=cdi_noz_lr+(cdi_noz_hr-cdi_noz_lr)
     &        *(reynolds-40000d0)/(50000d0-40000d0)
      else
!
!     formerly pk_cdi_noz_hr for high Reynolds numbers
!
         cdi_noz=0.9975d0-0.00653d0*dsqrt(1000000d0/reynolds)
      endif
      
      return
      end
