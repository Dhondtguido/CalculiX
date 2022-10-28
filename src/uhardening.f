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
      subroutine uhardening(amat,iel,iint,t1l,epini,ep,dtime,fiso,dfiso,
     &                      fkin,dfkin)
!
!     hardening user subroutine
!
!     INPUT:
!
!     amat:    material name (maximum 20 characters)
!     iel:     element number
!     iint:    integration point number
!     t1l:     temperature at the end of the increment
!     epini:   equivalent irreversible strain at the start 
!              of the increment
!     ep:      present equivalent irreversible strain
!     dtime:   time increment
!
!     OUTPUT:
!
!     fiso:    present isotropic hardening Von Mises stress
!     dfiso:   present isotropic hardening tangent (derivative
!              of the Von Mises stress with respect to the
!              equivalent irreversible strain)
!     fkin:    present kinematic hardening Von Mises stress
!     dfkin:   present kinematic hardening tangent (derivative
!              of the Von Mises stress with respect to the
!              equivalent irreversible strain)
!
      implicit none
!
      character*80 amat
      integer iel,iint
      real*8 t1l,epini,ep,dtime,fiso,dfiso,fkin,dfkin
!
      return
      end
