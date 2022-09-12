!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine checkimpacts(ne,neini,temax,sizemaxinc,
     & energyref,tmin,tmax,tper,idivergence,
     & iforceincsize,istab,dtheta,r_abs,energy,energyini,allwk,
     & allwkini,dampwk,dampwkini,emax,mortar,maxdecay,enetoll)
!     
!     # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
!     Routine that contains the implementation of the logic to
!       rule the increment size during contact conditions.
!     The values of tolerances have been tested for the ball
!       model, the two sliding blades (simplified model of 
!       blade from Mr. Wittig) and the real blade model.
!     Friction has not been tested deeply.
! 
!     Main variables and meaning
!
!     sizemaxinc    : maximum size of the current increment
!     iforceincsize : flag to set "dtheta=sizemaxinc" in 
!                       checkconvergence.
!     cstate        : vector containing contact data
!     temax            : max. natural period of oscillation 
!     icase         : flag to print debug informations
!     emax          : maximum energy of the system over the ti-
!                       me history
!     r         : energy residual before contact (or re-
!                       adapted)
!     delta         : eneres normalized
!     fact          : factor to set sizemaxinc according to the  
!                       contact formulation
!     stab_th       : \hat{r}_{e}(t_n) -> mod belytschko before
!                       contact (or initial value). This is the 
!                       stability threshold
!     delta_r_rel   : \hat{r}_{e}(t) - \hat{r}_{e}(t-1) -> varia-
!                       tion of the modified belitschko criterion 
!                       used to control jumps
!     r_rel         : \hat{r}_{e}(t) actual value of the modified  
!                       belitschko criterion
! 
!     Proposed by Matteo Pacher
! 
!     # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
!     
      implicit none
!     
      integer idivergence,
     & iforceincsize,ne,neini,istab,mortar
!     
      real*8 temax,energyref,sizemaxinc,tmin,tmax,tper,dtheta,
     & delta_r_rel,r_rel,delta,allwk,allwkini,energy(4),
     & energyini(4),dampwk,dampwkini,emax,fact,
     & r_rel_bc,maxdecay,r_abs,enetoll,delta_r_abs
!     
!     Initialization
!     
      iforceincsize=0
!     
!     Adaption of the energy residual (absolute/relative check)
!     
      if(dabs(r_abs).lt.(enetoll/4.d0))then
         delta=r_abs*emax
      else
         delta=r_abs
      endif
!     
      if(mortar.eq.0)then
         fact=10.d0
      elseif(mortar.eq.1) then
         fact=1.d0
      endif
!     
!     Compute thresholds and energy values
!     
      delta_r_abs=energy(1)+energy(2)+energy(3)+energy(4)-allwk-
     &     dampwk-(energyini(1)+energyini(2)+energyini(3)+
     &     energyini(4)-allwkini-dampwkini)
!
      if(emax.le.0.d0)then
!     
!     No kinetic energy in the structure: energyref is the internal energy
!     this happens at the beginning of the calculation
!     
         r_rel=(energy(1)+energy(2)+energy(3)+energy(4)-allwk-
     &        dampwk-energyref)/energyref
         delta_r_rel=delta_r_abs/energyref
         r_rel_bc=delta/energyref
      else
         r_rel=(energy(1)+energy(2)+energy(3)+energy(4)-allwk-
     &        dampwk-energyref)/emax
         delta_r_rel=delta_r_abs/emax
         r_rel_bc=delta/emax
      endif
!     
!     Logic to adapt the increment size
!     
      if(mortar.eq.0)then
!     
!     Energy conservation rules for NODE TO SURFACE penalty contact
!     
         if((delta_r_rel.lt.(-0.008d0)).and.(ne.ge.neini))then
!     
!     Impact (or too high variation during pers. contact)
!     delta_r_rel = r_rel-r_rel_ini
!     
            idivergence=1
            sizemaxinc=dtheta*0.25d0
            iforceincsize=1
         elseif((r_rel-r_rel_bc.gt.0.0025d0).and.(ne.le.neini))then
!     
!     Rebound (or too high variation during pers. contact)
!     r_rel_bc is r_rel before contact
!     
            idivergence=1
            sizemaxinc=dtheta*0.5d0
            iforceincsize=1
         else
!     
!     Persistent Contact
!     
            if(r_rel.gt.(-0.9d0*maxdecay))then
               sizemaxinc=max(fact*temax/tper,1.01d0*dtheta)
               sizemaxinc=min(sizemaxinc,100.d0*temax/tper)
            else
               sizemaxinc=max(temax/tper/10.d0,0.5d0*dtheta)
               istab=1
            endif
!            
         endif
!     
      elseif(mortar.eq.1)then
!     
!     Energy conservation rules for SURFACE TO SURFACE penalty contact
!     
         if((delta_r_rel.lt.(-0.008d0)).and.(ne.ge.neini))then
!     
!     Impact (or too high variation during pers. contact)
!     delta_r_rel = r_rel-r_rel_ini
!     
            idivergence=1
            sizemaxinc=dtheta*0.25d0
            iforceincsize=1
!     
         elseif((r_rel-r_rel_bc.gt.0.0025d0).and.(ne.le.neini))then     
!     
!     Rebound (or too high variation during pers. contact)
!     r_rel_bc is r_rel before contact
!     
            idivergence=1
            sizemaxinc=dtheta*0.5d0
            iforceincsize=1
!     
         else
!     
!     Persistent Contact
!     
            if(r_rel.gt.(-0.9d0*maxdecay))then
               sizemaxinc=min(fact*temax/tper,1.1d0*dtheta)
               sizemaxinc=min(sizemaxinc,100.d0*temax/tper)
            else
               sizemaxinc=max(temax/tper/10.d0,0.5d0*dtheta)
               istab=1
            endif
         endif
      endif                     !(mortar)
!     
      if(sizemaxinc.lt.tmin)then
        sizemaxinc=tmin
      elseif(sizemaxinc.gt.tmax) then
        sizemaxinc=tmax
      endif
!
      write(*,*) '*INFO in checkimpacts: due to impact rules the'
      write(*,*) '      maximum allowed time increment has been'
      write(*,*) '      changed to',sizemaxinc*tper
!     
      return
      end
