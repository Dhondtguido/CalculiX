!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine con2phys(vold,vcon,nk,ntmat_,shcon,nshcon,rhcon,nrhcon,
     &     physcon,ithermal,compressible,iturbulent,inomat,mi,
     &     ierr,ifreesurface,dgravity,depth,nka,nkb)
!     
!     calculates the physical variable from the conservative
!     variables:
!     for gases: temperature,velocity and pressure 
!     for thermal liquids: temperature, velocity and density
!     for athermal liquids: velocity and density
!     
      implicit none
!     
      integer compressible,ifreesurface,nrhcon(*),ntmat_,iturbulent,
     &     mi(*),nshcon(*),nk,ithermal(*),i,j,k,imat,ierr,inomat(*),
     &     nka,nkb
!     
      real*8 vold(0:mi(2),*),vcon(nk,0:mi(2)),rhcon(0:1,ntmat_,*),rho,
     &     c1,cp,r,temp,temp0,c2,shcon(0:3,ntmat_,*),physcon(*),
     &     dgravity,depth(*)
!     
      if(ithermal(1).gt.1) then
!     
        do i=nka,nkb
          imat=inomat(i)
          temp=vold(0,i)
!     
          if(compressible.eq.1) then
!     
!     gas : rho*epsilon_t, rho and rho*v known
!     
!     calculation of the static temperature
!     
            rho=vcon(i,4)
            r=shcon(3,1,imat)
            c1=(vcon(i,0)-(vcon(i,1)**2+vcon(i,2)**2+
     &           vcon(i,3)**2)/(2.d0*rho))/rho
!     
!     the energy density per volume unit cannot be
!     negative
!     
            if(c1.lt.1.d-10) c1=1.d-10
!     
            temp0=temp
            j=0
            do
              call materialdata_cp_sec(imat,ntmat_,temp,shcon,
     &             nshcon,cp,physcon)
              temp=c1/(cp-r)+physcon(1)
              j=j+1
              if((dabs(temp-temp0).lt.1.d-4*dabs(temp)).or.
     &             (dabs(temp-temp0).lt.1.d-10)) then
                vold(0,i)=temp
                exit
              endif
              if(j.gt.100) then
                write(*,*) 
     &               '*ERROR in con2phys: too many iterations'
                write(*,*) '       for node',i
                write(*,*) '       increment is recalculated'
                write(*,*) '       with a higher shock smoothing'
                ierr=1
                return
              endif
              temp0=temp
            enddo
!     
            if(ifreesurface.eq.0) then
!     
!     determining the pressure (gas equation)
!     
              vold(4,i)=rho*r*(temp-physcon(1))
            else
!     
!     determining the pressure (shallow water)
!     
              vold(4,i)=dgravity*(vcon(i,4)**2-depth(i)**2)/2.d0
            endif
!     
!     calculating the velocity from the conservative variables
!     
            do k=1,3
              vold(k,i)=vcon(i,k)/rho
            enddo
!     
!     calculating the turbulent quantities from the conservative
!     variables
!     
            if(iturbulent.ne.0) then
              do k=5,6
                vold(k,i)=vcon(i,k)/rho
              enddo
            endif
!     
            cycle
          else
!     
!     thermal liquid
!     
            c1=vcon(i,0)
            c2=(vcon(i,1)**2+vcon(i,2)**2+vcon(i,3)**2)/2.d0
            temp0=temp
            j=0
!     
!     iterating to find the temperature
!     
            do
              call materialdata_cp_sec(imat,ntmat_,temp,shcon,
     &             nshcon,cp,physcon)
              call materialdata_rho(rhcon,nrhcon,imat,rho,
     &             temp,ntmat_,ithermal)
              temp=(c1-c2/rho)/(rho*cp)+physcon(1)
              j=j+1
              if((dabs(temp-temp0).lt.1.d-4*dabs(temp)).or.
     &             (dabs(temp-temp0).lt.1.d-10)) then
                vold(0,i)=temp
                exit
              endif
              if(j.gt.100) then
                write(*,*) 
     &               '*ERROR in con2phys: too many iterations'
                write(*,*) '       for node',i
                write(*,*) '       actual temperature ',temp,' K'
                stop
              endif
              temp0=temp
            enddo
!     
!     storing the density
!     
            vcon(i,4)=rho
!     
!     calculating the velocity
!     
            do k=1,3
              vold(k,i)=vcon(i,k)/rho
            enddo
!     
!     calculating the turbulent quantities
!     
            if(iturbulent.ne.0) then
              do k=5,6
                vold(k,i)=vcon(i,k)/rho
              enddo
            endif
!     
          endif
        enddo
      else
!     
!     athermal calculations
!     
        do i=nka,nkb
!     
!     calculating the fictitious pressure for shallow water
!     calculations
!     
          if(ifreesurface.eq.1) then
            vold(4,i)=dgravity*(vcon(i,4)**2-depth(i)**2)/2.d0
          endif
!     
!     calculating the velocity
!     
          rho=vcon(i,4)
          do k=1,3
            vold(k,i)=vcon(i,k)/rho
          enddo
!     
!     calculating the turbulent quantities
!     
          if(iturbulent.ne.0) then
            do k=5,6
              vold(k,i)=vcon(i,k)/rho
            enddo
          endif
!     
        enddo
      endif
!     
      return
      end
