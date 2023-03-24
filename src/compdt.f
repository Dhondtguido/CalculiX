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
      subroutine compdt(nk,dt,nshcon,shcon,vold,ntmat_,
     &     iponoel,inoel,dtimef,ielmat,dh,cocon,
     &     ncocon,ithermal,mi,
     &     vcon,compressible,tincf,ierr,ifreesurface,dgravity,iit)
!     
!     - determine the time step for each node (stored in field dt
!     and the minimum value across all nodes (dtimef)
!     
      implicit none
!     
      integer nk,i,iponoel(*),inoel(2,*),index,nelem,ithermal(*),
     &     mi(*),compressible,ierr,ifreesurface,iit,
     &     nshcon(*),ntmat_,ielmat(mi(3),*),imat,ncocon(2,*),
     &     iflag
!     
      real*8 dtimef,dt(*),dvi,r,cp,rho,shcon(0:3,ntmat_,*),tincf,
     &     vold(0:mi(2),*),temp,vel,dtcon,dtmed,dtimefold,
     &     dh(*),cocon(0:6,ntmat_,*),dtthd,cond,c,
     &     vcon(nk,0:mi(2)),dgravity
!     
      data iflag /3/
!     
!     determining the time increment dt for each node.
!     
!     edge nodes (fields iponoel and inoel are determined in precfd.f)
!
      dtimefold=dtimef
      dtimef=1.d30
!     
      do i=1,nk
        index=iponoel(i)
        if(index.le.0) cycle
!     
!     look at an element belonging to the edge node
!     
        nelem=inoel(1,index)
!     
!     determining the time increment
!     
        imat=ielmat(1,nelem)
        temp=vold(0,i)
!     
!     convective time step (dtcon)
!     
        vel=dsqrt(vold(1,i)**2+vold(2,i)**2+vold(3,i)**2)
!     
        if(compressible.eq.1) then
          if(ifreesurface.eq.0) then
!     
!     gas: speed of sound = dsqrt(kappa*r*temp)
!     
            rho=vcon(i,4)
            call materialdata_cp(imat,ntmat_,temp,shcon,nshcon,cp)
            r=shcon(3,1,imat)
            dtcon=dh(i)/(dsqrt(cp*r*temp/(cp-r))+vel)
            c=cp*r*temp/(cp-r)
          else
!     
!     shallow water: speed of sound = dsqrt(g*depth)
!     
            dtcon=dh(i)/(dsqrt(dgravity*vcon(i,4))+vel)
          endif
        else
!     
!     liquid
!     
          rho=vcon(i,4)
          if(vel.lt.1.d-10) vel=1.d-10
          dtcon=dh(i)/vel
        endif
        dt(i)=dtcon
!     
!     mechanical diffusion time step (dtmed)
!     
        call materialdata_dvifem(imat,ntmat_,temp,shcon,nshcon,dvi)
        if(dvi.gt.1.d-20) then
          dtmed=dh(i)*dh(i)*rho/(2.d0*dvi)
          dt(i)=dtcon*dtmed/(dtcon+dtmed)
        endif
!     
!     thermal diffusion time step (dtthd)
!     
        if(ithermal(1).gt.1) then
          call materialdata_cond(imat,ntmat_,temp,cocon,ncocon,
     &         cond)
          call materialdata_cp(imat,ntmat_,temp,shcon,nshcon,cp)
          if(cond.gt.1.d-20) then
            dtthd=dh(i)*dh(i)*rho*cp/(2.d0*cond)
            dt(i)=(dt(i)*dtthd)/(dt(i)+dtthd)
          endif
        endif
!     
!     safety factor (default for compressible fluids: 1.25,
!     default for incompressible fluids: 1.00;
!     can be changed by user: fifth entry underneath *CFD)
!     
        dt(i)=dt(i)/tincf
c        dt(i)=dt(i)/10.d0
!     
        if(dt(i).lt.dtimef) then
          dtimef=dt(i)
        endif
!     
      enddo
!     
!     increased damping for incompressible fluids
!     the use of an internal (for the damping) and an external
!     (for the time derivative) time step stems from Zienkiewicz,
!     Taylor and Nithiarasu, The Finite Element Method for Fluid
!     Dynamics, 6th edition, p94 bottom and p 95 top.
!     
      write(*,*)
      write(*,*) 'dtimef= ',dtimef
      write(*,*)
      if(dtimef.le.0.d0) then
        write(*,*) '*ERROR in compdt;'
        write(*,*) '       negative time increment;'
        write(*,*) '       the solution diverged'
        ierr=1
      elseif(((iit.gt.1).and.(dtimef.lt.dtimefold/100.d0))) then
        write(*,*) '*ERROR in compdt;'
        write(*,*) '       strongly decreasing time increment;'
        write(*,*) '       the solution diverged'
        ierr=1
      endif
!     
      return
      end
