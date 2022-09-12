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
      subroutine phys2con(inomat,vold,ntmat_,shcon,nshcon,physcon,
     &     compressible,vcon,rhcon,nrhcon,ithermal,mi,ifreesurface,ierr,
     &     dgravity,depth,nk,nka,nkb)
!
!     calculates the conservative variables from the physical variables
!     only for vcon(0...4), i.e. rho*epsilon,rho*vx,rho*vy,rho*vz,rho
!
!     NOT for the turbulent parameters      
!      
      implicit none
!
      integer inomat(*),node,ntmat_,nshcon(*),compressible,nrhcon(*),
     &  ithermal(*),mi(*),imat,k,ifreesurface,ierr,nk,nka,nkb
!
      real*8 rhcon(0:1,ntmat_,*),shcon(0:3,ntmat_,*),vold(0:mi(2),*),
     &  vcon(nk,0:mi(2)),physcon(*),temp,cp,r,rho,dgravity,depth(*)
!
      do node=nka,nkb
        imat=inomat(node)
        temp=vold(0,node)
        call materialdata_cp_sec(imat,ntmat_,temp,shcon,nshcon,
     &       cp,physcon)
!     
        if(compressible.eq.1) then
!     
!     compressible calculations (gas or shallow water)
!     
          if(ifreesurface.eq.0) then
            r=shcon(3,1,imat)
            rho=vold(4,node)/(r*(vold(0,node)-physcon(1)))
            vcon(node,0)=rho*(cp*(temp-physcon(1))+
     &           (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &           /2.d0)-vold(4,node)
          else
!     
!     shallow water equations
!     
            rho=2.d0*vold(4,node)/dgravity+depth(node)**2
            if(rho.le.0.d0) then
              write(*,*) '*ERROR in phys2con: fluid depth cannot'
              write(*,*) '       be determined'
              ierr=1
              return
            else
              rho=dsqrt(rho)
            endif
            if(ithermal(1).gt.1) then
              vcon(node,0)=rho*(cp*(temp-physcon(1))+
     &             (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &             /2.d0)
            endif
          endif
          vcon(node,4)=rho
        else
!     
!     incompressible calculations (liquid) 
!     
          if(ithermal(1).gt.1) then
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
            vcon(node,0)=rho*(cp*(temp-physcon(1))+
     &           (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &           /2.d0)
            vcon(node,4)=rho
          else
            rho=vcon(node,4)
          endif
        endif
!     
        do k=1,3
          vcon(node,k)=rho*vold(k,node)
        enddo
      enddo
!     
      return
      end
      
