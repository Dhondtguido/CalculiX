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
!
!     this subroutine enables to calculate the discharge coefcieint of  
!     a labyrinth with one fin as a function of the ratio b/s and the 
!     pressure ratio Pdownstream/Pupstream
!     the results are interpolated 
!     
!     the relevant data can be found in:
!     "Air System Correlations Part 1: Labyrinth Seals"
!     H.Zimmermann and K.H. Wollf
!     ASME 98-GT-206
!     fig 11 p 7
!
!     author: Yannick Muller
!
      subroutine cd_lab_1spike(pt0zps1,s,b,cd_1spike)
!
      implicit none
! 
      integer nx,ny,idx,idy
!     
      real*8 pt0zps1,s,b,cd_1spike,z1,z2,z3,z4,xi,et,pdszpus,bzs
!
      real*8 Pdszpus_tab(7)
      data pdszpus_tab
     &    /0.400d0,0.500d0,0.555d0,0.625d0,0.714d0,0.833d0,1.000d0/
!
      real*8 bzs_tab(9)
      data bzs_tab
     &    /0.250d0,0.285d0,0.330d0,0.400d0,0.5000d0,0.660d0,1d0,2d0,4d0/
!     
      real*8 cd_1spike_tab(7,9)
      data cd_1spike_tab
     &     /0.930d0,0.875d0,0.830d0,0.790d0,0.750d0,0.700d0,0.650d0,
     &      0.930d0,0.875d0,0.830d0,0.800d0,0.750d0,0.710d0,0.660d0,
     &      0.930d0,0.875d0,0.830d0,0.800d0,0.750d0,0.710d0,0.660d0,
     &      0.918d0,0.875d0,0.830d0,0.800d0,0.750d0,0.710d0,0.670d0,
     &      0.912d0,0.875d0,0.830d0,0.800d0,0.750d0,0.710d0,0.675d0,
     &      0.900d0,0.875d0,0.830d0,0.800d0,0.750d0,0.710d0,0.687d0,
     &      0.900d0,0.875d0,0.830d0,0.800d0,0.750d0,0.725d0,0.687d0,
     &      0.912d0,0.875d0,0.862d0,0.837d0,0.800d0,0.785d0,0.743d0,
     &      0.912d0,0.880d0,0.870d0,0.860d0,0.860d0,0.855d0,0.850d0/
      bzs=b/s
      pdszpus=1/pt0zps1
      nx=7
      ny=9
!
      call ident(pdszpus_tab,pdszpus,nx,idx)
      call ident(bzs_tab,bzs,ny,idy)
!     
      if (idx.eq.0) then
         if(idy.eq.0) then
            cd_1spike=cd_1spike_tab(1,1)
         else
            if(idy.eq.ny) then
               cd_1spike=cd_1spike_tab(1,ny)
            else
               cd_1spike=cd_1spike_tab(1,idy)+(cd_1spike_tab(1,idy+1)
     &              -cd_1spike_tab(1,idy))
     &              *(bzs-bzs_tab(idy))/(bzs_tab(idy+1)-bzs_tab(idy))
            endif 
         endif
!     
      elseif(idx.ge.nx) then
         if(idy.le.0) then
            cd_1spike=cd_1spike_tab(nx,1)
         else
            if(idy.ge.ny) then
               cd_1spike=cd_1spike_tab(nx,ny)
            else
               cd_1spike=cd_1spike_tab(nx,idy)+
     &              (cd_1spike_tab(nx,idy+1)-cd_1spike_tab(nx,idy))
     &              *(bzs-bzs_tab(idy))/(bzs_tab(idy+1)-bzs_tab(idy))
            endif 
         endif
      else
         if(idy.le.0) then
!     
            cd_1spike=cd_1spike_tab(idx,1)+(cd_1spike_tab(idx+1,1)
     &           -cd_1spike_tab(idx,1))
     &           *(pdszpus-pdszpus_tab(idx))/(pdszpus_tab(idx+1)
     &           -pdszpus_tab(idx))
         elseif(idy.ge.ny) then
            cd_1spike=cd_1spike_tab(idx,ny)+(cd_1spike_tab(idx+1,ny)
     &           -cd_1spike_tab(idx,ny))
     &           *(pdszpus-pdszpus_tab(idx))/(pdszpus_tab(idx+1)
     &           -pdszpus_tab(idx))
         else
            xi=(pdszpus-pdszpus_tab(idx))/(pdszpus_tab(idx+1)
     &           -pdszpus_tab(idx))
            et=(bzs-bzs_tab(idy))/(bzs_tab(idy+1)-bzs_tab(idy))
            z1=cd_1spike_tab(idx,idy)
            z2=cd_1spike_tab(idx+1,idy)
            z3=cd_1spike_tab(idx,idy+1)
            z4=cd_1spike_tab(idx+1,idy+1)
            cd_1spike=(1-xi)*(1-et)*z1+(1-xi)*et*z3
     &           +xi*(1-et)*z2+xi*et*z4        
         endif
      endif
!     
      return
      end
