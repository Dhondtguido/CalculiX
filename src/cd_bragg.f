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
!     this subroutine enables to calculate a compressibility correction factor 
!     following the results that can be found in:
!     
!     S.L.Bragg
!     "Effect of conpressibility on the discharge coefficient of orifices 
!     and convergent nozzles"
!     Journal of Mechanical engineering vol 2 No 1 1960
!
!     author: Yannick Muller
!     
      subroutine cd_bragg(cd,p2p1,cdbragg,itype)
!     
      implicit none
!     
!     itype is used in the company version of cd_bragg
!
      integer nx,ny,idx,idy,i,j,itype
!     
      real*8 cd,p2p1,cdbragg,z1,z2,z3,z4,et,xi
!     
      real*8 cd_tab (12)
      data cd_tab
     &     /0.457d0,0.500d0,0.550d0,0.600d0,0.650d0,0.700d0,
     &     0.750d0,0.800d0,0.850d0,0.900d0,0.950d0,1.000d0/
!     
      real*8 p2p1_tab (19)
      data p2p1_tab
     &     /0.00d0,0.10d0,0.15d0,0.20d0,0.25d0,0.30d0,0.35d0,0.40d0,
     &     0.45d0,0.50d0,0.55d0,0.60d0,0.65d0,0.70d0,0.75d0,0.80d0,
     &     0.85d0,0.90d0,1.00d0/
!     
      real*8 cd_bragg_tab(19,12)
      data ((cd_bragg_tab(i,j),i=1,19),j=1,12)
     &     /0.754d0,0.735d0,0.724d0,0.712d0,0.701d0,0.688d0,0.672d0,
     &     0.655d0,0.633d0,0.610d0,0.590d0,0.570d0,0.549d0,0.530d0,
     &     0.514d0,0.500d0,0.488d0,0.477d0,0.454d0,
!     
     &     0.789d0,0.770d0,0.760d0,0.749d0,0.747d0,0.733d0,0.709d0,
     &     0.691d0,0.672d0,0.650d0,0.628d0,0.606d0,0.588d0,0.572d0,
     &     0.558d0,0.543d0,0.531d0,0.520d0,0.500d0,
!     
     &     0.833d0,0.815d0,0.805d0,0.796d0,0.783d0,0.771d0,0.758d0,
     &     0.740d0,0.720d0,0.700d0,0.675d0,0.655d0,0.638d0,0.621d0,
     &     0.607d0,0.592d0,0.580d0,0.569d0,0.550d0,
!     
     &     0.870d0,0.855d0,0.846d0,0.828d0,0.827d0,0.815d0,0.801d0,
     &     0.786d0,0.769d0,0.749d0,0.725d0,0.704d0,0.685d0,0.670d0,
     &     0.654d0,0.641d0,0.630d0,0.619d0,0.600d0,
!     
     &     0.902d0,0.890d0,0.882d0,0.875d0,0.867d0,0.855d0,0.842d0,
     &     0.830d0,0.811d0,0.792d0,0.773d0,0.751d0,0.732d0,0.718d0,
     &     0.700d0,0.689d0,0.678d0,0.668d0,0.650d0,
!     
     &     0.929d0,0.920d0,0.912d0,0.908d0,0.900d0,0.890d0,0.880d0,
     &     0.869d0,0.852d0,0.835d0,0.815d0,0.794d0,0.778d0,0.761d0,
     &     0.749d0,0.736d0,0.725d0,0.716d0,0.700d0,
!     
     &     0.952d0,0.946d0,0.940d0,0.936d0,0.930d0,0.921d0,0.913d0,
     &     0.903d0,0.889d0,0.873d0,0.854d0,0.836d0,0.820d0,0.808d0,
     &     0.796d0,0.785d0,0.775d0,0.766d0,0.750d0,
!     
     &     0.970d0,0.966d0,0.962d0,0.958d0,0.953d0,0.948d0,0.941d0,
     &     0.935d0,0.923d0,0.909d0,0.890d0,0.874d0,0.860d0,0.849d0,
     &     0.838d0,0.829d0,0.820d0,0.812d0,0.800d0,
!     
     &     0.983d0,0.9805d0,0.98d0,0.978d0,0.975d0,0.970d0,0.965d0,
     &     0.958d0,0.950d0,0.949d0,0.926d0,0.911d0,0.900d0,0.890d0,
     &     0.881d0,0.874d0,0.867d0,0.860d0,0.850d0,
!
     &     0.992d0,0.991d0,0.990d0,0.989d0,0.988d0,0.985d0,0.981d0,
     &     0.980d0,0.973d0,0.967d0,0.956d0,0.943d0,0.935d0,0.928d0,
     &     0.920d0,0.915d0,0.910d0,0.907d0,0.900d0,
!     
     &     0.999d0,0.999d0,0.998d0,0.998d0,0.998d0,0.997d0,0.995d0,
     &     0.992d0,0.990d0,0.988d0,0.981d0,0.975d0,0.970d0,0.964d0,
     &     0.960d0,0.958d0,0.954d0,0.952d0,0.950d0,
!     
     &     1.000d0,1.000d0,1.000d0,1.000d0,1.000d0,1.000d0,1.000d0,
     &     1.000d0,1.000d0,1.000d0,1.000d0,1.000d0,1.000d0,1.000d0,
     &     1.000d0,1.000d0,1.000d0,1.000d0,1.000d0/
!     
      nx=19
      ny=12
!     
      call ident(p2p1_tab,p2p1,nx,idx)
      call ident(cd_tab,cd,ny,idy)
!     
      if (idx.eq.0) then
         if(idy.eq.0) then
            cdbragg=cd_bragg_tab(1,1)
         else
            if(idy.eq.ny) then
               cdbragg=cd_bragg_tab(1,ny)
            else
               cdbragg=cd_bragg_tab(1,idy)+(cd_bragg_tab(1,idy+1)
     &              -cd_bragg_tab(1,idy))
     &              *(cd-cd_tab(idy))/(cd_tab(idy+1)-cd_tab(idy))
            endif 
         endif
!     
      elseif(idx.ge.nx) then
         if(idy.le.0) then
            cdbragg=cd_bragg_tab(nx,1)
         else
            if(idy.ge.ny) then
               cdbragg=cd_bragg_tab(nx,ny)
            else
               cdbragg=cd_bragg_tab(nx,idy)+
     &              (cd_bragg_tab(nx,idy+1)-cd_bragg_tab(nx,idy))
     &              *(cd-cd_tab(idy))/(cd_tab(idy+1)-cd_tab(idy))
            endif 
         endif
      else
         if(idy.le.0) then
!     
            cdbragg=cd_bragg_tab(idx,1)+(cd_bragg_tab(idx+1,1)
     &           -cd_bragg_tab(idx,1))
     &           *(p2p1-p2p1_tab(idx))/(p2p1_tab(idx+1)
     &           -p2p1_tab(idx))
         elseif(idy.ge.ny) then
            cdbragg=cd_bragg_tab(idx,ny)+(cd_bragg_tab(idx+1,ny)
     &           -cd_bragg_tab(idx,ny))
     &           *(p2p1-p2p1_tab(idx))/(p2p1_tab(idx+1)
     &           -p2p1_tab(idx))
         else
            xi=(p2p1-p2p1_tab(idx))/(p2p1_tab(idx+1)
     &           -p2p1_tab(idx))
            et=(cd-cd_tab(idy))/(cd_tab(idy+1)-cd_tab(idy))
            z1=cd_bragg_tab(idx,idy)
            z2=cd_bragg_tab(idx+1,idy)
            z3=cd_bragg_tab(idx,idy+1)
            z4=cd_bragg_tab(idx+1,idy+1)
            cdbragg=(1-xi)*(1-et)*z1+(1-xi)*et*z3
     &           +xi*(1-et)*z2+xi*et*z4        
         endif
      endif
!     
      return
      end
