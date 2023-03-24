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
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!     this subroutine enables to calculate the dicharge coefficient 
!     for a labyrinth with more than one spike
!     as a function of the number of spikes(n), the pressure ratio (p2p1), 
!     the ratio between the gap and the breadth of the spike (s/b),
!     the number of reynolds (reynolds)
!
!     H.Zimmermann and K.h. Wolff
!     "Air system correlations part 1 Labyrinth seals"
!     asme 98-GT-206
!
!     author: Yannick Muller
!
      subroutine cd_lab_straight (n,p2p1,s,b,reynolds,cd_lab)
!
      implicit none
!
      integer i,j,n,idx,idy,nx,ny
!
      real*8 szb,p2p1,p1p2,s,b,reynolds,cd_lab,z1,z2,z3,z4,
     &     et,xi
!
      real*8 szb1(3)
      data szb1
     &     /0.230000d0,0.440000d0,0.830000d0/
!
      real*8 reynlds1(21)
      data reynlds1
     &     /100.0d0,200.0d0,300.d0,400.0d0,500.00d0,1000.0d0,
     &      2000.d0,3000.d0,5000.d0,7000.d0,9000.d0,11000.d0,13000.d0,
     &      15000.d0,18000.d0,21000.d0,25000.d0,30000.d0,35000.d0,
     &       40000.d0,50000.d0/
!
      real*8 tcd1(3,21)
      data ((tcd1(i,j),i=1,3),j=1,21)
     &     /0.470d0,0.330d0,0.230d0,
     &      0.500d0,0.365d0,0.274d0,
     &      0.517d0,0.385d0,0.300d0,
     &      0.520d0,0.400d0,0.320d0,
     &      0.530d0,0.415d0,0.333d0,
     &      0.550d0,0.449d0,0.376d0,
     &      0.575d0,0.483d0,0.420d0,
     &      0.590d0,0.500d0,0.450d0,
     &      0.607d0,0.530d0,0.480d0,
     &      0.620d0,0.550d0,0.500d0,
     &      0.625d0,0.565d0,0.515d0,
     &      0.630d0,0.570d0,0.527d0,
     &      0.630d0,0.580d0,0.540d0,
     &      0.630d0,0.585d0,0.555d0,
     &      0.630d0,0.589d0,0.565d0,
     &      0.630d0,0.589d0,0.576d0,
     &      0.630d0,0.590d0,0.580d0,
     &      0.630d0,0.590d0,0.588d0,
     &      0.630d0,0.590d0,0.590d0,
     &      0.630d0,0.590d0,0.590d0,
     &      0.630d0,0.590d0,0.590d0/
!
      real*8 szb2(3)
      data szb2 
     &     /0.230000d0,0.440000d0,0.830000d0/
!
      real*8 reynlds2(21)
      data reynlds2
     &     /100.0d0,200.0d0,300.d0,400.0d0,500.00d0,1000.0d0,
     &      2000.d0,3000.d0,5000.d0,7000.d0,9000.d0,11000.d0,13000.d0,
     &      15000.d0,18000.d0,21000.d0,25000.d0,30000.d0,35000.d0,
     &       40000.d0,50000.d0/
!
      real*8 tcd2(3,21)
      data ((tcd2(i,j),i=1,3),j=1,21)
     &     /0.400d0,0.335d0,0.250d0,
     &     0.445d0,0.390d0,0.308d0,
     &     0.470d0,0.420d0,0.340d0,
     &     0.490d0,0.440d0,0.360d0,
     &     0.505d0,0.455d0,0.380d0,
     &     0.550d0,0.500d0,0.442d0,
     &     0.600d0,0.555d0,0.500d0,
     &     0.625d0,0.580d0,0.525d0,
     &     0.650d0,0.615d0,0.570d0,
     &     0.660d0,0.640d0,0.600d0,
     &     0.660d0,0.650d0,0.617d0,
     &     0.660d0,0.655d0,0.635d0,
     &     0.660d0,0.657d0,0.645d0,
     &     0.660d0,0.660d0,0.650d0,
     &     0.660d0,0.660d0,0.655d0,
     &     0.660d0,0.660d0,0.660d0,
     &     0.660d0,0.660d0,0.660d0,
     &     0.660d0,0.660d0,0.660d0,
     &     0.660d0,0.660d0,0.660d0,
     &     0.660d0,0.660d0,0.660d0,
     &     0.660d0,0.660d0,0.660d0/
!
      p1p2=1/p2p1
      szb=s/b
!
!     which table is to be used?
!
      if(n.le.2) then
!     cd is interpolated in tcd1
!
       nx=3
       ny=22
!     interpolation in the 2d table.
!
       call ident(szb1,szb,nx,idx)
       call ident(reynlds1,reynolds,ny,idy)
!     
       if (idx.eq.0) then
          if(idy.eq.0) then
             cd_lab=tcd1(1,1)
          else
             if(idy.eq.ny) then
                cd_lab=tcd1(1,ny)
             else
                cd_lab=tcd1(1,idy)+(tcd1(1,idy+1)-tcd1(1,idy))
     &               *(reynolds-reynlds1(idy))
     &               /(reynlds1(idy+1)-reynlds1(idy))
             endif 
          endif
!     
       elseif(idx.ge.nx) then
          if(idy.le.0) then
             cd_lab=tcd1(nx,1)
          else
             if(idy.ge.ny) then
                cd_lab=tcd1(nx,ny)
             else
                cd_lab=tcd1(nx,idy)+(tcd1(nx,idy+1)-tcd1(nx,idy))
     &               *(reynolds-reynlds1(idy))
     &               /(reynlds1(idy+1)-reynlds1(idy))
             endif 
          endif
       else
          if(idy.le.0) then
             
             cd_lab=tcd1(idx,1)+(tcd1(idx+1,1)-tcd1(idx,1))
     &            *(szb-szb1(idx))/(szb1(idx+1)-szb1(idx))
          elseif(idy.ge.ny) then
             cd_lab=tcd1(idx,ny)+(tcd1(idx+1,ny)-tcd1(idx,ny))
     &            *(szb-szb1(idx))/(szb1(idx+1)-szb1(idx))
          else
             xi=(szb-szb1(idx))/(szb1(idx+1)-szb1(idx))
             et=(reynolds-reynlds1(idy))/
     &            (reynlds1(idy+1)-reynlds1(idy))
             z1=tcd1(idx,idy)
             z2=tcd1(idx+1,idy)
             z3=tcd1(idx,idy+1)
             z4=tcd1(idx+1,idy+1)
             cd_lab=(1-xi)*(1-et)*z1+(1-xi)*et*z3
     &            +xi*(1-et)*z2+xi*et*z4        
          endif
       endif         
!     
      else
!     cd is interpolated in tcd2
!     
         nx=3
         ny=21
!     interpolation in the 2d table.
!     
         call ident(szb2,szb,nx,idx)
         call ident(reynlds2,reynolds,ny,idy)
!     
         if (idx.eq.0) then
            if(idy.eq.0) then
               cd_lab=tcd2(1,1)
            else
               if(idy.eq.ny) then
                  cd_lab=tcd2(1,ny)
               else
                  cd_lab=tcd2(1,idy)+(tcd2(1,idy+1)-tcd2(1,idy))
     &                 *(reynolds-reynlds2(idy))
     &                 /(reynlds2(idy+1)-reynlds2(idy))
               endif 
            endif
!     
         elseif(idx.ge.nx) then
            if(idy.le.0) then
               cd_lab=tcd2(nx,1)
            else
               if(idy.ge.ny) then
                  cd_lab=tcd2(nx,ny)
               else
                  cd_lab=tcd2(nx,idy)+(tcd2(nx,idy+1)-tcd2(nx,idy))
     &                 *(reynolds-reynlds2(idy))
     &                 /(reynlds2(idy+1)-reynlds2(idy))
               endif 
            endif
         else
            if(idy.le.0) then
               
               cd_lab=tcd2(idx,1)+(tcd2(idx+1,1)-tcd2(idx,1))
     &              *(szb-szb2(idx))/(szb2(idx+1)-szb2(idx))
            elseif(idy.ge.ny) then
               cd_lab=tcd2(idx,ny)+(tcd2(idx+1,ny)-tcd2(idx,ny))
     &              *(szb-szb2(idx))/(szb2(idx+1)-szb2(idx))
            else
               xi=(szb-szb2(idx))/(szb2(idx+1)-szb2(idx))
               et=(reynolds-reynlds2(idy))/
     &              (reynlds2(idy+1)-reynlds2(idy))
               z1=tcd2(idx,idy)
               z2=tcd2(idx+1,idy)
               z3=tcd2(idx,idy+1)
               z4=tcd2(idx+1,idy+1)
               cd_lab=(1-xi)*(1-et)*z1+(1-xi)*et*z3
     &              +xi*(1-et)*z2+xi*et*z4        
            endif
         endif      
!     
      endif
!     
      return
      end
