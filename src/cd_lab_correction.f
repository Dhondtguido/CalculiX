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
!     this subroutine enables to calculate thecorrection factor of the discharge  
!     coefficient of a labyrinth with one fin as a function of the ratio b/s and the 
!     pressure ratio Pdownstream/Pupstream
!     the results are interpolated 
!     
!     the relevant data can be found in:
!     "Air System Correlations Part 1: Labyrinth Seals"
!     H.Zimmermann and K.H. Wollf
!     ASME 98-GT-206
!     fig 12 p 7
!
!     author: Yannick Muller
!
      subroutine cd_lab_correction(p1p2,s,b,cd_correction)
!
      implicit none
! 
      integer nx,ny,idx,idy
!     
      real*8 s,b,cd_correction,z1,z2,z3,z4,xi,et,szb,p1p2
!
      real*8 puszpds_tab(7)
      data puszpds_tab
     &    /1.d0,1.2d0,1.4d0,1.6d0,1.8d0,2.d0,2.5d0/
!
      real*8 szb_tab(9)
      data szb_tab
     &    /0.25d0,0.5d0,1.d0,1.5d0,2d0,2.5d0,3d0,3.5d0,4d0/
!     
      real*8 cd_correction_tab(9,7)
      data cd_correction_tab
     &  /1.05d0,1.07d0,1.03d0,0.98d0,0.95d0,0.94d0,0.95d0,0.95d0,0.95d0,
     &   1.15d0,1.07d0,1.02d0,0.95d0,0.92d0,0.91d0,0.91d0,0.92d0,0.92d0,
     &   1.15d0,1.05d0,0.98d0,0.91d0,0.88d0,0.86d0,0.86d0,0.87d0,0.87d0,
     &   1.15d0,1.04d0,0.95d0,0.87d0,0.85d0,0.84d0,0.83d0,0.83d0,0.83d0,
     &   1.15d0,1.03d0,0.91d0,0.85d0,0.81d0,0.80d0,0.80d0,0.80d0,0.80d0, 
     &   1.15d0,1.01d0,0.90d0,0.82d0,0.79d0,0.79d0,0.77d0,0.77d0,0.77d0,
     &   1.10d0,1.00d0,0.88d0,0.79d0,0.75d0,0.74d0,0.73d0,0.72d0,0.70d0/
!
      szb=s/b
!      
      nx=9
      ny=7
!
!      p1p2=1/p2p1
!      if ((p1p2.ge.2.5d0).or.(szb.ge.4d0))then
!         write(*,*) '*WARNING in cd_lab_correction'
!         write(*,*) 'p1p2>2.5 or szb>4'
!         write(*,*) 'check input file'
!         write(*,*) 'calculation will proceed using cd_lab_correction=1'
!         cd_correction=1.d0
!        return
!      endif
!
      call ident(puszpds_tab,p1p2,ny,idy)
      call ident(szb_tab,szb,nx,idx)
!     
      if (idx.eq.0) then
         if(idy.eq.0) then
            cd_correction=cd_correction_tab(1,1)
         else
            if(idy.eq.ny) then
               cd_correction=cd_correction_tab(1,ny)
            else
               cd_correction=cd_correction_tab(1,idy)
     &           +(cd_correction_tab(1,idy+1)-cd_correction_tab(1,idy))
     &              *(szb-szb_tab(idx))/(szb_tab(idx+1)-szb_tab(idx))
            endif 
         endif
!     
      elseif(idx.ge.nx) then
         if(idy.le.0) then
            cd_correction=cd_correction_tab(nx,1)
         else
            if(idy.ge.ny) then
               cd_correction=cd_correction_tab(nx,ny)
            else
            cd_correction=cd_correction_tab(nx,idy)
     &     +(cd_correction_tab(nx,idy+1)-cd_correction_tab(nx,idy))
     &     *(szb-szb_tab(idx))/(szb_tab(idx+1)-szb_tab(idx))
            endif 
         endif
      else
         if(idy.le.0) then
!     
            cd_correction=cd_correction_tab(idx,1)
     &          +(cd_correction_tab(idx+1,1)-cd_correction_tab(idx,1))
     &           *(p1p2-puszpds_tab(idy))/(puszpds_tab(idy+1)
     &           -puszpds_tab(idy))
         elseif(idy.ge.ny) then
            cd_correction=cd_correction_tab(idx,ny)
     &         +(cd_correction_tab(idx+1,ny)-cd_correction_tab(idx,ny))
     &           *(p1p2-puszpds_tab(idy))/(puszpds_tab(idy+1)
     &           -puszpds_tab(idy))
         else
            et=(p1p2-puszpds_tab(idy))/(puszpds_tab(idy+1)
     &           -puszpds_tab(idy))
            xi=(szb-szb_tab(idx))/(szb_tab(idx+1)-szb_tab(idx))
            z1=cd_correction_tab(idx,idy)
            z2=cd_correction_tab(idx+1,idy)
            z3=cd_correction_tab(idx,idy+1)
            z4=cd_correction_tab(idx+1,idy+1)
            cd_correction=(1-xi)*(1-et)*z1+(1-xi)*et*z3
     &           +xi*(1-et)*z2+xi*et*z4        
         endif
      endif
!
!      if (cd_correction.ge.1.d0)then
!         cd_correction=1.d0
!      endif
!     
      return
      end
