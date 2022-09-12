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
!     this function enable to determine the discharge coefficient of bleed
!     tappings
!
      subroutine cd_bleedtapping(ps2,ps1,ps1pt1,nummer,curve,x_tab,y_tab
     &     ,cd)
!
!
!     in : SImultation of the secondary air system of aero engines
!          K.J.KUTZ T.M. SPEER
!          Transactions of the ASME vol.116 April 1994
!
!     author: Yannick Muller
!
      implicit none
!
      integer nummer,id,i,number,curve,index
      real*8 x_tab(15),y_tab(15)
!
!     Fig.7 tapping with lip
!
      real*8 cdx1(9)
      data  cdx1
     &     /0.24d0,0.52d0,0.8d0,1.14d0,1.42d0,1.9d0,2.5d0,3d0,3.4d0/
!

      real*8 cdy1(9)
      data cdy1
     &     /0.167d0,0.310d0,0.467d0,0.611d0,0.711d0,0.789d0,0.833d0,
     &     0.866d0,0.888d0/ 
!
!    Fig.7 tapping without lip
!
      real*8 cdx2(7)
      data  cdx2
     &     /1.0d0,1.14d0,1.42d0,1.9d0,2.5d0,3.0d0,3.4d0/

      real*8 cdy2(7)
      data cdy2
     &     /0.d0,0.122d0,0.377d0,0.7d0,0.766d0,0.769d0,0.772d0/

      real*8 ps2,ps1,dab,ps2pt1,ps1pt1,cdy(15),cd,cdx(20),
     &     dabmax
!
      ps2pt1=ps2/ps1
      dabmax=100.d0
! 
      if(nummer.eq.0) then
         if (curve.eq.1) then
            index=9
            write(*,*)
            write(*,*) 'Cd calculations will be performed using'
            write(*,*) 'Cd-Kurven HP3 Schlitz;Kurve Nr. 1'
            do i=1,index
               cdx(i)=cdx1(i)
               cdy(i)=cdy1(i)
            enddo
!     
         elseif(curve.eq.2) then
            index=7
            write(*,*)
            write(*,*) 'Cd calculations will be performed using'
            write(*,*) 'Cd-Kurven HP3 Schlitz;Kurve Nr. 2'
            do i=1,index
               cdx(i)=cdx2(i)
               cdy(i)=cdy2(i)
            enddo
!     
         elseif(curve.gt.2) then
            write(*,*)
            write(*,*) 'no characteristic available under this index'
            write(*,*) 'cd is implicitely assumed equal to 1'
            cd=1.d0
            return
         endif
!     
!     psvptv  ratio between the static pressure in the main canal 
!     and the total pressure in the main canal
!     
!     check whether ps1/pt1 less than 1 , if not then a warning is sent and 
!     the calculation will peroceed with an  "oversized" dab
!     
         if(abs(1.d0-ps2pt1).le.dabmax*(1.d0-ps1pt1)) then
            dab=(1.d0-ps2pt1)/(1.d0-ps1pt1)
         else 
            dab=dabmax
            write(*,*) 'in cd_bleedtapping.f: ps1/pt1=',ps1pt1
            write(*,*) 'the calculation will proceed with DAB=100.'
         endif
!     
!     determination of cd with the caracteristics
!     
         call ident(cdx,dab,index,id)
         if(id.eq.0) then
            cd=cdy(1)
         elseif(id.ge.index) then
            cd=cdy(index)
         else
            cd=cdy(id)+(cdy(id+1)-cdy(id))
     &           *(dab-cdx(id))/(cdx(id+1)-cdx(id))
         endif
!     
      else
         if(abs(1.d0-ps2pt1).le.dabmax*(1.d0-ps1pt1)) then
            dab=(1.d0-ps2pt1)/(1.d0-ps1pt1)
         else 
            dab=dabmax
            write(*,*) 'in cd_bleedtapping.f: ps1/pt1=',ps1pt1
            write(*,*) 'the calculation will proceed with DAB=100.'
         endif

         call ident(x_tab,dab,nummer,id)
         if(id.eq.0) then
            cd=y_tab(1)
         elseif(id.ge.nummer) then
            cd=y_tab(nummer)
         else
            cd=y_tab(id)+(y_tab(id+1)-y_tab(id))
     &           *(dab-x_tab(id))/(x_tab(id+1)-x_tab(id))
         endif   
      endif
      return
      end
