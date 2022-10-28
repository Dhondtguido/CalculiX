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
      subroutine cd_chamfer(l,d,p_up,p_down,angle,cd)
!
!     calculates the discharge coefficient of holes with chamfered inlets 
!     using N. Hay and A.Spencer 
!     "Disharge coefficient of Cooling holes with radiused and chamfered 
!     inlets" ASME 91-GT-269
!     
!     Nota:the radius correction is not used here due to the unreliability 
!     of the results proposed check first line of table 1
!
!     author: Yannick Muller
!     
      implicit none
!
      integer i,j,idx,idy,nx,ny
!
      real*8 l,d,p_up,p_down,angle,puzpd,lzd,xi,et,z1,z2,z3,z4,
     &     cd, tab_cd(3,4), tab30(3,4),tab45(3,4)
!
      real*8 xpuzpd(3)
      data xpuzpd /1.2d0,1.6d0,2.2d0/
!
      real*8 ylzd (4)
      data ylzd /0.25d0,0.50d0,1.00d0,2.00d0/
!
      data ((tab30(i,j),i=1,3),j=1,4)
     &           /1.45d0,1.31d0,1.24d0,
     &            1.35d0,1.28d0,1.21d0,
     &            1.23d0,1.19d0,1.13d0,
     &            1.20d0,1.18d0,1.10d0/
!
      data ((tab45(i,j),i=1,3),j=1,4) 
     &           /1.19d0,1.19d0,1.16d0,
     &            1.23d0,1.19d0,1.13d0,
     &            1.14d0,1.11d0,1.07d0,
     &            1.11d0,1.09d0,1.03d0/
!     
      nx=3
      ny=4
!     
      lzd=l/d
      puzpd=p_up/p_down
!
      call ident(xpuzpd,puzpd,nx,idx)
      call ident(ylzd,lzd,ny,idy)
!     
      if (abs(angle-30.d0).le.0.1d0) then
         do i=1,3
            do j=1,4
            tab_cd(i,j)=tab30(i,j)
         enddo
      enddo
!
      elseif(abs(angle-45.d0).le.0.1d0) then
         do i=1,3
            do j=1,4
               tab_cd(i,j)=tab45(i,j)
            enddo
         enddo
      else 
         write(*,*) '*WARNING in cd_chamfer.f :unacceptable angle'
     &,angle,'grad'
        write(*,*) 'Chamfer correction is assumed Cd_chamfer=1'
      endif
!
      if (idx.eq.0) then
         if(idy.eq.0) then
            cd=tab_cd(1,1)
         else
            if(idy.eq.ny) then
               cd=tab_cd(1,ny)
            else
               cd=tab_cd(1,idy)+(tab_cd(1,idy+1)-tab_cd(1,idy))
     &              *(lzd-ylzd(idy))/(ylzd(idy+1)-ylzd(idy))
            endif 
         endif
!     
      elseif(idx.ge.nx) then
         if(idy.le.0) then
            cd=tab_cd(nx,1)
         else
            if(idy.ge.ny) then
               cd=tab_cd(nx,ny)
            else
               cd=tab_cd(nx,idy)+(tab_cd(nx,idy+1)-tab_cd(nx,idy))
     &              *(lzd-ylzd(idy))/(ylzd(idy+1)-ylzd(idy))
            endif 
         endif
      else
         if(idy.le.0) then
            cd=tab_cd(idx,1)+(tab_cd(idx+1,1)-tab_cd(idx,1))
     &           *(puzpd-xpuzpd(idx))/(xpuzpd(idx+1)-xpuzpd(idx))
         elseif(idy.ge.ny) then
            cd=tab_cd(idx,ny)+(tab_cd(idx+1,ny)-tab_cd(idx,ny))
     &           *(puzpd-xpuzpd(idx))/(xpuzpd(idx+1)-xpuzpd(idx))
         else
            xi=(puzpd-xpuzpd(idx))/(xpuzpd(idx+1)-xpuzpd(idx))
            et=(lzd-ylzd(idy))/(ylzd(idy+1)-ylzd(idy))
            z1=tab_cd(idx,idy)
            z2=tab_cd(idx+1,idy)
            z3=tab_cd(idx,idy+1)
            z4=tab_cd(idx+1,idy+1)
            cd=(1-xi)*(1-et)*z1+(1-xi)*et*z3
     &           +xi*(1-et)*z2+xi*et*z4 
         endif
      endif
!     
!      write(*,*)'chamfer correction equals to',cd
!     
      return
      end
      
