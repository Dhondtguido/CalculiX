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
      subroutine printoutnodefluidfem(prlab,v,vold,physcon,ii,
     &  node,trab,inotr,ntrans,co,mi,vcon,nk,nknew)
!
!     stores results in the .dat file
!
      implicit none
!
      character*6 prlab(*)
!
      integer node,ii,j,inotr(2,*),ntrans,mi(*),nk,nodf,nknew(*)
!
      real*8 v(nk,0:mi(2)),trab(7,*),vcon(nk,0:mi(2)),
     &     co(3,*),a(3,3),physcon(*),vold(0:mi(2),*)
!
      nodf=nknew(node)
!
      if(prlab(ii)(1:4).eq.'VF  ') then
         if((ntrans.eq.0).or.(prlab(ii)(6:6).eq.'G')) then
            write(5,'(i10,1p,3(1x,e13.6))') node,
     &           (vold(j,nodf),j=1,3)
         elseif(inotr(1,nodf).eq.0) then
            write(5,'(i10,1p,3(1x,e13.6))') node,
     &           (vold(j,nodf),j=1,3)
         else
            call transformatrix(trab(1,inotr(1,nodf)),co(1,nodf),a)
            write(5,'(i10,1p,3(1x,e13.6))') node,
     &      vold(1,nodf)*a(1,1)+vold(2,nodf)*a(2,1)+vold(3,nodf)*a(3,1),
     &      vold(1,nodf)*a(1,2)+vold(2,nodf)*a(2,2)+vold(3,nodf)*a(3,2),
     &      vold(1,nodf)*a(1,3)+vold(2,nodf)*a(2,3)+vold(3,nodf)*a(3,3)
         endif
      elseif(prlab(ii)(1:4).eq.'PSF ') then
         write(5,'(i10,1x,1p,e13.6)') node,
     &           vold(4,nodf)
      elseif(prlab(ii)(1:4).eq.'TSF ') then
         write(5,'(i10,1x,1p,e13.6)') node,
     &           vold(0,nodf)
      elseif(prlab(ii)(1:4).eq.'PTF ') then
         write(5,'(i10,1x,1p,e13.6)') node,vold(4,nodf)*
     &        (1.d0+(v(nodf,0)-1.d0)/2*v(nodf,1)**2)**(v(nodf,0)/
     &        (v(nodf,0)-1.d0))
      elseif(prlab(ii)(1:4).eq.'TTF ') then
         write(5,'(i10,1x,1p,e13.6)') node,
     &     vold(0,nodf)*(1.d0+(v(nodf,0)-1.d0)/2*v(nodf,1)**2)
      elseif(prlab(ii)(1:4).eq.'CP  ') then
         write(5,'(i10,1x,1p,e13.6)') node,
     &            (vold(4,nodf)-physcon(6))*2.d0/
     &            (physcon(7)*physcon(5)**2)
      elseif(prlab(ii)(1:4).eq.'TURB') then
         write(5,'(i10,1x,1p,e13.6,1p,e13.6)') node,
     &            vold(5,nodf),vold(6,nodf)
      elseif(prlab(ii)(1:4).eq.'MACH') then
         write(5,'(i10,1x,1p,e13.6)') node,v(nodf,1)
      elseif(prlab(ii)(1:4).eq.'DEPF') then
        write(5,'(i10,1x,1p,e13.6)') node,vcon(nodf,4)
      endif
!
      flush(5)
!
      return
      end






