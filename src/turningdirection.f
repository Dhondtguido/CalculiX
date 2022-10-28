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
      subroutine turningdirection(v,e1,e2,xn,mi,nk,turdir,lakon,ipkon,
     &     kon,ne,co)
!     
!     determines a local axis system based on the rotation axis
!     defined on a CENTRIF card 
!     
      implicit none
!     
      character*1 turdir
      character*8 lakon(*),lakonl
!     
      integer mi(*),i,j,k,nk,iflag,indexe,ipkon(*),kon(*),konl(20),
     &     loopa,ne,nope
!     
      real*8 v(0:mi(2),*),e1(3),e2(3),
     &     t(3,3),ur(3),ui(3),psi,psin,r,phi,phin,col(3),cog(3),
     &     coln(3),cogn(3),a,dist,xi,et,ze,xn(3),gs(8,4),xl(3,20),
     &     shp(4,20),xsj,ratio(20),co(3,*),umag,ureal,
     &     uimag,umagmax
!     
      turdir=' '
      iflag=1
!     
!     transformationsmatrix transforming the prime coordinates p'
!     into the global ones p: p=T.p'
!     
      do i=1,3
        t(i,1)=e1(i)
        t(i,2)=e2(i)
        t(i,3)=xn(i)
      enddo
!     
      umagmax=0.d0
      do i=1,ne
!     
!     only existing, volumetric elements (elements may have
!     been deactivated using *MODEL CHANGE; in that case lakon(i)
!     may still contain 'C3D'
!     
        if((ipkon(i).lt.0).or.(lakon(i)(1:3).ne.'C3D')) cycle
        indexe=ipkon(i)
        lakonl=lakon(i)
!     
!     number of nodes belonging to the element
!     
        if(lakonl(1:4).eq.'C3D8') then
          nope=8
        elseif(lakonl(4:5).eq.'20') then
          nope=20
        elseif(lakonl(4:5).eq.'10') then
          nope=10
        elseif(lakonl(4:4).eq.'4') then
          nope=4
        elseif(lakonl(4:5).eq.'15') then
          nope=15
        elseif(lakonl(4:4).eq.'6') then
          nope=6
        endif
!     
!     storing the global coordinates of the nodes
!     
        do j=1,nope
          konl(j)=kon(indexe+j)
          do k=1,3
            xl(k,j)=co(k,konl(j))
          enddo
        enddo
!     
!     location of the center of the element in local
!     coordinates
!     
        if((nope.eq.8).or.(nope.eq.20)) then
          xi=0.d0
          et=0.d0
          ze=0.d0
        elseif((nope.eq.4).or.(nope.eq.10)) then
          xi=0.25d0
          et=0.25d0
          ze=0.25d0
        else
          xi=0.33d0
          et=0.33d0
          ze=0.d0
        endif
!     
!     calculation of the shape functions and
!     in the center
!     
        if(lakonl(1:5).eq.'C3D8R') then
          call shape8hr(xl,xsj,shp,gs,a)
        elseif(lakonl(1:5).eq.'C3D8I') then
          call shape8hu(xi,et,ze,xl,xsj,shp,iflag)
        elseif(nope.eq.20) then
          if(lakonl(7:7).eq.'A') then
            call shape20h_ax(xi,et,ze,xl,xsj,shp,iflag)
          elseif((lakonl(7:7).eq.'E').or.(lakonl(7:7).eq.'S')) then
            call shape20h_pl(xi,et,ze,xl,xsj,shp,iflag)
          else
            call shape20h(xi,et,ze,xl,xsj,shp,iflag)
          endif
        elseif(nope.eq.8) then
          call shape8h(xi,et,ze,xl,xsj,shp,iflag)
        elseif(nope.eq.10) then
          call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
        elseif(nope.eq.4) then
          call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
        elseif(nope.eq.15) then
          call shape15w(xi,et,ze,xl,xsj,shp,iflag)
        else
          call shape6w(xi,et,ze,xl,xsj,shp,iflag)
        endif
!     
!     calculating the global coordinates of the center
!     
        do k=1,3
          cog(k)=0.d0
          do j=1,nope
            cog(k)=cog(k)+shp(4,j)*xl(k,j)
          enddo
        enddo
!     
!     transforming the global coordinates into prime coordinates
!     p'=T^T.p
!     
        do k=1,3
          col(k)=t(1,k)*cog(1)+t(2,k)*cog(2)+t(3,k)*cog(3)
        enddo
!     
!     calculating radius and angle
!     
        r=dsqrt(col(1)*col(1)+col(2)*col(2))
        phi=datan2(col(2),col(1))
!     
!     calculating the real and imaginary part of the displacements
!     at the center
!     
        do k=1,3
          ur(k)=0.d0
          ui(k)=0.d0
          do j=1,nope
            ur(k)=ur(k)+shp(4,j)*v(k,konl(j))
            ui(k)=ui(k)+shp(4,j)*v(k,nk+konl(j))
          enddo
        enddo
!     
!     no traveling wave if the imaginary part of the displacements
!     is zero
!     
        if(dabs(ui(1)*ui(1)+ui(2)*ui(2)+ui(3)*ui(3)).eq.0.d0) cycle
!     
!     looking for the turning direction of         
!     (ur(1)+i*ui(1))**2+(ur(2)+i*ui(2))**2+(ur(3)+i*ui(3))**2 =
!     ur(1)**2+ur(2)**2+ur(3)**2-ui(1)**2-ui(2)**2-ui(3)**2         
!     +2*i*(ur(1)*ui(1)+ur(2)*ui(2)+ur(3)*ui(3))
!     The turning direction is given by the sign of the phase, or,
!     equivalently, the sign of tan(phase)/2 
!     
        ureal=ur(1)*ur(1)+ur(2)*ur(2)+ur(3)*ur(3)-
     &       ui(1)*ui(1)-ui(2)*ui(2)-ui(3)*ui(3)
        uimag=ur(1)*ui(1)+ur(2)*ui(2)+ur(3)*ui(3)
        psi=uimag/ureal
        umag=ureal*ureal+uimag*uimag
        if(umag.gt.umagmax) then
          umagmax=umag
        else
          cycle
        endif
!     
!     perturbing phi slightly (this is an experimental value)
!     
        phin=phi+0.000001d0
!     
!     prime coordinates at the perturbed position      
!     
        coln(1)=r*dcos(phin)
        coln(2)=r*dsin(phin)
        coln(3)=col(3)
!     
!     calculating the global coordinates at the perturbed position
!     
        do k=1,3
          cogn(k)=t(k,1)*coln(1)+t(k,2)*coln(2)+t(k,3)*coln(3)
        enddo
!     
!     determining the local coordinates at the perturbed position        
!     
        loopa=8
        call attach_3d(xl,cogn,nope,ratio,dist,xi,et,ze,loopa)
!     
!     calculation of the shape functions
!     at the perturbed position
!     
        if(lakonl(1:5).eq.'C3D8R') then
          call shape8hr(xl,xsj,shp,gs,a)
        elseif(lakonl(1:5).eq.'C3D8I') then
          call shape8hu(xi,et,ze,xl,xsj,shp,iflag)
        elseif(nope.eq.20) then
          if(lakonl(7:7).eq.'A') then
            call shape20h_ax(xi,et,ze,xl,xsj,shp,iflag)
          elseif((lakonl(7:7).eq.'E').or.(lakonl(7:7).eq.'S')) then
            call shape20h_pl(xi,et,ze,xl,xsj,shp,iflag)
          else
            call shape20h(xi,et,ze,xl,xsj,shp,iflag)
          endif
        elseif(nope.eq.8) then
          call shape8h(xi,et,ze,xl,xsj,shp,iflag)
        elseif(nope.eq.10) then
          call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
        elseif(nope.eq.4) then
          call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
        elseif(nope.eq.15) then
          call shape15w(xi,et,ze,xl,xsj,shp,iflag)
        else
          call shape6w(xi,et,ze,xl,xsj,shp,iflag)
        endif
!     
!     calculating the real and imaginary part of the displacements
!     at the perturbed position
!     
        do k=1,3
          ur(k)=0.d0
          ui(k)=0.d0
          do j=1,nope
            ur(k)=ur(k)+shp(4,j)*v(k,konl(j))
            ui(k)=ui(k)+shp(4,j)*v(k,nk+konl(j))
          enddo
        enddo
!     
!     calculating psi at the perturbed position
!     
        psin=(ur(1)*ui(1)+ur(2)*ui(2)+ur(3)*ui(3))/
     &       (ur(1)*ur(1)+ur(2)*ur(2)+ur(3)*ur(3)-
     &       ui(1)*ui(1)-ui(2)*ui(2)-ui(3)*ui(3))
!     
!     determining the turning direction
!     
        if(psin.gt.psi) then
          turdir='B'
        else
          turdir='F'
        endif
      enddo
!     
      return
      end

