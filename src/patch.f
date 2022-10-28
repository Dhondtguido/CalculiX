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
      subroutine patch(iterms,node,sti,scpav,mi,kon,ipkon,
     & ipoints,members,linpatch,co,lakon,iavflag)
!
!     computes the smoothed nodal stresses for an element patch
!
!     author: Sascha Merz
!
      implicit none
!
      integer i,j,nope,mint3d,indexe,ielem,kon(*),ipkon(*),
     & iterms,k,iflag,nrhs,info,node,ipnt,
     & mi(*),irow,ipoints,members(*),linpatch,ielidx,iavflag
!
      real*8 xl(3,20),co(3,*),shp(4,20),
     & pgauss(3),pol(20),pntdist,w,scpav(6,*),xsj,
     & sti(6,mi(1),*),xi,et,ze,rv1(ipoints),pp(ipoints,iterms),
     & pdat(ipoints,6),pfit(6),pwrk(iterms),pre(ipoints,iterms),
     & z(ipoints,ipoints)
!
      real*8 tmpstr(6),gauss3d5e(3,4)
!
      character*8 lakon(*)
!      
      logical matu,matv
!
      include 'gauss.f'
!
!     initialize
!
!     iavflag: if 1, then build average of integration point values
!
      if(iavflag.eq.0) then
         do i=1,6
            pfit(i)=0.d0
         enddo
         iflag=1
!
!        irow: row number of the rectangular matrix of the overdetermined
!        system of equations
!
         irow=0
!
!        loop over patch elements
!        linpatch: number of elements in patch
!     
         do ielidx=1,linpatch
!
            ielem=members(ielidx)
            if(lakon(ielem)(1:5).eq.'C3D20') then
!
!              nope: nodes per element
!
               nope=20
               if(lakon(ielem)(6:6).eq.'R') then
                  mint3d=8
               else
                  mint3d=27
               endif
            elseif(lakon(ielem)(1:4).eq.'C3D8') then
               nope=8
               if(lakon(ielem)(5:5).eq.'R') then
                  mint3d=1
               else
                  mint3d=8
               endif
            elseif(lakon(ielem)(1:5).eq.'C3D10') then
               nope=10
               mint3d=4
            endif
!
            indexe=ipkon(ielem)
!     
!           coordinates of the nodes belonging to the element
!               
            do j=1,nope
               do k=1,3
                  xl(k,j)=co(k,kon(indexe+j))
               enddo
            enddo
!     
!           loop over the integration points in one element
!
            do ipnt=1,mint3d
               irow=irow+1
               if(nope.eq.10) then
                  xi=gauss3d5(1,ipnt)
                  et=gauss3d5(2,ipnt)
                  ze=gauss3d5(3,ipnt)
                  call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.20) then
                  if(mint3d.eq.8) then
                     xi=gauss3d2(1,ipnt)
                     et=gauss3d2(2,ipnt)
                     ze=gauss3d2(3,ipnt)
                  elseif(mint3d.eq.27) then
                     xi=gauss3d3(1,ipnt)
                     et=gauss3d3(2,ipnt)
                     ze=gauss3d3(3,ipnt)
                  endif
                  call shape20h(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.8) then
                  if(mint3d.eq.1) then
                     xi=gauss3d1(1,ipnt)
                     et=gauss3d1(2,ipnt)
                     ze=gauss3d1(3,ipnt)
                  elseif(mint3d.eq.8) then
                     xi=gauss3d2(1,ipnt)
                     et=gauss3d2(2,ipnt)
                     ze=gauss3d2(3,ipnt)
                  endif
                  call shape8h(xi,et,ze,xl,xsj,shp,iflag)
               endif
!
!              in pgauss the global coordinates of the integration
!              point are saved.
!
               do j=1,3
                  pgauss(j)=0.d0
                  do k=1,nope
                     pgauss(j)=pgauss(j)+shp(4,k)*xl(j,k)
                  enddo
!
!                 the origin of the coordinate system is moved to the
!                 evaluated node for higher numerical stability
!
                  pgauss(j)=pgauss(j)-co(j,node)
               enddo
!
!              evaluate patch polynomial for the integration point
!  
               pol(1)=1.d0
               pol(2)=pgauss(1)
               pol(3)=pgauss(2)
               pol(4)=pgauss(3)
               pol(5)=pgauss(1)*pgauss(2)
               pol(6)=pgauss(1)*pgauss(3)
               pol(7)=pgauss(2)*pgauss(3)
               pol(8)=pgauss(1)*pgauss(1)
               pol(9)=pgauss(2)*pgauss(2)
               pol(10)=pgauss(3)*pgauss(3)
               pol(11)=pgauss(1)*pgauss(2)*pgauss(3)
               if(iterms.gt.11) then
                  pol(12)=pgauss(1)*pgauss(1)*pgauss(2)
                  pol(13)=pgauss(1)*pgauss(2)*pgauss(2)
                  pol(14)=pgauss(1)*pgauss(1)*pgauss(3)
                  pol(15)=pgauss(1)*pgauss(3)*pgauss(3)
                  pol(16)=pgauss(2)*pgauss(2)*pgauss(3)
                  pol(17)=pgauss(2)*pgauss(3)*pgauss(3)
                  if(iterms.gt.17) then
                     pol(18)=pgauss(1)*pgauss(1)*pgauss(1)
                     pol(19)=pgauss(2)*pgauss(2)*pgauss(2)
                     pol(20)=pgauss(3)*pgauss(3)*pgauss(3)
                  endif
               endif
!     
!              weighting for integration point
!
               pntdist=dsqrt(pgauss(1)*pgauss(1)
     &              +pgauss(2)*pgauss(2)
     &              +pgauss(3)*pgauss(3))
               w=pntdist**(-1.5d0)
!     
               do j=1,6
                  pdat(irow,j)=sti(j,ipnt,ielem)*w
               enddo
!     
               do j=1,iterms
                  pp(irow,j)=pol(j)*w
               enddo
            enddo
         enddo
!
!        using singular value decomposition for the least squares fit
!
         matu=.false.
         matv=.true.
         nrhs=6
!
         call hybsvd(ipoints,ipoints,ipoints,ipoints,ipoints,
     &        ipoints,iterms,pp,pwrk,matu,pp,matv,
     &        pre,z,pdat,nrhs,info,rv1)
         if(info.ne.0) then
            write(*,*) '*ERROR in patch: Bad conditioned matrix,',
     &           ' using average of sampling point values.'
            iavflag=1
         endif
      endif
!
!     matrix multiplication. only the first value of the
!     solution vector is needed. the singular values are manipulated
!     to increase the numerical stability
!     
      if(iavflag.eq.0) then
         do j=1,iterms
            if(pwrk(j).lt.1.d-22) then
               pwrk(j)=0.d0
            else 
               pwrk(j)=1.d0/pwrk(j)
            endif
         enddo
         do j=1,iterms
            pre(1,j)=pre(1,j)*pwrk(j)
         enddo
         do j=1,nrhs
            do k=1,iterms
               pfit(j)=pfit(j)+pre(1,k)*pdat(k,j)
            enddo
         enddo
!
!        pfit is an array containing the coefficients for the polynom
!        for the six stress components
!
!        solution in the node
!
         do j=1,6
            scpav(j,node)=scpav(j,node)+pfit(j)
         enddo
      endif
!
!     if there are not enough elements to fit a polynomial,
!     build average value of the sampling point values
!
      if(iavflag.eq.1) then
         do j=1,6
            tmpstr(j)=0.d0
         enddo
         irow=0
         do ielidx=1,linpatch
            ielem=members(ielidx)
            if(lakon(ielem)(1:5).eq.'C3D20') then
               if(lakon(ielem)(6:6).eq.'R') then
                  mint3d=8
               else
                  mint3d=27
               endif
            elseif(lakon(ielem)(1:5).eq.'C3D10') then
               mint3d=4
            elseif(lakon(ielem)(1:4).eq.'C3D8') then
               if(lakon(ielem)(5:5).eq.'R') then
                  mint3d=1
               else
                  mint3d=8
               endif
            endif
            do ipnt=1,mint3d
               irow=irow+1
               do j=1,6
                  tmpstr(j)=tmpstr(j)+sti(j,ipnt,ielem)
               enddo 
            enddo
         enddo
         do j=1,6
            scpav(j,node)=scpav(j,node)+tmpstr(j)/irow
         enddo
      endif
!
      return
      end
