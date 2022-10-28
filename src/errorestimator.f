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
      subroutine errorestimator(yi,yn,ipkon,kon,lakon,nk,
     &  ne,mi,ielmat,nterms,inum,co,vold,cflag,ielprop,prop)
!
!     the error in the node is calculated based on the maximum difference
!     between the max principal stress (mechanical calculations)
!     or heat flux size (thermal calculations) at the integration
!     points in the elements belonging to the node
!
!     for mechanical applications this error is converted into a
!     true stress error with heuristic equations
!
      implicit none
!
      logical force
!
      character*1 cflag
      character*8 lakon(*),lakonl
!
      integer ipkon(*),kon(*),mi(*),ne,indexe,null,nonei20(3,12),
     &  nonei10(3,6),nk,i,j,k,node,nonei15(3,9),nopev,nterms,
     &  mint3d,ielmat(mi(3),*),inum(*),ielprop(*)
!
      real*8 yi(nterms,mi(1),*),yn(nterms,*),size,wpsmin,wpsmax,
     &  absdiff,reldiff,sizemax,al(3),sizemin,c(3,3),prop(*),
     &  wpsmin1,wpsmax1,wpsmin3,wpsmax3,co(3,*),vold(0:mi(2),*)
!
      data nonei10 /5,1,2,6,2,3,7,3,1,8,1,4,9,2,4,10,3,4/
!
      data nonei15 /7,1,2,8,2,3,9,3,1,10,4,5,11,5,6,12,6,4,
     &  13,1,4,14,2,5,15,3,6/
!
      data nonei20 /9,1,2,10,2,3,11,3,4,12,4,1,
     &  13,5,6,14,6,7,15,7,8,16,8,5,
     &  17,1,5,18,2,6,19,3,7,20,4,8/
!
      null=0
!
      do i=1,nk
         do j=1,nterms
            yn(j,i)=0.d0
         enddo
      enddo
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         indexe=ipkon(i)
         lakonl=lakon(i)
!
         if(lakonl(7:8).eq.'LC') cycle
!
         if(lakonl(1:1).eq.'F') then
            cycle
         elseif(lakonl(4:4).eq.'2') then
            nopev=8
         elseif(lakonl(4:4).eq.'8') then
            nopev=8
         elseif(lakonl(4:5).eq.'10') then
            nopev=4
         elseif(lakonl(4:4).eq.'4') then
            nopev=4
         elseif(lakonl(4:5).eq.'15') then
            nopev=6
         elseif(lakonl(4:4).eq.'6') then
            nopev=6
         else
            cycle
         endif
!
         if(lakonl(4:5).eq.'8R') then
            mint3d=1
         elseif((lakonl(4:4).eq.'8').or.
     &           (lakonl(4:6).eq.'20R')) then
            mint3d=8
         elseif(lakonl(4:4).eq.'2') then
            mint3d=27
         elseif(lakonl(4:5).eq.'10') then
            mint3d=4
         elseif(lakonl(4:4).eq.'4') then
            mint3d=1
         elseif(lakonl(4:5).eq.'15') then
            mint3d=9
         elseif(lakonl(4:5).eq.'6') then
            mint3d=2
         elseif(lakonl(1:2).eq.'ES') then
            cycle
         endif
!
!        calculating the maximal differences of first principal
!        stress across the integration points or of the
!        temperatures across the nodes
!
         absdiff=0.d0
         reldiff=0.d0
!
         if(nterms.eq.6) then
!
!           mechanical calculation: max principal stress
!
            wpsmin1=1.d30
            wpsmax1=-1.e30
!
            wpsmin3=1.d30
            wpsmax3=-1.e30
!
            do j=1,mint3d
               c(1,1)=yi(1,j,i)
               c(2,2)=yi(2,j,i)
               c(3,3)=yi(3,j,i)
               c(1,2)=yi(4,j,i)
               c(1,3)=yi(5,j,i)
               c(2,3)=yi(6,j,i)
!     
!              calculate the eigenvalues
!
!              al(1): smallest eigenvalue
!              al(3): largest eigenvalue
!     
               call calceigenvalues(c,al)
!     
               wpsmin1=min(wpsmin1,al(1))
               wpsmax1=max(wpsmax1,al(1))
               wpsmin3=min(wpsmin3,al(3))
               wpsmax3=max(wpsmax3,al(3))
!
            enddo
!
!           check which eigenvalue is the largest in
!           absolute value       
!
            if(wpsmax3.ge.-wpsmin1) then
               wpsmin=wpsmin3
               wpsmax=wpsmax3
            else
               wpsmin=wpsmin1
               wpsmax=wpsmax1
            endif
!
            absdiff=wpsmax-wpsmin
            if(max(dabs(wpsmax),dabs(wpsmin)).lt.1.d-30) then
               reldiff=0.d0
            else
               reldiff=absdiff/(max(dabs(wpsmax),dabs(wpsmin)))
            endif
         else
!
!           thermal calculation: temperature at the vertex nodes
!            
            sizemin=1.e30
            sizemax=0.d0
!
c            do j=1,mint3d
c               c(1,1)=yi(1,j,i)
c               c(2,2)=yi(2,j,i)
c               c(3,3)=yi(3,j,i)
c!
c               size=dsqrt(c(1,1)**2+c(2,2)**2+c(3,3)**2)
c               sizemin=min(sizemin,size)
c               sizemax=max(sizemax,size)
c!
c            enddo
            do j=1,nopev
               size=vold(0,kon(indexe+j))
               sizemin=min(sizemin,size)
               sizemax=max(sizemax,size)
            enddo
            absdiff=sizemax-sizemin
            if(max(sizemax,sizemin).lt.1.d-30) then
               reldiff=0.d0
            else
               reldiff=absdiff/(max(dabs(sizemax),dabs(sizemin)))
            endif
         endif
!
!        transferring the maximum to the nodes belonging to the
!        element
!
         do j=1,nopev
            yn(2,kon(indexe+j))=max(yn(2,kon(indexe+j)),reldiff)
         enddo
!
      enddo
!
!     converting the error estimator into a stress error (%)
!     through heuristic relationships
!
!     not covered: C3D8R* and C3D6* elements
!
      if(nterms.eq.6) then
         do i=1,ne
!     
            if(ipkon(i).lt.0) cycle
            indexe=ipkon(i)
            lakonl=lakon(i)
!     
            if(lakonl(7:8).eq.'LC') cycle
!     
            if(lakonl(1:1).eq.'F') then
               cycle
            elseif(lakonl(4:4).eq.'2') then
               nopev=8
            elseif(lakonl(4:4).eq.'8') then
               nopev=8
            elseif(lakonl(4:5).eq.'10') then
               nopev=4
            elseif(lakonl(4:4).eq.'4') then
               nopev=4
            elseif(lakonl(4:5).eq.'15') then
               nopev=6
            elseif(lakonl(4:4).eq.'6') then
               nopev=6
            else
               cycle
            endif
!     
            if(lakonl(4:5).eq.'10') then
!
!              10-node tetrahedral element
!
               do j=1,nopev
                  node=kon(indexe+j)
                  yn(1,node)=max(yn(1,node),30.000d0*yn(2,node))
               enddo
            elseif((lakonl(4:7).eq.'20  ').or.
     &              (lakonl(4:7).eq.'20 L').or.
     &              (lakonl(4:7).eq.'20 B')) then
!     
!              true 20-node brick element or S8 or B32 (shell/beam)
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  if(yn(2,node).le.0.26d0) then
                     yn(1,node)=max(yn(1,node),5.115d0*yn(2,node))
                  else
                     yn(1,node)=max(yn(1,node),
     &                              27.792d0*yn(2,node)-5.895d0)
                  endif
               enddo
            elseif(lakonl(4:6).eq.'20 ') then
!     
!              expanded 20-node brick element (plane stress,
!              plane strain or axisymmetric)
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  if(yn(2,node).le.0.325d0) then
                     yn(1,node)=max(yn(1,node),9.538d0*yn(2,node))
                  else
                     yn(1,node)=
     &                    max(yn(1,node),53.695d0*yn(2,node)-14.351d0)
                  endif
               enddo
            elseif((lakonl(4:7).eq.'20R ').or.
     &              (lakonl(4:7).eq.'20RL').or.
     &              (lakonl(4:7).eq.'20RB')) then
!     
!              true 20-node brick element with reduced integration
!              or S8R or B32R (shell/beam)
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  if(yn(2,node).le.0.18d0) then
                     yn(1,node)=max(yn(1,node),20.278d0*yn(2,node))
                  else
                     yn(1,node)=max(yn(1,node),
     &                              74.318d0*yn(2,node)-9.727d0)
                  endif
               enddo
            elseif(lakonl(4:6).eq.'20R') then
!     
!              expanded 20-node brick element with reduced integration
!              (plane stress, plane strain or axisymmetric)
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  yn(1,node)=max(yn(1,node),54.054d0*yn(2,node))
               enddo
            elseif(lakonl(4:5).eq.'8I') then
!     
!              true C3D8I-element or S4 or BE31 (shell/beam)
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  if(yn(2,node).le.0.165d0) then
                     yn(1,node)=max(yn(1,node),30.303d0*yn(2,node))
                  else
                     yn(1,node)=
     &                    max(yn(1,node),139.535d0*yn(2,node)-18.023d0)
                  endif
               enddo
            elseif(lakonl(4:7).eq.'8   ') then
!     
!              true 8-node brick element
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  if(yn(2,node).le.0.157d0) then
                     yn(1,node)=max(yn(1,node),31.847d0*yn(2,node))
                  else
                     yn(1,node)=max(yn(1,node),
     &                              85.324d0*yn(2,node)-8.396d0)
                  endif
               enddo
            elseif(lakonl(4:5).eq.'8 ') then
!     
!              expanded 8-node brick element (plane stress, plane strain
!              or axisymmetric)
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  yn(1,node)=max(yn(1,node),74.074d0*yn(2,node))
               enddo
            elseif(lakonl(4:5).eq.'15') then
!     
!              true or expanded 15-node wedge element
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  yn(1,node)=max(yn(1,node),46.189d0*yn(2,node))
               enddo
            endif
         enddo
!
!     converting the error estimator into a temperature error (%)
!     through heuristic relationships
!
!     not covered: C3D8R, C3D8I and C3D15 elements
!
      elseif(nterms.eq.3) then
         do i=1,ne
!     
            if(ipkon(i).lt.0) cycle
            indexe=ipkon(i)
            lakonl=lakon(i)
!     
            if(lakonl(7:8).eq.'LC') cycle
!     
            if(lakonl(1:1).eq.'F') then
               cycle
            elseif(lakonl(4:4).eq.'2') then
               nopev=8
            elseif(lakonl(4:4).eq.'8') then
               nopev=8
            elseif(lakonl(4:5).eq.'10') then
               nopev=4
            elseif(lakonl(4:4).eq.'4') then
               nopev=4
            elseif(lakonl(4:5).eq.'15') then
               nopev=6
            elseif(lakonl(4:4).eq.'6') then
               nopev=6
            else
               cycle
            endif
!     
            if(lakonl(4:4).eq.'4') then
!
!              4-node tetrahedral element
!
               do j=1,nopev
                  node=kon(indexe+j)
                  yn(1,node)=max(yn(1,node),
     &                           12.5064d0*yn(2,node)+0.62694d0)
               enddo
            elseif(lakonl(4:5).eq.'10') then
!
!              10-node tetrahedral element
!
               do j=1,nopev
                  node=kon(indexe+j)
                  yn(1,node)=max(yn(1,node),
     &                           6.6739d0*yn(2,node)+0.0668d0)
               enddo
            elseif((lakonl(4:7).eq.'20  ').or.
     &              (lakonl(4:7).eq.'20 L').or.
     &              (lakonl(4:7).eq.'20 B')) then
!     
!              true 20-node brick element or S8 or B32 (shell/beam)
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  yn(1,node)=max(yn(1,node),
     &                           3.0980d0*yn(2,node)+0.0744d0)
               enddo
            elseif(lakonl(4:6).eq.'20 ') then
!     
!              expanded 20-node brick element (plane stress,
!              plane strain or axisymmetric)
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  yn(1,node)=max(yn(1,node),2.791d0*yn(2,node)+0.146d0)
               enddo
            elseif(lakonl(4:7).eq.'20R ') then
!     
!              true 20-node brick element with reduced integration
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  yn(1,node)=max(yn(1,node),
     &                           2.6085d0*yn(2,node)+0.03934d0)
               enddo
            elseif((lakonl(4:7).eq.'20RL').or.
     &              (lakonl(4:7).eq.'20RB')) then
!     
!              S8R or B32R (shell/beam)
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  yn(1,node)=max(yn(1,node),2.566d0*yn(2,node)+0.0485d0)
               enddo
            elseif(lakonl(4:6).eq.'20R') then
!     
!              expanded 20-node brick element with reduced integration
!              (plane stress, plane strain or axisymmetric)
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  yn(1,node)=max(yn(1,node),2.337d0*yn(2,node)+0.0985d0)
               enddo
            elseif(lakonl(4:5).eq.'8I') then
!     
!              true C3D8I element S4 or BE31 (shell/beam)
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  if(yn(2,node).le.0.0490d0) then
                     yn(1,node)=max(yn(1,node),
     &                              7.3703d0*yn(2,node)+0.413d0)
                  else
                     yn(1,node)=
     &                    max(yn(1,node),3.8261d0*yn(2,node)+0.587d0)
                  endif
               enddo
            elseif(lakonl(4:7).eq.'8   ') then
!     
!              true 8-node brick element
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  if(yn(2,node).le.0.0412d0) then
                     yn(1,node)=max(yn(1,node),
     &                              8.069d0*yn(2,node)+0.351d0)
                  else
                     yn(1,node)=max(yn(1,node),
     &                              3.858d0*yn(2,node)+0.525d0)
                  endif
               enddo
            elseif(lakonl(4:5).eq.'8 ') then
!     
!              expanded 8-node brick element (plane stress, plane strain
!              or axisymmetric)
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  if(yn(2,node).le.0.0412d0) then
                     yn(1,node)=max(yn(1,node),
     &                              8.069d0*yn(2,node)+0.351d0)
                  else
                     yn(1,node)=max(yn(1,node),
     &                              3.858d0*yn(2,node)+0.525d0)
                  endif
               enddo
            elseif(lakonl(4:5).eq.'6  L') then
!     
!              S3-element expanded into C3D6
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  if(yn(2,node).le.0.0203d0) then
                     yn(1,node)=max(yn(1,node),
     &                              26.951d0*yn(2,node)+0.514d0)
                  elseif(yn(2,node).le.0.0608d0) then
                     yn(1,node)=max(yn(1,node),
     &                              3.0485d0*yn(2,node)+0.999d0)
                  else
                     yn(1,node)=max(yn(1,node),
     &                              10.335d0*yn(2,node)+0.556d0)
                  endif
               enddo
            elseif(lakonl(4:4).eq.'6') then
!     
!              true or expanded 6-node wedge element
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  yn(1,node)=max(yn(1,node),
     &                           13.3096d0*yn(2,node)+0.179d0)
               enddo
            elseif((lakonl(4:5).eq.'15 A').or.
     &             (lakonl(4:5).eq.'15 E').or.
     &             (lakonl(4:5).eq.'15 S')) then
!     
!              CAX6, CPS6 or CPE6 element
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  yn(1,node)=max(yn(1,node),3.576d0*yn(2,node)+0.1552d0)
               enddo
            elseif(lakonl(4:5).eq.'15') then
!     
!              true or expanded 15-node wedge element
!     
               do j=1,nopev
                  node=kon(indexe+j)
                  yn(1,node)=max(yn(1,node),
     &                           3.5303d0*yn(2,node)+0.15927d0)
               enddo
            endif
         enddo
      endif
!     
!     determining the field values in the midside nodes
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         indexe=ipkon(i)
         lakonl=lakon(i)
!
         if(lakonl(7:8).eq.'LC') cycle
!
         if(lakonl(4:5).eq.'20') then
            do j=9,20
               do k=1,2
                  yn(k,kon(indexe+j))=(
     &                 yn(k,kon(indexe+nonei20(2,j-8)))+
     &                 yn(k,kon(indexe+nonei20(3,j-8))))/2.d0
               enddo
            enddo
         elseif(lakonl(4:5).eq.'10') then
            do j=5,10
               do k=1,2
                  yn(k,kon(indexe+j))=(yn(k,kon(indexe+nonei10(2,j-4)))+
     &                 yn(k,kon(indexe+nonei10(3,j-4))))/2.d0
               enddo
            enddo
         elseif(lakonl(4:5).eq.'15') then
            do j=7,15
               do k=1,2
                  yn(k,kon(indexe+j))=(
     &                 yn(k,kon(indexe+nonei15(2,j-6)))+
     &                 yn(k,kon(indexe+nonei15(3,j-6))))/2.d0
               enddo
            enddo
         endif
      enddo
!
!     mapping 3D on 1D/2D
!
      if(cflag.eq.'I') then
         force=.false.
         call map3dto1d2d(yn,ipkon,inum,kon,lakon,nterms,nk,
     &  ne,cflag,co,vold,force,mi,ielprop,prop)
      endif
!     
      return
      end
