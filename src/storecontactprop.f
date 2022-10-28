!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine storecontactprop(ne,ne0,lakon,kon,ipkon,mi,
     & ielmat,elcon,mortar,adb,nactdof,springarea,
     & ncmat_,ntmat_,stx,temax)
!     
!     # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
!     Routine that computes the natural period of oscillation 
!       of the contact elements. A simplified two-dofs model
!       has been used to compute the natural frequencies due 
!       to the relative motion between slave and master surf-
!       ace. 
!     
!     Main variables and meaning
!     
!     temax       : max. natural period of oscillation 
!     stx      : vector containig results (from results.c)
!     springmm : average mass of master surface
!     springms : average mass of slave surface
!     xk : spring stiffness between the surfaces
!     xmacont : mass of the actual node of the element
!     areaslav : area of the slave surface (stiffness comp)
!     
!     Proposed by Matteo Pacher
!     
!     # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
!     
      implicit none
!     
      character*8 lakonl,lakon(*)
!     
      integer mi(*),ne,ne0,kon(*),ipkon(*),nactdof(0:mi(2),*),
     &     indexe,j,nope,ielmat(mi(3),*),imat,mortar,nopem,
     &     indexn,nelem,ncmat_,ntmat_
!     
      real*8 adb(*),springarea(2,*),pi,te,stx(6,mi(1),ne),
     &     springmm,springms,elcon(0:ncmat_,ntmat_,*),springfac,
     &     xmacont,areaslav,xk,temax
!     
!     Initialization
!
      temax=0.d0
      pi=4.d0*datan(1.d0)
      
!     Loop over the contact elements      
      do nelem=ne0+1,ne
         indexe=ipkon(nelem)
         lakonl=lakon(nelem)
         imat=ielmat(1,nelem) 
         
         xk=elcon(2,1,imat)     ! element stiffness
         springmm=0.0d0         ! mass master surface
         springms=0.0d0         ! mass slave node/surface
         
         if(mortar.eq.0) then  
!     
!     Contact node to face
!   
!           nope = number of master nodes + 1 (1 slave node)
!  
            nope=ichar(lakonl(8:8))-47
            springfac=1.d0
!     
            do indexn=1,nope                   
               xmacont=0.0d0
               
               do j=1,3           
                  if(nactdof(j,kon(indexe+indexn)).gt.0)then
                     xmacont=max(xmacont,
     &                    adb(nactdof(j,kon(indexe+indexn))))
                  endif
               enddo
!     mass accumulation
               if(indexn.eq.nope)then
                  springms=springms+xmacont
!
!                 mean mass at master nodes
!
                  springmm=springmm/(nope-1.0d0)    
               else
                  springmm=springmm+xmacont    
               endif        
            enddo
!
            areaslav=springarea(1,kon(indexe+nope+1))
            xk=xk*areaslav
!     
!           checking whether mass or slave side is constrained
!     
            if((springmm.le.0.d0).and.(springms.le.0.d0)) then
               cycle
            elseif(springmm.le.0.d0) then
               springmm=springms
            elseif(springms.le.0.d0) then
               springms=springmm
            endif
!     
         elseif(mortar.eq.1) then
!     
!           Contact face to face
!     
!           nopem = number of master nodes
!
            nopem=ichar(lakonl(8:8))-48
            springfac=0.1d0
!     
            do indexn=1,kon(indexe)  
!                 
               xmacont=0.0d0
               do j=1,3           
                  if(nactdof(j,kon(indexe+indexn)).gt.0)then
                     xmacont=max(xmacont,
     &                    adb(nactdof(j,kon(indexe+indexn))))
                  endif
               enddo
!     
               if(indexn.gt.nopem)then !slave
                  springms=springms+xmacont    
               else             !master
                  springmm=springmm+xmacont    
               endif        
            enddo
!     
!           mean mass at slave and master nodes
!
            springms=springms/(kon(indexe)-nopem)
            springmm=springmm/nopem
!     
            areaslav=springarea(1,kon(1+indexe+kon(indexe)))
            xk=xk*areaslav*springfac
!     
!           checking whether mass or slave side is constrained
!     
            if((springmm.le.0.d0).and.(springms.le.0.d0)) then
               cycle
            elseif(springmm.le.0.d0) then
               springmm=springms
            elseif(springms.le.0.d0) then
               springms=springmm
            endif
         endif                  ! face-to-face
!     
!     Calculation of the natural period according to the 2-dof model
!     
         te=0.d0
         if(xk.ne.0.d0)then
            te=2.d0*pi*dsqrt((springmm*springms)/
     &           ((springmm+springms)*xk))
!
!           exponential pressure-overclosure behavior
!
            if(int(elcon(3,1,imat)).eq.1) then
               te=te/10.d0
            endif
         endif
!
         if(stx(4,1,nelem).gt.0.d0) temax=max(te,temax)
!     
      enddo                     ! loop over the contact elements
!     
      return
      end
