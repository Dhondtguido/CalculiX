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
      subroutine postinitialnet(ieg,lakon,v,ipkon,kon,nflow,prop,
     &     ielprop,ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,mi,iponoel,
     &     inoel,itg,ntg,nactdog)
!
!     this routine only applies to compressible networks
!
!     determination of initial values based on the boundary conditions
!     and the initial values given by the user by propagating these
!     through the network (using information on the mass flow direction
!     derived from unidirectional network elements or mass flow given
!     by the user (boundary conditions or initial conditions)) 
!
!     it is assumed that mass flows cannot
!     be identically zero (a zero mass flow leads to convergence problems).
!
!     This routine is used for elements in which the pressure gradient
!     does not allow to determine the mass flow, e.g. the free vortex,
!     the forced vortex and the rotating pipe
!     
      implicit none
!
      character*8 lakon(*)
!           
      integer mi(*),ieg(*),nflow,i,ielmat(mi(3),*),ntmat_,node1,node2,
     &     nelem,index,nshcon(*),ipkon(*),kon(*),nodem,imat,ielprop(*),
     &     nrhcon(*),neighbor,ichange,iponoel(*),inoel(2,*),indexe,
     &     itg(*),ntg,j,nactdog(0:3,*)
!     
      real*8 prop(*),shcon(0:3,ntmat_,*),xflow,v(0:mi(2),*),cp,r,
     &     dvi,rho,rhcon(0:1,ntmat_,*),kappa,cti,Ti,ri,ro,p1zp2,omega,
     &     p2zp1
!
c      write(*,*) 'postinitialnet '
c      do i=1,ntg
c         write(*,'(i10,3(1x,e11.4))') itg(i),(v(j,itg(i)),j=0,2)
c      enddo
!
      do
         ichange=0
!
!        propagation of the mass flow through the network
!
         do i=1,nflow
            nelem=ieg(i)
            indexe=ipkon(nelem)
            nodem=kon(indexe+2)
!
            if((dabs(v(1,nodem)).le.0.d0).and.
     &         (nactdog(1,nodem).ne.0)) then
!
!              no initial mass flow given yet
!              check neighbors for mass flow (only if not
!              branch nor joint)
!
!              first end node
!
               node1=kon(indexe+1)
!
               if(node1.ne.0) then
                  index=iponoel(node1)
!
                  if(inoel(2,inoel(2,index)).eq.0) then
!
!                 no branch nor joint; determine neighboring element
!
                     if(inoel(1,index).eq.nelem) then
                        neighbor=inoel(1,inoel(2,index))
                     else
                        neighbor=inoel(1,index)
                     endif
!
!                 initial mass flow in neighboring element
!
                     xflow=v(1,kon(ipkon(neighbor)+2))
!
                     if(dabs(v(1,nodem)).gt.0.d0) then
!
!                    propagate initial mass flow
!
                        if(dabs(xflow).gt.0.d0) then
                           v(1,nodem)=xflow
                           ichange=1
                           cycle
                        endif
                     else
!
!                    propagate only the sign of the mass flow
!
                        if(dabs(xflow).gt.0.d0) then
                           v(1,nodem)=xflow
                           ichange=1
                           cycle
                        endif
                     endif
                  endif
               endif
!
!              second end node
!
               node2=kon(indexe+3)
!
               if(node2.ne.0) then
                  index=iponoel(node2)
!
                  if(inoel(2,inoel(2,index)).eq.0) then
!
!                 no branch nor joint; determine neighboring element
!
                     if(inoel(1,index).eq.nelem) then
                        neighbor=inoel(1,inoel(2,index))
                     else
                        neighbor=inoel(1,index)
                     endif
!
!                 initial mass flow in neighboring element
!
                     xflow=v(1,kon(ipkon(neighbor)+2))
!
                     if(dabs(v(1,nodem)).gt.0.d0) then
!
!                    propagate initial mass flow
!
                        if(dabs(xflow).gt.0.d0) then
                           v(1,nodem)=xflow
                           ichange=1
                           cycle
                        endif
                     else
!
!                    propagate only the sign of the mass flow
!
                        if(dabs(xflow).gt.0.d0) then
                           v(1,nodem)=xflow
                           ichange=1
                           cycle
                        endif
                     endif
                  endif
               endif
            endif
         enddo
c         write(*,*) 'postinitialnet '
c         do i=1,ntg
c            write(*,'(i10,3(1x,e11.4))') itg(i),(v(j,itg(i)),j=0,2)
c         enddo
         if(ichange.eq.0) exit
      enddo
!     
      return
      end
      
      
