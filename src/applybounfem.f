!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
!                          !     
      subroutine applybounfem(nodeboun,ndirboun,xbounact,
     &     nk,vold,isolidsurf,xsolidsurf,ifreestream,iturbulent,
     &     vcon,shcon,nshcon,ntmat_,physcon,v,
     &     compressible,nodempc,ipompc,coefmpc,inomat,
     &     mi,ilboun,ilmpc,labmpc,coefmodmpc,iexplicit,nbouna,
     &     nbounb,nmpca,nmpcb,nfreestreama,nfreestreamb,
     &     nsolidsurfa,nsolidsurfb)
!     
!     1) applies temperature and velocity SPC's for
!     incompressible fluids (liquids)
!     materials (the pressure BC's have been applied in applybounp.f)
!     2) applies temperature, velocity and pressure SPC's for
!     compressible fluids (gases)
!     3) applies MPC's for the conservative variables    
!     
      implicit none
!     
      character*20 labmpc(*)
!     
      integer iturbulent,compressible,mi(*),ntmat_,
     &     nodeboun(*),isolidsurf(*),j,ilboun(*),ilmpc(*),
     &     ndirboun(*),nshcon(*),nk,i,node,imat,
     &     iexplicit,ifreestream(*),
     &     index,nodei,nodempc(3,*),ipompc(*),
     &     ist,ndir,ndiri,inomat(*),nref,nbouna,nbounb,
     &     nmpca,nmpcb,nfreestreama,nfreestreamb,nsolidsurfa,
     &     nsolidsurfb
!     
      real*8 rho,vold(0:mi(2),*),xbounact(*),
     &     shcon(0:3,ntmat_,*),xnorm,coefmodmpc(*),
     &     temp,xsolidsurf(*),sum,vcon(nk,0:mi(2)),physcon(*),
     &     coefmpc(*),residu,correction,xkin,xtu,sumk,sumt,
     &     dvi,v(nk,0:mi(2))
!     
      nref=0
!
!     SPC's: temperature, velocity and pressure (latter only for
!     compressible fluids)
!     
      do j=nbouna,nbounb
!     
!     monotonically increasing DOF-order
!     
        i=ilboun(j)
!     
!     pressure boundary conditions for incompressible materials
!     have already been treated
!     
        ndir=ndirboun(i)
        if((iexplicit.eq.0).and.(ndir.gt.3)) cycle
!     
!     check whether fluid node
!     
        node=nodeboun(i)
!     
!     calculating the physical variables for the node at stake
!     
        vold(ndir,node)=xbounact(i)
        if(node.le.nk) v(node,ndir)=0.d0
      enddo
!     
!     MPC's: temperature, velocity and pressure
!     
      if(compressible.eq.0) then
!
!     incompressible fluids
!
        nref=0
!     
        do j=nmpca,nmpcb
!     
!     monotonically increasing DOF-order
!     
          i=ilmpc(j)
          ist=ipompc(i)
!     
!     pressure multiple point constraints for incompressible materials
!     have already been treated
!     
          ndir=nodempc(2,ist)
!     
!     check whether fluid node
!     
          node=nodempc(1,ist)
!     
!     calculating the value of the dependent DOF of the MPC
!     
          index=nodempc(3,ist)
          if(index.eq.0) cycle
          sum=0.d0
          do
            nodei=nodempc(1,index)
            ndiri=nodempc(2,index)
            sum=sum+coefmpc(index)*vold(ndiri,nodei)
            index=nodempc(3,index)
            if(index.eq.0) exit
          enddo
!     
!     calculating the physical variables for the node at stake
!     
          vold(ndir,node)=-sum/coefmpc(ist)
          if(node.le.nk) v(node,ndir)=0.d0
        enddo
      else
!
!     compressible fluids
!     
!     MPC's are treated by distributing the residual proportional to
!     the coefficients
!     
!     Right now it is assumed that the MPC's are independent of each 
!     other, i.e. degrees of freedom used in one MPC are not used in 
!     any other MPC
!
        nref=0
!        
        do j=nmpca,nmpcb
          i=ilmpc(j)
          index=ipompc(i)
!     
!     pressure multiple point constraints for incompressible materials
!     have already been treated
!     
          ndir=nodempc(2,index)
!     
!     check whether fluid node
!     
          node=nodempc(1,index)
!     
!     calculating the value of the dependent DOF of the MPC
!     
          residu=coefmpc(index)*vold(ndir,node)
          xnorm=1.d0
          if(index.eq.0) cycle
          do
            index=nodempc(3,index)
            if(index.eq.0) exit
            nodei=nodempc(1,index)
            ndiri=nodempc(2,index)
            residu=residu+coefmpc(index)*vold(ndiri,nodei)
          enddo
!     
!     correcting all terms of the MPC
!
          index=ipompc(i)
          do
            nodei=nodempc(1,index)
            ndiri=nodempc(2,index)
!     
            correction=-residu*coefmodmpc(index)
            vold(ndiri,nodei)=vold(ndiri,nodei)+correction
            if(nodei.le.nk) v(nodei,ndiri)=0.d0
            index=nodempc(3,index)
            if(index.eq.0) exit
          enddo
        enddo
!     
      endif
!     
      if(iturbulent.ne.0) then
!     
!     freestream conditions for the iturbulent variables
!     
        xtu=10.d0*physcon(5)/physcon(8)
        xkin=10.d0**(-3.5d0)*xtu
        do j=nfreestreama,nfreestreamb
          node=ifreestream(j)
          imat=inomat(node)
          if(imat.eq.0) cycle
          temp=vold(0,node)
          call materialdata_dvifem(imat,ntmat_,temp,shcon,nshcon,dvi)
!     
!     density 
!     
          rho=vcon(node,4)
!     
          vcon(node,5)=xkin*dvi
          vcon(node,6)=xtu*rho
!          
          vold(5,node)=vcon(node,5)/rho
          vold(6,node)=xtu
!
          v(node,5)=0.d0
          v(node,6)=0.d0
        enddo
!     
!     solid boundary conditions for the turbulent variables 
!     
        do j=nsolidsurfa,nsolidsurfb
!
!         turbulent kinetic energy is applied at the wall
!
          node=isolidsurf(j)
          imat=inomat(node)
          if(imat.eq.0) cycle
          temp=vold(0,node)
          call materialdata_dvifem(imat,ntmat_,temp,shcon,nshcon,dvi)
          rho=vcon(node,4)
          vcon(node,5)=0.d0
          vcon(node,6)=800.d0*dvi/(xsolidsurf(j)**2)
!
          vold(5,node)=0.d0
          vold(6,node)=vcon(node,6)/rho
!
          v(node,5)=0.d0
          v(node,6)=0.d0
        enddo
!     
!     taking fluid pressure MPC's into account: it is assumed
!     that cyclic fluid pressure MPC's also apply to the iturbulent
!     conservative variables
!     
        do j=nmpca,nmpcb
          i=ilmpc(j)
          if(labmpc(i)(1:6).ne.'CYCLIC') cycle
          ist=ipompc(i)
          ndir=nodempc(2,ist)
          if(ndir.ne.4) cycle
          node=nodempc(1,ist)
!     
!     check whether fluid MPC
!     
          imat=inomat(node)
          if(imat.eq.0) cycle
!     
          index=nodempc(3,ist)
          if(index.eq.0) cycle
!          
          sumk=0.d0
          sumt=0.d0
!     
          do
            nodei=nodempc(1,index)
            sumk=sumk+coefmpc(index)*vcon(nodei,5)
            sumt=sumt+coefmpc(index)*vcon(nodei,6)
            index=nodempc(3,index)
            if(index.eq.0) exit
          enddo
          vcon(node,5)=-sumk/coefmpc(ist)
          vcon(node,6)=-sumt/coefmpc(ist)
        enddo
      endif
!     
      return
      end
      
