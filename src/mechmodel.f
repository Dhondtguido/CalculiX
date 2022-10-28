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
      subroutine mechmodel(elconloc,elas,emec,kode,emec0,ithermal,
     &     icmd,beta,stre,xkl,ckl,vj,xikl,vij,plconloc,xstate,xstateini,
     &     ielas,amat,t1l,dtime,time,ttime,iel,iint,nstate_,mi,
     &     iorien,pgauss,orab,eloc,mattyp,pnewdt,istep,iinc,ipkon,
     &     nmethod,iperturb,depvisc,nlgeom_undo,physcon)
!
!     kode=-1: Arruda-Boyce
!          -2: Mooney-Rivlin
!          -3: Neo-Hooke
!          -4: Ogden (N=1)
!          -5: Ogden (N=2)
!          -6: Ogden (N=3)
!          -7: Polynomial (N=1)
!          -8: Polynomial (N=2)
!          -9: Polynomial (N=3)
!          -10: Reduced Polynomial (N=1)
!          -11: Reduced Polynomial (N=2)
!          -12: Reduced Polynomial (N=3)
!          -13: Van der Waals (not implemented yet)
!          -14: Yeoh
!          -15: Hyperfoam (N=1)
!          -16: Hyperfoam (N=2)
!          -17: Hyperfoam (N=3)
!          -50: deformation plasticity
!          -51: incremental plasticity (no viscosity)
!          -52: viscoplasticity
!       < -100: user material routine with -kode-100 user
!               defined constants with keyword *USER MATERIAL
!
      implicit none
!
      character*80 amat
!
      integer kode,ithermal(*),icmd,ielas,iel,iint,nstate_,mi(*),iorien,
     &  mattyp,istep,iinc,ipkon(*),nmethod,iperturb(*),nlgeom_undo
!
      real*8 elconloc(*),elas(21),emec(*),emec0(*),beta(*),stre(*),
     &  ckl(*),vj,plconloc(*),t1l,xkl(*),xikl(*),vij,depvisc,
     &  dtime,didc(27),d2idc2(243),dibdc(27),d2ibdc2(243),
     &  dudc(9),d2udc2(81),dldc(27),d2ldc2(243),dlbdc(27),d2lbdc2(243),
     &  pgauss(3),orab(7,*),time,ttime,eloc(6),pnewdt,physcon(*)
!
      real*8 xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*)
!
      if(kode.gt.0) then
         call linel(kode,mattyp,beta,emec,stre,elas,elconloc,
     &  iorien,orab,pgauss)
      elseif(kode.gt.-50) then
         mattyp=3
         call rubber(elconloc,elas,emec,kode,didc,d2idc2,
     &     dibdc,d2ibdc2,dudc,d2udc2,dldc,d2ldc2,dlbdc,d2lbdc2,
     &     ithermal,icmd,beta,stre)
      elseif(kode.eq.-50) then
         mattyp=3
c         call defplas(elconloc,elas,emec,ithermal,icmd,beta,stre,
c     &     ckl,vj,xstate,nstate_,iel,iint,mi)
         call umat_abaqusnl_total(amat,iel,iint,kode,elconloc,emec,
     &        emec0,beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),nstate_,xstateini,xstate,stre,elas,
     &        iorien,pgauss,orab,istep,iinc,pnewdt,nmethod,iperturb)
      elseif(kode.gt.-100) then
         mattyp=3
         if(iperturb(2).eq.1) then
            call incplas(elconloc,plconloc,xstate,xstateini,elas,emec,
     &           ithermal,icmd,beta,stre,vj,kode,ielas,amat,t1l,dtime,
     &           time,ttime,iel,iint,nstate_,mi(1),eloc,pgauss,nmethod,
     &           pnewdt,depvisc)
         else
            call incplas_lin(elconloc,plconloc,xstate,xstateini,elas,
     &           emec,
     &           ithermal,icmd,beta,stre,vj,kode,ielas,amat,t1l,dtime,
     &           time,ttime,iel,iint,nstate_,mi(1),eloc,pgauss,nmethod,
     &           pnewdt,depvisc)
         endif
      else
         mattyp=3
         call umat_main(amat,iel,iint,kode,elconloc,emec,emec0,beta,
     &        xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,icmd,ielas,
     &        mi(1),nstate_,xstateini,xstate,stre,elas,iorien,pgauss,
     &        orab,pnewdt,istep,iinc,ipkon,nmethod,iperturb,depvisc,
     &        eloc,nlgeom_undo,physcon)
      endif
!
      return
      end
