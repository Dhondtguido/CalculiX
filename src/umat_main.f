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
      subroutine umat_main(amat,iel,iint,kode,elconloc,emec,emec0,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi,nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,pnewdt,istep,iinc,ipkon,nmethod,
     &        iperturb,depvisc,eloc,nlgeom_undo,physcon)
!
!     calculates stiffness and stresses for a user defined material
!     law
!
      implicit none
!
      character*80 amat,amatloc
!
      integer ithermal(*),icmd,kode,ielas,iel,iint,nstate_,mi(*),iorien,
     &  istep,iinc,ipkon(*),nmethod,iperturb(*),nlgeom_undo
!
      real*8 elconloc(*),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &  vj,t1l,dtime,xkl(3,3),xikl(3,3),vij,pgauss(3),orab(7,*),
     &  time,ttime,pnewdt,depvisc,eloc(6),physcon(*)
!
      real*8 xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*)
!
      if(amat(1:8).eq.'ABAQUSNL') then

         amatloc(1:72)=amat(9:80)
         amatloc(73:80)='        '
         call umat_abaqusnl(amatloc,iel,iint,kode,elconloc,emec,
     &        emec0,beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,istep,iinc,pnewdt,nmethod,iperturb)
!         
      elseif(amat(1:9).eq.'@ABAQUSNL') then
!
         call umat_abaqusnl(amat,iel,iint,kode,elconloc,emec,
     &        emec0,beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,istep,iinc,pnewdt,nmethod,iperturb)
!
      elseif(amat(1:6).eq.'ABAQUS') then
!
         amatloc(1:74)=amat(7:80)
         amatloc(75:80)='      '
         call umat_abaqus(amatloc,iel,iint,kode,elconloc,emec,
     &        emec0,beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,istep,iinc,pnewdt,nmethod,iperturb)
!
      elseif(amat(1:7).eq.'@ABAQUS') then

         call umat_abaqus(amat,iel,iint,kode,elconloc,emec,
     &        emec0,beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,istep,iinc,pnewdt,nmethod,iperturb)
!        
      elseif(amat(1:10).eq.'ANISO_PLAS') then
!
         amatloc(1:70)=amat(11:80)
         amatloc(71:80)='          '
         call umat_aniso_plas(amatloc,
     &        iel,iint,kode,elconloc,emec,emec0,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,nmethod,pnewdt)
!
      elseif(amat(1:11).eq.'ANISO_CREEP') then
!
         amatloc(1:69)=amat(12:80)
         amatloc(70:80)='           '
         call umat_aniso_creep(amatloc,
     &        iel,iint,kode,elconloc,emec,emec0,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,nmethod,pnewdt,depvisc)
!
      elseif(amat(1:10).eq.'CIARLET_EL') then
!
         amatloc(1:70)=amat(11:80)
         amatloc(71:80)='          '
         call umat_ciarlet_el(amatloc,
     &        iel,iint,kode,elconloc,emec,
     &        emec0,beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),
     &        nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,orab)
!
      elseif(amat(1:16).eq.'COMPRESSION_ONLY') then
!
         amatloc(1:64)=amat(17:80)
         amatloc(65:80)='                '
         call umat_compression_only(amatloc,
     &        iel,iint,kode,elconloc,emec,
     &        emec0,beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),
     &        nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,orab)
!
      elseif(amat(1:13).eq.'ELASTIC_FIBER') then
!
         amatloc(1:67)=amat(14:80)
         amatloc(68:80)='             '
         call umat_elastic_fiber(amat(14:80),
     &        iel,iint,kode,elconloc,emec,emec0,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),
     &        nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,orab)
!
      elseif(amat(1:11).eq.'JOHNSONCOOK') then
!
         call umat_abaqusnl_total(amat,iel,iint,kode,elconloc,emec,
     &        emec0,beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,istep,iinc,pnewdt,nmethod,iperturb)
!
      elseif(amat(1:12).eq.'LIN_EL_COROT') then
!
         amatloc(1:68)=amat(13:80)
         amatloc(69:80)='            '
         call umat_lin_el_corot(amatloc,iel,iint,kode,
     &        elconloc,emec,emec0,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,eloc,nlgeom_undo)
!
      elseif(amat(1:10).eq.'LIN_ISO_EL') then
!
         amatloc(1:70)=amat(11:80)
         amatloc(71:80)='          '
         call umat_lin_iso_el(amatloc,
     &        iel,iint,kode,elconloc,emec,emec0,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),
     &        nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,orab)
!
      elseif(amat(1:9).eq.'IDEAL_GAS') then
!
         amatloc(1:71)=amat(10:80)
         amatloc(72:80)='          '
         call umat_ideal_gas(amatloc,
     &        iel,iint,kode,elconloc,emec,emec0,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),
     &        nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,orab,
     &        physcon)
!
      elseif(amat(1:20).eq.'SINGLE_CRYSTAL_CREEP') then
!
         amatloc(1:60)=amat(21:80)
         amatloc(61:80)='                    '
         call umat_single_crystal_creep(amatloc,
     &        iel,iint,kode,elconloc,emec,
     &        emec0,beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),
     &        nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,orab,
     &        pnewdt)
!
      elseif(amat(1:14).eq.'SINGLE_CRYSTAL') then
!
         amatloc(1:66)=amat(15:80)
         amatloc(67:80)='              '
         call umat_single_crystal(amatloc,
     &        iel,iint,kode,elconloc,emec,
     &        emec0,beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),
     &        nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,orab)
!
      elseif(amat(1:12).eq.'TENSION_ONLY') then
!
         amatloc(1:68)=amat(13:80)
         amatloc(69:80)='            '
         call umat_tension_only(amatloc,
     &        iel,iint,kode,elconloc,emec,
     &        emec0,beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),
     &        nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,orab)
!
      elseif(amat(1:4).eq.'USER') then
!
         amatloc(1:76)=amat(5:80)
         amatloc(77:80)='    '
         call umat_user(amatloc,iel,iint,kode,elconloc,emec,emec0,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,pnewdt,ipkon)
!
      elseif(amat(1:18).eq.'UNDO_NLGEOM_LIN_EL') then
!
         amatloc(1:62)=amat(19:80)
         amatloc(63:80)='                  '
         call umat_undo_nlgeom_lin_el(amatloc,iel,iint,kode,
     &        elconloc,emec,emec0,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,eloc,nlgeom_undo)
!
      elseif(amat(1:22).eq.'UNDO_NLGEOM_LIN_ISO_EL') then
!
         amatloc(1:58)=amat(23:80)
         amatloc(59:80)='                      '
         call umat_undo_nlgeom_lin_iso_el(amatloc,iel,iint,kode,
     &        elconloc,emec,emec0,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi(1),nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,eloc,nlgeom_undo)
!
      elseif(amat(1:1).eq.'@') then
!
         call call_external_umat_user(amat,iel,iint,kode,elconloc,
     &        emec,emec0,beta,xikl,vij,xkl,vj,ithermal,t1l,
     &        dtime,time,ttime,icmd,ielas,mi(1),nstate_,xstateini,
     &        xstate,stre,stiff,iorien,pgauss,orab,pnewdt,ipkon)
      else
         write(*,*) '*ERROR in umat: no user material subroutine'
         write(*,*) '       defined for material ',amat
         call exit(201)
      endif
!
      return
      end
