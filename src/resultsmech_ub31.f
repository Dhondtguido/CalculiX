!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2025 Guido Dhondt
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
      subroutine resultsmech_ub31(co,kon,ipkon,lakon,ne,v,
     &  stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &  ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,
     &  iprestr,eme,iperturb,fn,iout,qa,vold,nmethod,
     &  veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
     &  xstateini,xstiff,xstate,npmat_,matname,mi,ielas,icmd,
     &  ncmat_,nstate_,stiini,vini,ener,eei,enerini,istep,iinc,
     &  reltime,calcul_fn,calcul_qa,calcul_cauchy,nener,
     &  ikin,nal,ne0,thicke,emeini,i,ielprop,prop,t0g,t1g)
!
!     calculates stresses and strains for user element UB31
!     Stresses stx are populated with:
!     - stx(1): combined max axial stress (axial + bending)
!     - stx(2): shear stress tau_xy
!     - stx(3): shear stress tau_xz
!     - stx(4): torsional stress tau_torsion
!     - stx(5): bending moment M_y
!     - stx(6): bending moment M_z
!
      use ub31_module
      use codecheck_module
      implicit none
!
      character*8 lakon(*)
      character*80 matname(*)
!
      integer kon(*),mi(*),
     &  nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),
     &  ielorien(mi(3),*),ntmat_,ipkon(*),ne0,
     &  istep,iinc,ne,ithermal(*),iprestr,
     &  nener,norien,iperturb(*),iout,
     &  nal,icmd,nmethod,ielas,
     &  ncmat_,nstate_,ikin,ielprop(*),
     &  nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,calcul_fn,
     &  calcul_cauchy,calcul_qa,nelem,i
!
      real*8 co(3,*),v(0:mi(2),*),stiini(6,mi(1),*),t0g(2,*),t1g(2,*),
     &     stx(6,mi(1),*),prop(*),elcon(0:ncmat_,ntmat_,*),
     &     rhcon(0:1,ntmat_,*),alcon(0:6,ntmat_,*),vini(0:mi(2),*),
     &     alzero(*),orab(7,*),fn(0:mi(2),*),t0(*),t1(*),
     &     prestr(6,mi(1),*),eme(6,mi(1),*),vold(0:mi(2),*),
     &     ener(2,mi(1),*),eei(6,mi(1),*),enerini(2,mi(1),*),
     &     veold(0:mi(2),*),qa(*),dtime,time,ttime,
     &     plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &     xstiff(27,mi(1),*),xstate(nstate_,mi(1),*),
     &     xstateini(nstate_,mi(1),*),reltime,thicke(mi(3),*),
     &     emeini(6,mi(1),*)
!
      integer k, j, imat, iorien, ihyper, index, indexe, nope, ndof
      integer kk, i1
      integer release_codes(2), code, r, node1, node2, inode, konl(2)
      integer node_id
      logical released(12)
      real*8 sect_type, dims(6), rot_angle, offsets(3,2),
     &  x_beam(3,2), e2_in(3), len_e2, j_tors, principal_angle,
     &  rot_angle_total, e, un, um, xl(3,2), e1(3), e2(3), e3(3),
     &  dl, c1, tm(3,3), t0l, t1l, elconloc(ncmat_), rho, stiff(21),
     &  coords(3), eth(6), plconloc(802),
     &  xstiff_pt(27), a, xi11, xi22, xk_y, xk_z,
     &  xl_def(3,2),e1_init(3),e1_curr(3),len_init,len_curr,
     &  cross_e1(3),dot_e1,u_rot(3),cross_u(3),dot_u_e2,len_cross,
     &  e2_init(3), e3_init(3), off_vec(3), rot_disp(3),
     &  s(12,12), tmg(12,12), u_glob(12), u_loc_node(12),
     &  u_loc_beam(12), f_local(12), axial_f, vy, vz, tx, my, mz,
     &  cy, cz, sigma_axial, tau_xy, tau_xz, tau_tor, sigma_max,
     &  diag_val, factor, x_loc, w1, w2, alpha, beta, mz_x, my_x,
     &  sigma_x, sigma_b_y, sigma_b_z, w_load_val, xi_station,
     &  m1_eq_y, m2_eq_y, m1_eq_z,
     &  m2_eq_z, M1_y, M2_y, M1_z, M2_z, M0_y, M0_z, R1, R2,
     &  phi_y, phi_z, vy_x, vz_x,
     &  spring_k(12), col_r(12), denom, ks, s_rr
      integer station_i, load_dir, load_type, station_k, num_sub
      integer station_pct, m_idx, itarg, iunit
      logical is_active_inc, is_open, has_target
      real*8 ux_val, uy_val, uz_val, rx_val, ry_val, rz_val
      real*8 h1, h2, h3, h4, dh1, dh2, dh3, dh4, vh, wh
      real*8 thetazh, thetayh, v0, w0, thetaz0, thetay0
      real*8 qx_val, qy_val, qz_val
      real*8 f_csv(3), m_csv(3), u_csv(3), rot_csv(3), q_csv(3)
      real*8 y_pts(4), z_pts(4), s_pts(4), s_max_tens, s_min_comp
      real*8 tau_vy_val, tau_vz_val, tau_tor_val, s_peak, s_vm_val
      real*8 s_bend_res, alpha_t, beta_c, alpha_l, beta_l
      real*8 s_norm_val, r_dist, tau_tor_pt, s_vm_pt
      character*16 pt_labels(4)
      integer ipt
      character*80 amat
!
!

      nelem = i
      indexe = ipkon(nelem)
      node1 = kon(indexe+1)
      node2 = kon(indexe+2)
!
      imat=ielmat(1,nelem)
      amat=matname(imat)
      if(norien.gt.0) then
         iorien=max(0,ielorien(1,nelem))
      else
         iorien=0
      endif
      if(nelcon(1,imat).lt.0) then
         ihyper=1
      else
         ihyper=0
      endif
      rho=rhcon(1,imat,1)
!
      index=ielprop(nelem)
      if (index .gt. 0) then
         sect_type=prop(index+1)
         dims(1)=prop(index+2)
         dims(2)=prop(index+3)
         dims(3)=prop(index+4)
         dims(4)=prop(index+5)
         dims(5)=prop(index+6)
         dims(6)=prop(index+7)
         rot_angle=prop(index+8)
         
         offsets(1,1)=prop(index+9)
         offsets(2,1)=prop(index+10)
         offsets(3,1)=prop(index+11)
         
         offsets(1,2)=prop(index+12)
         offsets(2,2)=prop(index+13)
         offsets(3,2)=prop(index+14)
         
         if (allocated(ielrelease)) then
            if (nelem .le. size(ielrelease, 2) .and.
     &          ielrelease(1, nelem) .ge. 0) then
               release_codes(1) = ielrelease(1, nelem)
            else
               release_codes(1) = nint(prop(index+15))
            endif
            if (nelem .le. size(ielrelease, 2) .and.
     &          ielrelease(2, nelem) .ge. 0) then
               release_codes(2) = ielrelease(2, nelem)
            else
               release_codes(2) = nint(prop(index+16))
            endif
         else
            release_codes(1) = nint(prop(index+15))
            release_codes(2) = nint(prop(index+16))
         endif
         
         e2_in(1)=prop(index+17)
         e2_in(2)=prop(index+18)
         e2_in(3)=prop(index+19)
      else
         write(*,*) '*ERROR in resultsmech_ub31: element ',nelem,
     &        ' has no section properties.'
         write(*,*) '  Use *USER SECTION or *BEAM SECTION'
     &        ,' to define cross-section.'
         call exit(201)
      endif
!
      nope=2
      ndof=6
!
      do k=1,nope
         konl(k)=kon(indexe+k)
         do j=1,3
            xl(j,k)=co(j,konl(k))
         enddo
      enddo
!
      t0l=0.d0
      t1l=0.d0
      if(ithermal(1).eq.1) then
         do k=1,nope
            t0l=t0l+t0(konl(k))/2.d0
            t1l=t1l+t1(konl(k))/2.d0
         enddo
      endif
!
      call materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &     imat,amat,iorien,coords,orab,ntmat_,stiff,rho,
     &     nelem,ithermal,alzero,0,t0l,t1l,
     &     ihyper,0,elconloc,eth,nelcon(1,imat),plicon,
     &     nplicon,plkcon,nplkcon,npmat_,
     &     plconloc,mi(1),dtime,kk,
     &     xstiff_pt,ncmat_,iperturb)
!
      e=elconloc(1)
      un=elconloc(2)
      um=e/(2.d0*(1.d0+un))
!
      if (allocated(elem_E)) then
         if (nelem .le. size(elem_E)) then
            elem_E(nelem) = e
            elem_nu(nelem) = un
            if (npmat_ .gt. 0 .and. ntmat_ .gt. 0) then
               if (plicon(1,1,imat) .gt. 0.d0) then
                  elem_fy(nelem) = plicon(1,1,imat)
               endif
            endif
         endif
      endif
!
      call compute_section_properties(sect_type, dims, e, un, 
     &  a, xi22, xi11, j_tors, xk_y, xk_z, principal_angle)
!
      do j=1,3
         e1(j)=xl(j,2)-xl(j,1)
      enddo
      dl=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
      do j=1,3
         e1(j)=e1(j)/dl
      enddo
!
      len_e2=dsqrt(e2_in(1)**2 + e2_in(2)**2 + e2_in(3)**2)
      if (len_e2 .gt. 1.d-6) then
         do j=1,3
            e2(j)=e2_in(j)/len_e2
         enddo
         c1 = e2(1)*e1(1) + e2(2)*e1(2) + e2(3)*e1(3)
         do j=1,3
            e2(j) = e2(j) - c1*e1(j)
         enddo
         len_e2=dsqrt(e2(1)**2 + e2(2)**2 + e2(3)**2)
         if (len_e2 .lt. 1.d-6) then
            if (dabs(e1(3)) .lt. 0.9999d0) then
               e2(1)=-e1(2)
               e2(2)= e1(1)
               e2(3)= 0.d0
            else
               e2(1)=0.d0
               e2(2)=1.d0
               e2(3)=0.d0
            endif
            len_e2=dsqrt(e2(1)**2 + e2(2)**2 + e2(3)**2)
         endif
         do j=1,3
            e2(j)=e2(j)/len_e2
         enddo
      else
         ! Auto-orientation: standard structural upright rule (Z up)
         if (dabs(e1(3)) .lt. 0.9999d0) then
            e2(1)=-e1(2)
            e2(2)= e1(1)
            e2(3)= 0.d0
            len_e2=dsqrt(e2(1)**2 + e2(2)**2)
            do j=1,3
               e2(j)=e2(j)/len_e2
            enddo
         else
            e2(1)=0.d0
            e2(2)=1.d0
            e2(3)=0.d0
         endif
      endif
!
      rot_angle_total = rot_angle + principal_angle
      if (dabs(rot_angle_total) .gt. 1.d-6) then
         call rotate_vector_local(e2, e1, rot_angle_total, e2)
      endif
!
      e3(1)=e1(2)*e2(3)-e1(3)*e2(2)
      e3(2)=e1(3)*e2(1)-e1(1)*e2(3)
      e3(3)=e1(1)*e2(2)-e1(2)*e2(1)
!
      e2(1)=e3(2)*e1(3)-e3(3)*e1(2)
      e2(2)=e3(3)*e1(1)-e3(1)*e1(3)
      e2(3)=e3(1)*e1(2)-e3(2)*e1(1)
!
      do k=1,nope
         do j=1,3
            x_beam(j,k) = xl(j,k) + offsets(1,k)*e1(j) + 
     &                    offsets(2,k)*e2(j) + offsets(3,k)*e3(j)
         enddo
      enddo
!
      do j=1,3
         e1(j)=x_beam(j,2)-x_beam(j,1)
      enddo
      dl=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
      do j=1,3
         e1(j)=e1(j)/dl
      enddo
!
      c1 = e2(1)*e1(1) + e2(2)*e1(2) + e2(3)*e1(3)
      do j=1,3
         e2(j) = e2(j) - c1*e1(j)
      enddo
      len_e2=dsqrt(e2(1)**2 + e2(2)**2 + e2(3)**2)
      do j=1,3
         e2(j)=e2(j)/len_e2
      enddo
      e3(1)=e1(2)*e2(3)-e1(3)*e2(2)
      e3(2)=e1(3)*e2(1)-e1(1)*e2(3)
      e3(3)=e1(1)*e2(2)-e1(2)*e2(1)
!
!     deformed physical beam coordinates including offset rotations
      do i1=1,nope
         node_id = konl(i1)
         off_vec(1) = x_beam(1,i1) - xl(1,i1)
         off_vec(2) = x_beam(2,i1) - xl(2,i1)
         off_vec(3) = x_beam(3,i1) - xl(3,i1)

         ! v(1:3) = translation, v(4:6) = rotation vector
         ! Cross product: Theta x Offset
         rot_disp(1) = v(5,node_id)*off_vec(3) -
     &                 v(6,node_id)*off_vec(2)
         rot_disp(2) = v(6,node_id)*off_vec(1) -
     &                 v(4,node_id)*off_vec(3)
         rot_disp(3) = v(4,node_id)*off_vec(2) -
     &                 v(5,node_id)*off_vec(1)

         do j=1,3
            xl_def(j,i1) = x_beam(j,i1) + v(j,node_id) + rot_disp(j)
         enddo
      enddo

!     initial unit chord vector (node 1 to last node)
      do j=1,3
         if (nope.eq.2) then
            e1_init(j)=x_beam(j,2)-x_beam(j,1)
         else
            e1_init(j)=x_beam(j,3)-x_beam(j,1)
         endif
      enddo
      len_init=dsqrt(e1_init(1)**2+e1_init(2)**2+e1_init(3)**2)
      do j=1,3
         e1_init(j)=e1_init(j)/len_init
      enddo

!     current unit chord vector
      do j=1,3
         if (nope.eq.2) then
            e1_curr(j)=xl_def(j,2)-xl_def(j,1)
         else
            e1_curr(j)=xl_def(j,3)-xl_def(j,1)
         endif
      enddo
      len_curr=dsqrt(e1_curr(1)**2+e1_curr(2)**2+e1_curr(3)**2)
      do j=1,3
         e1_curr(j)=e1_curr(j)/len_curr
      enddo

!     initial e2_init
      len_e2=dsqrt(e2_in(1)**2 + e2_in(2)**2 + e2_in(3)**2)
      if (len_e2 .gt. 1.d-6) then
         do j=1,3
            e2_init(j)=e2_in(j)/len_e2
         enddo
         c1 = e2_init(1)*e1_init(1)+e2_init(2)*e1_init(2)+
     &        e2_init(3)*e1_init(3)
         do j=1,3
            e2_init(j)=e2_init(j)-c1*e1_init(j)
         enddo
         len_e2=dsqrt(e2_init(1)**2+e2_init(2)**2+e2_init(3)**2)
         if (len_e2 .lt. 1.d-6) then
            if (dabs(e1_init(3)) .lt. 0.9999d0) then
               e2_init(1)=-e1_init(2)
               e2_init(2)= e1_init(1)
               e2_init(3)= 0.d0
            else
               e2_init(1)=0.d0
               e2_init(2)=1.d0
               e2_init(3)=0.d0
            endif
            len_e2=dsqrt(e2_init(1)**2+e2_init(2)**2+e2_init(3)**2)
         endif
         do j=1,3
            e2_init(j)=e2_init(j)/len_e2
         enddo
      else
         if (dabs(e1_init(3)) .lt. 0.9999d0) then
            e2_init(1)=-e1_init(2)
            e2_init(2)= e1_init(1)
            e2_init(3)= 0.d0
            len_e2=dsqrt(e2_init(1)**2+e2_init(2)**2)
            do j=1,3
               e2_init(j)=e2_init(j)/len_e2
            enddo
         else
            e2_init(1)=0.d0
            e2_init(2)=1.d0
            e2_init(3)=0.d0
         endif
      endif
      if (dabs(rot_angle_total) .gt. 1.d-6) then
         call rotate_vector_local(e2_init, e1_init, rot_angle_total,
     &                            e2_init)
      endif

!     rotate e2_init to e2_curr
      cross_e1(1)=e1_init(2)*e1_curr(3)-e1_init(3)*e1_curr(2)
      cross_e1(2)=e1_init(3)*e1_curr(1)-e1_init(1)*e1_curr(3)
      cross_e1(3)=e1_init(1)*e1_curr(2)-e1_init(2)*e1_curr(1)
      len_cross=dsqrt(cross_e1(1)**2+cross_e1(2)**2+cross_e1(3)**2)
      dot_e1=e1_init(1)*e1_curr(1)+e1_init(2)*e1_curr(2)+
     &       e1_init(3)*e1_curr(3)

      if (len_cross .gt. 1.d-6) then
         do j=1,3
            u_rot(j)=cross_e1(j)/len_cross
         enddo
         cross_u(1)=u_rot(2)*e2_init(3)-u_rot(3)*e2_init(2)
         cross_u(2)=u_rot(3)*e2_init(1)-u_rot(1)*e2_init(3)
         cross_u(3)=u_rot(1)*e2_init(2)-u_rot(2)*e2_init(1)
         dot_u_e2=u_rot(1)*e2_init(1)+u_rot(2)*e2_init(2)+
     &            u_rot(3)*e2_init(3)
         do j=1,3
            e2(j)=e2_init(j)*dot_e1 + cross_u(j)*len_cross +
     &            u_rot(j)*dot_u_e2*(1.d0 - dot_e1)
         enddo
      else
         do j=1,3
            e2(j)=e2_init(j)
         enddo
      endif

!     orthogonalize e2 with respect to e1_curr
      c1 = e2(1)*e1_curr(1)+e2(2)*e1_curr(2)+e2(3)*e1_curr(3)
      do j=1,3
         e2(j)=e2(j)-c1*e1_curr(j)
      enddo
      len_e2=dsqrt(e2(1)**2+e2(2)**2+e2(3)**2)
      do j=1,3
         e2(j)=e2(j)/len_e2
      enddo

!     deformed e3
      e3(1)=e1_curr(2)*e2(3)-e1_curr(3)*e2(2)
      e3(2)=e1_curr(3)*e2(1)-e1_curr(1)*e2(3)
      e3(3)=e1_curr(1)*e2(2)-e1_curr(2)*e2(1)

!     populate tm with axes (initial un-deformed for linear static)
      if (iperturb(1) .le. 1) then
         e3_init(1)=e1_init(2)*e2_init(3)-e1_init(3)*e2_init(2)
         e3_init(2)=e1_init(3)*e2_init(1)-e1_init(1)*e2_init(3)
         e3_init(3)=e1_init(1)*e2_init(2)-e1_init(2)*e2_init(1)
         do j=1,3
            tm(1,j)=e1_init(j)
            tm(2,j)=e2_init(j)
            tm(3,j)=e3_init(j)
         enddo
      else
         do j=1,3
            tm(1,j)=e1_curr(j)
            tm(2,j)=e2(j)
            tm(3,j)=e3(j)
         enddo
      endif
!
      do j=1,12
         do k=1,12
            s(j,k)=0.d0
         enddo
      enddo
!
      s(1,1) = e*a/dl
      s(1,7) = -s(1,1)
      s(7,7) = s(1,1)

      phi_z = 12.d0 * e * xi11 / (um * a * xk_y * dl**2)
      s(2,2) = 12.d0*e*xi11/((1.d0 + phi_z)*(dl**3))
      s(2,6) = 6.d0*e*xi11/((1.d0 + phi_z)*(dl**2))
      s(2,8) = -s(2,2)
      s(2,12) = s(2,6)

      s(6,6) = (4.d0 + phi_z)*e*xi11/((1.d0 + phi_z)*dl)
      s(6,8) = -6.d0*e*xi11/((1.d0 + phi_z)*(dl**2))
      s(6,12) = (2.d0 - phi_z)*e*xi11/((1.d0 + phi_z)*dl)

      s(8,8) = s(2,2)
      s(8,12) = -s(2,6)
      s(12,12) = s(6,6)

      phi_y = 12.d0 * e * xi22 / (um * a * xk_z * dl**2)
      s(3,3) = 12.d0*e*xi22/((1.d0 + phi_y)*(dl**3))
      s(3,5) = -6.d0*e*xi22/((1.d0 + phi_y)*(dl**2))
      s(3,9) = -s(3,3)
      s(3,11) = s(3,5)

      s(5,5) = (4.d0 + phi_y)*e*xi22/((1.d0 + phi_y)*dl)
      s(5,9) = -s(3,5)
      s(5,11) = (2.d0 - phi_y)*e*xi22/((1.d0 + phi_y)*dl)

      s(9,9) = s(3,3)
      s(9,11) = -s(3,5)
      s(11,11) = s(5,5)

      s(4,4) = um*j_tors/dl
      s(4,10) = -s(4,4)
      s(10,10) = s(4,4)
!
      do j=1,12
         do k=j,12
            s(k,j)=s(j,k)
         enddo
      enddo
!
!     extract global displacements and transform to local frame
!
      do k=1,6
         u_glob(k) = v(k, node1)
         u_glob(k+6) = v(k, node2)
      enddo
!
      do inode=1,2
         k = (inode-1)*6
         do j=1,3
            u_loc_node(k+j) = tm(j,1)*u_glob(k+1) + 
     &                        tm(j,2)*u_glob(k+2) + 
     &                        tm(j,3)*u_glob(k+3)
            u_loc_node(k+3+j) = tm(j,1)*u_glob(k+4) + 
     &                          tm(j,2)*u_glob(k+5) + 
     &                          tm(j,3)*u_glob(k+6)
         enddo
      enddo
!
!     apply local offsets
!
      do inode=1,2
         k = (inode-1)*6
         u_loc_beam(k+1) = u_loc_node(k+1) + 
     &     offsets(3,inode)*u_loc_node(k+5) - 
     &     offsets(2,inode)*u_loc_node(k+6)
         u_loc_beam(k+2) = u_loc_node(k+2) + 
     &     offsets(1,inode)*u_loc_node(k+6) - 
     &     offsets(3,inode)*u_loc_node(k+4)
         u_loc_beam(k+3) = u_loc_node(k+3) + 
     &     offsets(2,inode)*u_loc_node(k+4) - 
     &     offsets(1,inode)*u_loc_node(k+5)
         u_loc_beam(k+4) = u_loc_node(k+4)
         u_loc_beam(k+5) = u_loc_node(k+5)
         u_loc_beam(k+6) = u_loc_node(k+6)
      enddo
!
!     stiffness condensation for releases & semi-rigid springs
!
      do k=1,12
         released(k) = .false.
         spring_k(k) = 0.d0
      enddo
      do inode = 1, nope
         code = release_codes(inode)
         r = (inode-1)*6
         if (code .gt. 0) then
            if (iand(code, 1) .ne. 0) then
               released(r+1) = .true.
               if (allocated(relspring)) then
                  if (nelem .le. size(relspring,3)) then
                     spring_k(r+1) = relspring(1, inode, nelem)
                  endif
               endif
            endif
            if (iand(code, 2) .ne. 0) then
               released(r+2) = .true.
               if (allocated(relspring)) then
                  if (nelem .le. size(relspring,3)) then
                     spring_k(r+2) = relspring(2, inode, nelem)
                  endif
               endif
            endif
            if (iand(code, 4) .ne. 0) then
               released(r+3) = .true.
               if (allocated(relspring)) then
                  if (nelem .le. size(relspring,3)) then
                     spring_k(r+3) = relspring(3, inode, nelem)
                  endif
               endif
            endif
            if (iand(code, 8) .ne. 0) then
               released(r+4) = .true.
               if (allocated(relspring)) then
                  if (nelem .le. size(relspring,3)) then
                     spring_k(r+4) = relspring(4, inode, nelem)
                  endif
               endif
            endif
            if (iand(code, 16) .ne. 0) then
               released(r+5) = .true.
               if (allocated(relspring)) then
                  if (nelem .le. size(relspring,3)) then
                     spring_k(r+5) = relspring(5, inode, nelem)
                  endif
               endif
            endif
            if (iand(code, 32) .ne. 0) then
               released(r+6) = .true.
               if (allocated(relspring)) then
                  if (nelem .le. size(relspring,3)) then
                     spring_k(r+6) = relspring(6, inode, nelem)
                  endif
               endif
            endif
         endif
      enddo
!
      do r = 1, 12
         if (released(r)) then
            s_rr = s(r,r)
            ks = spring_k(r)
            if (ks .lt. 0.d0) ks = 0.d0
            denom = s_rr + ks
            if (dabs(denom) .gt. 1.d-12) then
               do k = 1, 12
                  col_r(k) = s(k,r)
               enddo
               do k = 1, 12
                  do j = 1, 12
                     s(k,j) = s(k,j) - (col_r(k) * col_r(j)) / denom
                  enddo
               enddo
               if (ks .le. 1.d-12) then
                  s(r,r) = 1.d-9 * s_rr
               endif
            endif
         endif
      enddo
!
!     compute internal local forces: f_local = s * u_loc_beam
!
       do k=1,12
          f_local(k) = 0.d0
          do j=1,12
             f_local(k) = f_local(k) + s(k,j) * u_loc_beam(j)
          enddo
          if (allocated(ff_saved) .and. iout.gt.0) then
             f_local(k) = f_local(k) - ff_saved(k, i)
          endif
       enddo
!
!     transform and accumulate global internal forces
      do j=1,3
         fn(j,node1) = fn(j,node1) + tm(1,j)*f_local(1) +
     &                               tm(2,j)*f_local(2) +
     &                               tm(3,j)*f_local(3)
         fn(j+3,node1) = fn(j+3,node1) + tm(1,j)*f_local(4) +
     &                                   tm(2,j)*f_local(5) +
     &                                   tm(3,j)*f_local(6)
         
         fn(j,node2) = fn(j,node2) + tm(1,j)*f_local(7) +
     &                               tm(2,j)*f_local(8) +
     &                               tm(3,j)*f_local(9)
         fn(j+3,node2) = fn(j+3,node2) + tm(1,j)*f_local(10) +
     &                                   tm(2,j)*f_local(11) +
     &                                   tm(3,j)*f_local(12)
      enddo
!
!     internal forces/moments at representative node (node 2)
!
      axial_f = f_local(7)
      vy      = f_local(8)
      vz      = f_local(9)
      tx      = f_local(10)
      my      = 0.5d0 * (f_local(5) - f_local(11))
      mz      = 0.5d0 * (-f_local(6) + f_local(12))
!
!     outer fiber coordinates for combined stress calculation
!
      if (nint(sect_type) .eq. 1) then
         cy = dims(1)/2.d0
         cz = dims(2)/2.d0
      elseif (nint(sect_type) .eq. 2) then
         cy = dims(1)
         cz = dims(1)
      elseif (nint(sect_type) .eq. 3) then
         cy = max(dims(2), dims(4))/2.d0
         cz = dims(1)/2.d0
      elseif (nint(sect_type) .eq. 4) then
         ! T-section (dims: h, b, t_f, t_w)
         cy = dims(2)/2.d0
         w1 = dims(2)*dims(3)
         w2 = dims(4)*(dims(1)-dims(3))
         alpha = (w1*(dims(1)-dims(3)/2.d0) + 
     &            w2*(dims(1)-dims(3))/2.d0)/(w1+w2)
         cz = max(alpha, dims(1) - alpha)
      elseif (nint(sect_type) .eq. 5) then
         cy = dims(2)/2.d0
         cz = dims(1)/2.d0
      elseif (nint(sect_type) .eq. 6) then
         ! L-section (dims: b, h, t)
         w1 = dims(3)*dims(2)
         w2 = (dims(1)-dims(3))*dims(3)
         alpha = (w1*(dims(3)/2.d0) + 
     &            w2*(dims(3)+(dims(1)-dims(3))/2.d0))/(w1+w2)
         beta = (w1*(dims(2)/2.d0) + w2*(dims(3)/2.d0))/(w1+w2)
         cy = max(alpha, dims(1) - alpha)
         cz = max(beta, dims(2) - beta)
      elseif (nint(sect_type) .eq. 7) then
         cy = dims(2)/2.d0
         cz = dims(1)/2.d0
      else
         cy = 0.05d0
         cz = 0.05d0
      endif
!
!     compute stress values
!
      sigma_axial = axial_f / a
      sigma_max = dabs(sigma_axial) + dabs(my)*cz/xi22 + 
     &            dabs(mz)*cy/xi11
      tau_xy = vy / (a * xk_y)
      tau_xz = vz / (a * xk_z)
      tau_tor = tx * max(cy,cz) / j_tors
!
!     write continuous diagram along element length to beam_stresses.txt
!
      w1        = prop(index+24)
      w2        = prop(index+25)
      alpha     = prop(index+26)
      beta      = prop(index+27)
      load_type = nint(prop(index+28))
      load_dir  = nint(prop(index+29))
      
      ! Calculate nodal end moments
      M1_y = f_local(5)
      M2_y = -f_local(11)
      M1_z = -f_local(6)
      M2_z = f_local(12)
      
!     allocate ub31_stx and ub31_for lazily on first call
      if (.not. allocated(ub31_stx)) then
         allocate(ub31_stx(6, 11, ne0))
         ub31_stx = 0.d0
      endif
      if (.not. allocated(ub31_for)) then
         allocate(ub31_for(6, 11, ne0))
         ub31_for = 0.d0
      endif

!     1. Populate standard 11 stations for .frd and CCX internal use
      do station_i = 0, 10
         xi_station = dble(station_i) / 10.d0
         x_loc = xi_station * dl
         
         ! Free bending moment M0(x)
         M0_y = 0.d0
         M0_z = 0.d0
         
         if (load_dir .eq. 1) then
            if (load_type .eq. 1) then
               M0_z = 0.5d0 * w1 * x_loc * (dl - x_loc)
            elseif (load_type .eq. 2) then
               if (w1 .eq. 0.d0) then
                  M0_z = (w2 * (dl**2) / 6.d0) * 
     &                   (xi_station - xi_station**3)
               else
                  M0_z = (w1 * (dl**2) / 6.d0) * 
     &                   (2.d0*xi_station - 3.d0*xi_station**2 + 
     &                    xi_station**3)
               endif
            elseif (load_type .eq. 3) then
               R1 = w1 * dl * (beta - alpha) * 
     &              (1.d0 - 0.5d0 * (alpha + beta))
               R2 = w1 * dl * (beta - alpha) * 
     &              0.5d0 * (alpha + beta)
               if (x_loc .lt. alpha * dl) then
                  M0_z = R1 * x_loc
               elseif (x_loc .gt. beta * dl) then
                  M0_z = R2 * (dl - x_loc)
               else
                  M0_z = R1 * x_loc - 0.5d0 * w1 * 
     &                   ((x_loc - alpha * dl)**2)
               endif
            endif
         elseif (load_dir .eq. 2) then
            if (load_type .eq. 1) then
               M0_y = -0.5d0 * w1 * x_loc * (dl - x_loc)
            elseif (load_type .eq. 2) then
               if (w1 .eq. 0.d0) then
                  M0_y = -(w2 * (dl**2) / 6.d0) * 
     &                   (xi_station - xi_station**3)
               else
                  M0_y = -(w1 * (dl**2) / 6.d0) * 
     &                   (2.d0*xi_station - 3.d0*xi_station**2 + 
     &                    xi_station**3)
               endif
            elseif (load_type .eq. 3) then
               R1 = w1 * dl * (beta - alpha) * 
     &              (1.d0 - 0.5d0 * (alpha + beta))
               R2 = w1 * dl * (beta - alpha) * 
     &              0.5d0 * (alpha + beta)
               if (x_loc .lt. alpha * dl) then
                  M0_y = -R1 * x_loc
               elseif (x_loc .gt. beta * dl) then
                  M0_y = -R2 * (dl - x_loc)
               else
                  M0_y = -R1 * x_loc + 0.5d0 * w1 * 
     &                   ((x_loc - alpha * dl)**2)
               endif
            endif
         endif
         
         ! Total bending moments at station x
         my_x = M1_y * (1.d0 - xi_station) + M2_y * xi_station + M0_y
         mz_x = M1_z * (1.d0 - xi_station) + M2_z * xi_station + M0_z
         
         ! Station shear forces
         vy_x = f_local(8)
         vz_x = f_local(9)
         if (load_dir .eq. 1) then
            if (load_type .eq. 1) then
               vy_x = -f_local(2) - w1 * x_loc
            endif
         elseif (load_dir .eq. 2) then
            if (load_type .eq. 1) then
               vz_x = -f_local(3) - w1 * x_loc
            endif
         endif
         
         ! Individual stress components at station x
         sigma_b_y = dabs(my_x) * cz / xi22
         sigma_b_z = dabs(mz_x) * cy / xi11
         sigma_x   = dabs(sigma_axial) + sigma_b_y + sigma_b_z

         ! Store in dedicated station arrays (j=1..11 for stations 0..10)
         ub31_stx(1, station_i+1, nelem) = sigma_x
         ub31_stx(2, station_i+1, nelem) = tau_xy
         ub31_stx(3, station_i+1, nelem) = tau_xz
         ub31_stx(4, station_i+1, nelem) = tau_tor
         ub31_stx(5, station_i+1, nelem) = sigma_b_y
         ub31_stx(6, station_i+1, nelem) = sigma_b_z

         ub31_for(1, station_i+1, nelem) = -f_local(1)
         ub31_for(2, station_i+1, nelem) = vy_x
         ub31_for(3, station_i+1, nelem) = vz_x
         ub31_for(4, station_i+1, nelem) = tx
         ub31_for(5, station_i+1, nelem) = my_x
         ub31_for(6, station_i+1, nelem) = mz_x
      enddo

!     2. Dynamic Multi-Station CSV Results Streaming (*USER BEAM OUTPUT)
      if (.not. out_active) goto 880

!     Check increment filter
      call check_ub31_inc_active(istep, iinc, iout, out_inc_mode, 
     &     out_inc_freq, out_inc_list, out_ninc_list, is_active_inc)
      if (.not. is_active_inc) goto 880

!     Check if element belongs to any active target
      has_target = .false.
      if (allocated(out_elem_active)) then
         if (nelem .le. size(out_elem_active, 1)) then
            do itarg = 1, out_num_targets
               if (out_elem_active(nelem, itarg)) then
                  has_target = .true.
                  exit
               endif
            enddo
         endif
      else
         has_target = .true.
      endif
      if (.not. has_target) goto 880

!     Open CSV unit(s) if not yet open
      do itarg = 1, out_num_targets
         if (allocated(out_elem_active)) then
            if (nelem .gt. size(out_elem_active, 1)) cycle
            if (.not. out_elem_active(nelem, itarg)) cycle
         endif
         iunit = 90 + itarg
         inquire(unit=iunit, opened=is_open)
         if (.not. is_open) then
            if (istep .eq. 1 .and. iinc .le. 1 .and. 
     &          .not. out_target_file_init(itarg)) then
               open(iunit, file=trim(out_target_filename(itarg)), 
     &              status="unknown", position="rewind")
               call write_ub31_csv_header(iunit, out_flag_f, 
     &              out_flag_u, out_flag_s, out_flag_pts, out_flag_q)
               out_target_file_init(itarg) = .true.
            else
               open(iunit, file=trim(out_target_filename(itarg)), 
     &              status="unknown", position="append")
            endif
         endif
      enddo

!     Evaluate along span with N subdivisions (N+1 stations)
      num_sub = max(1, out_subdivisions)
      do station_k = 0, num_sub
         xi_station = dble(station_k) / dble(num_sub)
         x_loc = xi_station * dl
         station_pct = nint(xi_station * 100.d0)

         ! Free bending moment M0(x)
         M0_y = 0.d0
         M0_z = 0.d0
         if (load_dir .eq. 1) then
            if (load_type .eq. 1) then
               M0_z = 0.5d0 * w1 * x_loc * (dl - x_loc)
            elseif (load_type .eq. 2) then
               if (w1 .eq. 0.d0) then
                  M0_z = (w2 * (dl**2) / 6.d0) * 
     &                   (xi_station - xi_station**3)
               else
                  M0_z = (w1 * (dl**2) / 6.d0) * 
     &                   (2.d0*xi_station - 3.d0*xi_station**2 + 
     &                    xi_station**3)
               endif
            elseif (load_type .eq. 3) then
               R1 = w1 * dl * (beta - alpha) * 
     &              (1.d0 - 0.5d0 * (alpha + beta))
               R2 = w1 * dl * (beta - alpha) * 
     &              0.5d0 * (alpha + beta)
               if (x_loc .lt. alpha * dl) then
                  M0_z = R1 * x_loc
               elseif (x_loc .gt. beta * dl) then
                  M0_z = R2 * (dl - x_loc)
               else
                  M0_z = R1 * x_loc - 0.5d0 * w1 * 
     &                   ((x_loc - alpha * dl)**2)
               endif
            endif
         elseif (load_dir .eq. 2) then
            if (load_type .eq. 1) then
               M0_y = -0.5d0 * w1 * x_loc * (dl - x_loc)
            elseif (load_type .eq. 2) then
               if (w1 .eq. 0.d0) then
                  M0_y = -(w2 * (dl**2) / 6.d0) * 
     &                   (xi_station - xi_station**3)
               else
                  M0_y = -(w1 * (dl**2) / 6.d0) * 
     &                   (2.d0*xi_station - 3.d0*xi_station**2 + 
     &                    xi_station**3)
               endif
            elseif (load_type .eq. 3) then
               R1 = w1 * dl * (beta - alpha) * 
     &              (1.d0 - 0.5d0 * (alpha + beta))
               R2 = w1 * dl * (beta - alpha) * 
     &              0.5d0 * (alpha + beta)
               if (x_loc .lt. alpha * dl) then
                  M0_y = -R1 * x_loc
               elseif (x_loc .gt. beta * dl) then
                  M0_y = -R2 * (dl - x_loc)
               else
                  M0_y = -R1 * x_loc + 0.5d0 * w1 * 
     &                   ((x_loc - alpha * dl)**2)
               endif
            endif
         endif
         
         my_x = M1_y * (1.d0 - xi_station) + M2_y * xi_station + M0_y
         mz_x = M1_z * (1.d0 - xi_station) + M2_z * xi_station + M0_z
         
         vy_x = f_local(8)
         vz_x = f_local(9)
         if (load_dir .eq. 1) then
            if (load_type .eq. 1) then
               vy_x = -f_local(2) - w1 * x_loc
            endif
         elseif (load_dir .eq. 2) then
            if (load_type .eq. 1) then
               vz_x = -f_local(3) - w1 * x_loc
            endif
         endif
         
         sigma_b_y = dabs(my_x) * cz / xi22
         sigma_b_z = dabs(mz_x) * cy / xi11
         sigma_x   = dabs(sigma_axial) + sigma_b_y + sigma_b_z

         ! Displacements along span via Hermite shape functions + sag
         ux_val = (1.d0 - xi_station) * u_loc_beam(1) + 
     &            xi_station * u_loc_beam(7)
         rx_val = (1.d0 - xi_station) * u_loc_beam(4) + 
     &            xi_station * u_loc_beam(10)

         h1 = 1.d0 - 3.d0*(xi_station**2) + 2.d0*(xi_station**3)
         h2 = xi_station - 2.d0*(xi_station**2) + (xi_station**3)
         h3 = 3.d0*(xi_station**2) - 2.d0*(xi_station**3)
         h4 = -(xi_station**2) + (xi_station**3)

         dh1 = -6.d0*xi_station + 6.d0*(xi_station**2)
         dh2 = 1.d0 - 4.d0*xi_station + 3.d0*(xi_station**2)
         dh3 = 6.d0*xi_station - 6.d0*(xi_station**2)
         dh4 = -2.d0*xi_station + 3.d0*(xi_station**2)

         vh = h1*u_loc_beam(2) + dl*h2*u_loc_beam(6) + 
     &        h3*u_loc_beam(8) + dl*h4*u_loc_beam(12)
         thetazh = (dh1*u_loc_beam(2) + dh3*u_loc_beam(8))/dl + 
     &             dh2*u_loc_beam(6) + dh4*u_loc_beam(12)

         wh = h1*u_loc_beam(3) - dl*h2*u_loc_beam(5) + 
     &        h3*u_loc_beam(9) - dl*h4*u_loc_beam(11)
         thetayh = -((dh1*u_loc_beam(3) + dh3*u_loc_beam(9))/dl - 
     &               dh2*u_loc_beam(5) - dh4*u_loc_beam(11))

         v0 = 0.d0
         w0 = 0.d0
         thetaz0 = 0.d0
         thetay0 = 0.d0
         if (load_dir .eq. 1 .and. load_type .eq. 1) then
            v0 = (w1 * x_loc * (dl - x_loc) * 
     &           (dl**2 + x_loc*(dl - x_loc))) / (24.d0 * e * xi11)
            thetaz0 = (w1 * (dl**3 - 6.d0*dl*(x_loc**2) + 
     &                4.d0*(x_loc**3))) / (24.d0 * e * xi11)
         elseif (load_dir .eq. 2 .and. load_type .eq. 1) then
            w0 = (w1 * x_loc * (dl - x_loc) * 
     &           (dl**2 + x_loc*(dl - x_loc))) / (24.d0 * e * xi22)
            thetay0 = -(w1 * (dl**3 - 6.d0*dl*(x_loc**2) + 
     &                4.d0*(x_loc**3))) / (24.d0 * e * xi22)
         endif

         uy_val = vh + v0
         uz_val = wh + w0
         ry_val = thetayh + thetay0
         rz_val = thetazh + thetaz0

         ! Active distributed loads at station
         qx_val = 0.d0
         qy_val = 0.d0
         qz_val = 0.d0
         if (load_dir .eq. 1) then
            if (load_type .eq. 1) then
               qy_val = w1
            elseif (load_type .eq. 2) then
               qy_val = w1*(1.d0 - xi_station) + w2*xi_station
            elseif (load_type .eq. 3) then
               if (x_loc .ge. alpha*dl .and. x_loc .le. beta*dl) then
                  qy_val = w1
               endif
            endif
         elseif (load_dir .eq. 2) then
            if (load_type .eq. 1) then
               qz_val = w1
            elseif (load_type .eq. 2) then
               qz_val = w1*(1.d0 - xi_station) + w2*xi_station
            elseif (load_type .eq. 3) then
               if (x_loc .ge. alpha*dl .and. x_loc .le. beta*dl) then
                  qz_val = w1
               endif
            endif
         endif

         ! Coordinate system transformation
         if (out_coords .eq. 2) then
            ! Global coordinates: T^T * v_local
            do m_idx = 1, 3
               f_csv(m_idx) = axial_f*e1(m_idx) + vy_x*e2(m_idx) + 
     &                        vz_x*e3(m_idx)
               m_csv(m_idx) = f_local(10)*e1(m_idx) + my_x*e2(m_idx) + 
     &                        mz_x*e3(m_idx)
               u_csv(m_idx) = ux_val*e1(m_idx) + uy_val*e2(m_idx) + 
     &                        uz_val*e3(m_idx)
               rot_csv(m_idx) = rx_val*e1(m_idx) + ry_val*e2(m_idx) + 
     &                          rz_val*e3(m_idx)
               q_csv(m_idx) = qx_val*e1(m_idx) + qy_val*e2(m_idx) + 
     &                        qz_val*e3(m_idx)
            enddo
         else
            ! Local beam coordinates
            f_csv(1) = axial_f
            f_csv(2) = vy_x
            f_csv(3) = vz_x
            m_csv(1) = f_local(10)
            m_csv(2) = my_x
            m_csv(3) = mz_x
            u_csv(1) = ux_val
            u_csv(2) = uy_val
            u_csv(3) = uz_val
            rot_csv(1) = rx_val
            rot_csv(2) = ry_val
            rot_csv(3) = rz_val
            q_csv(1) = qx_val
            q_csv(2) = qy_val
            q_csv(3) = qz_val
         endif

         ! Section Points Stress Recovery (Signed Normal Stresses at Corners/Fibers)
         if (nint(sect_type) .eq. 1) then
            ! RECT (dims: b, h) -> y in [-h/2, +h/2], z in [-b/2, +b/2]
            y_pts(1) = dims(2)/2.d0
            z_pts(1) = -dims(1)/2.d0
            pt_labels(1) = '1_TL'
            y_pts(2) = dims(2)/2.d0
            z_pts(2) = dims(1)/2.d0
            pt_labels(2) = '2_TR'
            y_pts(3) = -dims(2)/2.d0
            z_pts(3) = -dims(1)/2.d0
            pt_labels(3) = '3_BL'
            y_pts(4) = -dims(2)/2.d0
            z_pts(4) = dims(1)/2.d0
            pt_labels(4) = '4_BR'
         elseif (nint(sect_type) .eq. 2) then
            ! CIRC / PIPE (dims: r_o, t)
            y_pts(1) = dims(1)
            z_pts(1) = 0.d0
            pt_labels(1) = '1_Top'
            y_pts(2) = 0.d0
            z_pts(2) = dims(1)
            pt_labels(2) = '2_Right'
            y_pts(3) = -dims(1)
            z_pts(3) = 0.d0
            pt_labels(3) = '3_Bottom'
            y_pts(4) = 0.d0
            z_pts(4) = -dims(1)
            pt_labels(4) = '4_Left'
         elseif (nint(sect_type) .eq. 3) then
            ! I-BEAM (dims: h, b_top, t_f_top, b_bot, t_f_bot, t_w)
            y_pts(1) = dims(1)/2.d0
            z_pts(1) = -dims(2)/2.d0
            pt_labels(1) = '1_TL'
            y_pts(2) = dims(1)/2.d0
            z_pts(2) = dims(2)/2.d0
            pt_labels(2) = '2_TR'
            y_pts(3) = -dims(1)/2.d0
            z_pts(3) = -dims(4)/2.d0
            pt_labels(3) = '3_BL'
            y_pts(4) = -dims(1)/2.d0
            z_pts(4) = dims(4)/2.d0
            pt_labels(4) = '4_BR'
         elseif (nint(sect_type) .eq. 4) then
            ! T-SECTION (dims: h, b, t_f, t_w)
            w1 = dims(2)*dims(3)
            w2 = dims(4)*(dims(1)-dims(3))
            alpha_t = (w1*(dims(1)-dims(3)/2.d0) + 
     &                 w2*(dims(1)-dims(3))/2.d0)/(w1+w2)
            y_pts(1) = dims(1) - alpha_t
            z_pts(1) = -dims(2)/2.d0
            pt_labels(1) = '1_TopL'
            y_pts(2) = dims(1) - alpha_t
            z_pts(2) = dims(2)/2.d0
            pt_labels(2) = '2_TopR'
            y_pts(3) = -alpha_t
            z_pts(3) = 0.d0
            pt_labels(3) = '3_StemBot'
            y_pts(4) = dims(1) - alpha_t - dims(3)
            z_pts(4) = dims(4)/2.d0
            pt_labels(4) = '4_Junction'
         elseif (nint(sect_type) .eq. 5) then
            ! U-CHANNEL (dims: h, b, t_f, t_w)
            w1 = dims(4)*dims(1)
            w2 = 2.d0*(dims(2)-dims(4))*dims(3)
            beta_c = (w1*(dims(4)/2.d0) + 
     &                w2*(dims(4)+(dims(2)-dims(4))/2.d0))/(w1+w2)
            y_pts(1) = dims(1)/2.d0
            z_pts(1) = dims(2) - beta_c
            pt_labels(1) = '1_TopFlg'
            y_pts(2) = dims(1)/2.d0
            z_pts(2) = -beta_c
            pt_labels(2) = '2_TopWeb'
            y_pts(3) = -dims(1)/2.d0
            z_pts(3) = -beta_c
            pt_labels(3) = '3_BotWeb'
            y_pts(4) = -dims(1)/2.d0
            z_pts(4) = dims(2) - beta_c
            pt_labels(4) = '4_BotFlg'
         elseif (nint(sect_type) .eq. 6) then
            ! L-ANGLE (dims: b, h, t)
            w1 = dims(3)*dims(2)
            w2 = (dims(1)-dims(3))*dims(3)
            alpha_l = (w1*(dims(3)/2.d0) + 
     &                 w2*(dims(3)+(dims(1)-dims(3))/2.d0))/(w1+w2)
            beta_l = (w1*(dims(2)/2.d0) + w2*(dims(3)/2.d0))/(w1+w2)
            y_pts(1) = dims(2) - alpha_l
            z_pts(1) = -beta_l
            pt_labels(1) = '1_VertTop'
            y_pts(2) = -alpha_l
            z_pts(2) = -beta_l
            pt_labels(2) = '2_Heel'
            y_pts(3) = -alpha_l
            z_pts(3) = dims(1) - beta_l
            pt_labels(3) = '3_HorizTip'
            y_pts(4) = -alpha_l + dims(3)
            z_pts(4) = -beta_l + dims(3)
            pt_labels(4) = '4_Inner'
         elseif (nint(sect_type) .eq. 7) then
            ! BOX (dims: h, b, t1..t4)
            y_pts(1) = dims(1)/2.d0
            z_pts(1) = -dims(2)/2.d0
            pt_labels(1) = '1_TL'
            y_pts(2) = dims(1)/2.d0
            z_pts(2) = dims(2)/2.d0
            pt_labels(2) = '2_TR'
            y_pts(3) = -dims(1)/2.d0
            z_pts(3) = -dims(2)/2.d0
            pt_labels(3) = '3_BL'
            y_pts(4) = -dims(1)/2.d0
            z_pts(4) = dims(2)/2.d0
            pt_labels(4) = '4_BR'
         else
            y_pts(1) = cy
            z_pts(1) = -cz
            pt_labels(1) = '1_TL'
            y_pts(2) = cy
            z_pts(2) = cz
            pt_labels(2) = '2_TR'
            y_pts(3) = -cy
            z_pts(3) = -cz
            pt_labels(3) = '3_BL'
            y_pts(4) = -cy
            z_pts(4) = cz
            pt_labels(4) = '4_BR'
         endif

         tau_vy_val = vy_x / (a * max(xk_y, 0.01d0))
         tau_vz_val = vz_x / (a * max(xk_z, 0.01d0))

         if (out_flag_pts) then
            do ipt = 1, 4
               s_norm_val = (axial_f / a) - (mz_x * y_pts(ipt) / xi11) + 
     &                      (my_x * z_pts(ipt) / xi22)
               r_dist = dsqrt(y_pts(ipt)**2 + z_pts(ipt)**2)
               tau_tor_pt = dabs(f_local(10)) * r_dist / j_tors
               s_vm_pt = dsqrt(s_norm_val**2 + 3.d0*(tau_vy_val**2 + 
     &                   tau_vz_val**2 + tau_tor_pt**2))

               do itarg = 1, out_num_targets
                  if (allocated(out_elem_active)) then
                     if (nelem .gt. size(out_elem_active, 1)) cycle
                     if (.not. out_elem_active(nelem, itarg)) cycle
                  endif
                  iunit = 90 + itarg
                  call write_ub31_point_row(iunit, istep, iinc, time,
     &                 nelem, station_pct, x_loc, pt_labels(ipt),
     &                 y_pts(ipt), z_pts(ipt), out_flag_f, f_csv,
     &                 m_csv, out_flag_u, u_csv, rot_csv, 
     &                 (out_flag_s .or. out_flag_pts), sigma_axial,
     &                 s_norm_val, tau_vy_val, tau_vz_val,
     &                 tau_tor_pt, s_vm_pt, out_flag_q, q_csv)
               enddo
            enddo
         else
            ! Write single row per station
            do itarg = 1, out_num_targets
               if (allocated(out_elem_active)) then
                  if (nelem .gt. size(out_elem_active, 1)) cycle
                  if (.not. out_elem_active(nelem, itarg)) cycle
               endif
               iunit = 90 + itarg
               call write_ub31_station_row(iunit, istep, iinc, time,
     &              nelem, station_pct, x_loc, out_flag_f, f_csv,
     &              m_csv, out_flag_u, u_csv, rot_csv, out_flag_s, 
     &              sigma_axial, sigma_b_y, sigma_b_z, sigma_x, 
     &              tau_xy, tau_xz, tau_tor, out_flag_q, q_csv)
            enddo
         endif
      enddo

      do itarg = 1, out_num_targets
         if (allocated(out_elem_active)) then
            if (nelem .gt. size(out_elem_active, 1)) cycle
            if (.not. out_elem_active(nelem, itarg)) cycle
         endif
         flush(90 + itarg)
      enddo

 880  continue
!
!     populate stx array for .frd output:
!     stx(1) = Combined Max Normal Stress sigma_max (Axial + Bending)
!     stx(2) = Transverse Shear Stress tau_xy
!     stx(3) = Transverse Shear Stress tau_xz
!     stx(4) = Torsional Shear Stress tau_tor
!     stx(5) = Bending Moment My
!     stx(6) = Bending Moment Mz
!
!     also fill stx(1..mi(1)) with end-point values for CCX internal use
      do j=1,mi(1)
         ! For CCX internal stx: station 1 -> xi=0.0, last -> xi=1.0
         if (j .le. 11) then
            station_i = j - 1
            xi_station = dble(station_i) / 10.d0
            x_loc = xi_station * dl
            
            M0_y = 0.d0
            M0_z = 0.d0
            if (load_dir .eq. 1) then
               if (load_type .eq. 1) then
                  M0_z = 0.5d0 * w1 * x_loc * (dl - x_loc)
               elseif (load_type .eq. 2) then
                  if (w1 .eq. 0.d0) then
                     M0_z = (w2 * (dl**2) / 6.d0) * 
     &                      (xi_station - xi_station**3)
                  else
                     M0_z = (w1 * (dl**2) / 6.d0) * 
     &                      (2.d0*xi_station - 3.d0*xi_station**2 + 
     &                       xi_station**3)
                  endif
               elseif (load_type .eq. 3) then
                  R1 = w1 * dl * (beta - alpha) * 
     &                 (1.d0 - 0.5d0 * (alpha + beta))
                  R2 = w1 * dl * (beta - alpha) * 
     &                 0.5d0 * (alpha + beta)
                  if (x_loc .lt. alpha * dl) then
                     M0_z = R1 * x_loc
                  elseif (x_loc .gt. beta * dl) then
                     M0_z = R2 * (dl - x_loc)
                  else
                     M0_z = R1 * x_loc - 0.5d0 * w1 * 
     &                      ((x_loc - alpha * dl)**2)
                  endif
               endif
            elseif (load_dir .eq. 2) then
               if (load_type .eq. 1) then
                  M0_y = -0.5d0 * w1 * x_loc * (dl - x_loc)
               elseif (load_type .eq. 2) then
                  if (w1 .eq. 0.d0) then
                     M0_y = -(w2 * (dl**2) / 6.d0) * 
     &                      (xi_station - xi_station**3)
                  else
                     M0_y = -(w1 * (dl**2) / 6.d0) * 
     &                      (2.d0*xi_station - 3.d0*xi_station**2 + 
     &                       xi_station**3)
                  endif
               elseif (load_type .eq. 3) then
                  R1 = w1 * dl * (beta - alpha) * 
     &                 (1.d0 - 0.5d0 * (alpha + beta))
                  R2 = w1 * dl * (beta - alpha) * 
     &                 0.5d0 * (alpha + beta)
                  if (x_loc .lt. alpha * dl) then
                     M0_y = -R1 * x_loc
                  elseif (x_loc .gt. beta * dl) then
                     M0_y = -R2 * (dl - x_loc)
                  else
                     M0_y = -R1 * x_loc + 0.5d0 * w1 * 
     &                      ((x_loc - alpha * dl)**2)
                  endif
               endif
            endif
            my_x = M1_y * (1.d0 - xi_station) + M2_y * xi_station + M0_y
            mz_x = M1_z * (1.d0 - xi_station) + M2_z * xi_station + M0_z
            sigma_x = dabs(sigma_axial) + dabs(my_x)*cz/xi22 + 
     &                dabs(mz_x)*cy/xi11

            stx(1,j,i) = sigma_axial
            stx(2,j,i) = tau_xy
            stx(3,j,i) = tau_xz
            stx(4,j,i) = tau_tor
            stx(5,j,i) = my_x
            stx(6,j,i) = mz_x
         else
            stx(1,j,i) = sigma_axial
            stx(2,j,i) = tau_xy
            stx(3,j,i) = tau_xz
            stx(4,j,i) = tau_tor
            stx(5,j,i) = my
            stx(6,j,i) = mz
         endif
         
         ! strain array (restored elastic strains)
         eme(1,j,i) = axial_f / (e*a)
         eme(2,j,i) = vy / (um*a*xk_y)
         eme(3,j,i) = vz / (um*a*xk_z)
         eme(4,j,i) = tx / (um*j_tors)
         eme(5,j,i) = my / (e*xi22)
         eme(6,j,i) = mz / (e*xi11)
         
         if (loc(eei) .ne. 0) then
            do k=1,6
               eei(k,j,i) = eme(k,j,i)
            enddo
         endif
      enddo
!
      if (calcul_qa .eq. 1) then
         do inode = 1, 2
            do j = 1, 3
               qa(1) = qa(1) + dabs(fn(j, konl(inode)))
            enddo
         enddo
         nal = nal + 6
      endif
!
      return
      end
!
! !
! ======================================================================
!     Helper routine to write dynamic CSV header
! ======================================================================
      subroutine write_ub31_csv_header(iunit, flag_f, flag_u,
     &     flag_s, flag_pts, flag_q)
      implicit none
      integer iunit
      logical flag_f, flag_u, flag_s, flag_pts, flag_q
      character*800 h_line
!
      if (flag_pts) then
         h_line = 'Step,Increment,Time,Element,Station_Pct,X_local,'//
     &            'Section_Pt,Y_local,Z_local'
         if (flag_f) then
            h_line = trim(h_line) //
     &           ',Fx_Axial,Vy_Shear,Vz_Shear,Mx_Torsion,'//
     &           'My_Bending,Mz_Bending'
         endif
         if (flag_u) then
            h_line = trim(h_line) //
     &           ',Ux,Uy,Uz,Rot_X,Rot_Y,Rot_Z'
         endif
         if (flag_s .or. flag_pts) then
            h_line = trim(h_line) //
     &           ',Sxx_Axial,Sxx_Normal,Sxy_Shear,Sxz_Shear,'//
     &           'Stors_Torsion,S_vonMises'
         endif
         if (flag_q) then
            h_line = trim(h_line) //
     &           ',Qx_Load,Qy_Load,Qz_Load'
         endif
      else
         h_line = 'Step,Increment,Time,Element,Station_Pct,X_local'
         if (flag_f) then
            h_line = trim(h_line) //
     &           ',Fx_Axial,Vy_Shear,Vz_Shear,Mx_Torsion,'//
     &           'My_Bending,Mz_Bending'
         endif
         if (flag_u) then
            h_line = trim(h_line) //
     &           ',Ux,Uy,Uz,Rot_X,Rot_Y,Rot_Z'
         endif
         if (flag_s) then
            h_line = trim(h_line) //
     &           ',Sxx_Axial,Sxx_Bending_Y,Sxx_Bending_Z,'//
     &           'Sxx_Max_Combined,Sxy_Shear,Sxz_Shear,Stors_Torsion'
         endif
         if (flag_q) then
            h_line = trim(h_line) //
     &           ',Qx_Load,Qy_Load,Qz_Load'
         endif
      endif
      write(iunit, '(A)') trim(h_line)
      end subroutine write_ub31_csv_header
!
! ======================================================================
!     Helper routine to write single station row
! ======================================================================
      subroutine write_ub31_station_row(iunit, istep, iinc, time_val,
     &     nelem, st_pct, x_loc, flag_f, f_v, m_v, flag_u, u_v,
     &     r_v, flag_s, s_ax, s_by, s_bz, s_max, t_xy, t_xz,
     &     t_tor, flag_q, q_v)
      implicit none
      integer iunit, istep, iinc, nelem, st_pct
      logical flag_f, flag_u, flag_s, flag_q
      real*8 time_val, x_loc, f_v(3), m_v(3), u_v(3), r_v(3),
     &  s_ax, s_by, s_bz, s_max, t_xy, t_xz, t_tor, q_v(3)
      character*1200 l_buf
      character*300 n_buf
!
      write(l_buf, '(I4,A,I6,A,F14.6,A,I8,A,I4,A,F10.4)')
     &     istep, ',', iinc, ',', time_val, ',', nelem, ',',
     &     st_pct, ',', x_loc
!
      if (flag_f) then
         write(n_buf, 1005)
     &        ',', f_v(1), ',', f_v(2), ',', f_v(3),
     &        ',', m_v(1), ',', m_v(2), ',', m_v(3)
 1005    format(A,F16.4,A,F16.4,A,F16.4,A,F16.4,A,F16.4,A,F16.4)
         l_buf = trim(l_buf) // trim(n_buf)
      endif
!
      if (flag_u) then
         write(n_buf, 1006)
     &        ',', u_v(1), ',', u_v(2), ',', u_v(3),
     &        ',', r_v(1), ',', r_v(2), ',', r_v(3)
 1006    format(A,F16.6,A,F16.6,A,F16.6,A,F16.6,A,F16.6,A,F16.6)
         l_buf = trim(l_buf) // trim(n_buf)
      endif
!
      if (flag_s) then
         write(n_buf, 1010)
     &        ',', s_ax, ',', s_by, ',', s_bz, ',', s_max,
     &        ',', t_xy, ',', t_xz, ',', t_tor
 1010    format(A,F16.4,A,F16.4,A,F16.4,A,F16.4,A,F16.4,A,F16.4,A,F16.4)
         l_buf = trim(l_buf) // trim(n_buf)
      endif
!
      if (flag_q) then
         write(n_buf, '(A,F16.4,A,F16.4,A,F16.4)')
     &        ',', q_v(1), ',', q_v(2), ',', q_v(3)
         l_buf = trim(l_buf) // trim(n_buf)
      endif
!
      write(iunit, '(A)') trim(l_buf)
      end subroutine write_ub31_station_row
!
! ======================================================================
!     Helper routine to write single section-point row
! ======================================================================
      subroutine write_ub31_point_row(iunit, istep, iinc, time_val,
     &     nelem, st_pct, x_loc, pt_name, y_pt, z_pt, flag_f, f_v,
     &     m_v, flag_u, u_v, r_v, flag_s, s_ax, s_norm, t_vy, t_vz,
     &     t_tor, s_vm, flag_q, q_v)
      implicit none
      integer iunit, istep, iinc, nelem, st_pct
      character*(*) pt_name
      real*8 y_pt, z_pt
      logical flag_f, flag_u, flag_s, flag_q
      real*8 time_val, x_loc, f_v(3), m_v(3), u_v(3), r_v(3),
     &  s_ax, s_norm, t_vy, t_vz, t_tor, s_vm, q_v(3)
      character*1200 l_buf
      character*300 n_buf
!
      write(l_buf, 1004)
     &     istep, ',', iinc, ',', time_val, ',', nelem, ',',
     &     st_pct, ',', x_loc, ',', trim(pt_name), ',',
     &     y_pt, ',', z_pt
 1004 format(I4,A,I6,A,F14.6,A,I8,A,I4,A,F10.4,A,A,A,F12.4,A,F12.4)
!
      if (flag_f) then
         write(n_buf, 1005)
     &        ',', f_v(1), ',', f_v(2), ',', f_v(3),
     &        ',', m_v(1), ',', m_v(2), ',', m_v(3)
 1005    format(A,F16.4,A,F16.4,A,F16.4,A,F16.4,A,F16.4,A,F16.4)
         l_buf = trim(l_buf) // trim(n_buf)
      endif
!
      if (flag_u) then
         write(n_buf, 1006)
     &        ',', u_v(1), ',', u_v(2), ',', u_v(3),
     &        ',', r_v(1), ',', r_v(2), ',', r_v(3)
 1006    format(A,F16.6,A,F16.6,A,F16.6,A,F16.6,A,F16.6,A,F16.6)
         l_buf = trim(l_buf) // trim(n_buf)
      endif
!
      if (flag_s) then
         write(n_buf, 1010)
     &        ',', s_ax, ',', s_norm, ',', t_vy, ',', t_vz,
     &        ',', t_tor, ',', s_vm
 1010    format(A,F16.4,A,F16.4,A,F16.4,A,F16.4,A,F16.4,A,F16.4)
         l_buf = trim(l_buf) // trim(n_buf)
      endif
!
      if (flag_q) then
         write(n_buf, '(A,F16.4,A,F16.4,A,F16.4)')
     &        ',', q_v(1), ',', q_v(2), ',', q_v(3)
         l_buf = trim(l_buf) // trim(n_buf)
      endif
!
      write(iunit, '(A)') trim(l_buf)
      end subroutine write_ub31_point_row
!
! ======================================================================
!     Helper routine to check if current increment should produce output
! ======================================================================
      subroutine check_ub31_inc_active(istep, iinc, iout, inc_mode,
     &     inc_freq, inc_list, ninc_list, is_active)
      implicit none
      integer istep, iinc, iout, inc_freq, inc_list(50), ninc_list
      character*4 inc_mode
      logical is_active
      integer k
!
      is_active = .false.
!     Do not output during initial passes (iout <= 0, except -2)
      if (iout .le. 0 .and. iout .ne. -2) return
!
      if (inc_mode .eq. 'ALL ') then
         is_active = .true.
      elseif (inc_mode .eq. 'FREQ') then
         if (inc_freq .gt. 0) then
            if (iinc .gt. 0 .and. mod(iinc, inc_freq) .eq. 0) then
               is_active = .true.
            endif
         endif
         if (iout .eq. 2 .or. iout .eq. -2) is_active = .true.
      elseif (inc_mode .eq. 'LIST') then
         do k = 1, ninc_list
            if (iinc .eq. inc_list(k)) then
               is_active = .true.
               exit
            endif
         enddo
      else
!        Default: LAST
         is_active = (iout .eq. 2 .or. iout .eq. -2 .or. iout .eq. 1)
      endif
      end subroutine check_ub31_inc_active

