!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2025 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
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
      subroutine e_c3d_ub31(co,kon,lakonl,p1,p2,omx,bodyfx,nbody,s,sm,
     &  ff,nelem,nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &  ielmat,ielorien,norien,orab,ntmat_,
     &  t0,t1,ithermal,vold,iperturb,nelemload,
     &  sideload,xload,nload,idist,sti,stx,iexpl,plicon,
     &  nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &  matname,mi,ncmat_,mass,stiffness,buckling,rhsi,intscheme,
     &  ttime,time,istep,iinc,coriolis,xloadold,reltime,
     &  ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,veold,
     &  ne0,ipkon,thicke,
     &  integerglob,doubleglob,tieset,istartset,iendset,ialset,ntie,
     &  nasym,ielprop,prop)
!
!     computation of the element matrix and rhs for the user element
!     of type UB31 (2-node Euler-Bernoulli beam element)
!     Supports:
!     - Analytical section properties (Rectangular, Pipe, I-beam, Box)
!     - Local offsets at node 1 & node 2
!     - Section rotation
!     - Member end releases (hinges) at node 1 & node 2
!
      use ub31_module
      implicit none
!
      integer mass,stiffness,buckling,rhsi,coriolis
!
      character*8 lakonl
      character*20 sideload(*)
      character*80 matname(*),amat
      character*81 tieset(3,*)
!
      integer konl(26),nelemload(2,*),nbody,nelem,mi(*),kon(*),
     &  ielprop(*),index,mattyp,ithermal(*),iperturb(*),nload,idist,
     &  i,j,k,i1,nmethod,kk,nelcon(2,*),nrhcon(*),nalcon(2,*),
     &  ielmat(mi(3),*),ielorien(mi(3),*),ipkon(*),indexe,
     &  ntmat_,nope,norien,ihyper,iexpl,kode,imat,iorien,istiff,
     &  ncmat_,intscheme,istep,iinc,ipompc(*),nodempc(3,*),
     &  nmpc,ikmpc(*),ilmpc(*),ne0,ndof,istartset(*),iendset(*),
     &  ialset(*),ntie,integerglob(*),nasym,nplicon(0:ntmat_,*),
     &  nplkcon(0:ntmat_,*),npmat_
!
      real*8 co(3,*),xl(3,20),veold(0:mi(2),*),rho,s(60,60),bodyfx(3),
     &  ff(60),elconloc(ncmat_),coords(3),p1(3),c1,c2,c3,c4,
     &  p2(3),eth(6),rhcon(0:1,ntmat_,*),reltime,prop(*),tm(3,3),
     &  alcon(0:6,ntmat_,*),alzero(*),orab(7,*),t0(*),t1(*),
     &  xloadold(2,*),vold(0:mi(2),*),xload(2,*),omx,e,un,um,
     &  sm(60,60),sti(6,mi(1),*),stx(6,mi(1),*),t0l,t1l,coefmpc(*),
     &  stiff(21),thicke(mi(3),*),doubleglob(*),dl,e2(3),e3(3),
     &  plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &  xstiff(27,mi(1),*),plconloc(802),dtime,ttime,time,tmg(12,12),
     &  a,xi11,xi22,xk_y,xk_z,e1(3),phi_y,phi_z,c_rot_y,c_rot_z,
     &  c_mass,c_tors_mass,
     &  sg(12,12),elcon(0:ncmat_,ntmat_,*),smg(12,12),
     &  sect_type, dims(6), rot_angle, offsets(3,2),
     &  x_beam(3,2), e2_in(3), len_e2, j_tors, principal_angle,
     &  rot_angle_total, ff_loc(60), w_glob(3), w_loc(3), val, n_th,
     &  w1, w2, alpha, beta, f1, f2, m1, m2
      integer release_codes(2), code, inode, r,
     &  load_type, load_dir, a_pct, b_pct, node1, node2,
     &  lumped_flag, env_stat, node_id, iorient_error
      character(len=10) :: env_lump
      logical released(12)
      real*8 spring_k(12)
      real*8 px1, py1, pz1, px2, py2, pz2, c1_cnt, c2_cnt,
     &  ax1, ay1, az1, ax2, ay2, az2, q1_glob(3), q2_glob(3),
     &  q1_loc(3), q2_loc(3)
      real*8 xl_def(3,20),e1_init(3),e1_curr(3),len_init,len_curr,
     &  cross_e1(3),dot_e1,u_rot(3),cross_u(3),dot_u_e2,len_cross,
     &  e2_init(3),e3_init(3),n_axial,c1_g,c2_g,c3_g,c4_g,
     &  c_tors_g,eps_th,timo_val,off_vec(3),rot_disp(3)
!
!

      indexe=ipkon(nelem)
!
!     material and orientation
!
      imat=ielmat(1,nelem)
      amat=matname(imat)
      if(norien.gt.0) then
         iorien=max(0,ielorien(1,nelem))
      else
         iorien=0
      endif
!
      if(nelcon(1,imat).lt.0) then
         ihyper=1
      else
         ihyper=0
      endif
      rho=rhcon(1,imat,1)
!
!     properties of the cross section and offsets
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
         timo_val=prop(index+20)
!
!        Native direct-beam element names encode the formulation in
!        lakon character 5.  Plain UB31 retains the section-card setting.
         if(lakonl(5:5).eq.'T') then
!           Timoshenko: use automatically calculated section shear factors.
            timo_val=0.d0
         elseif(lakonl(5:5).eq.'E') then
!           Euler-Bernoulli: suppress transverse shear deformation.
            timo_val=-1.d0
         endif
      else
         write(*,*) '*ERROR in e_c3d_ub31: element ',nelem,
     &        ' has no section properties.'
         write(*,*) '  Use *USER SECTION or *BEAM SECTION'
     &        ,' to define cross-section.'
         call exit(201)
      endif
!
      nope=2
      ndof=6
!
!     computation of the coordinates of the local nodes
!
      do i=1,nope
         konl(i)=kon(indexe+i)
         do j=1,3
            xl(j,i)=co(j,konl(i))
         enddo
      enddo
!
!     initialisation for distributed forces
!
      if(rhsi.eq.1) then
        if(idist.ne.0) then
          do i=1,ndof*nope
            ff(i)=0.d0
          enddo
        endif
      endif
!
!     initialisation of sm and s
!
      if((mass.eq.1).or.(buckling.eq.1)) then
        do i=1,ndof*nope
          do j=1,ndof*nope
            sm(i,j)=0.d0
          enddo
        enddo
      endif
!
      do i=1,ndof*nope
        do j=1,ndof*nope
          s(i,j)=0.d0
        enddo
      enddo
!
!     material data
!
      t0l=0.d0
      t1l=0.d0
      if(ithermal(1).eq.1) then
         do i1=1,nope
            t0l=t0l+t0(konl(i1))/2.d0
            t1l=t1l+t1(konl(i1))/2.d0
         enddo
      endif
!
      kode=nelcon(1,imat)
      if(kode.eq.2) then
         mattyp=1
      else
         mattyp=0
      endif
!
      istiff=0
      call materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &     imat,amat,iorien,coords,orab,ntmat_,stiff,rho,
     &     nelem,ithermal,alzero,mattyp,t0l,t1l,
     &     ihyper,istiff,elconloc,eth,kode,plicon,
     &     nplicon,plkcon,nplkcon,npmat_,
     &     plconloc,mi(1),dtime,kk,
     &     xstiff,ncmat_,iperturb)
!
      if(elconloc(1).gt.0.d0) then
         e=elconloc(1)
         un=elconloc(2)
         um=e/(2.d0*(1.d0+un))
      else
         write(*,*) '*ERROR in e_c3d_ub31: invalid elastic properties'
         call exit(201)
      endif
!
!     analytical calculation of cross-section properties (including principal axes rotation)
!
!     NOTE on inertia naming convention:
!     compute_section_properties returns (A, Iyy, Izz, J, xk_y, xk_z)
!     where Iyy = moment of inertia about local y-axis (=> bending in x-z plane)
!           Izz = moment of inertia about local z-axis (=> bending in x-y plane)
!     These map here as: xi22 <- Iyy  (x-z bending, phi_y block)
!                        xi11 <- Izz  (x-y bending, phi_z block)
!     Shear factors:      xk_y acts on the x-y shear (v direction)
!                         xk_z acts on the x-z shear (w direction)
!
      call compute_section_properties(sect_type, dims, e, un,
     &  a, xi22, xi11, j_tors, xk_y, xk_z, principal_angle)
!
!     A single M1 or M2 release on a rotated principal system
!     is a linear combination of principal-axis rotational DOFs.
!     The current scalar condensation cannot represent it.
!     Reject this case rather than silently releasing a different
!     geometrical direction. Dual bending release and torsion
!     remain valid for asymmetric sections.
      if(dabs(principal_angle).gt.1.d-8) then
         do inode=1,2
            code=release_codes(inode)
            if (btest(code,4).neqv.btest(code,5)) then
               write(*,*) '*ERROR in B31T/B31E: single M1/M2'
               write(*,*) '  release on asymmetric section,'
               write(*,*) '  element ',nelem
               write(*,*) '  Use M1-M2 or no bending release.'
               call exit(201)
            endif
         enddo
      endif
!
!     tentative local axes (based on node line)
!
      do j=1,3
         e1(j)=xl(j,2)-xl(j,1)
      enddo
      dl=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
      do j=1,3
         e1(j)=e1(j)/dl
      enddo
!
!     Native CCX section normal; explicit input is projected.
      call ccx_beam_section_normal(e1,e2_in,e2,iorient_error)
      if(iorient_error.ne.0) then
         write(*,*) '*ERROR: beam section axis parallel'
         write(*,*) '  to element tangent (element ',nelem,')'
         call exit(201)
      endif
!
!     apply section rotation angle + principal axes rotation angle
!
      rot_angle_total = rot_angle + principal_angle
      if (dabs(rot_angle_total) .gt. 1.d-6) then
         call rotate_vector_local(e2, e1, rot_angle_total, e2)
      endif
!
!     calculate e3 = e1 x e2
!
      e3(1)=e1(2)*e2(3)-e1(3)*e2(2)
      e3(2)=e1(3)*e2(1)-e1(1)*e2(3)
      e3(3)=e1(1)*e2(2)-e1(2)*e2(1)
!
!     re-orthogonalize e2 = e3 x e1
!
      e2(1)=e3(2)*e1(3)-e3(3)*e1(2)
      e2(2)=e3(3)*e1(1)-e3(1)*e1(3)
      e2(3)=e3(1)*e1(2)-e3(2)*e1(1)
!
      do i=1,nope
         do j=1,3
            x_beam(j,i) = xl(j,i) + offsets(1,i)*e1(j) +
     &                    offsets(2,i)*e2(j) + offsets(3,i)*e3(j)
         enddo
      enddo
!
!     calculate physical beam axis e1, length dl
!
      do j=1,3
         e1(j)=x_beam(j,2)-x_beam(j,1)
      enddo
      dl=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
      do j=1,3
         e1(j)=e1(j)/dl
      enddo
!
!     re-compute e3 and e2 based on physical beam axis
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
!     transformation matrix from global to local system
!
!     deformed physical beam coordinates including offset rotations
      do i1=1,nope
         node_id = konl(i1)
         off_vec(1) = x_beam(1,i1) - xl(1,i1)
         off_vec(2) = x_beam(2,i1) - xl(2,i1)
         off_vec(3) = x_beam(3,i1) - xl(3,i1)

         ! vold(1:3) = translation, vold(4:6) = rotation vector
         ! Cross product: Theta x Offset
         rot_disp(1) = vold(5,node_id)*off_vec(3) -
     &                 vold(6,node_id)*off_vec(2)
         rot_disp(2) = vold(6,node_id)*off_vec(1) -
     &                 vold(4,node_id)*off_vec(3)
         rot_disp(3) = vold(4,node_id)*off_vec(2) -
     &                 vold(5,node_id)*off_vec(1)

         do j=1,3
            xl_def(j,i1) = x_beam(j,i1) + vold(j,node_id) + rot_disp(j)
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

!     Reconstruct undeformed axes using the same convention.
      call ccx_beam_section_normal(e1_init,e2_in,e2_init,
     &                             iorient_error)
      if(iorient_error.ne.0) then
         write(*,*) '*ERROR: invalid initial beam section axis'
         write(*,*) '  element ',nelem
         call exit(201)
      endif
      if (dabs(rot_angle_total) .gt. 1.d-6) then
         call rotate_vector_local(e2_init,e1_init,rot_angle_total,
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
!     computation of the local element stiffness & mass matrices
!
      if((stiffness.eq.1).or.(mass.eq.1)) then
!
!        stiffness matrix S' in local coordinates
!
         ! Axial stiffness
         s(1,1) = e*a/dl
         s(1,7) = -s(1,1)
         s(7,7) = s(1,1)

         ! Torsion stiffness
         s(4,4) = um*j_tors/dl
         s(4,10) = -s(4,4)
         s(10,10) = s(4,4)

          ! Bending in local x-y plane (resisted by xi11, involves v (2/8) and r_z (6/12))
          if (timo_val .lt. -1.d-6) then
             ! Pure Euler-Bernoulli mode (KAPPA=0, TIMOSHENKO=NO)
             phi_z = 0.d0
             phi_y = 0.d0
          else if (timo_val .gt. 1.d-6) then
             ! User-specified custom shear coefficient (e.g. KAPPA=0.85, KAPPA=0.8333)
             phi_z = 12.d0 * e * xi11 / (um * a * timo_val * dl**2)
             phi_y = 12.d0 * e * xi22 / (um * a * timo_val * dl**2)
          else
             ! Default: automatic analytical Cowper shear factor
             phi_z = 12.d0 * e * xi11 / (um * a * xk_y * dl**2)
             phi_y = 12.d0 * e * xi22 / (um * a * xk_z * dl**2)
          endif

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

          ! Bending in local x-z plane (resisted by xi22, involves w (3/9) and r_y (5/11))
          s(3,3) = 12.d0*e*xi22/((1.d0 + phi_y)*(dl**3))
          s(3,5) = -6.d0*e*xi22/((1.d0 + phi_y)*(dl**2))
          s(3,9) = -s(3,3)
          s(3,11) = s(3,5)

          s(5,5) = (4.d0 + phi_y)*e*xi22/((1.d0 + phi_y)*dl)
          s(5,9) = 6.d0*e*xi22/((1.d0 + phi_y)*(dl**2))
          s(5,11) = (2.d0 - phi_y)*e*xi22/((1.d0 + phi_y)*dl)

          s(9,9) = s(3,3)
          s(9,11) = -s(3,5)
          s(11,11) = s(5,5)

         ! complete symmetric parts of s
         do i=1,12
            do j=1,i
               s(i,j)=s(j,i)
            enddo
         enddo
!
!        geometric stiffness matrix contribution (P-delta member curvature)
!        n_axial < 0 indicates compression (P = -n_axial > 0)
!        Only computed in nonlinear geometric analysis (NLGEOM, iperturb >= 2)
!
         eps_th = 0.d0
         if (ithermal(1) .gt. 0) then
            eps_th = eth(1)
         endif
         n_axial = 0.d0
         if (iperturb(1) .ge. 2) then
            if (len_init .gt. 1.d-10) then
               n_axial = e * a *
     &              ((len_curr - len_init) / len_init - eps_th)
            endif
         endif

         if (dabs(n_axial) .gt. 1.d-12) then
            c1_g = 6.d0 / (5.d0 * dl)
            c2_g = 1.d0 / 10.d0
            c3_g = 2.d0 * dl / 15.d0
            c4_g = -dl / 30.d0
            c_tors_g = (xi11 + xi22) / (a * dl)

            ! Transverse v (DOFs 2, 6, 8, 12 in local x-y plane)
            s(2,2) = s(2,2) + n_axial * c1_g
            s(2,6) = s(2,6) + n_axial * c2_g
            s(2,8) = s(2,8) - n_axial * c1_g
            s(2,12) = s(2,12) + n_axial * c2_g

            s(6,6) = s(6,6) + n_axial * c3_g
            s(6,8) = s(6,8) - n_axial * c2_g
            s(6,12) = s(6,12) + n_axial * c4_g

            s(8,8) = s(8,8) + n_axial * c1_g
            s(8,12) = s(8,12) - n_axial * c2_g

            s(12,12) = s(12,12) + n_axial * c3_g

            ! Transverse w (DOFs 3, 5, 9, 11 in local x-z plane)
            s(3,3) = s(3,3) + n_axial * c1_g
            s(3,5) = s(3,5) - n_axial * c2_g
            s(3,9) = s(3,9) - n_axial * c1_g
            s(3,11) = s(3,11) - n_axial * c2_g

            s(5,5) = s(5,5) + n_axial * c3_g
            s(5,9) = s(5,9) + n_axial * c2_g
            s(5,11) = s(5,11) + n_axial * c4_g

            s(9,9) = s(9,9) + n_axial * c1_g
            s(9,11) = s(9,11) + n_axial * c2_g

            s(11,11) = s(11,11) + n_axial * c3_g

            ! Torsion (DOFs 4, 10)
            s(4,4) = s(4,4) + n_axial * c_tors_g
            s(4,10) = s(4,10) - n_axial * c_tors_g
            s(10,10) = s(10,10) + n_axial * c_tors_g

            ! complete symmetry
            do i = 1, 12
               do j = 1, i
                  s(i,j) = s(j,i)
               enddo
            enddo
         endif
!
!        consistent mass matrix SM' in local coordinates
!
         if(buckling.eq.1) then
            do i=1,12
               do j=1,12
                  sm(i,j)=0.d0
               enddo
            enddo
            c1 = 12.d0 / (5.d0 * dl)
            c2 = 0.2d0
            c3 = 4.d0 * dl / 15.d0
            c4 = -dl / 15.d0

            ! Transverse v (DOFs 2, 6, 8, 12)
            sm(2,2) = c1
            sm(2,6) = c2
            sm(2,8) = -c1
            sm(2,12) = c2

            sm(6,2) = c2
            sm(6,6) = c3
            sm(6,8) = -c2
            sm(6,12) = c4

            sm(8,2) = -c1
            sm(8,6) = -c2
            sm(8,8) = c1
            sm(8,12) = -c2

            sm(12,2) = c2
            sm(12,6) = c4
            sm(12,8) = -c2
            sm(12,12) = c3

            ! Transverse w (DOFs 3, 5, 9, 11)
            sm(3,3) = c1
            sm(3,5) = -c2
            sm(3,9) = -c1
            sm(3,11) = -c2

            sm(5,3) = -c2
            sm(5,5) = c3
            sm(5,9) = c2
            sm(5,11) = c4

            sm(9,3) = -c1
            sm(9,5) = c2
            sm(9,9) = c1
            sm(9,11) = c2

            sm(11,3) = -c2
            sm(11,5) = c4
            sm(11,9) = c2
            sm(11,11) = c3
         else if(mass.eq.1) then
            do i=1,12
               do j=1,12
                  sm(i,j)=0.d0
               enddo
            enddo

            lumped_flag = 0
            if (iexpl .gt. 1) lumped_flag = 1
            call get_environment_variable(
     &           "CCX_LUMPED_MASS",
     &           env_lump,
     &           status=env_stat)
            if (env_stat .eq. 0) then
               if (env_lump(1:1) .eq. '1' .or.
     &             env_lump(1:1) .eq. 'Y' .or.
     &             env_lump(1:1) .eq. 'y') then
                  lumped_flag = 1
               elseif (env_lump(1:1) .eq. '0' .or.
     &                 env_lump(1:1) .eq. 'N' .or.
     &                 env_lump(1:1) .eq. 'n') then
                  lumped_flag = 0
               endif
            endif

            if (lumped_flag .eq. 1) then
               c_mass = rho * a * dl / 2.d0
               c_tors_mass = rho * (xi11 + xi22) * dl / 2.d0
!              FIX: torsional mass uses Ip = xi11+xi22 (polar moment of area),
!              NOT j_tors (torsional stiffness constant J).  Using J in the
!              mass caused the J to cancel in omega^2=K/M, making the torsion
!              period section-independent.  Ip varies with section shape.
!              NOTE on rotational lumped mass:
!              c_rot_y/z uses the physical cross-section second moment of area
!              (rho*I*L/2), NOT the FEM row-sum of the consistent mass matrix.
!              This is the rigid-body rotational inertia of the cross-section.
!              For slender beams (L >> r_gyr) this understates modal rotational
!              inertia versus the consistent mass, making bending frequencies
!              slightly higher — conservative for explicit CFL stability.
               c_rot_y = rho * xi22 * dl / 2.d0
               c_rot_z = rho * xi11 * dl / 2.d0

               ! Node 1
               sm(1,1) = c_mass
               sm(2,2) = c_mass
               sm(3,3) = c_mass
               sm(4,4) = c_tors_mass
               sm(5,5) = c_rot_y
               sm(6,6) = c_rot_z

               ! Node 2
               sm(7,7) = c_mass
               sm(8,8) = c_mass
               sm(9,9) = c_mass
               sm(10,10) = c_tors_mass
               sm(11,11) = c_rot_y
               sm(12,12) = c_rot_z
            else
               ! Axial mass
               c1 = rho*a*dl/6.d0
               sm(1,1) = 2.d0*c1
               sm(7,7) = 2.d0*c1
               sm(1,7) = c1

               ! Torsional mass — uses Ip = xi11+xi22, NOT j_tors
               c2 = rho*(xi11 + xi22)*dl/6.d0
               sm(4,4) = 2.d0*c2
               sm(10,10) = 2.d0*c2
               sm(4,10) = c2

               ! Bending mass in x-y plane (v and r_z)
               c1 = rho*a*dl/420.d0
               sm(2,2) = 156.d0*c1
               sm(2,6) = 22.d0*dl*c1
               sm(2,8) = 54.d0*c1
               sm(2,12) = -13.d0*dl*c1

               sm(6,6) = 4.d0*dl*dl*c1
               sm(6,8) = 13.d0*dl*c1
               sm(6,12) = -3.d0*dl*dl*c1

               sm(8,8) = 156.d0*c1
               sm(8,12) = -22.d0*dl*c1

               sm(12,12) = 4.d0*dl*dl*c1

               ! Rotary inertia of cross-section for z-rotation (r_z)
               c_rot_z = rho*xi11*dl/6.d0
               sm(6,6) = sm(6,6) + 2.d0*c_rot_z
               sm(12,12) = sm(12,12) + 2.d0*c_rot_z
               sm(6,12) = sm(6,12) + c_rot_z

               ! Bending mass in x-z plane (w and r_y)
               sm(3,3) = 156.d0*c1
               sm(3,5) = -22.d0*dl*c1
               sm(3,9) = 54.d0*c1
               sm(3,11) = 13.d0*dl*c1

               sm(5,5) = 4.d0*dl*dl*c1
               sm(5,9) = -13.d0*dl*c1
               sm(5,11) = -3.d0*dl*dl*c1

               sm(9,9) = 156.d0*c1
               sm(9,11) = 22.d0*dl*c1

               sm(11,11) = 4.d0*dl*dl*c1

               ! Rotary inertia of cross-section for y-rotation (r_y)
               c_rot_y = rho*xi22*dl/6.d0
               sm(5,5) = sm(5,5) + 2.d0*c_rot_y
               sm(11,11) = sm(11,11) + 2.d0*c_rot_y
               sm(5,11) = sm(5,11) + c_rot_y
            endif

            ! complete symmetric parts of sm
            do i=1,12
               do j=1,i
                  sm(i,j)=sm(j,i)
               enddo
            enddo
         endif
!
!        Compute local distributed load vector
!
         do i=1,12
            ff_loc(i)=0.d0
         enddo
          if(rhsi.eq.1) then
             if(ithermal(1).eq.1) then
                n_th = elconloc(1) * a * eth(1)
                ff_loc(1) = ff_loc(1) - n_th
                ff_loc(7) = ff_loc(7) + n_th
             endif
            if(nbody.gt.0) then
               do j=1,3
                  w_glob(j) = bodyfx(j) * rho * a
               enddo
               do j=1,3
                  w_loc(j) = tm(j,1)*w_glob(1) +
     &                       tm(j,2)*w_glob(2) +
     &                       tm(j,3)*w_glob(3)
               enddo
                ff_loc(1) = ff_loc(1) + w_loc(1) * dl / 2.d0
                ff_loc(2) = ff_loc(2) + w_loc(2) * dl / 2.d0
                ff_loc(3) = ff_loc(3) + w_loc(3) * dl / 2.d0
                ff_loc(4) = ff_loc(4) + 0.d0
                ff_loc(5) = ff_loc(5) - w_loc(3) * (dl**2) / 12.d0
                ff_loc(6) = ff_loc(6) + w_loc(2) * (dl**2) / 12.d0
                ff_loc(7) = ff_loc(7) + w_loc(1) * dl / 2.d0
                ff_loc(8) = ff_loc(8) + w_loc(2) * dl / 2.d0
                ff_loc(9) = ff_loc(9) + w_loc(3) * dl / 2.d0
                ff_loc(10) = ff_loc(10) + 0.d0
                ff_loc(11) = ff_loc(11) + w_loc(3) * (dl**2) / 12.d0
                ff_loc(12) = ff_loc(12) - w_loc(2) * (dl**2) / 12.d0
             endif
             if(omx.gt.1.d-10) then
                node1 = kon(indexe+1)
                node2 = kon(indexe+2)

                px1 = co(1, node1) - p1(1)
                py1 = co(2, node1) - p1(2)
                pz1 = co(3, node1) - p1(3)
                c1_cnt = px1*p2(1) + py1*p2(2) + pz1*p2(3)
                ax1 = (px1 - c1_cnt*p2(1)) * omx
                ay1 = (py1 - c1_cnt*p2(2)) * omx
                az1 = (pz1 - c1_cnt*p2(3)) * omx

                px2 = co(1, node2) - p1(1)
                py2 = co(2, node2) - p1(2)
                pz2 = co(3, node2) - p1(3)
                c2_cnt = px2*p2(1) + py2*p2(2) + pz2*p2(3)
                ax2 = (px2 - c2_cnt*p2(1)) * omx
                ay2 = (py2 - c2_cnt*p2(2)) * omx
                az2 = (pz2 - c2_cnt*p2(3)) * omx

                q1_glob(1) = rho * a * ax1
                q1_glob(2) = rho * a * ay1
                q1_glob(3) = rho * a * az1

                q2_glob(1) = rho * a * ax2
                q2_glob(2) = rho * a * ay2
                q2_glob(3) = rho * a * az2

                do j=1,3
                   q1_loc(j) = tm(j,1)*q1_glob(1) +
     &                        tm(j,2)*q1_glob(2) +
     &                        tm(j,3)*q1_glob(3)
                   q2_loc(j) = tm(j,1)*q2_glob(1) +
     &                        tm(j,2)*q2_glob(2) +
     &                        tm(j,3)*q2_glob(3)
                enddo

                ff_loc(1) = ff_loc(1) + dl / 6.d0 *
     &            (2.d0 * q1_loc(1) + q2_loc(1))
                ff_loc(2) = ff_loc(2) + dl / 6.d0 *
     &            (2.d0 * q1_loc(2) + q2_loc(2))
                ff_loc(3) = ff_loc(3) + dl / 6.d0 *
     &            (2.d0 * q1_loc(3) + q2_loc(3))

                ff_loc(7) = ff_loc(7) + dl / 6.d0 *
     &            (q1_loc(1) + 2.d0 * q2_loc(1))
                ff_loc(8) = ff_loc(8) + dl / 6.d0 *
     &            (q1_loc(2) + 2.d0 * q2_loc(2))
                ff_loc(9) = ff_loc(9) + dl / 6.d0 *
     &            (q1_loc(3) + 2.d0 * q2_loc(3))
             endif
             if(idist.ne.0) then
                do k=1,nload
                   if(nelemload(1,k).eq.nelem) then
                      val = xload(1,k)
                       ! Initialize local load parameters
                       load_type = 1
                       load_dir = 0
                       if(sideload(k)(1:2).eq.'PX') then
                           ff_loc(1) = ff_loc(1) + val * dl / 2.d0
                           ff_loc(7) = ff_loc(7) + val * dl / 2.d0
                           load_dir = 3
                        elseif(sideload(k)(1:2).eq.'P1') then
                           load_dir = 1
                           if(sideload(k)(5:7).eq.'_T1') then
                              load_type = 2
                              w1 = 0.d0
                              w2 = val
                              alpha = 0.d0
                              beta = 1.d0
                              f1 = (3.d0/20.d0) * val * dl
                              m1 = (1.d0/30.d0) * val * (dl**2)
                              f2 = (7.d0/20.d0) * val * dl
                              m2 = -(1.d0/20.d0) * val * (dl**2)
                           elseif(sideload(k)(5:7).eq.'_T2') then
                              load_type = 2
                              w1 = val
                              w2 = 0.d0
                              alpha = 0.d0
                              beta = 1.d0
                              f1 = (7.d0/20.d0) * val * dl
                              m1 = (1.d0/20.d0) * val * (dl**2)
                              f2 = (3.d0/20.d0) * val * dl
                              m2 = -(1.d0/30.d0) * val * (dl**2)
                           elseif(sideload(k)(5:7).eq.'_P_') then
                              load_type = 3
                              read(sideload(k)(8:9), '(i2)') a_pct
                              read(sideload(k)(11:12), '(i2)') b_pct
                              alpha = dfloat(a_pct) / 100.d0
                              beta = dfloat(b_pct) / 100.d0
                              w1 = val
                              w2 = val
                              f1 = val * dl * ( (beta - alpha) -
     &                            (beta**3 - alpha**3) +
     &                            0.5d0*(beta**4 - alpha**4) )
                              m1 = val * (dl**2) * (
     &                            0.5d0*(beta**2 - alpha**2) -
     &                            (2.d0/3.d0)*(beta**3 - alpha**3) +
     &                            0.25d0*(beta**4 - alpha**4) )
                              f2 = val * dl * ( (beta**3 - alpha**3) -
     &                            0.5d0*(beta**4 - alpha**4) )
                              m2 = val * (dl**2) * (
     &                            -(1.d0/3.d0)*(beta**3 - alpha**3) +
     &                            0.25d0*(beta**4 - alpha**4) )
                           else
                              load_type = 1
                              w1 = val
                              w2 = val
                              alpha = 0.d0
                              beta = 1.d0
                              f1 = val * dl / 2.d0
                              m1 = val * (dl**2) / 12.d0
                              f2 = val * dl / 2.d0
                              m2 = -val * (dl**2) / 12.d0
                           endif
                           ff_loc(2) = ff_loc(2) + f1
                           ff_loc(6) = ff_loc(6) + m1
                           ff_loc(8) = ff_loc(8) + f2
                           ff_loc(12) = ff_loc(12) + m2
                        elseif(sideload(k)(1:2).eq.'P2') then
                           load_dir = 2
                           if(sideload(k)(5:7).eq.'_T1') then
                              load_type = 2
                              w1 = 0.d0
                              w2 = val
                              alpha = 0.d0
                              beta = 1.d0
                              f1 = (3.d0/20.d0) * val * dl
                              m1 = -(1.d0/30.d0) * val * (dl**2)
                              f2 = (7.d0/20.d0) * val * dl
                              m2 = (1.d0/20.d0) * val * (dl**2)
                           elseif(sideload(k)(5:7).eq.'_T2') then
                              load_type = 2
                              w1 = val
                              w2 = 0.d0
                              alpha = 0.d0
                              beta = 1.d0
                              f1 = (7.d0/20.d0) * val * dl
                              m1 = -(1.d0/20.d0) * val * (dl**2)
                              f2 = (3.d0/20.d0) * val * dl
                              m2 = (1.d0/30.d0) * val * (dl**2)
                           elseif(sideload(k)(5:7).eq.'_P_') then
                              load_type = 3
                              read(sideload(k)(8:9), '(i2)') a_pct
                              read(sideload(k)(11:12), '(i2)') b_pct
                              alpha = dfloat(a_pct) / 100.d0
                              beta = dfloat(b_pct) / 100.d0
                              w1 = val
                              w2 = val
                              f1 = val * dl * ( (beta - alpha) -
     &                            (beta**3 - alpha**3) +
     &                            0.5d0*(beta**4 - alpha**4) )
                              m1 = -val * (dl**2) * (
     &                            0.5d0*(beta**2 - alpha**2) -
     &                            (2.d0/3.d0)*(beta**3 - alpha**3) +
     &                            0.25d0*(beta**4 - alpha**4) )
                              f2 = val * dl * ( (beta**3 - alpha**3) -
     &                            0.5d0*(beta**4 - alpha**4) )
                              m2 = -val * (dl**2) * (
     &                            -(1.d0/3.d0)*(beta**3 - alpha**3) +
     &                            0.25d0*(beta**4 - alpha**4) )
                           else
                              load_type = 1
                              w1 = val
                              w2 = val
                              alpha = 0.d0
                              beta = 1.d0
                              f1 = val * dl / 2.d0
                              m1 = -val * (dl**2) / 12.d0
                              f2 = val * dl / 2.d0
                              m2 = val * (dl**2) / 12.d0
                           endif
                           ff_loc(3) = ff_loc(3) + f1
                           ff_loc(5) = ff_loc(5) + m1
                           ff_loc(9) = ff_loc(9) + f2
                           ff_loc(11) = ff_loc(11) + m2
                        else
                           write(*,*) 'ERROR: Unknown load type'
                        endif
                        ! Save load parameters for resultsmech
                        if (load_dir .gt. 0) then
                           prop(index+28) = dble(load_type)
                           prop(index+29) = dble(load_dir)
                           prop(index+24) = w1
                           prop(index+25) = w2
                           prop(index+26) = alpha
                           prop(index+27) = beta
                        endif
                    endif
                 enddo
              endif
         endif
!
!        Apply local offsets: transforms s and sm from physical beam axes to node offsets
!
         call apply_offsets_local(s, sm, ff_loc, mass, offsets, nope)
!
!        Apply member end releases (hinges/springs) at local DOFs via static condensation
!
         do k = 1, 12
            released(k) = .false.
            spring_k(k) = 0.d0
         enddo
         do inode = 1, nope
            r = (inode-1)*6
            code = release_codes(inode)
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

         call condense_element_local(s, ff_loc, 12, released, spring_k)
         if (rhsi.eq.1) then
            if (.not. allocated(ff_saved)) then
               allocate(ff_saved(60, 100000))
               ff_saved = 0.d0
            endif
            if (nelem.le.100000) then
               ff_saved(1:12, nelem) = ff_loc(1:12)
            endif
         endif
!
!        Transform local load vector to global coordinates
!
         do i=1,12
            ff(i)=0.d0
         enddo
         if(rhsi.eq.1) then
            do inode=1,2
               r = (inode-1)*6
               do j=1,3
                  ff(r+j) = tm(1,j)*ff_loc(r+1) +
     &                      tm(2,j)*ff_loc(r+2) +
     &                      tm(3,j)*ff_loc(r+3)
                  ff(r+3+j) = tm(1,j)*ff_loc(r+4) +
     &                        tm(2,j)*ff_loc(r+5) +
     &                        tm(3,j)*ff_loc(r+6)
               enddo
            enddo
         endif
!
!        Transform resulting local stiffness and mass matrices to global coordinates
!
         do i=1,12
            do j=1,12
               tmg(i,j)=0.d0
            enddo
         enddo
         do i=1,3
            do j=1,3
               tmg(i,j)=tm(i,j)
               tmg(i+3,j+3)=tm(i,j)
               tmg(i+6,j+6)=tm(i,j)
               tmg(i+9,j+9)=tm(i,j)
            enddo
         enddo
!
         if((mass.eq.1).or.(buckling.eq.1)) then
            do i=1,12
               do j=1,12
                  sg(i,j)=0.d0
                  do k=1,12
                     sg(i,j)=sg(i,j)+sm(i,k)*tmg(k,j)
                  enddo
               enddo
            enddo
            do i=1,12
               do j=1,12
                  smg(i,j)=0.d0
                  do k=1,12
                     smg(i,j)=smg(i,j)+tmg(k,i)*sg(k,j)
                  enddo
               enddo
            enddo
            do i=1,12
               do j=1,12
                  sm(i,j)=smg(i,j)
               enddo
            enddo
         endif
!
         do i=1,12
            do j=1,12
               sg(i,j)=0.d0
               do k=1,12
                  sg(i,j)=sg(i,j)+s(i,k)*tmg(k,j)
               enddo
            enddo
         enddo
         do i=1,12
            do j=1,12
               smg(i,j)=0.d0
               do k=1,12
                  smg(i,j)=smg(i,j)+tmg(k,i)*sg(k,j)
               enddo
            enddo
         enddo
         do i=1,12
            do j=1,12
               s(i,j)=smg(i,j)
            enddo
         enddo
       endif
 !
      return
      end
