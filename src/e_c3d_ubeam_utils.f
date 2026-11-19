!
!     Utility subroutines for custom user beam element (UB31)
!     Includes:
!     - Rotation of orientation vector
!     - Analytical calculation of cross-section properties (Rect, Pipe, I, T, U, L)
!     - Principal axes calculation for asymmetric sections to avoid coupling errors
!     - Transformation for local offset offsets
!     - Static condensation for member end releases
!
      subroutine rotate_vector_local(v, axis, angle_deg, v_rot)
      implicit none
      real*8 v(3), axis(3), angle_deg, v_rot(3)
      real*8 theta, cos_t, sin_t, dot_prod, cross_prod(3), pi
      
      pi = 4.d0 * datan(1.d0)
      theta = angle_deg * pi / 180.d0
      cos_t = dcos(theta)
      sin_t = dsin(theta)
      
      dot_prod = v(1)*axis(1) + v(2)*axis(2) + v(3)*axis(3)
      
      cross_prod(1) = axis(2)*v(3) - axis(3)*v(2)
      cross_prod(2) = axis(3)*v(1) - axis(1)*v(3)
      cross_prod(3) = axis(1)*v(2) - axis(2)*v(1)
      
      v_rot(1) = v(1)*cos_t + cross_prod(1)*sin_t + 
     &           axis(1)*dot_prod*(1.d0 - cos_t)
      v_rot(2) = v(2)*cos_t + cross_prod(2)*sin_t + 
     &           axis(2)*dot_prod*(1.d0 - cos_t)
      v_rot(3) = v(3)*cos_t + cross_prod(3)*sin_t + 
     &           axis(3)*dot_prod*(1.d0 - cos_t)
      end subroutine

      subroutine compute_composite_inertia(n_rect, rect_b, rect_h, 
     &  rect_y, rect_z, A, Iyy, Izz, Iyz)
      implicit none
      integer n_rect, i
      real*8 rect_b(n_rect), rect_h(n_rect), rect_y(n_rect), 
     &  rect_z(n_rect)
      real*8 A, Iyy, Izz, Iyz, bar_y, bar_z, sum_A, sum_Ay, sum_Az,
     &  Ai, Iyy_i, Izz_i
      
      sum_A = 0.d0
      sum_Ay = 0.d0
      sum_Az = 0.d0
      
      do i = 1, n_rect
         Ai = rect_b(i) * rect_h(i)
         sum_A = sum_A + Ai
         sum_Ay = sum_Ay + Ai * rect_y(i)
         sum_Az = sum_Az + Ai * rect_z(i)
      enddo
      
      A = sum_A
      bar_y = sum_Ay / A
      bar_z = sum_Az / A
      
      Iyy = 0.d0
      Izz = 0.d0
      Iyz = 0.d0
      
      do i = 1, n_rect
         Ai = rect_b(i) * rect_h(i)
         Iyy_i = (1.d0/12.d0) * rect_b(i) * (rect_h(i)**3)
         Izz_i = (1.d0/12.d0) * rect_h(i) * (rect_b(i)**3)
         
         Iyy = Iyy + Iyy_i + Ai * ((rect_z(i) - bar_z)**2)
         Izz = Izz + Izz_i + Ai * ((rect_y(i) - bar_y)**2)
         Iyz = Iyz + Ai * (rect_y(i) - bar_y) * (rect_z(i) - bar_z)
      enddo
      
      end subroutine

      subroutine compute_section_properties(sect_type, dims, E, nu, 
     &  A, Iyy, Izz, J, xk_y, xk_z, principal_angle)
      implicit none
      real*8 sect_type, dims(6), E, nu, A, Iyy, Izz, J, xk_y, xk_z,
     &  principal_angle
      real*8 b, h, d_o, d_i, r_o, r_i, t_f, t_w, pi,
     &  b_mean, h_mean, t, Iyz, theta_p, Iyy_p, Izz_p
      real*8 rb(6), rh(6), ry(6), rz(6), b_t, t_t, b_b, t_b, b_max
      integer st
      
      pi = 4.d0 * datan(1.d0)
      st = nint(sect_type)
      Iyz = 0.d0
      principal_angle = 0.d0

      if (st .eq. 1) then
         ! Rectangular section (dims: b, h)
         b = dims(1)
         h = dims(2)
         A = b * h
         Iyy = (1.d0/12.d0) * b * (h**3)
         Izz = (1.d0/12.d0) * h * (b**3)
         
         if (b .le. h) then
            J = (1.d0/3.d0) * (b**3) * h * (1.d0 - 0.63d0 * (b/h) * 
     &          (1.d0 - (b**4)/(12.d0 * (h**4))))
         else
            J = (1.d0/3.d0) * (h**3) * b * (1.d0 - 0.63d0 * (h/b) * 
     &          (1.d0 - (h**4)/(12.d0 * (b**4))))
         endif
         
         xk_y = 10.d0 * (1.d0 + nu) / (12.d0 + 11.d0 * nu)
         xk_z = xk_y

      elseif (st .eq. 2) then
         ! Circular/Pipe section (dims: r_outer, thickness)
         ! If thickness is 0, it is a solid circular section of radius r_outer.
         r_o = dims(1)
         t = dims(2)
         if (t .eq. 0.d0) then
            r_i = 0.d0
         else
            r_i = r_o - t
         endif
         
         A = pi * (r_o**2 - r_i**2)
         Izz = (pi / 4.d0) * (r_o**4 - r_i**4)
         Iyy = Izz
         J = (pi / 2.d0) * (r_o**4 - r_i**4)
         
         if (t .eq. 0.d0) then
            xk_y = 6.d0 * (1.d0 + nu) / (7.d0 + 6.d0 * nu)
         else
            xk_y = 2.d0 * (1.d0 + nu) / (4.d0 + 3.d0 * nu)
         endif
         xk_z = xk_y

      elseif (st .eq. 3) then
         ! I-beam section (dims: h, b_top, t_f_top, b_bot, t_f_bot, t_w)
         h = dims(1)
         b_t = dims(2)
         t_t = dims(3)
         b_b = dims(4)
         t_b = dims(5)
         t_w = dims(6)
         
         b_max = max(b_t, b_b)
         
         ! Top flange
         rb(1) = b_t
         rh(1) = t_t
         ry(1) = b_max/2.d0
         rz(1) = h - t_t/2.d0
         
         ! Web
         rb(2) = t_w
         rh(2) = h - t_t - t_b
         ry(2) = b_max/2.d0
         rz(2) = t_b + rh(2)/2.d0
         
         ! Bottom flange
         rb(3) = b_b
         rh(3) = t_b
         ry(3) = b_max/2.d0
         rz(3) = t_b/2.d0
         
         call compute_composite_inertia(3, rb, rh, ry, rz, 
     &        A, Iyy, Izz, Iyz)
         
         J = (1.d0/3.d0) * b_t * (t_t**3) + 
     &       (1.d0/3.d0) * b_b * (t_b**3) + 
     &       (1.d0/3.d0) * (h - t_t - t_b) * (t_w**3)
     
         xk_y = (rb(1)*rh(1) + rb(3)*rh(3)) / A
         xk_z = (rh(2)*rb(2)) / A

      elseif (st .eq. 4) then
         ! T-profile (dims: h, b, t_f, t_w)
         h = dims(1)
         b = dims(2)
         t_f = dims(3)
         t_w = dims(4)
         
         ! Flange (top)
         rb(1) = b
         rh(1) = t_f
         ry(1) = b/2.d0
         rz(1) = h - t_f/2.d0
         
         ! Web (bottom)
         rb(2) = t_w
         rh(2) = h - t_f
         ry(2) = b/2.d0
         rz(2) = rh(2)/2.d0
         
         call compute_composite_inertia(2, rb, rh, ry, rz, 
     &        A, Iyy, Izz, Iyz)
         
         J = (1.d0/3.d0) * b * (t_f**3) + 
     &       (1.d0/3.d0) * (h - t_f) * (t_w**3)
         
         xk_y = (rb(1)*rh(1)) / A
         xk_z = (rh(2)*rb(2)) / A

      elseif (st .eq. 5) then
         ! U-channel section (dims: h, b, t_f, t_w)
         h = dims(1)
         b = dims(2)
         t_f = dims(3)
         t_w = dims(4)
         
         ! Web (vertical)
         rb(1) = t_w
         rh(1) = h
         ry(1) = t_w/2.d0
         rz(1) = h/2.d0
         
         ! Top flange
         rb(2) = b - t_w
         rh(2) = t_f
         ry(2) = t_w + rb(2)/2.d0
         rz(2) = h - t_f/2.d0
         
         ! Bottom flange
         rb(3) = b - t_w
         rh(3) = t_f
         ry(3) = t_w + rb(3)/2.d0
         rz(3) = t_f/2.d0
         
         call compute_composite_inertia(3, rb, rh, ry, rz, 
     &        A, Iyy, Izz, Iyz)
         
         J = (2.d0/3.d0) * (b - t_w/2.d0) * (t_f**3) + 
     &       (1.d0/3.d0) * (h - t_f) * (t_w**3)
         
         xk_y = (rb(2)*rh(2) + rb(3)*rh(3)) / A
         xk_z = (rh(1)*rb(1)) / A

      elseif (st .eq. 6) then
         ! L-angle section (dims: b, h, t)
         b = dims(1)
         h = dims(2)
         t = dims(3)
         
         ! Vertical leg
         rb(1) = t
         rh(1) = h
         ry(1) = t/2.d0
         rz(1) = h/2.d0
         
         ! Horizontal leg
         rb(2) = b - t
         rh(2) = t
         ry(2) = t + rb(2)/2.d0
         rz(2) = t/2.d0
         
         call compute_composite_inertia(2, rb, rh, ry, rz, 
     &        A, Iyy, Izz, Iyz)
         
         J = (1.d0/3.d0) * (b - t/2.d0) * (t**3) + 
     &       (1.d0/3.d0) * (h - t/2.d0) * (t**3)
         
         xk_y = (rb(2)*rh(2)) / A
         xk_z = (rh(1)*rb(1)) / A

      elseif (st .eq. 9) then
!        Abaqus unequal-thickness L section, a,b,t1,t2.
!        Horizontal leg is a*t1, vertical remainder t2*(b-t1).
         b=dims(1)
         h=dims(2)
         t_b=dims(3)
         t_w=dims(4)
         if(b.le.0.d0.or.h.le.0.d0.or.t_b.le.0.d0.or.
     &      t_w.le.0.d0.or.t_b.ge.h.or.t_w.ge.b) then
            write(*,*) '*ERROR: invalid Abaqus L geometry'
            call exit(201)
         endif
         rb(1)=b
         rh(1)=t_b
         ry(1)=b/2.d0
         rz(1)=t_b/2.d0
         rb(2)=t_w
         rh(2)=h-t_b
         ry(2)=t_w/2.d0
         rz(2)=t_b+(h-t_b)/2.d0
         call compute_composite_inertia(2,rb,rh,ry,rz,
     &        A,Iyy,Izz,Iyz)
         J=(b*t_b**3+(h-t_b)*t_w**3)/3.d0
!        Shear factors for thin angles are provisional estimates.
         xk_y=rb(1)*rh(1)/A
         xk_z=rb(2)*rh(2)/A
!
      elseif (st .eq. 7) then
         ! Box section (dims: h, b, t1, t2, t3, t4)
         ! h = overall height, b = overall width
         ! t1=bot, t2=left web, t3=top, t4=right web
         h = dims(1)
         b = dims(2)
         t_b = dims(3) ! t1
         t_w = dims(4) ! t2
         t_f = dims(5) ! t3
         
         ! Top flange
         rb(1) = b
         rh(1) = t_f
         ry(1) = b/2.d0
         rz(1) = h - t_f/2.d0
         
         ! Bottom flange
         rb(2) = b
         rh(2) = t_b
         ry(2) = b/2.d0
         rz(2) = t_b/2.d0
         
         ! Left web
         rb(3) = t_w
         rh(3) = h - t_f - t_b
         ry(3) = t_w/2.d0
         rz(3) = t_b + rh(3)/2.d0
         
         ! Right web
         rb(4) = dims(6) ! t4
         rh(4) = h - t_f - t_b
         ry(4) = b - dims(6)/2.d0
         rz(4) = t_b + rh(4)/2.d0
         
         call compute_composite_inertia(4, rb, rh, ry, rz, 
     &        A, Iyy, Izz, Iyz)
         
         h_mean = h - 0.5d0 * (t_f + t_b)
         b_mean = b - 0.5d0 * (t_w + dims(6))
!
!        Closed boxes use the Bredt thin-wall torsion expression.
!        Mecway may encode open sections as BOX with one or more wall
!        thicknesses equal to zero.  In that case use the open thin-wall
!        sum J = sum(b_i*t_i^3/3) rather than dividing by a zero wall.
!
         if((t_f.gt.1.d-14).and.(t_b.gt.1.d-14).and.
     &      (t_w.gt.1.d-14).and.(dims(6).gt.1.d-14)) then
            J = 4.d0 * (h_mean * b_mean)**2 / 
     &          (b_mean*(1.d0/t_f + 1.d0/t_b) + 
     &           h_mean*(1.d0/t_w + 1.d0/dims(6)))
         else
            J = (1.d0/3.d0) * b * (t_f**3) +
     &          (1.d0/3.d0) * b * (t_b**3) +
     &          (1.d0/3.d0) * max(0.d0,h-t_f-t_b) * (t_w**3) +
     &          (1.d0/3.d0) * max(0.d0,h-t_f-t_b) *
     &          (dims(6)**3)
         endif
         
         xk_y = (rh(3)*rb(3) + rh(4)*rb(4)) / A
         xk_z = (rb(1)*rh(1) + rb(2)*rh(2)) / A

      elseif (st .eq. 8) then
!        CCX CIRC elliptical section, semiaxes along section n1/n2.
         b=dims(1)
         h=dims(2)
         if(b.le.0.d0.or.h.le.0.d0) then
            write(*,*) '*ERROR: invalid elliptical section'
            call exit(201)
         endif
         A=pi*b*h
         Iyy=pi*b*h**3/4.d0
         Izz=pi*h*b**3/4.d0
         J=pi*b**3*h**3/(b*b+h*h)
!        Approximate Cowper shear factor; ellipse needs verification.
         xk_y=6.d0*(1.d0+nu)/(7.d0+6.d0*nu)
         xk_z=xk_y
      else
         write(*,*) '*ERROR: unknown UB31 section ',st
         call exit(201)
      endif

      ! Rotate to principal coordinates if asymmetrical (Iyz is non-zero)
      if (dabs(Iyz) .gt. 1.d-12) then
         theta_p = 0.5d0 * datan2(-2.d0 * Iyz, Iyy - Izz)
         principal_angle = theta_p * 180.d0 / pi
         
         ! Principal moments of inertia
         Iyy_p = 0.5d0*(Iyy + Izz) +
     &           0.5d0*(Iyy - Izz)*dcos(2.d0*theta_p) -
     &           Iyz*dsin(2.d0*theta_p)
         Izz_p = 0.5d0*(Iyy + Izz) -
     &           0.5d0*(Iyy - Izz)*dcos(2.d0*theta_p) +
     &           Iyz*dsin(2.d0*theta_p)

         Iyy = Iyy_p
         Izz = Izz_p
      endif

      end subroutine

      subroutine apply_offsets_local(s, sm, ff, mass, offsets, nope)
      implicit none
      integer mass, nope, i, j, k, n, inode
      real*8 s(60,60), sm(60,60), ff(60), offsets(3,*), 
     &  T(60,60), Temp(60,60)

      n = nope * 6

      ! Initialize T as identity matrix
      do i = 1, n
         do j = 1, n
            T(i,j) = 0.d0
         enddo
         T(i,i) = 1.d0
      enddo

      ! Fill T with offset terms
      do inode = 1, nope
         k = (inode-1)*6
         ! u_beam = u_node + theta_y * oz - theta_z * oy
         T(k+1, k+5) = offsets(3, inode)  ! oz
         T(k+1, k+6) = -offsets(2, inode) ! -oy
         
         ! v_beam = v_node + theta_z * ox - theta_x * oz
         T(k+2, k+4) = -offsets(3, inode) ! -oz
         T(k+2, k+6) = offsets(1, inode)  ! ox
         
         ! w_beam = w_node + theta_x * oy - theta_y * ox
         T(k+3, k+4) = offsets(2, inode)  ! oy
         T(k+3, k+5) = -offsets(1, inode) ! -ox
      enddo

      ! Transform s: S_node = T^T * S_beam * T
      ! Temp = S_beam * T
      do i = 1, n
         do j = 1, n
            Temp(i,j) = 0.d0
            do k = 1, n
               Temp(i,j) = Temp(i,j) + s(i,k) * T(k,j)
            enddo
         enddo
      enddo
      ! s = T^T * Temp
      do i = 1, n
         do j = 1, n
            s(i,j) = 0.d0
            do k = 1, n
               s(i,j) = s(i,j) + T(k,i) * Temp(k,j)
            enddo
         enddo
      enddo

      ! Transform sm if mass = 1: SM_node = T^T * SM_beam * T
      if (mass .eq. 1) then
         do i = 1, n
            do j = 1, n
               Temp(i,j) = 0.d0
               do k = 1, n
                  Temp(i,j) = Temp(i,j) + sm(i,k) * T(k,j)
               enddo
            enddo
         enddo
         do i = 1, n
            do j = 1, n
               sm(i,j) = 0.d0
               do k = 1, n
                  sm(i,j) = sm(i,j) + T(k,i) * Temp(k,j)
               enddo
            enddo
         enddo
      endif

      ! Transform ff: ff_node = T^T * ff_beam
      do i = 1, n
         Temp(i,1) = ff(i)
      enddo
      do i = 1, n
         ff(i) = 0.d0
         do k = 1, n
            ff(i) = ff(i) + T(k,i) * Temp(k,1)
         enddo
      enddo

      end subroutine

      subroutine condense_element_local(s, ff, n, released, spring_k)
      implicit none
      integer n, i, j, r
      real*8 s(60,60), ff(60), spring_k(60), col_r(60)
      real*8 factor, denom, ks, s_rr
      logical released(60)

      do r = 1, n
         if (released(r)) then
            s_rr = s(r,r)
            ks = spring_k(r)
            if (ks .lt. 0.d0) ks = 0.d0
            denom = s_rr + ks
            if (dabs(denom) .gt. 1.d-12) then
               do i = 1, n
                  col_r(i) = s(i,r)
               enddo
               ! update load vector
               do i = 1, n
                  if (i .ne. r) then
                     ff(i) = ff(i) - (col_r(i) / denom) * ff(r)
                  endif
               enddo
               ff(r) = (ks / denom) * ff(r)
               ! update stiffness matrix: S*_ij = S_ij - col_r(i)*col_r(j)/denom
               do i = 1, n
                  do j = 1, n
                     s(i,j) = s(i,j) - (col_r(i) * col_r(j)) / denom
                  enddo
               enddo
               ! add small stabilization for pure pin (Ks == 0) to avoid zero pivot
               if (ks .le. 1.d-12) then
                  s(r,r) = 1.d-9 * s_rr
               endif
            endif
         endif
      enddo
      end subroutine
!     Orient the first section axis using the CalculiX convention.
!     A missing direction is treated as (0,0,-1); a parallel
!     direction is invalid (no arbitrary "upright" fallback).
      subroutine ccx_beam_section_normal(t,provided,n1,ierr)
      implicit none
      real*8 t(3),provided(3),n1(3),p(3),m,dotval
      integer ierr,j
      ierr=0
      m=dsqrt(provided(1)**2+provided(2)**2+provided(3)**2)
      if(m.le.1.d-12) then
         p(1)=0.d0
         p(2)=0.d0
         p(3)=-1.d0
      else
         do j=1,3
            p(j)=provided(j)/m
         enddo
      endif
      dotval=p(1)*t(1)+p(2)*t(2)+p(3)*t(3)
      do j=1,3
         n1(j)=p(j)-dotval*t(j)
      enddo
      m=dsqrt(n1(1)**2+n1(2)**2+n1(3)**2)
      if(m.lt.1.d-8) then
         ierr=1
         return
      endif
      do j=1,3
         n1(j)=n1(j)/m
      enddo
      end subroutine ccx_beam_section_normal
