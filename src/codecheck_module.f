!     ==================================================================
!     codecheck_module.f
!     Holds data structures, parameters, member overrides, and results
!     for CalculiX Beam Code-Checking Framework (Eurocode 3 EN 1993-1-1)
!     ==================================================================
      module codecheck_module
        implicit none

        logical, save :: codecheck_active = .false.
        character*16, save :: codecheck_code = 'EC3'
        character*132, save :: codecheck_filename = ' '
        character*81, save :: codecheck_elset = ' '
        integer, save :: codecheck_evaluated_step = -1
        integer, save :: codecheck_evaluated_inc = -1
        integer, save :: codecheck_output_level = 1
!       1: ALL_STATIONS (default, 11 stations per member per step)
!       2: GOVERNING (1 governing station row per member per step)
!       3: ENVELOPE (governing envelope across all steps)
!       4: FAILS_ONLY (only stations with UC > 1.0)
        integer, save :: codecheck_subdivisions = 10

!       Member design parameters (allocated dynamically)
        real*8, allocatable, save :: lcr_y(:)
        real*8, allocatable, save :: lcr_z(:)
        real*8, allocatable, save :: l_lt(:)
        real*8, allocatable, save :: k_y(:)
        real*8, allocatable, save :: k_z(:)
        real*8, allocatable, save :: k_lt(:)
        real*8, allocatable, save :: c1_lt(:)
        real*8, allocatable, save :: fy_user(:)
        real*8, allocatable, save :: gamma_m0(:)
        real*8, allocatable, save :: gamma_m1(:)
        real*8, allocatable, save :: phi_c(:)
        real*8, allocatable, save :: phi_b(:)
        real*8, allocatable, save :: phi_v(:)
        real*8, allocatable, save :: phi_t(:)

!       Cached material properties
        real*8, allocatable, save :: elem_E(:)
        real*8, allocatable, save :: elem_nu(:)
        real*8, allocatable, save :: elem_fy(:)

!       Governing and station results per member (Current Step)
        real*8, allocatable, save :: uc_tot(:)
        real*8, allocatable, save :: uc_ax(:)
        real*8, allocatable, save :: uc_sh(:)
        real*8, allocatable, save :: uc_bnd(:)
        real*8, allocatable, save :: uc_stab(:)
        real*8, allocatable, save :: uc_ltb(:)
        integer, allocatable, save :: uc_gov_station(:)
        character*20, allocatable, save :: uc_gov_mode(:)
        real*8, allocatable, save :: station_uc(:,:,:)

!       Envelope tracking across all analysis steps
        real*8, allocatable, save :: env_max_uc(:)
        integer, allocatable, save :: env_gov_step(:)
        integer, allocatable, save :: env_gov_inc(:)
        integer, allocatable, save :: env_gov_station(:)
        real*8, allocatable, save :: env_gov_xi(:)
        real*8, allocatable, save :: env_gov_xloc(:)
        real*8, allocatable, save :: env_gov_forces(:,:)
        real*8, allocatable, save :: env_gov_ucs(:,:)
        character*20, allocatable, save :: env_gov_mode(:)

      contains

!       ================================================================
!       ensure_codecheck_alloc: Ensures dynamic arrays are allocated
!       ================================================================
        subroutine ensure_codecheck_alloc(ne_in)
          implicit none
          integer, intent(in) :: ne_in
          integer :: cur_sz, req_sz

          req_sz = max(ne_in, 1000)

          if (.not. allocated(lcr_y)) then
             allocate(lcr_y(req_sz))
             allocate(lcr_z(req_sz))
             allocate(l_lt(req_sz))
             allocate(k_y(req_sz))
             allocate(k_z(req_sz))
             allocate(k_lt(req_sz))
             allocate(c1_lt(req_sz))
             allocate(fy_user(req_sz))
             allocate(gamma_m0(req_sz))
             allocate(gamma_m1(req_sz))
             allocate(phi_c(req_sz))
             allocate(phi_b(req_sz))
             allocate(phi_v(req_sz))
             allocate(phi_t(req_sz))

             allocate(elem_E(req_sz))
             allocate(elem_nu(req_sz))
             allocate(elem_fy(req_sz))

             allocate(uc_tot(req_sz))
             allocate(uc_ax(req_sz))
             allocate(uc_sh(req_sz))
             allocate(uc_bnd(req_sz))
             allocate(uc_stab(req_sz))
             allocate(uc_ltb(req_sz))
             allocate(uc_gov_station(req_sz))
             allocate(uc_gov_mode(req_sz))
             allocate(station_uc(6, 11, req_sz))

             allocate(env_max_uc(req_sz))
             allocate(env_gov_step(req_sz))
             allocate(env_gov_inc(req_sz))
             allocate(env_gov_station(req_sz))
             allocate(env_gov_xi(req_sz))
             allocate(env_gov_xloc(req_sz))
             allocate(env_gov_forces(6, req_sz))
             allocate(env_gov_ucs(5, req_sz))
             allocate(env_gov_mode(req_sz))

             lcr_y = -1.d0
             lcr_z = -1.d0
             l_lt  = -1.d0
             k_y   = 1.0d0
             k_z   = 1.0d0
             k_lt  = 1.0d0
             c1_lt = 1.0d0
             fy_user = 0.d0
             gamma_m0 = 1.0d0
             gamma_m1 = 1.0d0
             phi_c = 0.90d0
             phi_b = 0.90d0
             phi_v = 0.90d0
             phi_t = 0.90d0

             elem_E = 2.1d11
             elem_nu = 0.3d0
             elem_fy = 235.0d6

             uc_tot = 0.d0
             uc_ax = 0.d0
             uc_sh = 0.d0
             uc_bnd = 0.d0
             uc_stab = 0.d0
             uc_ltb = 0.d0
             uc_gov_station = 0
             uc_gov_mode = ' '
             station_uc = 0.d0

             env_max_uc = -1.d0
             env_gov_step = 0
             env_gov_inc = 0
             env_gov_station = 0
             env_gov_xi = 0.d0
             env_gov_xloc = 0.d0
             env_gov_forces = 0.d0
             env_gov_ucs = 0.d0
             env_gov_mode = ' '
          else
             cur_sz = size(lcr_y)
             if (req_sz .gt. cur_sz) then
                call reallocate_codecheck_arrays(cur_sz, req_sz)
             endif
          endif
        end subroutine ensure_codecheck_alloc

        subroutine reallocate_codecheck_arrays(cur_sz, req_sz)
          implicit none
          integer, intent(in) :: cur_sz, req_sz
          real*8, allocatable :: tmp_r(:), tmp_r2(:,:), tmp_r3(:,:,:)
          integer, allocatable :: tmp_i(:)
          character*20, allocatable :: tmp_c(:)

          allocate(tmp_r(req_sz))
          tmp_r = -1.d0
          tmp_r(1:cur_sz) = lcr_y(1:cur_sz)
          call move_alloc(tmp_r, lcr_y)

          allocate(tmp_r(req_sz))
          tmp_r = -1.d0
          tmp_r(1:cur_sz) = lcr_z(1:cur_sz)
          call move_alloc(tmp_r, lcr_z)

          allocate(tmp_r(req_sz))
          tmp_r = -1.d0
          tmp_r(1:cur_sz) = l_lt(1:cur_sz)
          call move_alloc(tmp_r, l_lt)

          allocate(tmp_r(req_sz))
          tmp_r = 1.0d0
          tmp_r(1:cur_sz) = k_y(1:cur_sz)
          call move_alloc(tmp_r, k_y)

          allocate(tmp_r(req_sz))
          tmp_r = 1.0d0
          tmp_r(1:cur_sz) = k_z(1:cur_sz)
          call move_alloc(tmp_r, k_z)

          allocate(tmp_r(req_sz))
          tmp_r = 1.0d0
          tmp_r(1:cur_sz) = k_lt(1:cur_sz)
          call move_alloc(tmp_r, k_lt)

          allocate(tmp_r(req_sz))
          tmp_r = 1.0d0
          tmp_r(1:cur_sz) = c1_lt(1:cur_sz)
          call move_alloc(tmp_r, c1_lt)

          allocate(tmp_r(req_sz))
          tmp_r = 0.0d0
          tmp_r(1:cur_sz) = fy_user(1:cur_sz)
          call move_alloc(tmp_r, fy_user)

          allocate(tmp_r(req_sz))
          tmp_r = 1.0d0
          tmp_r(1:cur_sz) = gamma_m0(1:cur_sz)
          call move_alloc(tmp_r, gamma_m0)

          allocate(tmp_r(req_sz))
          tmp_r = 1.0d0
          tmp_r(1:cur_sz) = gamma_m1(1:cur_sz)
          call move_alloc(tmp_r, gamma_m1)

          allocate(tmp_r(req_sz))
          tmp_r = 0.90d0
          tmp_r(1:cur_sz) = phi_c(1:cur_sz)
          call move_alloc(tmp_r, phi_c)

          allocate(tmp_r(req_sz))
          tmp_r = 0.90d0
          tmp_r(1:cur_sz) = phi_b(1:cur_sz)
          call move_alloc(tmp_r, phi_b)

          allocate(tmp_r(req_sz))
          tmp_r = 0.90d0
          tmp_r(1:cur_sz) = phi_v(1:cur_sz)
          call move_alloc(tmp_r, phi_v)

          allocate(tmp_r(req_sz))
          tmp_r = 0.90d0
          tmp_r(1:cur_sz) = phi_t(1:cur_sz)
          call move_alloc(tmp_r, phi_t)

          allocate(tmp_r(req_sz))
          tmp_r = 2.1d11
          tmp_r(1:cur_sz) = elem_E(1:cur_sz)
          call move_alloc(tmp_r, elem_E)

          allocate(tmp_r(req_sz))
          tmp_r = 0.3d0
          tmp_r(1:cur_sz) = elem_nu(1:cur_sz)
          call move_alloc(tmp_r, elem_nu)

          allocate(tmp_r(req_sz))
          tmp_r = 235.0d6
          tmp_r(1:cur_sz) = elem_fy(1:cur_sz)
          call move_alloc(tmp_r, elem_fy)

          allocate(tmp_r(req_sz))
          tmp_r = 0.0d0
          tmp_r(1:cur_sz) = uc_tot(1:cur_sz)
          call move_alloc(tmp_r, uc_tot)

          allocate(tmp_r(req_sz))
          tmp_r = 0.0d0
          tmp_r(1:cur_sz) = uc_ax(1:cur_sz)
          call move_alloc(tmp_r, uc_ax)

          allocate(tmp_r(req_sz))
          tmp_r = 0.0d0
          tmp_r(1:cur_sz) = uc_sh(1:cur_sz)
          call move_alloc(tmp_r, uc_sh)

          allocate(tmp_r(req_sz))
          tmp_r = 0.0d0
          tmp_r(1:cur_sz) = uc_bnd(1:cur_sz)
          call move_alloc(tmp_r, uc_bnd)

          allocate(tmp_r(req_sz))
          tmp_r = 0.0d0
          tmp_r(1:cur_sz) = uc_stab(1:cur_sz)
          call move_alloc(tmp_r, uc_stab)

          allocate(tmp_r(req_sz))
          tmp_r = 0.0d0
          tmp_r(1:cur_sz) = uc_ltb(1:cur_sz)
          call move_alloc(tmp_r, uc_ltb)

          allocate(tmp_i(req_sz))
          tmp_i = 0
          tmp_i(1:cur_sz) = uc_gov_station(1:cur_sz)
          call move_alloc(tmp_i, uc_gov_station)

          allocate(tmp_c(req_sz))
          tmp_c = ' '
          tmp_c(1:cur_sz) = uc_gov_mode(1:cur_sz)
          call move_alloc(tmp_c, uc_gov_mode)

          allocate(tmp_r3(6, 11, req_sz))
          tmp_r3 = 0.0d0
          tmp_r3(1:6, 1:11, 1:cur_sz) = station_uc(1:6, 1:11, 1:cur_sz)
          call move_alloc(tmp_r3, station_uc)

          allocate(tmp_r(req_sz))
          tmp_r = -1.d0
          tmp_r(1:cur_sz) = env_max_uc(1:cur_sz)
          call move_alloc(tmp_r, env_max_uc)

          allocate(tmp_i(req_sz))
          tmp_i = 0
          tmp_i(1:cur_sz) = env_gov_step(1:cur_sz)
          call move_alloc(tmp_i, env_gov_step)

          allocate(tmp_i(req_sz))
          tmp_i = 0
          tmp_i(1:cur_sz) = env_gov_inc(1:cur_sz)
          call move_alloc(tmp_i, env_gov_inc)

          allocate(tmp_i(req_sz))
          tmp_i = 0
          tmp_i(1:cur_sz) = env_gov_station(1:cur_sz)
          call move_alloc(tmp_i, env_gov_station)

          allocate(tmp_r(req_sz))
          tmp_r = 0.d0
          tmp_r(1:cur_sz) = env_gov_xi(1:cur_sz)
          call move_alloc(tmp_r, env_gov_xi)

          allocate(tmp_r(req_sz))
          tmp_r = 0.d0
          tmp_r(1:cur_sz) = env_gov_xloc(1:cur_sz)
          call move_alloc(tmp_r, env_gov_xloc)

          allocate(tmp_r2(6, req_sz))
          tmp_r2 = 0.d0
          tmp_r2(1:6, 1:cur_sz) = env_gov_forces(1:6, 1:cur_sz)
          call move_alloc(tmp_r2, env_gov_forces)

          allocate(tmp_r2(5, req_sz))
          tmp_r2 = 0.d0
          tmp_r2(1:5, 1:cur_sz) = env_gov_ucs(1:5, 1:cur_sz)
          call move_alloc(tmp_r2, env_gov_ucs)

          allocate(tmp_c(req_sz))
          tmp_c = ' '
          tmp_c(1:cur_sz) = env_gov_mode(1:cur_sz)
          call move_alloc(tmp_c, env_gov_mode)
        end subroutine reallocate_codecheck_arrays

      end module codecheck_module
