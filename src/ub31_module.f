!     Module holding shared data for UB31 custom beam element.
!
!     ff_saved(60, ne):   saved element force vectors from assembly.
!
!     ub31_stx(6, 11, ne): 11-station stress/force data per UB31.
!                          Indexed as ub31_stx(component, station, elem)
!                          Stations 1..11 correspond to xi = 0.0..1.0.
!
!     ielrelease(3, ne):  per-member release codes for beam elements.
!                          ielrelease(1, ielem): node 1 release code
!                          ielrelease(2, ielem): node 2 release code
!                          ielrelease(3, ielem): node 3 release code
!
!     relspring(6, 3, ne): per-member spring stiffnesses for releases.
!                          relspring(dof, inode, ielem): stiffness
!
!     *USER BEAM OUTPUT control variables:
!     out_active:          master switch for custom beam output
!     out_num_targets:     number of target (elset, file) pairs
!     out_target_elset:    ELSET name for target j (1..MAX_TARGETS)
!     out_target_filename: CSV filename for target j (1..MAX_TARGETS)
!     out_target_file_init:file initialization flag for target j
!     out_elem_active:     2D logical mask (elem, target_j)
!     out_subdivisions:    number of subdivisions along span (N)
!     out_inc_mode:        increment filtering ('LAST', 'ALL ', etc.)
!     out_inc_freq:        frequency N for FREQ mode
!     out_inc_list(50):    list of specific increment numbers
!     out_ninc_list:       number of increments in out_inc_list
!     out_coords:          coordinate system (1 = LOCAL, 2 = GLOBAL)
!     out_flag_f:          output forces & moments (Fx..Mz)
!     out_flag_u:          output displacements & rotations (ux..rz)
!     out_flag_s:          output stresses (Sxx_ax..Stors)
!     out_flag_q:          output active distributed loads (qx..qz)
!
      module ub31_module
        implicit none
        real*8, allocatable, save :: ff_saved(:,:)
        real*8, allocatable, save :: ub31_stx(:,:,:)
        real*8, allocatable, save :: ub31_for(:,:,:)
        integer, allocatable, save :: ielrelease(:,:)
        real*8, allocatable, save :: relspring(:,:,:)

        integer, parameter :: MAX_TARGETS = 50
        logical, save :: out_active = .true.
        integer, save :: out_num_targets = 1
        character*80, save :: out_target_filename(MAX_TARGETS) =
     &       'ub31_beam_forces.csv'
        character*81, save :: out_target_elset(MAX_TARGETS) = ' '
        logical, save :: out_target_file_init(MAX_TARGETS) = .false.
        character*80, save :: out_filename = 'ub31_beam_forces.csv'
        character*81, save :: out_elset = ' '
        integer, save :: out_subdivisions = 10
        character*4, save :: out_inc_mode = 'LAST'
        integer, save :: out_inc_freq = 1
        integer, save :: out_inc_list(50) = 0
        integer, save :: out_ninc_list = 0
        integer, save :: out_coords = 1
        logical, save :: out_flag_f = .true.
        logical, save :: out_flag_u = .false.
        logical, save :: out_flag_s = .true.
        logical, save :: out_flag_pts = .false.
        logical, save :: out_flag_q = .false.
        logical, save :: out_file_initialized = .false.
        logical, allocatable, save :: out_elem_active(:,:)
      end module ub31_module
