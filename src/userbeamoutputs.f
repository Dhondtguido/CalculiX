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
      subroutine userbeamoutputs(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,lakon,ne,irstrt,istep,istat,n,key,iline,ipol,inl,
     &  ipoinp,inp,ipoinpc,jobnamec,ier)
!
!     reading the input deck: *USER BEAM OUTPUT
!
!     Configures dynamic multi-station CSV output for UB31 beam elements
!
!     Syntax:
!       *USER BEAM OUTPUT [, ELSET=<set_name>] [, FILE=<filename.csv>]
!                         [, SUBDIVISIONS=<N>]
!                         [, INCREMENT=<LAST|ALL|N|(list)>]
!                         [, COORDINATES=<LOCAL|GLOBAL>]
!       <VARIABLE_1>, <VARIABLE_2>, ...
!     Parser for *USER BEAM OUTPUT keyword card and data lines
!
      use ub31_module
      implicit none
!
      character*132 textpart(16)
      character*81 set(*)
      character*8 lakon(*)
      character*132 jobnamec(*)
      character*1 inpc(*)
      character*1000 header_str
      character*80 token, param_val
      character*81 elset
      character*80 elset_list(MAX_TARGETS), file_list(MAX_TARGETS)
      character*80 single_file
!
      integer istartset(*), iendset(*), ialset(*), irstrt(*),
     &  ipoinp(2,*), inp(3,*), ipoinpc(0:*)
      integer nset, ne, istep, istat, n, key, iline, ipol, inl, ier
      integer i, m, ipos, len_t, istat_val, id_set, len_job, j, k, elem
      integer n_elsets, n_files, itarg, ipos_slash
      logical has_keys, has_file_param
!
      out_active = .true.
      has_file_param = .false.
!
!     Reassemble header parameters to handle tuple syntax
!
      header_str = ' '
      do i = 2, n
         if (i .eq. 2) then
            header_str = trim(textpart(i))
         else
            header_str = trim(header_str) // ',' // trim(textpart(i))
         endif
      enddo
!
!     Extract ELSET list and FILE list from header_str
!
      call get_param_tuple_or_val(header_str, 'ELSET=',
     &     elset_list, n_elsets, MAX_TARGETS)
      call get_param_tuple_or_val(header_str, 'FILE=',
     &     file_list, n_files, MAX_TARGETS)
!
      if (n_files .gt. 0) has_file_param = .true.
!
!     Determine number of target outputs
!
      if (n_elsets .le. 0) then
         out_num_targets = 1
         out_target_elset(1) = ' '
         if (n_files .ge. 1) then
            out_target_filename(1) = file_list(1)
         else
            len_job = len_trim(jobnamec(1))
            if (len_job .gt. 0) then
               out_target_filename(1) = jobnamec(1)(1:len_job)//
     &              '_beam.csv'
            else
               out_target_filename(1) = 'ub31_beam_forces.csv'
            endif
         endif
      else
         out_num_targets = n_elsets
         do itarg = 1, out_num_targets
            out_target_elset(itarg) = elset_list(itarg)
            if (itarg .le. n_files) then
               out_target_filename(itarg) = file_list(itarg)
            else
               len_job = len_trim(jobnamec(1))
               if (len_job .gt. 0) then
                  out_target_filename(itarg) = jobnamec(1)(1:len_job)//
     &                 '_' // trim(elset_list(itarg)) // '.csv'
               else
                  out_target_filename(itarg) =
     &                 trim(elset_list(itarg)) // '.csv'
               endif
            endif
         enddo
      endif
!
!     Prepend directory from jobnamec(1) if relative filename without /
!
      len_job = len_trim(jobnamec(1))
      ipos_slash = 0
      do m = len_job, 1, -1
         if (jobnamec(1)(m:m) .eq. '/') then
            ipos_slash = m
            exit
         endif
      enddo
      if (ipos_slash .gt. 0) then
         do itarg = 1, out_num_targets
            if (index(out_target_filename(itarg), '/') .eq. 0) then
               out_target_filename(itarg) =
     &              jobnamec(1)(1:ipos_slash) //
     &              trim(out_target_filename(itarg))
            endif
         enddo
      endif
!
!     Parse remaining options from textpart tokens
!
      do i = 2, n
         token = adjustl(textpart(i)(1:80))
         len_t = len_trim(token)
         if (len_t .eq. 0) cycle
!
!        Convert keyword parameter name up to '=' to uppercase
         ipos = index(token, '=')
         if (ipos .gt. 0) then
            do m = 1, ipos
               if (token(m:m) .ge. 'a' .and. token(m:m) .le. 'z') then
                  token(m:m) = char(ichar(token(m:m)) - 32)
               endif
            enddo
         else
            do m = 1, len_t
               if (token(m:m) .ge. 'a' .and. token(m:m) .le. 'z') then
                  token(m:m) = char(ichar(token(m:m)) - 32)
               endif
            enddo
         endif
!
         if (token(1:6) .eq. 'ELSET=' .or.
     &       token(1:5) .eq. 'FILE=') then
            ! Handled via get_param_tuple_or_val
            cycle
         elseif (token(1:13) .eq. 'SUBDIVISIONS=') then
            read(token(14:len_t), *, iostat=istat_val) out_subdivisions
            if (istat_val .ne. 0 .or. out_subdivisions .lt. 1) then
               write(*,*) '*WARNING reading *USER BEAM OUTPUT: '
               write(*,*) '         invalid SUBDIVISIONS value; '//
     &                    'defaulting to 10.'
               call inputwarning(inpc,ipoinpc,iline,
     &              "*USER BEAM OUTPUT%")
               out_subdivisions = 10
            endif
         elseif (token(1:9) .eq. 'SEGMENTS=') then
            read(token(10:len_t), *, iostat=istat_val) out_subdivisions
            if (istat_val .ne. 0 .or. out_subdivisions .lt. 1) then
               write(*,*) '*WARNING reading *USER BEAM OUTPUT: '
               write(*,*) '         invalid SEGMENTS value; '//
     &                    'defaulting to 10.'
               call inputwarning(inpc,ipoinpc,iline,
     &              "*USER BEAM OUTPUT%")
               out_subdivisions = 10
            endif
         elseif (token(1:9) .eq. 'STATIONS=') then
            read(token(10:len_t), *, iostat=istat_val) j
            if (istat_val .ne. 0 .or. j .lt. 2) then
               write(*,*) '*WARNING reading *USER BEAM OUTPUT: '
               write(*,*) '         invalid STATIONS value; '//
     &                    'defaulting to 11.'
               call inputwarning(inpc,ipoinpc,iline,
     &              "*USER BEAM OUTPUT%")
               out_subdivisions = 10
            else
               out_subdivisions = j - 1
            endif
         elseif (token(1:10) .eq. 'INCREMENT=') then
            param_val = token(11:len_t)
            call parse_increment_param(header_str, param_val,
     &           out_inc_mode, out_inc_freq, out_inc_list,
     &           out_ninc_list, inpc, ipoinpc, iline)
         elseif (token(1:12) .eq. 'COORDINATES=' .or.
     &           token(1:6) .eq. 'COORD=') then
            if (token(1:12) .eq. 'COORDINATES=') then
               param_val = token(13:len_t)
            else
               param_val = token(7:len_t)
            endif
            do m = 1, len_trim(param_val)
               if (param_val(m:m).ge.'a' .and.
     &             param_val(m:m).le.'z') then
                  param_val(m:m) = char(ichar(param_val(m:m)) - 32)
               endif
            enddo
            if (param_val(1:6) .eq. 'GLOBAL' .or.
     &          param_val(1:1) .eq. 'G') then
               out_coords = 2
            else
               out_coords = 1
            endif
         elseif (token(1:14) .eq. 'SECTION_POINTS' .or.
     &           token(1:13) .eq. 'SECTIONPOINTS' .or.
     &           token(1:6) .eq. 'POINTS' .or.
     &           token(1:7) .eq. 'CORNERS' .or.
     &           token(1:3) .eq. 'PTS' .or.
     &           token(1:4) .eq. 'SPTS') then
            if (index(token, '=NO').gt.0 .or. 
     &          index(token, '=FALSE').gt.0 .or.
     &          index(token, '=OFF').gt.0) then
               out_flag_pts = .false.
            else
               out_flag_pts = .true.
            endif
         else
            if (index(header_str, 'INCREMENT=(') .le. 0 .and.
     &          index(header_str, 'INCREMENT= (') .le. 0 .and.
     &          index(header_str, 'ELSET=(') .le. 0 .and.
     &          index(header_str, 'FILE=(') .le. 0) then
               write(*,*) '*WARNING reading *USER BEAM OUTPUT: '
               write(*,*) '         parameter not recognized: ',
     &                    trim(token)
               call inputwarning(inpc,ipoinpc,iline,
     &              "*USER BEAM OUTPUT%")
            endif
         endif
      enddo
!
!     Allocate and populate out_elem_active mask for each target
!
      if (allocated(out_elem_active)) deallocate(out_elem_active)
      allocate(out_elem_active(max(ne, 100000), MAX_TARGETS))
      out_elem_active = .false.
!
      do itarg = 1, out_num_targets
         out_target_file_init(itarg) = .false.
         if (out_target_elset(itarg) .eq. ' ') then
            do elem = 1, max(ne, 1)
               if (elem .le. size(out_elem_active, 1)) then
                  out_elem_active(elem, itarg) = .true.
               endif
            enddo
         else
            elset = out_target_elset(itarg)
            ipos = index(elset, ' ')
            if (ipos .le. 1) ipos = len_trim(elset) + 1
            elset(ipos:ipos) = 'E'
            call cident81(set, elset, nset, id_set)
            if (id_set .le. 0 .or. elset .ne. set(id_set)) then
               if (out_target_elset(itarg)(1:4) .eq. 'ALL ' .or.
     &             out_target_elset(itarg)(1:5) .eq. 'EALL ' .or.
     &             out_target_elset(itarg)(1:1) .eq. '*') then
                  do elem = 1, max(ne, 1)
                     if (elem .le. size(out_elem_active, 1)) then
                        out_elem_active(elem, itarg) = .true.
                     endif
                  enddo
               else
                  elset(ipos:ipos) = ' '
                  write(*,*) '*WARNING reading *USER BEAM OUTPUT: '
                  write(*,*) '         element set ', trim(elset),
     &                       ' has not yet been defined.'
                  call inputwarning(inpc,ipoinpc,iline,
     &                 "*USER BEAM OUTPUT%")
               endif
            else
               do j = istartset(id_set), iendset(id_set)
                  if (ialset(j) .gt. 0) then
                     elem = ialset(j)
                     if (elem .le. size(out_elem_active, 1)) then
                        out_elem_active(elem, itarg) = .true.
                     endif
                  else
                     k = ialset(j-2)
                     do
                        k = k - ialset(j)
                        if (k .ge. ialset(j-1)) exit
                        elem = k
                        if (elem .le. size(out_elem_active, 1)) then
                           out_elem_active(elem, itarg) = .true.
                        endif
                     enddo
                  endif
               enddo
            endif
         endif
      enddo
!
!     Read data lines for output variable keys
!
      out_flag_f = .false.
      out_flag_u = .false.
      out_flag_s = .false.
      out_flag_pts = .false.
      out_flag_q = .false.
      has_keys = .false.
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if ((istat .lt. 0) .or. (key .eq. 1)) exit
         if (n .le. 0) cycle
!
         do i = 1, n
            token = adjustl(textpart(i)(1:80))
            len_t = len_trim(token)
            if (len_t .eq. 0) cycle
            do m = 1, len_t
               if (token(m:m) .ge. 'a' .and. token(m:m) .le. 'z') then
                  token(m:m) = char(ichar(token(m:m)) - 32)
               endif
            enddo
!
            if (token(1:1) .eq. 'F' .or.
     &          token(1:5) .eq. 'FORCE' .or.
     &          token(1:6) .eq. 'FORCES' .or.
     &          token(1:4) .eq. 'BFOR') then
               out_flag_f = .true.
               has_keys = .true.
            elseif (token(1:1) .eq. 'U' .or.
     &              token(1:4) .eq. 'DISP' .or.
     &              token(1:12) .eq. 'DISPLACEMENT' .or.
     &              token(1:13) .eq. 'DISPLACEMENTS') then
               out_flag_u = .true.
               has_keys = .true.
            elseif (token(1:14) .eq. 'SECTION_POINTS' .or.
     &              token(1:13) .eq. 'SECTIONPOINTS' .or.
     &              token(1:12) .eq. 'SECTION_POINT' .or.
     &              token(1:6) .eq. 'POINTS' .or.
     &              token(1:5) .eq. 'POINT' .or.
     &              token(1:7) .eq. 'CORNERS' .or.
     &              token(1:6) .eq. 'CORNER' .or.
     &              token(1:13) .eq. 'STRESS_POINTS' .or.
     &              token(1:12) .eq. 'STRESS_POINT' .or.
     &              token(1:3) .eq. 'PTS' .or.
     &              token(1:4) .eq. 'SPTS') then
               out_flag_pts = .true.
               has_keys = .true.
            elseif (token(1:1) .eq. 'S' .or.
     &              token(1:6) .eq. 'STRESS' .or.
     &              token(1:8) .eq. 'STRESSES') then
               out_flag_s = .true.
               has_keys = .true.
            elseif (token(1:1) .eq. 'Q' .or.
     &              token(1:4) .eq. 'LOAD' .or.
     &              token(1:5) .eq. 'LOADS' .or.
     &              token(1:5) .eq. 'DLOAD') then
               out_flag_q = .true.
               has_keys = .true.
            elseif (token(1:3) .eq. 'ALL' .or.
     &              token(1:1) .eq. '*') then
               out_flag_f = .true.
               out_flag_u = .true.
               out_flag_s = .true.
               out_flag_q = .true.
               has_keys = .true.
            else
               write(*,*) '*WARNING reading *USER BEAM OUTPUT: '
               write(*,*) '         unrecognized output variable: ',
     &                    trim(token)
               call inputwarning(inpc,ipoinpc,iline,
     &              "*USER BEAM OUTPUT%")
            endif
         enddo
      enddo
!
!     Default to Forces + Stresses (matching standard beam output)
!
      if (.not. has_keys) then
         out_flag_f = .true.
         out_flag_s = .true.
      endif
!
      return
      end subroutine userbeamoutputs
!
! ======================================================================
!     Helper routine to strip quotes from a string
! ======================================================================
      subroutine strip_quotes(str)
      implicit none
      character*(*) str
      character*80 temp
      integer i, len_s, j
!
      temp = adjustl(str)
      len_s = len_trim(temp)
      j = 0
      str = ' '
      do i = 1, len_s
         if (temp(i:i) .ne. '''' .and. temp(i:i) .ne. '"') then
            j = j + 1
            str(j:j) = temp(i:i)
         endif
      enddo
      end subroutine strip_quotes
!
! ======================================================================
!     Helper routine to parse INCREMENT parameter
! ======================================================================
      subroutine parse_increment_param(header_str, param_val,
     &     inc_mode, inc_freq, inc_list, ninc_list, inpc, ipoinpc,
     &     iline)
      implicit none
      character*(*) header_str, param_val
      character*4 inc_mode
      character*1 inpc(*)
      integer inc_freq, inc_list(50), ninc_list, ipoinpc(0:*), iline
      character*1000 s_hdr
      character*80 s_val, sub_s
      integer p1, p2, istat_k, ival, m, len_v, idx, k
!
      s_val = adjustl(param_val)
      len_v = len_trim(s_val)
      do m = 1, len_v
         if (s_val(m:m) .ge. 'a' .and. s_val(m:m) .le. 'z') then
            s_val(m:m) = char(ichar(s_val(m:m)) - 32)
         endif
      enddo
!
!     Check for list format: INCREMENT=(...)
!
      s_hdr = header_str
      p1 = index(s_hdr, 'INCREMENT=(')
      if (p1 .le. 0) p1 = index(s_hdr, 'INCREMENT= (')
      if (p1 .gt. 0) then
         p1 = index(s_hdr(p1:), '(') + p1 - 1
         p2 = index(s_hdr(p1:), ')') + p1 - 1
         if (p2 .gt. p1) then
            sub_s = s_hdr(p1+1 : p2-1)
            inc_mode = 'LIST'
            ninc_list = 0
            do m = 1, len_trim(sub_s)
               if (sub_s(m:m) .eq. ',') sub_s(m:m) = ' '
            enddo
            do k = 1, 50
               sub_s = adjustl(sub_s)
               if (len_trim(sub_s) .eq. 0) exit
               read(sub_s, *, iostat=istat_k) ival
               if (istat_k .eq. 0) then
                  ninc_list = ninc_list + 1
                  inc_list(ninc_list) = ival
                  idx = index(sub_s, ' ')
                  if (idx .gt. 0) then
                     sub_s = sub_s(idx+1:)
                  else
                     sub_s = ' '
                  endif
               else
                  exit
               endif
            enddo
            return
         endif
      endif
!
!     Check standard keywords or integer
!
      if (s_val(1:4) .eq. 'LAST') then
         inc_mode = 'LAST'
         inc_freq = 1
      elseif (s_val(1:3) .eq. 'ALL') then
         inc_mode = 'ALL '
         inc_freq = 1
      else
         read(s_val, *, iostat=istat_k) ival
         if (istat_k .eq. 0 .and. ival .gt. 0) then
            inc_mode = 'FREQ'
            inc_freq = ival
         else
            write(*,*) '*WARNING reading *USER BEAM OUTPUT: '
            write(*,*) '         invalid INCREMENT option; '//
     &                 'defaulting to LAST.'
            call inputwarning(inpc,ipoinpc,iline,
     &           "*USER BEAM OUTPUT%")
            inc_mode = 'LAST'
            inc_freq = 1
         endif
      endif
      end subroutine parse_increment_param
!
! ======================================================================
!     Helper routine to extract tuple or single parameter value
! ======================================================================
      subroutine get_param_tuple_or_val(header_str, param_name,
     &     items, n_items, max_items)
      implicit none
      character*(*) header_str, param_name
      character*80 items(*)
      integer n_items, max_items
      integer p_len, h_len, i_pos, i, j, start_idx, end_idx
      character*1000 h_upper, p_upper
      character*500 val_str
!
      n_items = 0
      h_len = len_trim(header_str)
      p_len = len_trim(param_name)
      if (h_len .eq. 0 .or. p_len .eq. 0) return
!
      h_upper = header_str
      do i = 1, h_len
         if (h_upper(i:i) .ge. 'a' .and. h_upper(i:i) .le. 'z') then
            h_upper(i:i) = char(ichar(h_upper(i:i)) - 32)
         endif
      enddo
      p_upper = param_name
      do i = 1, p_len
         if (p_upper(i:i) .ge. 'a' .and. p_upper(i:i) .le. 'z') then
            p_upper(i:i) = char(ichar(p_upper(i:i)) - 32)
         endif
      enddo
!
      i_pos = index(h_upper(1:h_len), p_upper(1:p_len))
      if (i_pos .le. 0) return
!
      start_idx = i_pos + p_len
      do while (start_idx .le. h_len .and.
     &          header_str(start_idx:start_idx) .eq. ' ')
         start_idx = start_idx + 1
      enddo
      if (start_idx .gt. h_len) return
!
      if (header_str(start_idx:start_idx) .eq. '(') then
         start_idx = start_idx + 1
         end_idx = index(header_str(start_idx:h_len), ')')
         if (end_idx .gt. 0) then
            end_idx = start_idx + end_idx - 2
         else
            end_idx = h_len
         endif
         val_str = header_str(start_idx:end_idx)
      else
         end_idx = start_idx
         do while (end_idx .le. h_len .and.
     &             header_str(end_idx:end_idx) .ne. ',')
            end_idx = end_idx + 1
         enddo
         end_idx = end_idx - 1
         val_str = header_str(start_idx:end_idx)
      endif
!
      i = 1
      do while (i .le. len_trim(val_str) .and. n_items .lt. max_items)
         do while (i .le. len_trim(val_str) .and.
     &             val_str(i:i) .eq. ' ')
            i = i + 1
         enddo
         if (i .gt. len_trim(val_str)) exit
         j = i
         do while (j .le. len_trim(val_str) .and.
     &             val_str(j:j) .ne. ',')
            j = j + 1
         enddo
         n_items = n_items + 1
         items(n_items) = adjustl(val_str(i:j-1))
         call strip_quotes(items(n_items))
         i = j + 1
      enddo
      end subroutine get_param_tuple_or_val
