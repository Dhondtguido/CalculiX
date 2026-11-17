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
      subroutine userbeamsections(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,ielmat,matname,nmat,ielorien,orname,norien,
     &  thicke,ipkon,iponor,xnor,ixfree,
     &  offset,lakon,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &  ipoinpc,mi,ielprop,nprop,nprop_,prop,nelcon,ier)
!
!     reading the input deck: *USER BEAM SECTION
!
!     For user beam element (UB31), section data is stored
!     into the prop(*)/ielprop(*) arrays following the same CCX
!     methodology used by usersections.f.
!
!     prop layout written per element (27 slots, 0-padded):
!      1     : sect_type  (1=RECT,2=CIRC/PIPE,3=I,4=T,5=CHAN,6=L)
!      2-7   : dims(1..6)
!      8     : rot_angle = 0.0  (not settable via *USER BEAM SECTION)
!      9-11  : off_x1, off_y1, off_z1  (node 1 offsets)
!      12-14 : off_x2, off_y2, off_z2  (node 2 offsets)
!      15    : rel_1 = 0  (release code node 1)
!      16    : rel_2 = 0  (release code node 2)
!      17-19 : e2_x, e2_y, e2_z  (orientation vector)
!      20-23 : reserved
!      24-27 : LOAD CACHE — written at runtime by e_c3d_ub31 kernel:
!                24 = w1  (load intensity at node 1 / uniform value)
!                25 = w2  (load intensity at node 2)
!                26 = alpha  (patch load start fraction)
!                27 = beta   (patch load end fraction)
!     Total: 27 slots, all permanent geometry + 4 runtime cache.
!
      implicit none
!
      logical nodalthickness, native_section
!
      character*1 inpc(*)
      character*4 section
      character*8 lakon(*)
      character*80 matname(*),orname(*),material,orientation
      character*81 set(*),elset
      character*132 textpart(16)
      character*1000 header_str
      character*100 vals_str
!
      integer istartset(*),iendset(*),ialset(*),mi(*),ielmat(mi(3),*),
     &  ipoinpc(0:*),numnod,id,
     &  ielorien(mi(3),*),ipkon(*),iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),nset,nmat,norien,istep,istat,n,key,i,j,k,l,imaterial,
     &  iorientation,ipos,m,iponor(2,*),ixfree,
     &  indexx,indexe,irstrt(*),nelcon(2,*),ier,
     &  ielprop(*),nprop,nprop_,npropstart,iset,ndprop,
     &  commas, idx, idx_end
!
      real*8 thicke(mi(3),*),thickness1,thickness2,p(3),xnor(*),
     &  offset(2,*),offset1,offset2,dd,dims_temp(6),sect_type_val,
     &  original_dims(6), abaqus_i_l, geom_y,geom_z,
     &  ccx_offset1,ccx_offset2,ag,ac,as,ay,az,
     &  aa,ab,aw,zc,ctmp,theta_geom,wt1,wt2,
     &  Iyy_tmp,Izz_tmp,J_tmp,k1_tmp,k2_tmp,a_tmp,
     &  prop(*), off1(3), off2(3), off3(3), rel1, rel2, rel3,
     &  rot_angle_val, timo_val, custom_k, open_tol, flange_tol
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) 
     &       '*ERROR reading *USER BEAM SECTION: '
     &       //'*USER BEAM SECTION should'
         write(*,*) '  be placed before all step definitions'
         ier=1
         return
      endif
!
      nodalthickness=.false.
      offset1=0.d0
      offset2=0.d0
      rel1=0.d0
      rel2=0.d0
      rel3=0.d0
      rot_angle_val=0.d0
      timo_val=0.d0
      orientation='                                               '
      section='    '
      native_section=(textpart(1)(1:12).eq.'*BEAMSECTION')
      ipos=1
!
      header_str = ' '
      do i=2,n
         header_str = trim(header_str) // ',' // trim(textpart(i))
      enddo
!
      do i=2,n
         if(textpart(i)(1:9).eq.'MATERIAL=') then
            material=textpart(i)(10:89)
         elseif(textpart(i)(1:12).eq.'ORIENTATION=') then
            orientation=textpart(i)(13:92)
         elseif(textpart(i)(1:9).eq.'ROTATION=') then
            read(textpart(i)(10:29),'(f20.0)',iostat=istat)
     &           rot_angle_val
            if (istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*USER BEAM SECTION%",ier)
               return
            endif
         elseif(textpart(i)(1:6).eq.'ELSET=') then
            elset=textpart(i)(7:86)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
         elseif(textpart(i)(1:14).eq.'NODALTHICKNESS') then
            nodalthickness=.true.
         elseif(textpart(i)(1:8).eq.'SECTION=') then
            if(textpart(i)(9:12).eq.'CIRC') then
               section='CIRC'
            elseif(textpart(i)(9:12).eq.'RECT') then
               section='RECT'
            elseif(textpart(i)(9:10).eq.'I ') then
               section='I   '
            elseif(textpart(i)(9:10).eq.'T ') then
               section='T   '
            elseif(textpart(i)(9:12).eq.'PIPE') then
               section='PIPE'
            elseif(textpart(i)(9:12).eq.'CHAN') then
               section='CHAN'
            elseif(textpart(i)(9:10).eq.'L ') then
               section='L   '
            else
               section=textpart(i)(9:12)
            endif
         elseif(textpart(i)(1:9).eq.'RELEASE1=') then
            call parse_section_release_code(textpart(i)(10:), rel1)
         elseif(textpart(i)(1:9).eq.'RELEASE2=') then
            call parse_section_release_code(textpart(i)(10:), rel2)
         elseif(textpart(i)(1:9).eq.'RELEASE3=') then
            call parse_section_release_code(textpart(i)(10:), rel3)
         elseif(textpart(i)(1:6).eq.'KAPPA=') then
            if(textpart(i)(7:8).eq.'NO' .or.
     &         textpart(i)(7:9).eq.'OFF' .or.
     &         textpart(i)(7:11).eq.'FALSE') then
               timo_val = -1.d0
            else
               read(textpart(i)(7:), *, iostat=istat) custom_k
               if (istat .eq. 0) then
                  if (custom_k .le. 1.d-6) then
                     timo_val = -1.d0
                  else
                     timo_val = custom_k
                  endif
               else
                  timo_val = 0.d0
               endif
            endif
         elseif(textpart(i)(1:11).eq.'TIMOSHENKO=') then
            if(textpart(i)(12:13).eq.'NO' .or.
     &         textpart(i)(12:14).eq.'OFF' .or.
     &         textpart(i)(12:12).eq.'0' .or.
     &         textpart(i)(12:16).eq.'FALSE') then
               timo_val = -1.d0
            else
               read(textpart(i)(12:), *, iostat=istat) custom_k
               if (istat .eq. 0 .and. custom_k .gt. 1.d-6 .and.
     &             custom_k .lt. 10.d0) then
                  timo_val = custom_k
               else
                  timo_val = 0.d0
               endif
            endif
         else
            if (index(textpart(i), 'OFFSET').gt.0 .or.
     &          (index(header_str, 'OFFSET').gt.0 .and.
     &           index(textpart(i), '=').eq.0)) then
               ! Skip warning, this is part of the offset tuple
            else
               write(*,*) 
     &           '*WARNING reading *USER BEAM SECTION: '
     &           //'parameter not recognized:'
               write(*,*) '         ',
     &                     textpart(i)(1:index(textpart(i),' ')-1)
               call inputwarning(inpc,ipoinpc,iline,
     &"*USER BEAM SECTION%")
            endif
         endif
      enddo
!
!     check whether a section was defined
!
      if(section.eq.'    ') then
         write(*,*) '*ERROR in *USER BEAM SECTION: no section defined'
         ier=1
         return
      endif
!
!     check for the existence of the set, the material and orientation
!
      do i=1,nmat
         if(matname(i).eq.material) exit
      enddo
      if(i.gt.nmat) then
         write(*,*) '*ERROR in *USER BEAM SECTION: nonexistent material'
         write(*,*) '  '
         call inputerror(inpc,ipoinpc,iline,
     &         "*USER BEAM SECTION%",ier)
         return
      endif
      imaterial=i
!
      if(orientation.eq.' ') then
         iorientation=0
      else
         do i=1,norien
            if(orname(i).eq.orientation) exit
         enddo
         if(i.gt.norien) then
            write(*,*)
     &     '*ERROR in *USER BEAM SECTION: nonexistent orientation'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &            "*USER BEAM SECTION%",ier)
            return
         endif
         iorientation=i
      endif
!
      call cident81(set,elset,nset,id)
      i=nset+1
      if(id.gt.0) then
        if(elset.eq.set(id)) then
          i=id
        endif
      endif
      if(i.gt.nset) then
         elset(ipos:ipos)=' '
         write(*,*) '*ERROR reading *USER BEAM SECTION: element set ',
     &      elset(1:ipos)
         write(*,*) '  has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &         "*USER BEAM SECTION%",ier)
         return
      endif
      iset=i
!
!     Section type is finalized after the dimension line is read.
!     This permits native Mecway *BEAM SECTION cards which encode some
!     open sections using SECTION=BOX with one zero wall thickness.
!
!     --- read the first data line: up to 6 dimension values ---
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      do m = 1, 6
         dims_temp(m) = 0.d0
         if(.not.nodalthickness .and. m .le. n) then
            read(textpart(m)(1:20),'(f20.0)',iostat=istat) dims_temp(m)
            if(istat.gt.0) dims_temp(m) = 0.d0
         endif
      enddo
!
!     Native *BEAM SECTION geometry matches CCX rather than a
!     bespoke Mecway BOX arrangement. USER BEAM SECTION stays legacy.
      abaqus_i_l=0.d0
      geom_y=0.d0
      geom_z=0.d0
      if(native_section.and.section.eq.'I   ') then
!        Abaqus 2016 I shape: l,h,b1,b2,t1,t2,t3.
!        The first field l positions the local section origin above
!        the bottom flange. Convert geometry to the UB31 six-field
!        property definition (h,b_top,t_top,b_bot,t_bot,t_web).
         if(n.ne.7) then
            write(*,*) '*ERROR: native SECTION=I expects 7 values'
            ier=1
            return
         endif
         read(textpart(7)(1:20),'(f20.0)',iostat=istat) wt2
         if(istat.ne.0) then
            write(*,*) '*ERROR: I section web thickness'
            ier=1
            return
         endif
         abaqus_i_l=dims_temp(1)
         aa=dims_temp(2)
         ab=dims_temp(3)
         aw=dims_temp(4)
         wt1=dims_temp(5)
         ctmp=dims_temp(6)
         if(aa.le.0.d0.or.ab.le.0.d0.or.aw.le.0.d0.or.
     &      wt1.le.0.d0.or.ctmp.le.0.d0.or.wt2.le.0.d0.or.
     &      wt1+ctmp.ge.aa.or.wt2.ge.min(ab,aw)) then
            write(*,*) '*ERROR: invalid native I section geometry'
            ier=1
            return
         endif
         dims_temp(1)=aa
         dims_temp(2)=aw
         dims_temp(3)=ctmp
         dims_temp(4)=ab
         dims_temp(5)=wt1
         dims_temp(6)=wt2
         ag=ab*wt1+wt2*(aa-wt1-ctmp)+aw*ctmp
         geom_z=(ab*wt1*(wt1/2.d0)+
     &      wt2*(aa-wt1-ctmp)*
     &      (wt1+(aa-wt1-ctmp)/2.d0)+
     &      aw*ctmp*(aa-ctmp/2.d0))/ag-abaqus_i_l
      elseif(native_section.and.section.eq.'L   ') then
!        Abaqus L: a,b,t1(horizontal leg),t2(vertical leg).
!        Reference origin is the outside section corner.
         if(n.ne.4) then
            write(*,*) '*ERROR: native SECTION=L expects 4 values'
            ier=1
            return
         endif
         aa=dims_temp(1)
         ab=dims_temp(2)
         wt1=dims_temp(3)
         wt2=dims_temp(4)
         if(aa.le.0.d0.or.ab.le.0.d0.or.
     &      wt1.le.0.d0.or.wt2.le.0.d0.or.
     &      wt1.ge.ab.or.wt2.ge.aa) then
            write(*,*) '*ERROR: invalid native L section geometry'
            ier=1
            return
         endif
         ag=aa*wt1+wt2*(ab-wt1)
         geom_y=(aa*wt1*aa/2.d0+
     &      wt2*(ab-wt1)*wt2/2.d0)/ag
         geom_z=(aa*wt1*wt1/2.d0+
     &      wt2*(ab-wt1)*(ab+wt1)/2.d0)/ag
      endif
      if(native_section.and.section.eq.'RECT'.and.n.eq.1) then
         dims_temp(2)=dims_temp(1)
      endif
      if(native_section.and.section.eq.'BOX ') then
!        CCX: a(n1), b(n2), t1(right), t2(top), t3(left), t4(bottom)
!        UB31 internal: h, width, bottom, left, top, right.
         do m=1,6
            original_dims(m)=dims_temp(m)
         enddo
         dims_temp(1)=original_dims(2)
         dims_temp(2)=original_dims(1)
         dims_temp(3)=original_dims(6)
         dims_temp(4)=original_dims(5)
         dims_temp(5)=original_dims(4)
         dims_temp(6)=original_dims(3)
      endif
      open_tol=1.d-12*max(1.d0,dabs(dims_temp(1)),
     &                        dabs(dims_temp(2)))
      flange_tol=1.d-8*max(1.d0,dabs(dims_temp(3)),
     &                          dabs(dims_temp(5)))
      if(section.eq.'RECT') then
         sect_type_val=1.d0
      elseif(section.eq.'CIRC') then
         if(native_section) then
!           Native CIRC uses diameters of ellipse, not radius.
            if(n.eq.1) dims_temp(2)=dims_temp(1)
            dims_temp(1)=0.5d0*dims_temp(1)
            dims_temp(2)=0.5d0*dims_temp(2)
            sect_type_val=8.d0
         else
            sect_type_val=2.d0
         endif
      elseif(section.eq.'PIPE') then
         sect_type_val=2.d0
      elseif(section.eq.'I   ') then
         sect_type_val=3.d0
      elseif(section.eq.'T   ') then
         sect_type_val=4.d0
      elseif(section.eq.'CHAN') then
         sect_type_val=5.d0
      elseif(section.eq.'L   ') then
         if(native_section) then
            sect_type_val=9.d0
         else
            sect_type_val=6.d0
         endif
      elseif(section.eq.'BOX ') then
         if((.not.native_section).and.
     &      (dims_temp(3).gt.open_tol).and.
     &      (dims_temp(5).gt.open_tol).and.
     &      (dabs(dims_temp(3)-dims_temp(5)).le.flange_tol).and.
     &      (((dims_temp(4).gt.open_tol).and.
     &        (dabs(dims_temp(6)).le.open_tol)).or.
     &       ((dabs(dims_temp(4)).le.open_tol).and.
     &        (dims_temp(6).gt.open_tol)))) then
!           Preserve legacy UB31 USER BEAM SECTION channel shortcut.
            if(dims_temp(4).le.open_tol) then
               dims_temp(4)=dims_temp(6)
            endif
            dims_temp(3)=0.5d0*(dims_temp(3)+dims_temp(5))
            dims_temp(5)=0.d0
            dims_temp(6)=0.d0
            sect_type_val=5.d0
         else
            sect_type_val=7.d0
         endif
      else
         write(*,*) '*ERROR: unsupported UB31 section ',section
         ier=1
         return
      endif
!     Validate canonical section dimensions.
      if((section.eq.'RECT').or.(section.eq.'BOX ').or.
     &   (section.eq.'CIRC'.and.native_section)) then
         if(dims_temp(1).le.0.d0.or.dims_temp(2).le.0.d0) then
            write(*,*) '*ERROR: invalid beam section dimension'
            ier=1
            return
         endif
      endif
      if(section.eq.'PIPE') then
         if(dims_temp(1).le.0.d0.or.dims_temp(2).le.0.d0.or.
     &      dims_temp(2).ge.dims_temp(1)) then
            write(*,*) '*ERROR: invalid PIPE geometry'
            ier=1
            return
         endif
      endif
      if(section.eq.'BOX '.and.sect_type_val.eq.7.d0) then
         if(minval(dims_temp(3:6)).lt.0.d0.or.
     &      dims_temp(3)+dims_temp(5).ge.dims_temp(1).or.
     &      dims_temp(4)+dims_temp(6).ge.dims_temp(2).or.
     &      sum(dims_temp(3:6)).le.0.d0) then
            write(*,*) '*ERROR: invalid BOX geometry'
            ier=1
            return
         endif
      endif
!
      thickness1=dims_temp(1)
      thickness2=dims_temp(2)
!
!     --- read the second data line: orientation normal vector e2 ---
!     If omitted, native CCX section 1-axis is (0,0,-1).
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      p(1)=0.d0
      p(2)=0.d0
      p(3)=-1.d0
      if(.not.((istat.lt.0).or.(key.eq.1))) then
         indexx=-1
         read(textpart(1)(1:20),'(f20.0)',iostat=istat) p(1)
         if(istat.gt.0) p(1)=0.d0
         read(textpart(2)(1:20),'(f20.0)',iostat=istat) p(2)
         if(istat.gt.0) p(2)=0.d0
         read(textpart(3)(1:20),'(f20.0)',iostat=istat) p(3)
         if(istat.gt.0) p(3)=0.d0
         dd=dsqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3))
         if(dd.gt.1.d-10) then
            do j=1,3
               p(j)=p(j)/dd
            enddo
         else
            write(*,*) '*ERROR: explicitly zero beam normal'
            ier=1
            return
         endif
      endif
!
!     --- parse Abaqus-style offset parameters ---
!     Parses OFFSET=(y,z) (uniform local y/z offset) or
!     OFFSET1/2/3=(x,y,z) (independent per-node local x/y/z offsets)
!     Initialise offsets
      do m = 1, 3
         off1(m) = 0.d0
         off2(m) = 0.d0
         off3(m) = 0.d0
      enddo

!     1. Check uniform OFFSET=(y, z)
      idx = index(header_str, 'OFFSET=(')
      if (idx .gt. 0) then
         idx_end = index(header_str(idx+8:), ')') + idx + 8 - 1
         if (idx_end .gt. idx + 8) then
            vals_str = header_str(idx+8 : idx_end-1)
            commas = 0
            do m = 1, len_trim(vals_str)
               if (vals_str(m:m) .eq. ',') commas = commas + 1
            enddo
            if (commas .eq. 2) then
               read(vals_str, *, iostat=istat)
     &              off1(1), off1(2), off1(3)
            else if (commas .eq. 1) then
               read(vals_str, *, iostat=istat) off1(2), off1(3)
               off1(1) = 0.d0
            else
               read(vals_str, *, iostat=istat) off1(2)
               off1(1) = 0.d0
               off1(3) = 0.d0
            endif
            do m = 1, 3
               off2(m) = off1(m)
               off3(m) = off1(m)
            enddo
         endif
      endif

!     2. Check OFFSET1=(x, y, z)
      idx = index(header_str, 'OFFSET1=(')
      if (idx .gt. 0) then
         idx_end = index(header_str(idx+9:), ')') + idx + 9 - 1
         if (idx_end .gt. idx + 9) then
            vals_str = header_str(idx+9 : idx_end-1)
            commas = 0
            do m = 1, len_trim(vals_str)
               if (vals_str(m:m) .eq. ',') commas = commas + 1
            enddo
            if (commas .eq. 2) then
               read(vals_str, *, iostat=istat)
     &              off1(1), off1(2), off1(3)
            else if (commas .eq. 1) then
               read(vals_str, *, iostat=istat) off1(2), off1(3)
               off1(1) = 0.d0
            else
               read(vals_str, *, iostat=istat) off1(2)
               off1(1) = 0.d0
               off1(3) = 0.d0
            endif
         endif
      endif

!     3. Check OFFSET2=(x, y, z)
      idx = index(header_str, 'OFFSET2=(')
      if (idx .gt. 0) then
         idx_end = index(header_str(idx+9:), ')') + idx + 9 - 1
         if (idx_end .gt. idx + 9) then
            vals_str = header_str(idx+9 : idx_end-1)
            commas = 0
            do m = 1, len_trim(vals_str)
               if (vals_str(m:m) .eq. ',') commas = commas + 1
            enddo
            if (commas .eq. 2) then
               read(vals_str, *, iostat=istat)
     &              off2(1), off2(2), off2(3)
            else if (commas .eq. 1) then
               read(vals_str, *, iostat=istat) off2(2), off2(3)
               off2(1) = 0.d0
            else
               read(vals_str, *, iostat=istat) off2(2)
               off2(1) = 0.d0
               off2(3) = 0.d0
            endif
         endif
      endif

!     4. Check OFFSET3=(x, y, z)
      idx = index(header_str, 'OFFSET3=(')
      if (idx .gt. 0) then
         idx_end = index(header_str(idx+9:), ')') + idx + 9 - 1
         if (idx_end .gt. idx + 9) then
            vals_str = header_str(idx+9 : idx_end-1)
            commas = 0
            do m = 1, len_trim(vals_str)
               if (vals_str(m:m) .eq. ',') commas = commas + 1
            enddo
            if (commas .eq. 2) then
               read(vals_str, *, iostat=istat)
     &              off3(1), off3(2), off3(3)
            else if (commas .eq. 1) then
               read(vals_str, *, iostat=istat) off3(2), off3(3)
               off3(1) = 0.d0
            else
               read(vals_str, *, iostat=istat) off3(2)
               off3(1) = 0.d0
               off3(3) = 0.d0
            endif
         endif
      endif

!     Native CCX OFFSET1 and OFFSET2 are section-relative numbers.
!     They locate the reference line relative to the section axis,
!     not displacement vectors or physical lengths.
      if(native_section) then
         ccx_offset1=0.d0
         ccx_offset2=0.d0
         idx=index(header_str,'OFFSET1=')
         if(idx.gt.0) then
            read(header_str(idx+8:),*,iostat=istat) ccx_offset1
            if(istat.ne.0) then
               write(*,*) '*ERROR: invalid CCX OFFSET1'
               ier=1
               return
            endif
         endif
         idx=index(header_str,'OFFSET2=')
         if(idx.gt.0) then
            read(header_str(idx+8:),*,iostat=istat) ccx_offset2
            if(istat.ne.0) then
               write(*,*) '*ERROR: invalid CCX OFFSET2'
               ier=1
               return
            endif
         endif
         if(nodalthickness) then
            write(*,*) '*ERROR: direct beams do not support'
            write(*,*) '  *NODAL THICKNESS in *BEAM SECTION'
            ier=1
            return
         endif
!        Add the centroid-to-reference translation for I and L,
!        and the CCX thickness-normalised offsets. Positive CCX
!        OFFSET locates the reference line in the positive direction.
         geom_y=geom_y-ccx_offset1*thickness1
         geom_z=geom_z-ccx_offset2*thickness2
!
!        The stiffness kernel rotates asymmetric sections to
!        principal axes. Convert physical geometric offsets to
!        components of that principal system before storing them.
         call compute_section_properties(sect_type_val,dims_temp,
     &      1.d0,0.3d0,a_tmp,Iyy_tmp,Izz_tmp,J_tmp,
     &      k1_tmp,k2_tmp,theta_geom)
         theta_geom=(theta_geom+rot_angle_val)*
     &      4.d0*datan(1.d0)/180.d0
         ac=dcos(theta_geom)
         as=dsin(theta_geom)
         ay=geom_y*ac+geom_z*as
         az=-geom_y*as+geom_z*ac
         off1(2)=off1(2)+ay
         off1(3)=off1(3)+az
         off2(2)=off2(2)+ay
         off2(3)=off2(3)+az
      endif
!
!     Fallbacks for standard CCX arrays
      offset1 = off1(1)
      offset2 = off1(2)

!     --- write section data once into prop(*) at npropstart ---
!     A single 30-slot record covers UB31.
!
      ndprop = 30
      npropstart = max(1, nprop)
!
!     check capacity
!
      if(npropstart + ndprop .gt. nprop_) then
         write(*,*) '*ERROR in *USER BEAM SECTION: increase nprop_'
         ier=1
         return
      endif
!
!     initialise all slots to zero
!
      do m = 1, ndprop
         prop(npropstart + m) = 0.d0
      enddo
!
!     fill in the section data
!
      prop(npropstart + 1)  = sect_type_val
      prop(npropstart + 2)  = dims_temp(1)
      prop(npropstart + 3)  = dims_temp(2)
      prop(npropstart + 4)  = dims_temp(3)
      prop(npropstart + 5)  = dims_temp(4)
      prop(npropstart + 6)  = dims_temp(5)
      prop(npropstart + 7)  = dims_temp(6)
      prop(npropstart + 8)  = rot_angle_val  ! rot_angle
      prop(npropstart + 9)  = off1(1)        ! node 1 off_x
      prop(npropstart + 10) = off1(2)        ! node 1 off_y
      prop(npropstart + 11) = off1(3)        ! node 1 off_z
      prop(npropstart + 12) = off2(1)        ! node 2 off_x
      prop(npropstart + 13) = off2(2)        ! node 2 off_y
      prop(npropstart + 14) = off2(3)        ! node 2 off_z
      prop(npropstart + 15) = rel1           ! rel_1
      prop(npropstart + 16) = rel2           ! rel_2
      prop(npropstart + 17) = p(1)           ! e2_x
      prop(npropstart + 18) = p(2)           ! e2_y
      prop(npropstart + 19) = p(3)           ! e2_z
      prop(npropstart + 20) = timo_val       ! timoshenko flag (1=ON, 2=OFF)
      prop(npropstart + 21) = off3(2)        ! node 3 off_y
      prop(npropstart + 22) = off3(3)        ! node 3 off_z
      prop(npropstart + 23) = rel3           ! rel_3
!     slots 24-27: load cache, initialised to 0, written at runtime
      prop(npropstart + 24) = 0.d0           ! w1
      prop(npropstart + 25) = 0.d0           ! w2
      prop(npropstart + 26) = 0.d0           ! alpha
      prop(npropstart + 27) = 0.d0           ! beta
!
      nprop = npropstart + ndprop
!
!     --- assign material, orientation, offsets, and prop pointer
!         to every element in the set ---
!
      do j=istartset(iset),iendset(iset)
         if(ialset(j).gt.0) then
            if(lakon(ialset(j))(1:2).ne.'UB') then
               write(*,*) 
     &     '*ERROR: *USER BEAM SECTION only for UB elements.'
               write(*,*) '       Element ',ialset(j),
     &                    ' is not a user beam element.'
               ier=1
               return
            endif
            ielmat(1,ialset(j))=imaterial
            ielorien(1,ialset(j))=iorientation
            offset(1,ialset(j))=offset1
            offset(2,ialset(j))=offset2
!           assign prop pointer for user beam elements
            ielprop(ialset(j))=npropstart
!           store thicke for standard CCX beam elements
            indexe=ipkon(ialset(j))
            if ((lakon(ialset(j))(3:4) .eq. '31') .or.
     &          (lakon(ialset(j))(3:4) .eq. '21')) then
               numnod = 2
            else
               numnod = 3
            endif
            do l=1,numnod
               thicke(1,indexe+l)=thickness1
               thicke(2,indexe+l)=thickness2
            enddo
!           store xnor for UB elements (orientation normal)
            do l=1,numnod
               if(indexx.eq.-1) then
                  indexx=ixfree
                  do m=1,3
                     xnor(indexx+m)=p(m)
                  enddo
                  ixfree=ixfree+6
               endif
               iponor(1,indexe+l)=indexx
            enddo
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               if(lakon(k)(1:2).ne.'UB') then
                  write(*,*) 
     &     '*ERROR: *USER BEAM SECTION only for UB elements.'
                  write(*,*) '       Element ',k,
     &                       ' is not a user beam element.'
                  ier=1
                  return
               endif
               ielmat(1,k)=imaterial
               ielorien(1,k)=iorientation
               offset(1,k)=offset1
               offset(2,k)=offset2
               ielprop(k)=npropstart
               indexe=ipkon(k)
               if ((lakon(k)(3:4) .eq. '31') .or.
     &             (lakon(k)(3:4) .eq. '21')) then
                  numnod = 2
               else
                  numnod = 3
               endif
               do l=1,numnod
                  thicke(1,indexe+l)=thickness1
                  thicke(2,indexe+l)=thickness2
               enddo
               do l=1,numnod
                  if(indexx.eq.-1) then
                     indexx=ixfree
                     do m=1,3
                        xnor(indexx+m)=p(m)
                     enddo
                     ixfree=ixfree+6
                  endif
                  iponor(1,indexe+l)=indexx
               enddo
            enddo
         endif
      enddo
!
      if (key .ne. 1) then
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
      endif
!
      return
      end
!
! ==============================================================================
!     Helper routine to parse release code in *USER BEAM SECTION
! ==============================================================================
      subroutine parse_section_release_code(str, rel_val)
      implicit none
      character*(*) str
      real*8 rel_val
      character*80 s
      integer istat, len_s, m, ival
!
      s = adjustl(str)
      len_s = len_trim(s)
      rel_val = 0.d0
      if (len_s .eq. 0) return
!
      do m = 1, len_s
         if (s(m:m) .ge. 'a' .and. s(m:m) .le. 'z') then
            s(m:m) = char(ichar(s(m:m)) - 32)
         endif
      enddo
!
      if (s(1:len_s) .eq. 'M1' .or. s(1:len_s) .eq. 'RY' .or.
     &    s(1:len_s) .eq. 'ROTY' .or. s(1:len_s) .eq. 'MY') then
         rel_val = 16.d0
      elseif (s(1:len_s) .eq. 'M2' .or. s(1:len_s) .eq. 'RZ' .or.
     &        s(1:len_s) .eq. 'ROTZ' .or. s(1:len_s) .eq. 'MZ') then
         rel_val = 32.d0
      elseif (s(1:len_s) .eq. 'T' .or. s(1:len_s) .eq. 'RX' .or.
     &        s(1:len_s) .eq. 'ROTX' .or. s(1:len_s) .eq. 'MX' .or.
     &        s(1:len_s) .eq. 'TOR' .or.
     &        s(1:len_s) .eq. 'TORSION') then
         rel_val = 8.d0
      elseif (s(1:len_s) .eq. 'UX' .or. s(1:len_s) .eq. 'U1' .or.
     &        s(1:len_s) .eq. 'AXIAL' .or. s(1:len_s) .eq. 'N' .or.
     &        s(1:len_s) .eq. 'P0') then
         rel_val = 1.d0
      elseif (s(1:len_s) .eq. 'UY' .or. s(1:len_s) .eq. 'U2' .or.
     &        s(1:len_s) .eq. 'V1' .or. s(1:len_s) .eq. 'VY' .or.
     &        s(1:len_s) .eq. 'SHEAR1') then
         rel_val = 2.d0
      elseif (s(1:len_s) .eq. 'UZ' .or. s(1:len_s) .eq. 'U3' .or.
     &        s(1:len_s) .eq. 'V2' .or. s(1:len_s) .eq. 'VZ' .or.
     &        s(1:len_s) .eq. 'SHEAR2') then
         rel_val = 4.d0
      elseif (s(1:len_s) .eq. 'M1-M2' .or. s(1:len_s) .eq. 'M2-M1' .or.
     &        s(1:len_s) .eq. 'M1_M2' .or. s(1:len_s) .eq. 'M1+M2' .or.
     &        s(1:len_s) .eq. 'M2_M1' .or. s(1:len_s) .eq. 'M2+M1' .or.
     &        s(1:len_s) .eq. 'RY-RZ' .or. s(1:len_s) .eq. 'MY-MZ') then
         rel_val = 48.d0
      elseif (s(1:len_s) .eq. 'ALLM' .or. s(1:len_s) .eq. 'ALL_M' .or.
     &        s(1:len_s) .eq. 'ALL-M' .or. s(1:len_s) .eq. 'BALL' .or.
     &        s(1:len_s) .eq. 'SPHERICAL' .or.
     &        s(1:len_s) .eq. 'M1-M2-T' .or.
     &        s(1:len_s) .eq. 'T-M1-M2') then
         rel_val = 56.d0
      elseif (s(1:len_s) .eq. 'ALL' .or. s(1:len_s) .eq. 'FREE' .or.
     &        s(1:len_s) .eq. 'DISCONNECTED') then
         rel_val = 63.d0
      elseif (s(1:len_s) .eq. 'NONE' .or. s(1:len_s) .eq. 'FIXED' .or.
     &        s(1:len_s) .eq. 'RIGID' .or. s(1:len_s) .eq. '0') then
         rel_val = 0.d0
      else
         read(s(1:len_s), *, iostat=istat) ival
         if (istat .eq. 0 .and. ival .ge. 0 .and. ival .le. 63) then
            rel_val = dble(ival)
         else
            rel_val = 0.d0
         endif
      endif
      end subroutine parse_section_release_code