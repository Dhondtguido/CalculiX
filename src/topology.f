!
!     reading the input deck: *TOPOLOGY, INPUT=file
!     the file holds one scaling factor per element, in the order
!     of the *ELEMENT cards
!
      subroutine topology(inpc,textpart,istep,istat,n,iline,ipol,
     &  inl,ipoinp,inp,ipoinpc,ne,ier)
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16),fname
      character*256 line
!
      integer istep,istat,n,ipoinpc(0:*),key,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),i,ier,ne,ios,nval,ielem
!
      real*8 xval
!
      real*8, dimension(:), allocatable, save :: xscale
      integer, save :: nscalestore=0
!
      if(istep.gt.0) then
         write(*,*)
     &   '*ERROR reading *TOPOLOGY: *TOPOLOGY'
         write(*,*) '  should be placed before all step definitions'
         ier=1
         return
      endif
!
      if(ne.eq.0) then
         write(*,*)
     &   '*ERROR reading *TOPOLOGY: *TOPOLOGY'
         write(*,*) '  should be placed after the *ELEMENT cards'
         ier=1
         return
      endif
!
      if(nscalestore.ne.0) then
         write(*,*)
     &   '*ERROR reading *TOPOLOGY: only one card allowed'
         ier=1
         return
      endif
!
      fname=' '
      do i=2,n
         if(textpart(i)(1:6).eq.'INPUT=') then
            fname=textpart(i)(7:132)
         else
            write(*,*)
     &   '*WARNING in topology: parameter not recognized:'
            write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &           "*TOPOLOGY%")
         endif
      enddo
!
      if(fname.eq.' ') then
         write(*,*) '*ERROR reading *TOPOLOGY: *TOPOLOGY'
         write(*,*) '  requires the parameter INPUT=file'
         ier=1
         return
      endif
!
      open(77,file=fname,status='old',iostat=ios)
      if(ios.ne.0) then
         write(*,*) '*ERROR reading *TOPOLOGY: cannot open file'
         write(*,*) '  ',trim(fname)
         ier=1
         return
      endif
!
      allocate(xscale(ne))
      nval=0
      do
         read(77,'(a)',iostat=ios) line
         if(ios.ne.0) exit
         if(len_trim(line).eq.0) cycle
         nval=nval+1
         if(nval.le.ne) then
            read(line,*,iostat=ios) xscale(nval)
            if(ios.ne.0) then
               write(*,*)
     &         '*ERROR reading *TOPOLOGY: invalid value'
               write(*,*) '  ',trim(line)
               close(77)
               ier=1
               return
            endif
         endif
      enddo
      close(77)
!
      if(nval.ne.ne) then
         write(*,*) '*ERROR reading *TOPOLOGY: file ',trim(fname)
         write(*,*) '  contains',nval,' values, but the model has',
     &        ne,' elements'
         ier=1
         return
      endif
!
      nscalestore=ne
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
         write(*,*)
     &   '*ERROR reading *TOPOLOGY: no data lines allowed'
         write(*,*) '  the scaling factors are read from ',trim(fname)
         call inputerror(inpc,ipoinpc,iline,"*TOPOLOGY%",ier)
         return
      enddo
!
      return
!
      entry topologyget(ielem,xval)
!
      if(nscalestore.eq.0) then
         xval=1.d0
      else
         xval=xscale(ielem)
      endif
!
      return
      end
