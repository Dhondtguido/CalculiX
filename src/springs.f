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
      subroutine springs(inpc,textpart,nelcon,nmat,ntmat_,npmat_,
     &        plicon,nplicon,
     &        ncmat_,elcon,matname,irstrt,istep,istat,n,iline,ipol,
     &        inl,ipoinp,inp,nmat_,set,istartset,iendset,ialset,
     &        nset,ielmat,ielorien,ipoinpc,mi,norien,orname,ier)
!
!     reading the input deck: *SPRING
!
      implicit none
!
      logical linear
!
      character*1 inpc(*)
      character*80 matname(*),orientation,orname(*)
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer mi(*),nelcon(2,*),nmat,ntmat_,ntmat,npmat_,npmat,istep,
     &  n,key,i,nplicon(0:ntmat_,*),ncmat_,istat,istartset(*),id,
     &  iendset(*),irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),nmat_,
     &  ialset(*),ipos,nset,j,k,ielmat(mi(3),*),ielorien(mi(3),*),
     &  ipoinpc(0:*),idof,iorientation,norien,idof2,ier  
!
      real*8 plicon(0:2*npmat_,ntmat_,*),temperature,
     &  elcon(0:ncmat_,ntmat_,*)
!
      linear=.true.
!
      ntmat=0
      npmat=0
!
      orientation='
     &                           '
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *SPRING: *SPRING should be placed'
         write(*,*) '  before all step definitions'
         ier=1
         return
      endif
!
      nmat=nmat+1
      if(nmat.gt.nmat_) then
         write(*,*) '*ERROR reading *SPRING: increase nmat_'
         ier=1
         return
      endif
      matname(nmat)(1:6)='SPRING'
      do i=7,80
         matname(nmat)(i:i)=' '
      enddo
!
      do i=2,n
         if(textpart(i)(1:9).eq.'NONLINEAR') then
            linear=.false.
         elseif(textpart(i)(1:12).eq.'ORIENTATION=') then
            orientation=textpart(i)(13:92)
         elseif(textpart(i)(1:6).eq.'ELSET=') then
            elset=textpart(i)(7:86)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
         else
            write(*,*) 
     &        '*WARNING reading *SPRING: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*SPRING%")
         endif
      enddo
!
      if(orientation.eq.'                    ') then
         iorientation=0
      else
         do i=1,norien
            if(orname(i).eq.orientation) exit
         enddo
         if(i.gt.norien) then
            write(*,*)
     &       '*ERROR reading *SPRING: nonexistent orientation'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &           "*SPRING%",ier)
            return
         endif
         iorientation=i
      endif
!
      if(linear) then
         nelcon(1,nmat)=2
!
!        linear spring
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
!
!           check whether the first field in the first line 
!           underneath *SPRING contains a decimal point. If so,
!           this line is considered to be the start of material
!           data for SPRINGA elements. If not, it is considered
!           to contain degrees of freedom for SPRING1 or SPRING2 elements. 
!
            if(ntmat.eq.0) then
               idof=1
               do i=1,132
                  if(textpart(1)(i:i).eq.'.') then
                     idof=0
                     exit
                  endif
               enddo
               if(idof.eq.1) then
                  read(textpart(2)(1:10),'(i10)',iostat=istat) idof2
                  if(istat.gt.0) then
                     call inputerror(inpc,ipoinpc,iline,
     &                    "*SPRING%",ier)
                     return
                  endif
                  if(idof2.eq.0) then
                     if(ncmat_.lt.3) then
                        write(*,*) '*ERROR reading *SPRING: one degree'
                        write(*,*) '       of freedom was specified'
                        write(*,*) '       (no decimal point in entry),'
                        write(*,*) '       however, there are no'
                        write(*,*) '       SPRING1 elements'
                        write(*,*) '       in the input deck'
                        call inputerror(inpc,ipoinpc,iline,
     &                                              "*SPRING%",ier)
                        return
                     endif
                     read(textpart(1)(1:20),'(f20.0)',iostat=istat) 
     &                    elcon(3,1,nmat)
                  else
                     if(ncmat_.lt.4) then
                        write(*,*) '*ERROR reading *SPRING: two degrees'
                        write(*,*) '       of freedom were specified'
                        write(*,*) '       (no decimal point in entry),'
                        write(*,*) '       however, there are no'
                        write(*,*) '       SPRING2 elements'
                        write(*,*) '       in the input deck'
                        call inputerror(inpc,ipoinpc,iline,
     &                                              "*SPRING%",ier)
                        return
                     endif
                     read(textpart(1)(1:20),'(f20.0)',iostat=istat) 
     &                    elcon(3,1,nmat)
                     read(textpart(2)(1:20),'(f20.0)',iostat=istat) 
     &                    elcon(4,1,nmat)
                  endif
                  cycle
               endif
            endif
!
            ntmat=ntmat+1
            nelcon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*) '*ERROR reading *SPRING: increase ntmat_'
               ier=1
               return
            endif
            do i=1,2
               read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &                 elcon(i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*SPRING%",ier)
                  return
               endif
            enddo
            if(textpart(3)(1:1).ne.' ') then
               read(textpart(3)(1:20),'(f20.0)',iostat=istat)
     &                   elcon(0,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*SPRING%",ier)
                  return
               endif
            else
               elcon(0,ntmat,nmat)=0.d0
            endif
         enddo
      else
         nelcon(1,nmat)=-51
!
!        nonlinear spring behavior
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
!
!           check whether the first field in the first line 
!           underneath *SPRING contains a decimal point. If so,
!           this line is considered to be the start of material
!           data for SPRINGA elements. If not, it is considered
!           to contain degrees of freedom for SPRING1 or SPRING2 elements. 
!
            if(ntmat.eq.0) then
               idof=1
               do i=1,132
                  if(textpart(1)(i:i).eq.'.') then
                     idof=0
                     exit
                  endif
               enddo
               if(idof.eq.1) then
                  if(ncmat_.lt.4) then
                     write(*,*) '*ERROR reading *SPRING: a degree'
                     write(*,*) '       of freedom was specified'
                     write(*,*) '       (no decimal point in entry),'
                     write(*,*) '       however, there are neither'
                     write(*,*) '       SPRING1 nor SPRING2 elements'
                     write(*,*) '       in the input deck'
                     call inputerror(inpc,ipoinpc,iline,
     &                    "*SPRING%",ier)
                     return
                  endif
                  read(textpart(1)(1:20),'(f20.0)',iostat=istat) 
     &                 elcon(3,1,nmat)
                  read(textpart(2)(1:20),'(f20.0)',iostat=istat) 
     &                 elcon(4,1,nmat)
                  cycle
               endif
            endif
!
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) temperature
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*SPRING%",ier)
               return
            endif
!
!           first temperature
!
            if(ntmat.eq.0) then
               npmat=0
               ntmat=ntmat+1
               if(ntmat.gt.ntmat_) then
                  write(*,*) '*ERROR reading *SPRING: increase ntmat_'
                  ier=1
                  return
               endif
               nplicon(0,nmat)=ntmat
               plicon(0,ntmat,nmat)=temperature
!
!           new temperature
!
            elseif(plicon(0,ntmat,nmat).ne.temperature) then
               npmat=0
               ntmat=ntmat+1
               if(ntmat.gt.ntmat_) then
                  write(*,*) '*ERROR reading *SPRING: increase ntmat_'
                  ier=1
                  return
               endif
               nplicon(0,nmat)=ntmat
               plicon(0,ntmat,nmat)=temperature
            endif
            do i=1,2
               read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &              plicon(2*npmat+i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*SPRING%",ier)
                  return
               endif
            enddo
            npmat=npmat+1
            if(npmat.gt.npmat_) then
               write(*,*) '*ERROR reading *SPRING: increase npmat_'
               ier=1
               return
            endif
            nplicon(ntmat,nmat)=npmat
         enddo
      endif
!
      if(ntmat.eq.0) then
         write(*,*) '*ERROR reading *SPRING: *SPRING card without data'
         ier=1
         return
      endif
c      do i=1,nset
c         if(set(i).eq.elset) exit
c      enddo
      call cident81(set,elset,nset,id)
      i=nset+1
      if(id.gt.0) then
        if(elset.eq.set(id)) then
          i=id
        endif
      endif
      if(i.gt.nset) then
         elset(ipos:ipos)=' '
         write(*,*) '*ERROR reading *SPRING: element set ',elset
         write(*,*) '       has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &        "*SPRING%",ier)
         return
      endif
!
!     assigning the elements of the set the appropriate material
!
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            ielmat(1,ialset(j))=nmat
            ielorien(1,ialset(j))=iorientation
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               ielmat(1,k)=nmat
               ielorien(1,k)=iorientation
            enddo
         endif
      enddo
!
      return
      end

