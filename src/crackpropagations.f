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
      subroutine crackpropagations(inpc,textpart,nmethod,iperturb,
     &  isolver,istep,
     &  istat,n,tinc,tper,tmin,tmax,idrct,iline,ipol,inl,ipoinp,inp,
     &  ithermal,cs,ics,tieset,istartset,
     &  iendset,ialset,ipompc,nodempc,coefmpc,nmpc,nmpc_,ikmpc,
     &  ilmpc,mpcfree,mcs,set,nset,labmpc,ipoinpc,iexpl,nef,ttime,
     &  iaxial,nelcon,nmat,ier,jobnamec,matname)
!
!     reading the input deck: *CRACKPROPAGATION
!
!     isolver=0: SPOOLES
!             2: iterative solver with diagonal scaling
!             3: iterative solver with Cholesky preconditioning
!             4: sgi solver
!             5: TAUCS
!             7: pardiso
!             8: pastix
!
      implicit none
!
      logical timereset,input
!
      character*1 inpc(*)
      character*20 labmpc(*),solver
      character*80 material,matname(*)
      character*81 set(*),tieset(3,*)
      character*132 textpart(16),jobnamec(*),fnfrd
!
      integer nmethod,iperturb(*),isolver,istep,istat,n,key,i,idrct,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),ithermal(*),ics(*),iexpl,
     &  istartset(*),iendset(*),ialset(*),ipompc(*),nodempc(3,*),
     &  nmpc,nmpc_,ikmpc(*),ilmpc(*),mpcfree,nset,mcs,ipoinpc(0:*),
     &  nef,iaxial,nelcon(2,*),nmat,ier,j,k,l,ilen
!
      real*8 tinc,tper,tmin,tmax,cs(17,*),coefmpc(*),ttime
!
      idrct=0
      tinc=0.d0
      tper=0.d0
      tmin=0.d0
      tmax=3.5d0
      timereset=.false.
      input=.false.
!
      if(istep.ne.1) then
        write(*,*) '*WARNING reading *CRACK PROPAGATION:'
        write(*,*) '         this is step ',istep,'.'
        write(*,*) '         The frd-file of previous steps'
        write(*,*) '         will be deleted.'
!     
!       delete the .frd file (it is reopened in frd.c)
!
        ilen=index(jobnamec(1),char(0))-1
        fnfrd=jobnamec(1)(1:ilen)//'.frd'
        ilen=ilen+4
        open(7,file=fnfrd(1:ilen),status='unknown',err=71)
        close(7,status='delete',err=72)
!        
c         write(*,*) '*ERROR reading *CRACK PROPAGATION:'
c         write(*,*) '       *CRACK PROPAGATION can only be used'
c         write(*,*) '       within the first STEP'
c         ier=1
c         return
      endif
!
      do i=2,n
         if(textpart(i)(1:9).eq.'MATERIAL=') then
            material=textpart(i)(10:89)
          elseif(textpart(i)(1:7).eq.'LENGTH=') then
            if(textpart(i)(8:17).eq.'CUMULATIVE') then
              tmax=1.5d0
            elseif(textpart(i)(8:19).eq.'INTERSECTION') then
              tmax=2.5d0
            else
              write(*,*)
     &             '*ERROR reading *CRACK PROPAGATION: nonexistent'
              write(*,*) '       crack length determination method'
              write(*,*) '  '
              call inputerror(inpc,ipoinpc,iline,
     &             "*CRACK PROPAGATION%",ier)
            endif
          elseif(textpart(i)(1:6).eq.'INPUT=') then
            input=.true.
            jobnamec(4)(1:126)=textpart(i)(7:132)
            jobnamec(4)(127:132)='      '
            loop1: do j=1,126
               if(jobnamec(4)(j:j).eq.'"') then
                  do k=j+1,126
                     if(jobnamec(4)(k:k).eq.'"') then
                        do l=k-1,126
                           jobnamec(4)(l:l)=' '
                           exit loop1
                        enddo
                     endif
                     jobnamec(4)(k-1:k-1)=jobnamec(4)(k:k)
                  enddo
                  jobnamec(4)(126:126)=' '
               endif
            enddo loop1
         else
            write(*,*) 
     & '*WARNING reading *CRACK PROPAGATION: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*CRACKPROPAGATION%")
         endif
      enddo
!     
!     check for the INPUT parameter
!
      if(.not.input) then
        write(*,*) 
     &     '*ERROR reading *CRACK PROPAGATION: no input file specified:'
        write(*,*) '         ',
     &       textpart(i)(1:index(textpart(i),' ')-1)
        call inputerror(inpc,ipoinpc,iline,
     &       "*CRACK PROPAGATION%",ier)
      endif
!
!     check for the existence of the material
!
      do i=1,nmat
         if(matname(i).eq.material) exit
      enddo
      if(i.gt.nmat) then
         write(*,*) 
     &      '*ERROR reading *CRACK PROPAGATION: nonexistent material'
         write(*,*) '  '
         call inputerror(inpc,ipoinpc,iline,
     &        "*CRACK PROPAGATION%",ier)
         return
       endif
!
!     material name is stored in tmin
!
      tmin=i+0.5d0
!
      nmethod=15
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
            write(*,*) '*ERROR reading *CRACK PROPAGATION:'
            write(*,*)
     &           '         a crack propagation analysis is requested'
            write(*,*)
     &           '         but no maximum crack increment is specified'
            ier=1
            return
      endif
!
!     tinc: maximum crack increment
!      
      if(n.gt.0) then
        read(textpart(1)(1:20),'(f20.0)',iostat=istat) tinc
        if(istat.gt.0) then
          call inputerror(inpc,ipoinpc,iline,
     &         "*CRACK PROPAGATION%",ier)
          return
        endif
      endif
!
!     default: max. increment = min(a/5,rcur/5)
!
      if(tinc.le.0.d0) then
        tinc=1.d30
      endif
!
!     tper: maximum deflection angle (in degrees)
!      
      if(n.gt.1) then
        read(textpart(2)(1:20),'(f20.0)',iostat=istat) tper
        if(istat.gt.0) then
          call inputerror(inpc,ipoinpc,iline,
     &         "*CRACK PROPAGATION%",ier)
          return
        endif
      endif
!
!     default: no maximum deflection angle
!
      if(tper.le.0.d0) then
        tper=90.d0
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
 71   write(*,*) '*ERROR in openfile: could not open file ',fnfrd(1:i+4)
      call exit(201)
 72   write(*,*) '*ERROR in openfile: could not delete file ',
     &  fnfrd(1:i+4)
      call exit(201)
      end

