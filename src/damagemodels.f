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
      subroutine damagemodels(inpc,textpart,matname,nmat,nmat_,
     &  irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &  imat,ier,elcon,nelcon,ntmat_,ncmat_)
!
!     reading the input deck: *DAMAGE MODEL
!
      implicit none
!
      character*1 inpc(*)
      character*80 matname(*)
      character*132 textpart(16)
!
      integer nmat,nmat_,istep,istat,n,key,i,irstrt(*),iline,ipol,inl,
     &     ipoinp(2,*),inp(3,*),ipoinpc(0:*),imat,ier,nelcon(2,*),
     &     ntmat_,ncmat_,j,nconstants,isum,imax,ntmat
!
      real*8 type,elcon(0:ncmat_,ntmat_,*)
!
      ntmat=0
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) 
     &    '*ERROR reading *DAMAGE MODEL: *DAMAGE MODEL should be placed'
         write(*,*) '  before all step definitions'
         ier=1
         return
      endif
!
      nmat=nmat+1
      if(nmat.gt.nmat_) then
         write(*,*) '*ERROR reading *DAMAGE MODEL: increase nmat_'
         ier=1
         return
      endif
!
      imat=nmat
!
      do i=2,n
         if(textpart(i)(1:5).eq.'NAME=') then
            matname(nmat)=textpart(i)(6:85)
            if(textpart(i)(86:86).ne.' ') then
               write(*,*) 
     &             '*ERROR reading *DAMAGE MODEL: damage model name too'
               write(*,*) '       long (more than 80 characters)'
               write(*,*) '       damage model name:',textpart(i)(1:132)
               ier=1
               return
            endif
            exit
          else if(textpart(i)(1:5).eq.'TYPE=') then
            if(textpart(i)(6:15).eq.'RICETRACEY') then
              type=1.5
              nconstants=3
            elseif(textpart(i)(6:16).eq.'JOHNSONCOOK') then
              type=2.5
              nconstants=10
            else
              write(*,*) 
     &             '*ERROR reading *DAMAGE MODEL: damage model type'
              write(*,*) '       is not known.'
              write(*,*) '       damage model type:',textpart(i)(6:85)
              ier=1
              return
            endif
         else
            write(*,*) 
     &       '*WARNING reading *DAMAGE MODEL: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*DAMAGE MODEL%")
         endif
      enddo
!
!     the damage model is stored as a mechanical user material
!
      nelcon(1,nmat)=-100-nconstants
!
      do
        do j=1,(nconstants)/8+1
          if(j.eq.1) then
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            ntmat=ntmat+1
            nelcon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
              write(*,*) 
     &             '*ERROR reading *DAMAGE MODEL: increase ntmat_'
              ier=1
              return
            endif
            isum=1
            elcon(1,ntmat,nmat)=type
          else
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
              write(*,*) 
     &          '*ERROR reading *DAMAGE MODEL: anisotropic definition'
              write(*,*) '  is not complete. '
              call inputerror(inpc,ipoinpc,iline,
     &             "*DAMAGE MODEL%",ier)
              return
            endif
          endif
          imax=8
          if(isum.gt.nconstants+1) then
            imax=nconstants-isum+1
          endif
          do i=1,imax
            if(isum+i.le.nconstants) then
              read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &             elcon(isum+i,ntmat,nmat)
            else
              read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &             elcon(0,ntmat,nmat)
            endif
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*DAMAGE MODEL%",ier)
              return
            endif
          enddo
          isum=isum+imax
!     
        enddo
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

