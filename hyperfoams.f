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
      subroutine hyperfoams(inpc,textpart,elcon,nelcon,
     &  nmat,ntmat_,ncmat_,irstrt,istep,istat,n,iperturb,iline,ipol,
     &  inl,ipoinp,inp,ipoinpc,ier)
!
!     reading the input deck: *HYPERFOAM
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nelcon(2,*),nmat,ntmat,ntmat_,istep,istat,ipoinpc(0:*),
     &  n,key,i,ityp,iperturb(*),iend,ncmat_,irstrt(*),iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ier
!
      real*8 elcon(0:ncmat_,ntmat_,*)
!
      ntmat=0
      iperturb(1)=3
      iperturb(2)=1
      write(*,*) '*INFO reading *HYPERFOAM: nonlinear geometric'
      write(*,*) '      effects are turned on'
      write(*,*)
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *HYPERFOAM: *HYPERFOAM should be'
         write(*,*) '  placed before all step definitions'
         ier=1
         return
      endif
!
      if(nmat.eq.0) then
         write(*,*) '*ERROR reading *HYPERFOAM: *HYPERFOAM should be'
         write(*,*) '  preceded by a *MATERIAL card'
         ier=1
         return
      endif
!
      ityp=-15
!
      do i=2,n
         if(textpart(i)(1:2).eq.'N=') then
            if(textpart(i)(3:3).eq.'1') then
            elseif(textpart(i)(3:3).eq.'2') then
               ityp=-16
            elseif(textpart(i)(3:3).eq.'3') then
               ityp=-17
            else
               write(*,*) '*WARNING reading *HYPERFOAM: only N=1, N=2, o 
     &r N=3 are allowed; '
               call inputerror(inpc,ipoinpc,iline,
     &              "*HYPERFOAM%",ier)
               return
            endif
         else
            write(*,*) 
     &        '*WARNING reading *HYPERFOAM: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*HYPERFOAM%")
         endif
      enddo
!
      nelcon(1,nmat)=ityp
!
      if(ityp.ne.-17) then
         if(ityp.eq.-15) then
            iend=3
         elseif(ityp.eq.-16) then
            iend=6
         endif
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            ntmat=ntmat+1
            nelcon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*) '*ERROR reading *HYPERFOAM: increase ntmat_'
               ier=1
               return
            endif
            do i=1,iend
               read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &                    elcon(i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*HYPERFOAM%",ier)
                  return
               endif
            enddo
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) 
     &              elcon(0,ntmat,nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*HYPERFOAM%",ier)
               return
            endif
         enddo
      else
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            ntmat=ntmat+1
            nelcon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*) '*ERROR reading *HYPERFOAM: increase ntmat_'
               ier=1
               return
            endif
            do i=1,8
               read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &                       elcon(i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*HYPERFOAM%",ier)
                  return
               endif
            enddo
!
            iend=1
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
               write(*,*) 
     &           '*ERROR reading *HYPERFOAM: orthotropic definition'
               write(*,*) '  is not complete. '
               call inputerror(inpc,ipoinpc,iline,
     &              "*HYPERFOAM%",ier)
               return
            endif
            do i=1,iend
               read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &                 elcon(8+i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*HYPERFOAM%",ier)
                  return
               endif
            enddo
            read(textpart(2)(1:20),'(f20.0)',iostat=istat) 
     &                  elcon(0,ntmat,nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*HYPERFOAM%",ier)
               return
            endif
         enddo
      endif
!
      return
      end

