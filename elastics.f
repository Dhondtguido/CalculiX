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
      subroutine elastics(inpc,textpart,elcon,nelcon,
     &  nmat,ntmat_,ncmat_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,ipoinpc,ier)
!
!     reading the input deck: *ELASTIC
!
      implicit none
!
      logical engineering
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nelcon(2,*),nmat,ntmat,ntmat_,istep,istat,ipoinpc(0:*),
     &  n,key,i,ityp,ncmat_,irstrt(*),iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ier
!
      real*8 elcon(0:ncmat_,ntmat_,*),e1,e2,e3,un12,un21,un13,un31,
     &  un23,un32,gam
!
      ntmat=0
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *ELASTIC: *ELASTIC should be placed'
         write(*,*) '       before all step definitions'
         ier=1
         return
      endif
!
      if(nmat.eq.0) then
         write(*,*) 
     &       '*ERROR reading *ELASTIC: *ELASTIC should be preceded'
         write(*,*) '  by a *MATERIAL card'
         ier=1
         return
      endif
!
      ityp=2
!
      do i=2,n
         if(textpart(i)(1:5).eq.'TYPE=') then
            if(textpart(i)(6:8).eq.'ISO') then
               ityp=2
            elseif(textpart(i)(6:10).eq.'ORTHO') then
               ityp=9
               engineering=.false.
            elseif(textpart(i)(6:25).eq.'ENGINEERINGCONSTANTS') then
               ityp=9
               engineering=.true.
            elseif(textpart(i)(6:10).eq.'ANISO') then
               ityp=21
            endif
            exit
         else
            write(*,*) 
     &        '*WARNING reading *ELASTIC: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*ELASTIC%")
         endif
      enddo
!
      nelcon(1,nmat)=ityp
!
      if(ityp.eq.2) then
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            ntmat=ntmat+1
            nelcon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*) '*ERROR reading *ELASTIC: increase ntmat_'
               ier=1
               return
            endif
            if(n.lt.2) then
               write(*,*) '*ERROR reading *ELASTIC: not enough'
               write(*,*) '       constants on the input line'
               call inputerror(inpc,ipoinpc,iline,
     &              "*ELASTIC%",ier)
               return
            endif
            do i=1,2
               read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &                 elcon(i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*ELASTIC%",ier)
                  return
               endif
            enddo
!
!           check on the Young's modulus and the Poisson
!           coefficient
!
            if(elcon(1,ntmat,nmat).le.0.d0) then
               write(*,*) '*ERROR reading *ELASTIC: the Young'
               write(*,*) '       modulus should exceed 0.0'
               call inputerror(inpc,ipoinpc,iline,
     &              "*ELASTIC%",ier)
               return
            endif
            if(elcon(2,ntmat,nmat).ge.0.5d0) then
               write(*,*) '*ERROR reading *ELASTIC: Poisson'
               write(*,*) '       coefficient should be less than 0.5'
               call inputerror(inpc,ipoinpc,iline,
     &              "*ELASTIC%",ier)
               return
            endif
!
            if(textpart(3)(1:1).ne.' ') then
               read(textpart(3)(1:20),'(f20.0)',iostat=istat)
     &                   elcon(0,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*ELASTIC%",ier)
                  return
               endif
            else
               elcon(0,ntmat,nmat)=0.d0
            endif
         enddo
      elseif(ityp.eq.9) then
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            ntmat=ntmat+1
            nelcon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*) '*ERROR reading *ELASTIC: increase ntmat_'
               ier=1
               return
            endif
            if(n.lt.8) then
               write(*,*) '*ERROR reading *ELASTIC: not enough'
               write(*,*) '       constants on the input line'
               call inputerror(inpc,ipoinpc,iline,
     &              "*ELASTIC%",ier)
               return
            endif
            do i=1,8
               read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &                 elcon(i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*ELASTIC%",ier)
                  return
               endif
            enddo
!
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
               write(*,*) 
     &           '*ERROR reading *ELASTIC: orthotropic definition'
               write(*,*) '  is not complete. '
               call inputerror(inpc,ipoinpc,iline,
     &              "*ELASTIC%",ier)
               return
            endif
            do i=1,1
               read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &                 elcon(8+i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*ELASTIC%",ier)
                  return
               endif
            enddo
            if(textpart(2)(1:1).ne.' ') then
               read(textpart(2)(1:20),'(f20.0)',iostat=istat)
     &                       elcon(0,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*ELASTIC%",ier)
                  return
               endif
            else
               elcon(0,ntmat,nmat)=0.d0
            endif
            if(engineering) then
               e1=elcon(1,ntmat,nmat)
               e2=elcon(2,ntmat,nmat)
               e3=elcon(3,ntmat,nmat)
               un12=elcon(4,ntmat,nmat)
               un13=elcon(5,ntmat,nmat)
               un23=elcon(6,ntmat,nmat)
               un21=un12*e2/e1
               un31=un13*e3/e1
               un32=un23*e3/e2
               gam=1.d0/(1.d0-un12*un21-un23*un32-un31*un13
     &               -2.d0*un21*un32*un13)
               elcon(1,ntmat,nmat)=e1*(1.d0-un23*un32)*gam
               elcon(2,ntmat,nmat)=e1*(un21+un31*un23)*gam
               elcon(3,ntmat,nmat)=e2*(1.d0-un13*un31)*gam
               elcon(4,ntmat,nmat)=e1*(un31+un21*un32)*gam
               elcon(5,ntmat,nmat)=e2*(un32+un12*un31)*gam
               elcon(6,ntmat,nmat)=e3*(1.d0-un12*un21)*gam
            endif
         enddo
      elseif(ityp.eq.21) then
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            ntmat=ntmat+1
            nelcon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*) '*ERROR reading *ELASTIC: increase ntmat_'
               ier=1
               return
            endif
            if(n.lt.8) then
               write(*,*) '*ERROR reading *ELASTIC: not enough'
               write(*,*) '       constants on the input line'
               call inputerror(inpc,ipoinpc,iline,
     &              "*ELASTIC%",ier)
               return
            endif
            do i=1,8
               read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &                   elcon(i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*ELASTIC%",ier)
                  return
               endif
            enddo
!
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
               write(*,*) 
     &            '*ERROR reading *ELASTIC: anisotropic definition'
               write(*,*) '  is not complete. '
               call inputerror(inpc,ipoinpc,iline,
     &              "*ELASTIC%",ier)
               return
            endif
            if(n.lt.2) then
               write(*,*) '*ERROR reading *ELASTIC: not enough'
               write(*,*) '       constants on the input line'
               call inputerror(inpc,ipoinpc,iline,
     &              "*ELASTIC%",ier)
               return
            endif
            do i=1,8
               read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &            elcon(8+i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*ELASTIC%",ier)
                  return
               endif
            enddo
!
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
               write(*,*) 
     &           '*ERROR reading *ELASTIC: anisotropic definition'
               write(*,*) '  is not complete. '
               call inputerror(inpc,ipoinpc,iline,
     &              "*ELASTIC%",ier)
               return
            endif
            if(n.lt.5) then
               write(*,*) '*ERROR reading *ELASTIC: not enough'
               write(*,*) '       constants on the input line'
               call inputerror(inpc,ipoinpc,iline,
     &              "*ELASTIC%",ier)
               return
            endif
            do i=1,5
               read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &                  elcon(16+i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*ELASTIC%",ier)
                  return
               endif
            enddo
            if(textpart(6)(1:1).ne.' ') then
               read(textpart(6)(1:20),'(f20.0)',iostat=istat)
     &                  elcon(0,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*ELASTIC%",ier)
                  return
               endif
            else
               elcon(0,ntmat,nmat)=0.d0
            endif
         enddo
      endif
!
      if(ntmat.eq.0) nelcon(1,nmat)=0
!
      return
      end

