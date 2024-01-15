!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2023 Guido Dhondt
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
      subroutine mohrcoulombs(inpc,textpart,elcon,nelcon,
     &     nmat,ntmat_,ncmat_,irstrt,istep,istat,n,iperturb,iline,ipol,
     &     inl,ipoinp,inp,ipoinpc,ier,iplas,matname,nstate_)
!     
!     reading the input deck: *MOHR COULOMB
!     
      implicit none
!     
      character*1 inpc(*)
      character*80 matname(*)
      character*132 textpart(16)
!     
      integer nelcon(2,*),nmat,ntmat,ntmat_,istep,istat,ier,
     &     n,key,i,iperturb(*),iend,ncmat_,irstrt(*),iline,ipol,inl,
     &     ipoinp(2,*),inp(3,*),ipoinpc(0:*),iplas,nstate_
!     
      real*8 elcon(0:ncmat_,ntmat_,*),conversion
!     
      ntmat=0
      iperturb(1)=3
      iplas=1
!     
!     conversion factor from radians into degrees
!     
      conversion=4.d0*datan(1.d0)/180.d0
!     
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
        write(*,*) '*ERROR reading *MOHR COULOMB:'
        write(*,*) '       *MOHR COULOMB should be placed'
        write(*,*) '       before all step definitions'
        ier=1
        return
      endif
!     
      if(nmat.eq.0) then
        write(*,*) '*ERROR reading *MOHR COULOMB:'
        write(*,*) '       *MOHR COULOMB should be'
        write(*,*) '       preceded by a *MATERIAL card'
        ier=1
        return
      endif
!     
      do i=2,n
        write(*,*) 
     &       '*WARNING reading *MOHR COULOMB:'
        write(*,*) '         parameter not recognized:'
        write(*,*) '         ',
     &       textpart(i)(1:index(textpart(i),' ')-1)
        call inputwarning(inpc,ipoinpc,iline,
     &       "MOHR COULOMB%")
      enddo
!     
      nelcon(1,nmat)=-53
      nstate_=max(nstate_,7)
c     !
c     !     putting MOHRCOULOMB in front of the material name
c     !
c     if(matname(nmat)(70:80).ne.'           ') then
c     write(*,*) '*ERROR reading *MOHR COULOMB: the material name'
c     write(*,*) '       for a Mohr-Coulomb material must'
c     write(*,*) '       not exceed 69 characters'
c     ier=1
c     return
c     else
c     do i=80,12,-1
c     matname(nmat)(i:i)=matname(nmat)(i-11:i-11)
c     enddo
c     matname(nmat)(1:11)='MOHRCOULOMB'
c     endif
!     
      iend=2
      do
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
        if((istat.lt.0).or.(key.eq.1)) return
        ntmat=ntmat+1
        nelcon(2,nmat)=ntmat
        if(ntmat.gt.ntmat_) then
          write(*,*) '*ERROR reading *MOHR COULOMB:'
          write(*,*) '       increase ntmat_'
          ier=1
          return
        endif
        do i=1,iend
          read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &         elcon(i+2,ntmat,nmat)
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "MOHR COULOMB%",ier)
            return
          endif
        enddo
!
!       conversion of phi and psi into the parameters k and m
!
        elcon(3,ntmat,nmat)=conversion*elcon(3,ntmat,nmat)
        elcon(4,ntmat,nmat)=conversion*elcon(4,ntmat,nmat)
        elcon(3,ntmat,nmat)=(1.d0+dsin(elcon(3,ntmat,nmat)))/
     &       (1.d0-dsin(elcon(3,ntmat,nmat)))
        elcon(4,ntmat,nmat)=(1.d0+dsin(elcon(4,ntmat,nmat)))/
     &       (1.d0-dsin(elcon(4,ntmat,nmat)))
!        
        read(textpart(3)(1:20),'(f20.0)',iostat=istat) 
     &       elcon(0,ntmat,nmat)
        if(istat.gt.0) then
          call inputerror(inpc,ipoinpc,iline,
     &         "MOHR COULOMB%",ier)
          return
        endif
      enddo
!     
      return
      end

