!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine ratedependents(inpc,textpart,nelcon,nmat,ntmat_,
     &     iplas,iperturb,nstate_,ncmat_,elcon,matname,irstrt,istep,
     &     istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
!     reading the input deck: *RATE DEPENDENT
!     
      implicit none
!     
      logical johnsoncook
!     
      character*1 inpc(*)
      character*80 matname(*)
      character*132 textpart(16)
!     
      integer nelcon(2,*),nmat,ntmat_,ntmat,npmat,istep,
     &     n,key,i,ncmat_,
     &     iplas,iperturb(*),istat,nstate_,
     &     irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),ipoinpc(0:*),
     &     ier
!     
      real*8 elcon(0:ncmat_,ntmat_,*)
!     
      johnsoncook=.false.
!     
      npmat=0
!     
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
        write(*,*) '*ERROR reading *RATE DEPENDENT:'
        write(*,*) '       *RATE DEPENDENT should be placed'
        write(*,*) '       before all step definitions'
        ier=1
        return
      endif
!     
      if(nmat.eq.0) then
        write(*,*) '*ERROR reading *RATE DEPENDENT:'
        write(*,*) '       *RATE DEPENDENT should be preceded'
        write(*,*) '       by a *MATERIAL card'
        ier=1
        return
      endif
!     
      if((nelcon(1,nmat).ne.2).and.
     &     (matname(nmat)(1:11).ne.'JOHNSONCOOK')) then
        write(*,*) '*ERROR reading *RATE DEPENDENT:'
        write(*,*) '       *RATE DEPENDENT should be preceded'
        write(*,*) '       by an *ELASTIC,TYPE=ISO card'
        ier=1
        return
      endif
!     
      iperturb(1)=3
!     
      do i=2,n
        if(textpart(i)(1:16).eq.'TYPE=JOHNSONCOOK') then
          johnsoncook=.true.
        else
          write(*,*) '*WARNING reading *RATE DEPENDENT:'
          write(*,*) '         parameter not recognized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*RATE DEPENDENT%")
        endif
      enddo
!     
!     Johnson-Cook
!     user material; npmat=0; ntmat=1;
!     
      if(johnsoncook) then
        iplas=1
        if(matname(nmat)(1:11).ne.'JOHNSONCOOK') then
          write(*,*) '*ERROR reading *RATE DEPENDENT'
          write(*,*) '       the name of a Johnson Cook material'
          write(*,*) '       must start with JOHNSONCOOK'
          write(*,*) '       (blanks are allowed at any location'
          write(*,*) '        and the string is not case sensitive)'
          call inputerror(inpc,ipoinpc,iline,
     &         "*RATE DEPENDENT%",ier)
          return
        endif
        nelcon(1,nmat)=-111
        nstate_=max(nstate_,9)
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
        if((istat.lt.0).or.(key.eq.1).or.(n.lt.2)) then
          write(*,*) '*ERROR reading *RATE DEPENDENT'
          write(*,*) '       for the Johnson-Cook model at least'
          write(*,*) '       C and the reference strain rate must'
          write(*,*) '       be given'
          call inputerror(inpc,ipoinpc,iline,
     &         "*RATE DEPENDENT%",ier)
          return
        endif
        ntmat=1
        if(ntmat.gt.ntmat_) then
          write(*,*) 
     &         '*ERROR reading *RATE DEPENDENT: increase ntmat_'
          ier=1
          return
        endif
!     
!     reading C and the reference strain rate (elcon(7,8))
!     
        do i=1,2
          read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &         elcon(6+i,ntmat,nmat)
        enddo
        if(elcon(8,ntmat,nmat).le.0.d0) then
          write(*,*) '*ERROR reading *RATE DEPENDENT'
          write(*,*) '       the reference strain rate must be'
          write(*,*) '       strictly positive'
          call inputerror(inpc,ipoinpc,iline,
     &         "*RATE DEPENDENT%",ier)
          return
        endif
      else
        write(*,*) '*ERROR reading *RATE DEPENDENT'
        write(*,*) '       TYPE=JOHNSON COOK is lacking'
        call inputerror(inpc,ipoinpc,iline,
     &       "*RATE DEPENDENT%",ier)
        return
      endif
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      return
      end

