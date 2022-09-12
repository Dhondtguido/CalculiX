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
      subroutine masss(inpc,textpart,nrhcon,nmat,ntmat_,
     &        rhcon,matname,irstrt,istep,istat,n,iline,ipol,
     &        inl,ipoinp,inp,nmat_,set,istartset,iendset,ialset,
     &        nset,ielmat,ielorien,ipoinpc,mi,iaxial,ier)
!
!     reading the input deck: *MASS
!
      implicit none
!
      character*1 inpc(*)
      character*80 matname(*)
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer mi(*),nrhcon(*),nmat,ntmat_,ntmat,istep,id,
     &  n,key,i,istat,istartset(*),iaxial,ier,
     &  iendset(*),irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),nmat_,
     &  ialset(*),ipos,nset,j,k,ielmat(mi(3),*),ielorien(mi(3),*),
     &  ipoinpc(0:*)  
!
      real*8 rhcon(0:1,ntmat_,*)
!
      ntmat=0
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *MASS: *MASS should be placed'
         write(*,*) '  before all step definitions'
         ier=1
         return
      endif
!
      nmat=nmat+1
      if(nmat.gt.nmat_) then
         write(*,*) '*ERROR reading *MASS: increase nmat_'
         ier=1
         return
      endif
      matname(nmat)(1:4)='MASS'
      do i=5,80
         matname(nmat)(i:i)=' '
      enddo
!
      do i=2,n
         if(textpart(i)(1:6).eq.'ELSET=') then
            elset=textpart(i)(7:86)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
         else
            write(*,*) 
     &           '*WARNING reading *MASS: parameter not recognized:'
            write(*,*) '         ',
     &           textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*MASS%")
         endif
      enddo
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*) '*ERROR reading *MASS: definition of the'
         write(*,*) '       mass is not complete'
         call inputerror(inpc,ipoinpc,iline,
     &        "*MASS%",ier)
         return
      endif
!
!     reading the mass and storing it in rhcon
!
      ntmat=ntmat+1
      nrhcon(nmat)=ntmat
      if(ntmat.gt.ntmat_) then
         write(*,*) '*ERROR reading *MASS: increase ntmat_'
         ier=1
         return
      endif
      read(textpart(1)(1:20),'(f20.0)',iostat=istat)
     &     rhcon(1,ntmat,nmat)
!
!     dividing by 180 for axisymmetric structures
!
      if(iaxial.eq.180) rhcon(1,ntmat,nmat)=rhcon(1,ntmat,nmat)/iaxial
!
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*MASS%",ier)
         return
      endif
      rhcon(0,ntmat,nmat)=0.d0
!     
      if(ntmat.eq.0) then
         write(*,*) '*ERROR reading *MASS: *MASS card without data'
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
         write(*,*) '*ERROR reading *MASS: element set ',elset
         write(*,*) '       has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &        "*MASS%",ier)
         return
      endif
!     
!     assigning the elements of the set the appropriate material
!     
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            ielmat(1,ialset(j))=nmat
            ielorien(1,ialset(j))=0
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               ielmat(1,k)=nmat
               ielorien(1,k)=0
            enddo
         endif
      enddo
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      return
      end

