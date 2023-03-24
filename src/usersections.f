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
      subroutine usersections(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,ielmat,matname,nmat,
     &  irstrt,istep,istat,n,iline,ipol,
     &  inl,ipoinp,inp,lakon,ielprop,nprop,nprop_,prop,
     &  ipoinpc,mi,ier)
!
!     reading the input deck: *USER SECTION
!
      implicit none
!
      character*1 inpc(*)
      character*8 lakon(*)
      character*80 matname(*),material
      character*81 set(*),elset
      character*132 textpart(16)
!     
      integer mi(*),istartset(*),iendset(*),ialset(*),id,
     &     ielmat(mi(3),*),irstrt(*),nset,nmat,ndprop,npropstart,
     &     istep,istat,n,key,i,j,k,imaterial,ipos,lprop,ipoinpc(0:*),
     &     iline,ipol,inl,ipoinp(2,*),inp(3,*),ielprop(*),nprop,nprop_,
     &     iset,ier
!     
      real*8 prop(*)
!     
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
        write(*,*) 
     &       '*ERROR reading *USER SECTION: *USER SECTION should'
        write(*,*) '  be placed before all step definitions'
        ier=1
        return
      endif
!     
      do i=2,n
        if(textpart(i)(1:9).eq.'MATERIAL=') then
          material=textpart(i)(10:89)
        elseif(textpart(i)(1:6).eq.'ELSET=') then
          elset=textpart(i)(7:86)
          elset(81:81)=' '
          ipos=index(elset,' ')
          elset(ipos:ipos)='E'
!     
        elseif(textpart(i)(1:10).eq.'CONSTANTS=') then
          read(textpart(i)(11:20),'(i10)',iostat=istat) ndprop
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*USER SECTION%",ier)
            return
          endif
        else
          write(*,*) 
     &       '*WARNING reading *USER SECTION: parameter not recognized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*USER SECTION%")
        endif
      enddo
!     
!     check for the existence of the set and the material
!     
      do i=1,nmat
        if(matname(i).eq.material) exit
      enddo
      if(i.gt.nmat) then
        write(*,*) 
     &       '*ERROR reading *USER SECTION: nonexistent material'
        write(*,*) '  '
        call inputerror(inpc,ipoinpc,iline,
     &       "*USER SECTION%",ier)
        return
      endif
      imaterial=i
!     
c      do i=1,nset
c        if(set(i).eq.elset) exit
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
        write(*,*) '*ERROR reading *USER SECTION: element set ',elset
        write(*,*) '  has not yet been defined. '
        call inputerror(inpc,ipoinpc,iline,
     &       "*USER SECTION%",ier)
        return
      endif
      iset=i
!     
      npropstart=nprop
!     
!     reading the element properties
!     
      if(ndprop.gt.0) then
!     
!     general case
!     
        lprop=0
!     
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) exit
          do k=1,n
            lprop=lprop+1
            if(lprop.gt.ndprop) exit
            read(textpart(k),'(f20.0)',iostat=istat) 
     &           prop(nprop+lprop)
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*USER SECTION%",ier)
              return
            endif
          enddo
        enddo
        nprop=nprop+ndprop
      else
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
      endif
!     
      if(nprop.gt.nprop_) then
        write(*,*) '*ERROR reading *USER SECTION: increase nprop_'
        ier=1
        return
      endif
!     
!     assigning the elements of the set the appropriate material
!     and property pointer
!     
      do j=istartset(iset),iendset(iset)
        if(ialset(j).gt.0) then
          if(lakon(ialset(j))(1:1).ne.'U') then
            write(*,*) '*ERROR reading *USER SECTION: element ',
     &           ialset(j),' is no user element.'
            ier=1
            return
          endif
!     
          ielmat(1,ialset(j))=imaterial
          if(ndprop.gt.0) ielprop(ialset(j))=npropstart
        else
          k=ialset(j-2)
          do
            k=k-ialset(j)
            if(k.ge.ialset(j-1)) exit
            if(lakon(k)(1:1).ne.'U') then
              write(*,*) '*ERROR reading *USER SECTION: element ',
     &             k,' is no user element.'
              ier=1
              return
            endif
            ielmat(1,k)=imaterial
            if(ndprop.gt.0) ielprop(k)=npropstart
          enddo
        endif
      enddo
!     
      return
      end
