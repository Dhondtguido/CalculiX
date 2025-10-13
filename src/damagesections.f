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
      subroutine damagesections(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,ieldam,matname,nmat,ndam,
     &  lakon,kon,ipkon,irstrt,istep,istat,n,iline,ipol,inl,
     &  ipoinp,inp,ipoinpc,mi,co,ier)
!
!     reading the input deck: *DAMAGE SECTION
!
      implicit none
!
      character*1 inpc(*)
      character*8 lakon(*)
      character*80 matname(*),damagemodel
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer mi(*),istartset(*),iendset(*),ialset(*),ieldam(mi(3),*),
     &  kon(*),ipkon(*),irstrt(*),nset,nmat,ielem,node1,node2,m,id,
     &  istep,istat,n,key,i,j,k,l,idamagemodel,ipos,ndam,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),ipoinpc(0:*),ier,numnod
!
      real*8 xn(3),co(3,*),dd
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *DAMAGE SECTION: *DAMAGE SECTION'
         write(*,*)'       should be placed before all step definitions'
         ier=1
         return
      endif
!
      elset='
     &                      '
      ipos=0
!
      do i=2,n
         if(textpart(i)(1:12).eq.'DAMAGEMODEL=') then
            damagemodel=textpart(i)(13:92)
         elseif(textpart(i)(1:6).eq.'ELSET=') then
            elset=textpart(i)(7:86)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
         else
            write(*,*) 
     &     '*WARNING reading *DAMAGE SECTION: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*DAMAGE SECTION%")
         endif
      enddo
!
!     check for the existence of the damage model
!
      do i=1,nmat
         if(matname(i).eq.damagemodel) exit
      enddo
      if(i.gt.nmat) then
         write(*,*) 
     &      '*ERROR reading *DAMAGE SECTION: nonexistent damage model'
         write(*,*) '  '
         call inputerror(inpc,ipoinpc,iline,
     &        "*DAMAGE SECTION%",ier)
         return
      endif
      idamagemodel=i
!
!     check for the existence of the set
!
      if(ipos.eq.0) then
         write(*,*) '*ERROR reading *DAMAGE SECTION: no element set ',
     &        elset
         write(*,*) '       was been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &        "*DAMAGE SECTION%",ier)
         return
      endif
      call cident81(set,elset,nset,id)
      i=nset+1
      if(id.gt.0) then
        if(elset.eq.set(id)) then
          i=id
        endif
      endif
      if(i.gt.nset) then
         elset(ipos:ipos)=' '
         write(*,*) '*ERROR reading *DAMAGE SECTION: element set ',elset
         write(*,*) '       has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &        "*DAMAGE SECTION%",ier)
         return
      endif
!
!     assigning the elements of the set the appropriate damage model
!
      do j=istartset(i),iendset(i)
        if(ialset(j).gt.0) then
          k=ialset(j)
          ndam=1
          ieldam(1,k)=idamagemodel
        else
          k=ialset(j-2)
          do
            k=k-ialset(j)
            if(k.ge.ialset(j-1)) exit
            ndam=1
            ieldam(1,k)=idamagemodel
          enddo
        endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

