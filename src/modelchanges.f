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
      subroutine modelchanges(inpc,textpart,tieset,istat,n,iline,
     &           ipol,inl,ipoinp,inp,ntie,ipoinpc,istep,ipkon,nset,
     &           istartset,iendset,set,ialset,ne,mi,ielmat,iprestr,
     &           iperturb,ier,tietol)
!
!     reading the input deck: *MODEL CHANGE
!
      implicit none
!
      logical contactpair,add,remove,element,strainfree 
!
      character*1 inpc(*)
      character*81 tieset(3,*),noelset,set(*)
      character*132 textpart(16)
!
      integer istat,n,i,key,ipos,iline,ipol,inl,ipoinp(2,*),ipkon(*),
     &  inp(3,*),ntie,ipoinpc(0:*),iposslave,iposmaster,itie,istep,
     &  j,k,m,nset,ne,istartset(*),iendset(*),ialset(*),nelem,
     &     mi(*),ielmat(mi(3),*),iprestr,iperturb(*),ier
!
      real*8 tietol(4,*)
!
      if(istep.eq.0) then
         write(*,*) '*ERROR reading *MODEL CHANGE: *MODEL CHANGE'
         write(*,*) '       cannot be used before the first step'
         ier=1
         return
      endif
!
      contactpair=.false.
      add=.false.
      remove=.false.
      element=.false.
      strainfree=.false.
!
      do i=2,n
         if(textpart(i)(1:16).eq.'TYPE=CONTACTPAIR') then
            contactpair=.true.
         elseif(textpart(i)(1:12).eq.'TYPE=ELEMENT') then
            element=.true.
         elseif(textpart(i)(1:14).eq.'ADD=STRAINFREE') then
            add=.true.
            strainfree=.true.
         elseif(textpart(i)(1:14).eq.'ADD=WITHSTRAIN') then
            add=.true.
         elseif(textpart(i)(1:3).eq.'ADD') then
            add=.true.
            strainfree=.true.
         elseif(textpart(i)(1:20).eq.'MECHSTRAINTORESIDUAL')
     &     then
            strainfree=.true.
         elseif(textpart(i)(1:6).eq.'REMOVE') then
            remove=.true.
         else
            write(*,*) 
     &       '*WARNING reading *MODEL CHANGE: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*MODEL CHANGE%")
         endif
      enddo
!
!     checking the validity of the input
!
      if((.not.contactpair).and.(.not.element)) then
         write(*,*) '*ERROR reading *MODEL CHANGE: model change can'
         write(*,*) '       only be used for contact pairs or elements'
         ier=1
         return
      endif
!
      if((contactpair).and.(element)) then
         write(*,*) '*ERROR reading *MODEL CHANGE: model change cannot'
         write(*,*) '       be used for contact pairs and elements at'
         write(*,*) '       the same time'
         ier=1
         return
      endif
!
      if(.not.strainfree) then
         if((.not.add).and.(.not.remove)) then
            write(*,*) '*ERROR reading *MODEL CHANGE: at least ADD or'
            write(*,*) '        REMOVE has to be selected'
            ier=1
            return
         endif
      endif
!
      if(add.and.remove) then
         write(*,*) '*ERROR reading *MODEL CHANGE: ADD and REMOVE'
         write(*,*) '       cannot both be selected'
         ier=1
         return
      endif
!
      if((element).and.(strainfree)) then
         if(iperturb(1).gt.1) then
            iprestr=2
         else
            write(*,*) '*ERROR reading *MODEL CHANGE'
            write(*,*) '       a strain-free addition of elements'
            write(*,*) '       or stress removal in elements'
            write(*,*) '       is only possible for nonlinear'
            write(*,*) '       calculations'
            ier=1
            return
         endif
      endif
!
!     reading the slave and the master surface
!
      if(contactpair) then
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) then
            write(*,*)'*ERROR reading *MODEL CHANGE: definition of the '
            write(*,*) '      contact pair is not complete.'
            ier=1
            return
         endif
!
!     selecting the appropriate action
!
         iposslave=index(textpart(1)(1:80),' ')
         iposmaster=index(textpart(2)(1:80),' ')
         do i=1,ntie
            if((tieset(1,i)(81:81).ne.'C').and.
     &           (tieset(1,i)(81:81).ne.'-')) cycle
            ipos=index(tieset(2,i),' ')-1
            if(ipos.ne.iposslave) cycle
            if(tieset(2,i)(1:ipos-1).ne.textpart(1)(1:ipos-1)) cycle
            ipos=index(tieset(3,i),' ')-1
            if(ipos.ne.iposmaster) cycle
            if(tieset(3,i)(1:ipos-1).ne.textpart(2)(1:ipos-1)) cycle
            itie=i
            exit
         enddo
!
         if(add) then
           tieset(1,itie)(81:81)='C'
           tietol(4,itie)=1.d0
         else
            tieset(1,itie)(81:81)='-'
         endif
!     
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
      else
!
!       element change
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            do i=1,n
               read(textpart(i)(1:10),'(i10)',iostat=istat) 
     &              nelem
               if(istat.gt.0) then
!     
!                 set name
!     
                  noelset=textpart(i)(1:80)
                  noelset(81:81)=' '
                  ipos=index(noelset,' ')
                  noelset(ipos:ipos)='E'
c                  do j=1,nset
c                     if(noelset.eq.set(j)) then
                  call cident81(set,noelset,nset,j)
                  if(j.gt.0) then
                    if(noelset.eq.set(j)) then
                        m=iendset(j)-istartset(j)+1
                        do k=1,m
                           nelem=ialset(istartset(j)+k-1)
                           if((add).or.(remove)) then 
                              ipkon(nelem)=-2-ipkon(nelem)
                           endif
                           if(strainfree) ielmat(1,nelem)=
     &                                   -ielmat(1,nelem)
                        enddo
c                        exit
                      endif
                    endif
c                     endif
c                  enddo
                  if(noelset.ne.set(j)) then
                     noelset(ipos:ipos)=' '
                     write(*,*) '*ERROR reading *MODEL CHANGE:',
     &                  ' element set ',noelset
                     write(*,*) '       has not been defined yet'
                     ier=1
                     return
                  endif
               else
!     
!                 node or element number
!     
                  if(nelem.gt.ne) then
                     write(*,*) '*WARNING in noelsets: element ',
     &                    nelem
                     write(*,*) '          > ne;'
                  else
                     if((add).or.(remove)) then
                        ipkon(nelem)=-2-ipkon(nelem)
                     endif
                     if(strainfree) ielmat(1,nelem)=-ielmat(1,nelem)
                  endif
               endif
            enddo
         enddo
      endif
!     
      return
      end



