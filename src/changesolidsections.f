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
      subroutine changesolidsections(inpc,textpart,set,istartset,
     &  iendset,ialset,nset,ielmat,matname,nmat,ielorien,orname,norien,
     &  lakon,thicke,kon,ipkon,irstrt,istep,istat,n,iline,ipol,inl,
     &  ipoinp,inp,cs,mcs,iaxial,ipoinpc,mi,nelcon,ier)
!
!     reading the input deck: *CHANGE SOLID SECTION
!
      implicit none
!
      character*1 inpc(*)
      character*8 lakon(*)
      character*80 matname(*),orname(*),material,orientation
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer mi(*),istartset(*),iendset(*),ialset(*),ielmat(mi(3),*),
     &  ielorien(mi(3),*),kon(*),ipkon(*),indexe,irstrt(*),nset,nmat,
     &  norien,nelcon(2,*),ier,id,
     &  istep,istat,n,key,i,j,k,l,imaterial,iorientation,ipos,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),mcs,iaxial,ipoinpc(0:*)
!
      real*8 thicke(mi(3),*),thickness,pi,cs(17,*)
!
      if(istep.eq.0) then
         write(*,*) '*ERROR reading *CHANGE SOLID SECTION:'
         write(*,*) '       *CHANGE SOLID SECTION should'
         write(*,*) '       be placed within a step definition'
         ier=1
         return
      endif
!
      pi=4.d0*datan(1.d0)
!
      orientation='
     &                           '
      elset='
     &                      '
      ipos=0
!
      do i=2,n
         if(textpart(i)(1:9).eq.'MATERIAL=') then
            material=textpart(i)(10:89)
         elseif(textpart(i)(1:12).eq.'ORIENTATION=') then
            orientation=textpart(i)(13:92)
         elseif(textpart(i)(1:6).eq.'ELSET=') then
            elset=textpart(i)(7:86)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
         else
            write(*,*) '*WARNING reading *CHANGE SOLID SECTION:'
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*CHANGE SOLID SECTIONQ%")
         endif
      enddo
!
!     check for the existence of the set,the material and orientation
!
      do i=1,nmat
         if(matname(i).eq.material) exit
      enddo
      if(i.gt.nmat) then
         do i=1,nmat
            if(matname(i)(1:11).eq.'ANISO_CREEP') then
               if(matname(i)(12:20).eq.material(1:9)) exit
            elseif(matname(i)(1:10).eq.'ANISO_PLAS') then
               if(matname(i)(11:20).eq.material(1:10)) exit
            endif
         enddo
      endif
      if(i.gt.nmat) then
         write(*,*) 
     &     '*ERROR reading *CHANGE SOLID SECTION: nonexistent material'
         write(*,*) '  '
         call inputerror(inpc,ipoinpc,iline,
     &        "*CHANGE SOLID SECTION%",ier)
         return
      endif
      imaterial=i
!
      if(orientation.eq.'                    ') then
         iorientation=0
c      elseif(nelcon(1,i).eq.2) then
c         write(*,*) 
c     &      '*INFO reading *CHANGE SOLID SECTION: an orientation'
c         write(*,*) '      is for isotropic materials irrelevant'
c         call inputinfo(inpc,ipoinpc,iline,
c     &"*CHANGE SOLID SECTION%")
c         iorientation=0
      else
         do i=1,norien
            if(orname(i).eq.orientation) exit
         enddo
         if(i.gt.norien) then
            write(*,*) '*ERROR reading *CHANGE SOLID SECTION:'
            write(*,*) '       nonexistent orientation'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &           "*CHANGE SOLID SECTION%",ier)
            return
         endif
         iorientation=i
      endif
!
      if(ipos.eq.0) then
         write(*,*) 
     &    '*ERROR reading *CHANGE SOLID SECTION: no element set ',elset
         write(*,*) '       was been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &        "*CHANGE SOLID SECTION%",ier)
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
         write(*,*) '*ERROR reading *CHANGE SOLID SECTION:'
         write(*,*) '       element set ',elset
         write(*,*) '  has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &        "*CHANGE SOLID SECTION%",ier)
         return
      endif
!
!     assigning the elements of the set the appropriate material
!     and orientation number
!
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            if((lakon(ialset(j))(1:1).eq.'B').or.
     &         (lakon(ialset(j))(1:1).eq.'S')) then
               write(*,*) '*ERROR reading *CHANGE SOLID SECTION:'
               write(*,*) '       *CHANGE SOLID SECTION can'
               write(*,*) '       not be used for beam or shell elements
     &'
               write(*,*) '       Faulty element: ',ialset(j)
               ier=1
               return
            endif
            ielmat(1,ialset(j))=imaterial
            if(norien.gt.0) ielorien(1,ialset(j))=iorientation
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               if((lakon(k)(1:1).eq.'B').or.
     &              (lakon(k)(1:1).eq.'S')) then
                  write(*,*) '*ERROR reading *CHANGE SOLID SECTION:'
                  write(*,*) '       *CHANGE SOLID SECTION can'
                  write(*,*) '       not be used for beam or shell eleme
     &nts'
                  write(*,*) '       Faulty element: ',k
                  ier=1
                  return
               endif
               ielmat(1,k)=imaterial
               if(norien.gt.0) ielorien(1,k)=iorientation
            enddo
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
!     assigning a thickness to plane stress elements and an angle to
!     axisymmetric elements
!
      if(key.eq.0) then
         write(*,*) '*ERROR reading *CHANGE SOLID SECTION'
         write(*,*) '       no second line allowed'
         call inputerror(inpc,ipoinpc,iline,
     &        "*CHANGE SOLID SECTION%",ier)
         return
      endif
!
      return
      end

