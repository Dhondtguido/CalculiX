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
      subroutine membranesections(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,ielmat,matname,nmat,ielorien,orname,norien,
     &  thicke,kon,ipkon,offset,irstrt,istep,istat,n,iline,ipol,
     &  inl,ipoinp,inp,lakon,iaxial,ipoinpc,mi,icomposite,nelcon,ier)
!
!     reading the input deck: *MEMBRANE SECTION
!
      implicit none
!
      logical nodalthickness,composite
!
      character*1 inpc(*)
      character*8 lakon(*)
      character*80 matname(*),orname(*),material,orientation
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer mi(*),istartset(*),iendset(*),ialset(*),ielmat(mi(3),*),
     &  ielorien(mi(3),*),kon(*),ipkon(*),indexe,irstrt(*),nset,nmat,
     &  norien,nlayer,iset,icomposite,nelcon(2,*),ier,id,
     &  istep,istat,n,key,i,j,k,l,imaterial,iorientation,ipos,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),iaxial,ipoinpc(0:*)
!
      real*8 thicke(mi(3),*),thickness,offset(2,*),offset1
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) 
     &      '*ERROR reading *MEMBRANE SECTION: *MEMBRANE SECTION should'
         write(*,*) '  be placed before all step definitions'
         ier=1
         return
      endif
!
      nodalthickness=.false.
      composite=.false.
      offset1=0.d0
      material(1:1)=' '
      orientation(1:1)=' '
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
         elseif(textpart(i)(1:14).eq.'NODALTHICKNESS') then
            nodalthickness=.true.
         elseif(textpart(i)(1:7).eq.'OFFSET=') then
            read(textpart(i)(8:27),'(f20.0)',iostat=istat) offset1
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*MEMBRANE SECTION%",ier)
               return
            endif
         elseif(textpart(i)(1:9).eq.'COMPOSITE') then
            composite=.true.
         else
            write(*,*) 
     &   '*WARNING reading *MEMBRANE SECTION: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*MEMBRANE SECTION%")
         endif
      enddo
!
!     check for the existence of the material (not for composites)
!
      if(.not.composite) then
         do i=1,nmat
            if(matname(i).eq.material) exit
         enddo
         if(i.gt.nmat) then
            write(*,*) 
     &        '*ERROR reading *MEMBRANE SECTION: nonexistent material'
            call inputerror(inpc,ipoinpc,iline,
     &           "*MEMBRANE SECTION%",ier)
            return
         endif
         imaterial=i
      elseif(material(1:1).ne.' ') then
         write(*,*) '*ERROR reading *MEMBRANE SECTION: COMPOSITE and'
         write(*,*) '       MATERIAL are mutually exclusive parameters'
         ier=1
         return
      endif
!
!     check for the existence of the orientation, if any
!
      if(orientation(1:1).eq.' ') then
         iorientation=0
c      elseif(nelcon(1,i).eq.2) then
c         write(*,*) '*INFO reading *MEMBRANE SECTION: an orientation'
c         write(*,*) '      is for isotropic materials irrelevant'
c         call inputinfo(inpc,ipoinpc,iline,
c     &"*MEMBRANE SECTION%")
c         iorientation=0
      else
         do i=1,norien
            if(orname(i).eq.orientation) exit
         enddo
         if(i.gt.norien) then
            write(*,*)
     &       '*ERROR reading *MEMBRANE SECTION: nonexistent orientation'
            call inputerror(inpc,ipoinpc,iline,
     &           "*MEMBRANE SECTION%",ier)
            return
         endif
         iorientation=i
      endif
!
!     check for the existence of the element set
!
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
         write(*,*) '*ERROR reading *MEMBRANE SECTION: element set ',
     &        elset
         write(*,*) '       has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &        "*MEMBRANE SECTION%",ier)
         return
      endif
      iset=i
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
!     assigning a thickness to the elements
!
      if(.not.composite) then
         if(.not.nodalthickness) then
            read(textpart(1)(1:20),'(f20.0)',iostat=istat) thickness
            if(istat.gt.0) then
               write(*,*) 
     & '*ERROR reading *MEMBRANE SECTION: membrane thickness is lacking'
               call inputerror(inpc,ipoinpc,iline,
     &              "*MEMBRANE SECTION%",ier)
               return
            endif
            if(iaxial.eq.180) thickness=thickness/iaxial
            do j=istartset(iset),iendset(iset)
               if(ialset(j).gt.0) then
                  if(lakon(ialset(j))(1:1).ne.'M') then
                     write(*,*) 
     &         '*ERROR reading *MEMBRANE SECTION: *MEMBRANE SECTION can'
                     write(*,*) 
     &                 '       only be used for membrane elements.'
                     write(*,*) '       Element ',ialset(j),
     &                 ' is not a membrane element.'
                     ier=1
                     return
                  endif
                  indexe=ipkon(ialset(j))
                  do l=1,8
                     thicke(1,indexe+l)=thickness
                  enddo
                  ielmat(1,ialset(j))=imaterial
                  ielorien(1,ialset(j))=iorientation
                  offset(1,ialset(j))=offset1
               else
                  k=ialset(j-2)
                  do
                     k=k-ialset(j)
                     if(k.ge.ialset(j-1)) exit
                     if(lakon(k)(1:1).ne.'M') then
                        write(*,*) 
     &         '*ERROR reading *MEMBRANE SECTION: *MEMBRANE SECTION can'
                        write(*,*) 
     &                    '       only be used for membrane elements.'
                        write(*,*) '       Element ',k,
     &                    ' is not a membrane element.'
                        ier=1
                        return
                     endif
                     indexe=ipkon(k)
                     do l=1,8
                        thicke(1,indexe+l)=thickness
                     enddo
                     ielmat(1,k)=imaterial
                     ielorien(1,k)=iorientation
                     offset(1,k)=offset1
                  enddo
               endif
            enddo
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         endif
!     
      else
         if(nodalthickness) then
            write(*,*) '*ERROR reading *MEMBRANE SECTION: for composite'
            write(*,*) '       materials is the parameter NODAL'
            write(*,*) '       THICKNESS not allowed'
            ier=1
            return
         endif
!
c         icomposite=1
         nlayer=0
         do
            read(textpart(1)(1:20),'(f20.0)',iostat=istat) thickness
            if(istat.gt.0) then
               write(*,*) 
     & '*ERROR reading *MEMBRANE SECTION: membrane thickness is lacking'
               call inputerror(inpc,ipoinpc,iline,
     &              "*MEMBRANE SECTION%",ier)
               return
            endif
            if(iaxial.eq.180) thickness=thickness/iaxial
!     
!     reading the material name
!     
            read(textpart(3)(1:80),'(a80)',iostat=istat) material
            if(istat.gt.0) then
               write(*,*) 
     &           '*ERROR reading *MEMBRANE SECTION: no material defined'
               call inputerror(inpc,ipoinpc,iline,
     &              "*MEMBRANE SECTION%",ier)
               return
            endif
!     
!     check for the existence of the material
!     
            do i=1,nmat
               if(matname(i).eq.material) exit
            enddo
            if(i.gt.nmat) then
               write(*,*) 
     &      '*ERROR reading *MEMBRANE SECTION: nonexistent material'
               call inputerror(inpc,ipoinpc,iline,
     &              "*MEMBRANE SECTION%",ier)
               return
            endif
            imaterial=i
!     
!     reading the orientation, if any
!     if no orientation is specified, the global orientation defined
!     by the ORIENTATION parameter, if any, will be used
!     
            read(textpart(4)(1:80),'(a80)',iostat=istat) orientation
!     
c            if(orientation(1:1).eq.' ') then
c               iorientation=0
c            else
            if(orientation(1:1).ne.' ') then
               do i=1,norien
                  if(orname(i).eq.orientation) exit
               enddo
               if(i.gt.norien) then
                  write(*,*)
     &      '*ERROR reading *MEMBRANE SECTION: nonexistent orientation'
                  write(*,*) '  '
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*MEMBRANE SECTION%",ier)
                  return
               endif
               iorientation=i
            endif
!     
            nlayer=nlayer+1
!     
            do j=istartset(iset),iendset(iset)
               if(ialset(j).gt.0) then
                  if((lakon(ialset(j))(1:3).ne.'M3D8R').and.
     &               (lakon(ialset(j))(1:2).ne.'M3D6')) then
                     write(*,*) 
     &           '*ERROR reading *MEMBRANE SECTION: *MEMBRANE SECTION'
                     write(*,*) 
     &                 '       with the option COMPOSITE can'
                     write(*,*) 
     &        '       only be used for M3D8R or M3D6 membrane elements.'
                     write(*,*) '       Element ',ialset(j),
     &                 ' is not a M3D8R nor a M3D6 membrane element.'
                     ier=1
                     return
                  endif
                  indexe=ipkon(ialset(j))
                  do l=1,8
                     thicke(nlayer,indexe+l)=thickness
                  enddo
                  ielmat(nlayer,ialset(j))=imaterial
                  ielorien(nlayer,ialset(j))=iorientation
                  offset(1,ialset(j))=offset1
                  if(nlayer.gt.1) lakon(ialset(j))(8:8)='C'
               else
                  k=ialset(j-2)
                  do
                     k=k-ialset(j)
                     if(k.ge.ialset(j-1)) exit
                     if((lakon(k)(1:3).ne.'M3D8R').and.
     &                  (lakon(k)(1:2).ne.'M3D6')) then
                        write(*,*) 
     &            '*ERROR reading *MEMBRANE SECTION: *MEMBRANE SECTION'
                        write(*,*) 
     &                       '    with the option COMPOSITE can'
                        write(*,*) 
     &        '       only be used for M3D8R or M3D6 membrane elements.'
                        write(*,*) '       Element ',k,
     &                    ' is not a M3D8R nor a M3D6 membrane element.'
                        ier=1
                        return
                     endif
                     indexe=ipkon(k)
                     do l=1,8
                        thicke(nlayer,indexe+l)=thickness
                     enddo
                     ielmat(nlayer,k)=imaterial
                     ielorien(nlayer,k)=iorientation
                     offset(1,k)=offset1
                     if(nlayer.gt.1) lakon(k)(8:8)='C'
                  enddo
               endif
            enddo
!     
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
               if(nlayer.gt.1) icomposite=1
               return
            endif
         enddo
      endif
!     
      return
      end
      
