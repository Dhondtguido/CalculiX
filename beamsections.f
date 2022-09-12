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
      subroutine beamsections(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,ielmat,matname,nmat,ielorien,orname,norien,
     &  thicke,ipkon,iponor,xnor,ixfree,
     &  offset,lakon,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &  ipoinpc,mi,nelcon,ier)
!
!     reading the input deck: *BEAM SECTION
!
      implicit none
!
      logical nodalthickness
!
      character*1 inpc(*)
      character*4 section
      character*8 lakon(*)
      character*80 matname(*),orname(*),material,orientation
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer istartset(*),iendset(*),ialset(*),mi(*),ielmat(mi(3),*),
     &  ipoinpc(0:*),numnod,id,
     &  ielorien(mi(3),*),ipkon(*),iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),nset,nmat,norien,istep,istat,n,key,i,j,k,l,imaterial,
     &  iorientation,ipos,m,iponor(2,*),ixfree,
     &  indexx,indexe,irstrt(*),nelcon(2,*),ier
!
      real*8 thicke(mi(3),*),thickness1,thickness2,p(3),xnor(*),
     &  offset(2,*),offset1,offset2,dd
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) 
     &       '*ERROR reading *BEAM SECTION: *BEAM SECTION should'
         write(*,*) '  be placed before all step definitions'
         ier=1
         return
      endif
!
      nodalthickness=.false.
      offset1=0.d0
      offset2=0.d0
      orientation='                    
     &                           '
      section='    '
      ipos=1
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
         elseif(textpart(i)(1:8).eq.'SECTION=') then
            if(textpart(i)(9:12).eq.'CIRC') then
               section='CIRC'
            elseif(textpart(i)(9:12).eq.'RECT') then
               section='RECT'
            else
               write(*,*) 
     &           '*ERROR reading *BEAM SECTION: unknown section'
               ier=1
               return
            endif
         elseif(textpart(i)(1:8).eq.'OFFSET1=') then
            read(textpart(i)(9:28),'(f20.0)',iostat=istat) offset1
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*BEAM SECTION%",ier)
               return
            endif
         elseif(textpart(i)(1:8).eq.'OFFSET2=') then
            read(textpart(i)(9:28),'(f20.0)',iostat=istat) offset2
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*BEAM SECTION%",ier)
               return
            endif
         else
            write(*,*) 
     &       '*WARNING reading *BEAM SECTION: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*BEAM SECTION%")
         endif
      enddo
!
!     check whether a sections was defined
!
      if(section.eq.'    ') then
         write(*,*) '*ERROR reading *BEAM SECTION: no section defined'
         ier=1
         return
      endif
!
!     check for the existence of the set,the material and orientation
!
      do i=1,nmat
         if(matname(i).eq.material) exit
      enddo
      if(i.gt.nmat) then
         write(*,*) '*ERROR reading *BEAM SECTION: nonexistent material'
         write(*,*) '  '
         call inputerror(inpc,ipoinpc,iline,
     &        "*BEAM SECTION%",ier)
         return
      endif
      imaterial=i
!
      if(orientation.eq.'                    
     &                                 ') then
         iorientation=0
      else
         do i=1,norien
            if(orname(i).eq.orientation) exit
         enddo
         if(i.gt.norien) then
            write(*,*)
     &         '*ERROR reading *BEAM SECTION: nonexistent orientation'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &           "*BEAM SECTION%",ier)
            return
         endif
         iorientation=i
      endif
!
      call cident81(set,elset,nset,id)
      i=nset+1
      if(id.gt.0) then
        if(elset.eq.set(id)) then
          i=id
        endif
      endif
      if(i.gt.nset) then
         elset(ipos:ipos)=' '
         write(*,*) '*ERROR reading *BEAM SECTION: element set ',
     &      elset(1:ipos)
         write(*,*) '  has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &        "*BEAM SECTION%",ier)
         return
      endif
!
!     assigning the elements of the set the appropriate material,
!     orientation number, section and offset(s)
!
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            if(lakon(ialset(j))(1:1).ne.'B') then
               write(*,*) 
     &           '*ERROR reading *BEAM SECTION: *BEAM SECTION can'
               write(*,*) '       only be used for beam elements.'
               write(*,*) '       Element ',ialset(j),' is not a beam el
     &ement.'
               ier=1
               return
            endif
            ielmat(1,ialset(j))=imaterial
            ielorien(1,ialset(j))=iorientation
            offset(1,ialset(j))=offset1
            offset(2,ialset(j))=offset2
            if(section.eq.'RECT') then
               lakon(ialset(j))(8:8)='R'
            else
               lakon(ialset(j))(8:8)='C'
            endif
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               if(lakon(k)(1:1).ne.'B') then
                  write(*,*) 
     &              '*ERROR reading *BEAM SECTION: *BEAM SECTION can'
                  write(*,*) '       only be used for beam elements.'
                  write(*,*) '       Element ',k,' is not a beam element
     &.'
                  ier=1
                  return
               endif
               ielmat(1,k)=imaterial
               ielorien(1,k)=iorientation
               offset(1,k)=offset1
               offset(2,k)=offset2
               if(section.eq.'RECT') then
                  lakon(k)(8:8)='R'
               else
                  lakon(k)(8:8)='C'
               endif
            enddo
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
!     assigning a thickness to the elements
!
      if(.not.nodalthickness) then
         read(textpart(1)(1:20),'(f20.0)',iostat=istat) thickness1
         if(istat.gt.0) then
            write(*,*) 
     &   '*ERROR reading *BEAM SECTION: first beam thickness is lacking'
            call inputerror(inpc,ipoinpc,iline,
     &        "*BEAM SECTION%",ier)
            return
         endif
         if(n.gt.1) then
            read(textpart(2)(1:20),'(f20.0)',iostat=istat) thickness2
            if(istat.gt.0) then
               write(*,*) 
     &              '*ERROR reading *BEAM SECTION: ',
     &              'second beam thickness is lacking'
               call inputerror(inpc,ipoinpc,iline,
     &              "*BEAM SECTION%",ier)
               return
            endif
         else
            thickness2=thickness1
         endif
      else
!
!        for those elements for which nodal thickness is activated
!        the thickness is set to -1.d0
!     
         thickness1=-1.d0
         thickness2=thickness1
      endif
!
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            indexe=ipkon(ialset(j))
            read(lakon(ialset(j))(3:3),'(i1)') numnod
            numnod=numnod+1
            do l=1,numnod
               thicke(1,indexe+l)=thickness1
               thicke(2,indexe+l)=thickness2
            enddo
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               indexe=ipkon(k)
               read(lakon(k)(3:3),'(i1)') numnod
               numnod=numnod+1
               do l=1,numnod
                  thicke(1,indexe+l)=thickness1
                  thicke(2,indexe+l)=thickness2
               enddo
            enddo
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) return
!
!     assigning normal direction 1 for the beam
!
      indexx=-1
      read(textpart(1)(1:20),'(f20.0)',iostat=istat) p(1)
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*BEAM SECTION%",ier)
         return
      endif
      read(textpart(2)(1:20),'(f20.0)',iostat=istat) p(2)
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*BEAM SECTION%",ier)
         return
      endif
      read(textpart(3)(1:20),'(f20.0)',iostat=istat) p(3)
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*BEAM SECTION%",ier)
         return
      endif
      dd=dsqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3))
      if(dd.lt.1.d-10) then
         write(*,*) 
     &       '*ERROR reading *BEAM SECTION: normal in direction 1'
         write(*,*) '       has zero size'
         ier=1
         return
      endif
      do j=1,3
         p(j)=p(j)/dd
      enddo
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            indexe=ipkon(ialset(j))
            read(lakon(ialset(j))(3:3),'(i1)') numnod
            numnod=numnod+1
            do l=1,numnod
               if(indexx.eq.-1) then
                  indexx=ixfree
                  do m=1,3
                     xnor(indexx+m)=p(m)
                  enddo
                  ixfree=ixfree+6
               endif
               iponor(1,indexe+l)=indexx
            enddo
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               indexe=ipkon(k)
               read(lakon(k)(3:3),'(i1)') numnod
               numnod=numnod+1
               do l=1,numnod
               if(indexx.eq.-1) then
                  indexx=ixfree
                  do m=1,3
                     xnor(indexx+m)=p(m)
                  enddo
                  ixfree=ixfree+6
               endif
               iponor(1,indexe+l)=indexx
               enddo
            enddo
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

