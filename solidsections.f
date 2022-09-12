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
      subroutine solidsections(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,ielmat,matname,nmat,ielorien,orname,norien,
     &  lakon,thicke,kon,ipkon,irstrt,istep,istat,n,iline,ipol,inl,
     &  ipoinp,inp,cs,mcs,iaxial,ipoinpc,mi,co,ixfree,xnor,iponor,
     &  ier,orab)
!
!     reading the input deck: *SOLID SECTION
!
      implicit none
!
      logical nodalthickness
!
      character*1 inpc(*)
      character*8 lakon(*)
      character*80 matname(*),orname(*),material,orientation
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer mi(*),istartset(*),iendset(*),ialset(*),ielmat(mi(3),*),
     &  ielorien(mi(3),*),kon(*),ipkon(*),indexe,irstrt(*),nset,nmat,
     &  norien,ielem,node1,node2,m,indexx,ixfree,iponor(2,*),id,
     &  istep,istat,n,key,i,j,k,l,imaterial,iorientation,ipos,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),mcs,iaxial,ipoinpc(0:*),
     &  ier,numnod
!
      real*8 thicke(mi(3),*),thickness,pi,cs(17,*),xn(3),co(3,*),p(3),
     &     dd,xnor(*),orab(7,*)
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *SOLID SECTION: *SOLID SECTION'
         write(*,*)'       should be placed before all step definitions'
         ier=1
         return
      endif
!
      nodalthickness=.false.
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
         elseif(textpart(i)(1:14).eq.'NODALTHICKNESS') then
            nodalthickness=.true.
         else
            write(*,*) 
     &      '*WARNING reading *SOLID SECTION: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*SOLID SECTION%")
         endif
      enddo
!
!     check for the existence of the material
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
     &      '*ERROR reading *SOLID SECTION: nonexistent material'
         write(*,*) '  '
         call inputerror(inpc,ipoinpc,iline,
     &        "*SOLID SECTION%",ier)
         return
      endif
      imaterial=i
!
!     check for the existence of the orientation
!
      if(orientation.eq.'                    ') then
         iorientation=0
      else
         do i=1,norien
            if(orname(i).eq.orientation) exit
         enddo
         if(i.gt.norien) then
            write(*,*)
     &       '*ERROR reading *SOLID SECTION: nonexistent orientation'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &           "*SOLID SECTION%",ier)
            return
         endif
         iorientation=i
      endif
!
!     check for the existence of the set
!
      if(ipos.eq.0) then
         write(*,*) '*ERROR reading *SOLID SECTION: no element set ',
     &        elset
         write(*,*) '       was been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &        "*SOLID SECTION%",ier)
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
         write(*,*) '*ERROR reading *SOLID SECTION: element set ',elset
         write(*,*) '  has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &        "*SOLID SECTION%",ier)
         return
      endif
!
!     assigning the elements of the set the appropriate material
!     and orientation number
!
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            k=ialset(j)
            if((lakon(k)(1:1).eq.'B').or.
     &         (lakon(k)(1:1).eq.'S')) then
               write(*,*) 
     &          '*ERROR reading *SOLID SECTION: *SOLID SECTION can'
               write(*,*) '       not be used for beam or shell elements
     &'
               write(*,*) '       Faulty element: ',k
               ier=1
               return
            endif
            ielmat(1,k)=imaterial
            if(ielorien(1,k).lt.0) then
!
!              an orientation based on a distribution has the
!              same name as the distribution. However, to a
!              distribution which is not used in any orientation
!              definition no local coordinate system has been
!              assigned (i.e. orab(7,..)=0.d0)
!
               if((orname(-ielorien(1,k)).eq.orientation).and.
     &              (orab(7,-ielorien(1,k)).ne.0.d0)) then
                  ielorien(1,k)=-ielorien(1,k)
               else
                  ielorien(1,k)=iorientation
               endif
            else
               ielorien(1,k)=iorientation
            endif
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               if((lakon(k)(1:1).eq.'B').or.
     &              (lakon(k)(1:1).eq.'S')) then
                  write(*,*) '*ERROR reading *SOLID SECTION: *SOLID SECT
     &ION can'
                  write(*,*) '       not be used for beam or shell eleme
     &nts'
                  write(*,*) '       Faulty element: ',k
                  ier=1
                  return
               endif
               ielmat(1,k)=imaterial
               if(ielorien(1,k).lt.0) then
!
!                 an orientation based on a distribution has the
!                 same name as the distribution. However, to a
!                 distribution which is not used in any orientation
!                 definition no local coordinate system has been
!                 assigned (i.e. orab(7,..)=0.d0)
!     
                  if((orname(-ielorien(1,k)).eq.orientation).and.
     &                 (orab(7,-ielorien(1,k)).ne.0.d0)) then
                     ielorien(1,k)=-ielorien(1,k)
                  else
                     ielorien(1,k)=iorientation
                  endif
               else
                  ielorien(1,k)=iorientation
               endif
            enddo
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if(istat.lt.0) return
!
!     assigning a thickness to plane stress/strain elements and an angle to
!     axisymmetric elements
!
      if(key.eq.0) then
!
!        second line with thickness given
!        
         read(textpart(1)(1:20),'(f20.0)',iostat=istat) thickness
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*SOLID SECTION%",ier)
            return
         endif
!
!        for axial symmetric structures:
!           thickness for axial symmetric elements: 2 degrees
!           thickness for plane stress elements: reduced by 180
!           thickness for plane strain elements: reduced by 180
!
         if(.not.nodalthickness) then
            if(iaxial.eq.180) then
               if(lakon(ialset(istartset(i)))(1:2).eq.'CA') then
                  thickness=datan(1.d0)*8.d0/iaxial
               elseif(lakon(ialset(istartset(i)))(1:3).eq.'CPS') then
                  thickness=thickness/iaxial
               elseif(lakon(ialset(istartset(i)))(1:3).eq.'CPE') then
                  thickness=thickness/iaxial
               endif
            endif
         else
!
!           for those elements for which nodal thickness is activated
!           the thickness is set to -1.d0
!
            thickness=-1.d0
         endif
!
!        assigning the thickness to each node of the corresponding
!        elements (thickness specified)
!
         do j=istartset(i),iendset(i)
            if(ialset(j).gt.0) then
               if((lakon(ialset(j))(1:2).eq.'CP').or.
     &              (lakon(ialset(j))(1:2).eq.'CA')) then
!
!                 plane stress/strain or axisymmetric elements
!
                  indexe=ipkon(ialset(j))
                  read(lakon(ialset(j))(4:4),'(i1)') numnod
                  do l=1,numnod
                     thicke(1,indexe+l)=thickness
                  enddo
               elseif(lakon(ialset(j))(1:1).eq.'T') then
                  ielem=ialset(j)
!
!                 default cross section for trusses is the
!                 rectangular cross section
!
                  lakon(ielem)(8:8)='R'
                  indexe=ipkon(ielem)
                  node1=kon(indexe+1)
                  if(lakon(ielem)(4:4).eq.'2') then
                     node2=kon(indexe+2)
                  else
                     node2=kon(indexe+3)
                  endif
!
!                 determining a vector orthogonal to the truss
!                 element
!
                  do l=1,3
                     xn(l)=co(l,node2)-co(l,node1)
                  enddo
                  if(dabs(xn(1)).gt.0.d0) then
                     p(1)=-xn(3)
                     p(2)=0.d0
                     p(3)=xn(1)
                  elseif(dabs(xn(2)).gt.0.d0) then
                     p(1)=xn(2)
                     p(2)=-xn(1)
                     p(3)=0.d0
                  else
                     p(1)=0.d0
                     p(2)=xn(3)
                     p(3)=-xn(2)
                  endif
                  dd=dsqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3))
                  if(dd.lt.1.d-10) then
                     write(*,*) 
     &                    '*ERROR reading *SOLID SECTION: normal'
                     write(*,*) '       in direction 1 has zero size'
                     ier=1
                     return
                  endif
                  do l=1,3
                     p(l)=p(l)/dd
                  enddo
                  do l=1,3
                     thicke(1,indexe+l)=dsqrt(thickness)
                     thicke(2,indexe+l)=dsqrt(thickness)
                     indexx=ixfree
                     do m=1,3
                        xnor(indexx+m)=p(m)
                     enddo
                     ixfree=ixfree+6
                     iponor(1,indexe+l)=indexx
                  enddo
               endif
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  if((lakon(k)(1:2).eq.'CP').or.
     &                 (lakon(k)(1:2).eq.'CA')) then
                     indexe=ipkon(k)
                     read(lakon(k)(4:4),'(i1)') numnod
                     do l=1,numnod
                        thicke(1,indexe+l)=thickness
                     enddo
                  elseif(lakon(k)(1:1).eq.'T') then
!     
!                    default cross section for trusses is the
!                    rectangular cross section
!     
                     lakon(k)(8:8)='R'
                     indexe=ipkon(k)
                     node1=kon(indexe+1)
                     if(lakon(k)(4:4).eq.'2') then
                        node2=kon(indexe+2)
                     else
                        node2=kon(indexe+3)
                     endif
!
!                    determining a vector orthogonal to the truss
!                    element
!
                     do l=1,3
                        xn(l)=co(l,node2)-co(l,node1)
                     enddo
                     if(dabs(xn(1)).gt.0.d0) then
                        p(1)=-xn(3)
                        p(2)=0.d0
                        p(3)=xn(1)
                     elseif(dabs(xn(2)).gt.0.d0) then
                        p(1)=xn(2)
                        p(2)=-xn(1)
                        p(3)=0.d0
                     else
                        p(1)=0.d0
                        p(2)=xn(3)
                        p(3)=-xn(2)
                     endif
                     dd=dsqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3))
                     if(dd.lt.1.d-10) then
                        write(*,*) 
     &                       '*ERROR reading *SOLID SECTION: normal'
                        write(*,*) '       in direction 1 has zero size'
                        ier=1
                        return
                     endif
                     do l=1,3
                        p(l)=p(l)/dd
                     enddo
                     do l=1,3
                        thicke(1,indexe+l)=dsqrt(thickness)
                        thicke(2,indexe+l)=dsqrt(thickness)
                        indexx=ixfree
                        do m=1,3
                           xnor(indexx+m)=p(m)
                        enddo
                        ixfree=ixfree+6
                        iponor(1,indexe+l)=indexx
                     enddo
                  endif
               enddo
            endif
         enddo
      else
!
!        no second line (no thickness given)
!         
!        assigning the thickness to each node of the corresponding
!        elements (thickness not specified: only axisymmetric elements
!        or plane stress/strain elements with nodal thickness)
!
         thickness=datan(1.d0)*8.d0/iaxial
         do j=istartset(i),iendset(i)
            if(ialset(j).gt.0) then
               if(lakon(ialset(j))(1:2).eq.'CA') then
!
!                 axisymmetric elements
!
                  indexe=ipkon(ialset(j))
                  read(lakon(ialset(j))(4:4),'(i1)') numnod
                  do l=1,numnod
                     thicke(1,indexe+l)=thickness
                  enddo
               elseif((lakon(ialset(j))(1:3).eq.'CPS').or.
     &                 (lakon(ialset(j))(1:3).eq.'CPE')) then
                  if(nodalthickness) then
                     indexe=ipkon(ialset(j))
                     read(lakon(ialset(j))(4:4),'(i1)') numnod
                     do l=1,numnod
                        thicke(1,indexe+l)=-1.d0
                     enddo
                  endif
               endif
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  if(lakon(k)(1:2).eq.'CA') then
!
!                 axisymmetric elements
!
                     indexe=ipkon(k)
                     read(lakon(k)(4:4),'(i1)') numnod
                     do l=1,numnod
                        thicke(1,indexe+l)=thickness
                     enddo
                  elseif((lakon(k)(1:3).eq.'CPS').or.
     &                    (lakon(k)(1:3).eq.'CPE')) then
                     if(nodalthickness) then
                        indexe=ipkon(k)
                        read(lakon(k)(4:4),'(i1)') numnod
                        do l=1,numnod
                           thicke(1,indexe+l)=-1.d0
                        enddo
                     endif
                  endif
               enddo
            endif
         enddo
      endif
!     
!        defining cyclic symmetric conditions for axisymmetric
!        elements (needed for cavity radiation)
!
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            if(lakon(ialset(j))(1:2).eq.'CA') then
               if(mcs.gt.1) then
                  write(*,*) '*ERROR reading *SOLID SECTION: '
                  write(*,*) '       axisymmetric elements cannot be
     &combined with cyclic symmetry'
                  ier=1
                  return
               elseif(mcs.eq.1) then
                  if(int(cs(1,1)).ne.int(2.d0*pi/thickness+0.5d0)) 
     &                 then
                     write(*,*) '*ERROR reading *SOLID SECTION: '
                     write(*,*) '       it is not allowed to define t
     &wo different'
                     write(*,*) '       angles for an axisymmetric st
     &ructure'
                     ier=1
                     return
                  else
                     exit
                  endif
               endif
               mcs=1
               cs(1,1)=2.d0*pi/thickness+0.5d0
               cs(2,1)=-0.5d0
               cs(3,1)=-0.5d0
               cs(5,1)=1.5d0
               do k=6,9
                  cs(k,1)=0.d0
               enddo
               cs(10,1)=1.d0
               cs(11,1)=0.d0
               cs(12,1)=-1.d0
               cs(14,1)=0.5d0
               cs(15,1)=dcos(thickness)
               cs(16,1)=dsin(thickness)
               exit
            endif
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               if(lakon(k)(1:2).eq.'CA') then
                  if(mcs.gt.1) then
                     write(*,*) '*ERROR reading *SOLID SECTION: '
                     write(*,*) '       axisymmetric elements cannot 
     &be combined with cyclic symmetry'
                     ier=1
                     return
                  elseif(mcs.eq.1) then
                     if(int(cs(1,1)).ne.int(2.d0*pi/thickness+0.5d0)) 
     &                    then
                        write(*,*) '*ERROR reading *SOLID SECTION: '
                        write(*,*) '       it is not allowed to defin
     &e two different'
                        write(*,*) '       angles for an axisymmetric
     &structure'
                        ier=1
                        return
                     else
                        exit
                     endif
                  endif
                  mcs=1
                  cs(1,1)=2.d0*pi/thickness+0.5d0
                  cs(2,1)=-0.5d0
                  cs(3,1)=-0.5d0
                  cs(5,1)=1.5d0
                  do k=6,9
                     cs(k,1)=0.d0
                  enddo
                  cs(10,1)=1.d0
                  cs(11,1)=0.d0
                  cs(12,1)=-1.d0
                  cs(14,1)=0.5d0
                  cs(15,1)=dcos(thickness)
                  cs(16,1)=dsin(thickness)
                  exit
               endif
            enddo
         endif
      enddo
!     
      if(key.eq.0) then
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
      endif
!     
      return
      end

