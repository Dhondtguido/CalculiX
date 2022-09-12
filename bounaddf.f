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
      subroutine bounaddf(iface,is,ie,val,nodeboun,ndirboun,xboun,
     &  nboun,nboun_,iamboun,iamplitude,nam,ipompc,nodempc,
     &  coefmpc,nmpc,nmpc_,mpcfree,trab,
     &  ntrans,ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,type,
     &  typeboun,nmethod,iperturb,vold,mi,
     &  nelemload,sideload,xload,nload,lakon,ipkon,kon)
!
!     adds a facial boundary condition to the data base
!
!     for facial boundary conditions: dof=-[8*(face-1)+idof]
!        (where face = 10*element+local_face_number)
!     for nodal boundary conditions: dof=8*(node-1)+idof
!
      implicit none
!
      character*1 type,typeboun(*)
      character*8 lakon(*)
      character*20 labmpc(*),label,sideload(*)
!
      integer nodeboun(*),ndirboun(*),is,ie,nboun,nboun_,i,j,
     &  iamboun(*),iamplitude,nam,ipompc(*),nodempc(3,*),nmpc,nmpc_,
     &  mpcfree,ntrans,ikboun(*),ilboun(*),ikmpc(*),ipkon(*),indexe,
     &  ilmpc(*),itr,idof,newnode,number,id,idofnew,idnew,nk,nk_,
     &  mpcfreenew,nmethod,iperturb(*),ii,mi(*),three,kflag,
     &  iy(3),inumber,iface,nload,nelemload(2,*),nopes,kon(*),nope,
     &  nelem,loadid,ifacel,ifaceq(8,6),ifacet(6,4),ifacew(8,5)
!
      real*8 xboun(*),val,coefmpc(*),trab(7,*),a(3,3),co(3,*),cg(3),
     &  vold(0:mi(2),*),dx(3),xload(2,*)
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
!
      if(ntrans.le.0) then
         itr=0
      else
         nelem=int(iface/10.d0)
         ifacel=iface-10*nelem
         label(1:20)='T                   '
         write(label(2:2),'(i1)') ifacel
         call identifytransform(nelem,label,nelemload,sideload,
     &        nload,loadid)
         if(loadid.eq.0) then
            itr=0
         else
            itr=nelemload(2,loadid)
         endif
      endif
!
      loop: do ii=is,ie
         if((itr.eq.0).or.(ii.eq.0).or.(ii.eq.11).or.(ii.eq.8)) then
!
!        no transformation applies: simple SPC
!
            if(ii.le.3) then
               i=ii
            elseif(ii.eq.4) then
               write(*,*) '*ERROR in bounaddf: a boundary condition'
               write(*,*) '       on DOF 4 is not allowed'
               call exit(201)
            elseif(ii.eq.5) then
               write(*,*) '*ERROR in bounaddf: a boundary condition'
               write(*,*) '       on DOF 5 is not allowed'
               call exit(201)
            elseif(ii.eq.6) then
               write(*,*) '*ERROR in bounaddf: a boundary condition'
               write(*,*) '       on DOF 6 is not allowed'
               call exit(201)
            elseif(ii.eq.8) then
               i=4
            elseif(ii.eq.11) then
               i=0
            else
               write(*,*) '*ERROR in bounadd: unknown DOF: ',
     &              ii
               call exit(201)
            endif
            idof=-(8*(iface-1)+i)
            call nident(ikboun,idof,nboun,id)
            if(id.gt.0) then
               if(ikboun(id).eq.idof) then
                  j=ilboun(id)
                  if(typeboun(j).ne.type) cycle loop
                  xboun(j)=val
                  if(nam.gt.0) iamboun(j)=iamplitude
                  cycle loop
               endif
            endif
            nboun=nboun+1
c            write(*,*) 'bounaddf boun',iface,i,val,idof
            if(nboun.gt.nboun_) then
               write(*,*) '*ERROR in bounadd: increase nboun_'
               call exit(201)
            endif
            if((nmethod.eq.4).and.(iperturb(1).le.1)) then
               write(*,*) '*ERROR in bounadd: in a modal dynamic step'
               write(*,*) '       new SPCs are not allowed'
               call exit(201)
            endif
            nodeboun(nboun)=iface
            ndirboun(nboun)=i
            xboun(nboun)=val
            typeboun(nboun)=type
            if(nam.gt.0) iamboun(nboun)=iamplitude
!
!           updating ikboun and ilboun
!            
            do j=nboun,id+2,-1
               ikboun(j)=ikboun(j-1)
               ilboun(j)=ilboun(j-1)
            enddo
            ikboun(id+1)=idof
            ilboun(id+1)=nboun
         else
!
!        transformation applies: SPC is MPC in global carthesian
!        coordinates
!
!        number of nodes belonging to the face
!
            if(lakon(nelem)(4:4).eq.'8') then
               nope=8
               nopes=4
            elseif(lakon(nelem)(4:4).eq.'4') then
               nope=4
               nopes=3
            elseif(lakon(nelem)(4:4).eq.'6') then
               nope=6
               if(ifacel.le.2) then
                  nopes=3
               else
                  nopes=4
               endif
            endif
!
!        determining the center of gravity 
!
            do j=1,3
               cg(j)=0.d0
            enddo
!
            indexe=ipkon(nelem)
            if(nope.eq.8) then
               do i=1,nopes
                  do j=1,3
                     cg(j)=cg(j)+co(j,kon(indexe+ifaceq(i,ifacel)))
                  enddo
               enddo
            elseif(nope.eq.4) then
               do i=1,nopes
                  do j=1,3
                     cg(j)=cg(j)+co(j,kon(indexe+ifacet(i,ifacel)))
                  enddo
               enddo
            else
               do i=1,nopes
                  do j=1,3
                     cg(j)=cg(j)+co(j,kon(indexe+ifacew(i,ifacel)))
                  enddo
               enddo
            endif
            do j=1,3
               cg(j)=cg(j)/nopes
            enddo
!
!        determining the transformation coefficients at the center
!        of gravity
!
            call transformatrix(trab(1,itr),cg,a)
            if(ii.le.3) then
               i=ii
            elseif(ii.eq.4) then
               write(*,*) '*ERROR in bounaddf: a boundary condition'
               write(*,*) '       on DOF 4 is not allowed'
               call exit(201)
            elseif(ii.eq.5) then
               write(*,*) '*ERROR in bounaddf: a boundary condition'
               write(*,*) '       on DOF 5 is not allowed'
               call exit(201)
            elseif(ii.eq.6) then
               write(*,*) '*ERROR in bounaddf: a boundary condition'
               write(*,*) '       on DOF 6 is not allowed'
               call exit(201)
            elseif(ii.eq.8) then
               i=4
            elseif(ii.eq.11) then
               i=0
            else
               write(*,*) '*ERROR in bounadd: unknown DOF: ',
     &              ii
               call exit(201)
            endif
            if(int(xload(1,loadid)).ne.0) then
               newnode=int(xload(1,loadid))
               idofnew=8*(newnode-1)+i
               call nident(ikboun,idofnew,nboun,idnew)
               if(idnew.gt.0) then
                  if(ikboun(idnew).eq.idofnew) then
                     j=ilboun(idnew)
                     if(typeboun(j).ne.type) cycle
                     xboun(j)=val
                     if(nam.gt.0) iamboun(j)=iamplitude
                     cycle
                  endif
               endif
            else
!
!              new node is generated for the inhomogeneous MPC term
!
               if((nmethod.eq.4).and.(iperturb(1).le.1)) then
                  write(*,*)'*ERROR in bounadd: in a modal dynamic step'
                  write(*,*) '       new SPCs are not allowed'
                  call exit(201)
               endif
               nk=nk+1
               if(nk.gt.nk_) then
                  write(*,*) '*ERROR in bounadd: increase nk_'
                  call exit(201)
               endif
               newnode=nk
               xload(1,loadid)=newnode+0.5d0
               idofnew=8*(newnode-1)+i
               idnew=nboun
            endif
!
!           new mpc
!
            iy(1)=1
            iy(2)=2
            iy(3)=3
            dx(1)=dabs(a(1,i))
            dx(2)=dabs(a(2,i))
            dx(3)=dabs(a(3,i))
            three=3
            kflag=-2
            call dsort(dx,iy,three,kflag)
            do inumber=1,3
               number=iy(inumber)
               idof=-(8*(iface-1)+number)
               call nident(ikmpc,idof,nmpc,id)
               if(id.ne.0) then
                  if(ikmpc(id).eq.idof) cycle
               endif
               if(dabs(a(number,i)).lt.1.d-5) cycle
               nmpc=nmpc+1
               if(nmpc.gt.nmpc_) then
                  write(*,*) '*ERROR in bounadd: increase nmpc_'
                  call exit(201)
               endif
               labmpc(nmpc)='FLUIDSPC            '
c               write(*,*) nmpc,labmpc(nmpc),'bounaddf'
               ipompc(nmpc)=mpcfree
               do j=nmpc,id+2,-1
                  ikmpc(j)=ikmpc(j-1)
                  ilmpc(j)=ilmpc(j-1)
               enddo
               ikmpc(id+1)=idof
               ilmpc(id+1)=nmpc
               exit
            enddo
!
            inumber=inumber-1
            do j=1,3
               inumber=inumber+1
               if(inumber.gt.3) inumber=1
               number=iy(inumber)
               if(dabs(a(number,i)).lt.1.d-30) cycle
               nodempc(1,mpcfree)=iface
               nodempc(2,mpcfree)=number
               coefmpc(mpcfree)=a(number,i)
               mpcfree=nodempc(3,mpcfree)
               if(mpcfree.eq.0) then
                  write(*,*) '*ERROR in bounadd: increase memmpc_'
                  call exit(201)
               endif
            enddo
!
!           storage of the boundary condition (faster than
!           storage of the new node); the negative of the 
!           condition number is stored as tag for a SPC
!
            nodempc(1,mpcfree)=-(nboun+1)
            nodempc(2,mpcfree)=i
            coefmpc(mpcfree)=-1.d0
            mpcfreenew=nodempc(3,mpcfree)
            if(mpcfreenew.eq.0) then
               write(*,*) '*ERROR in bounadd: increase nmpc_'
               call exit(201)
            endif
            nodempc(3,mpcfree)=0
            mpcfree=mpcfreenew
!
!           nonhomogeneous term
!
            nboun=nboun+1
c            write(*,*) 'bounaddf inhom',newnode,i,val,idofnew
            if(nboun.gt.nboun_) then
               write(*,*) '*ERROR in bounadd: increase nboun_'
               call exit(201)
            endif
            nodeboun(nboun)=newnode
            ndirboun(nboun)=i
            xboun(nboun)=val
            typeboun(nboun)=type
            if(nam.gt.0) iamboun(nboun)=iamplitude
!
!           updating ikboun and ilboun
!            
            do j=nboun,idnew+2,-1
               ikboun(j)=ikboun(j-1)
               ilboun(j)=ilboun(j-1)
            enddo
            ikboun(idnew+1)=idofnew
            ilboun(idnew+1)=nboun
!            
         endif
      enddo loop
!
      return
      end

