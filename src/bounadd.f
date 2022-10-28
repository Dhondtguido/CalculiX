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
      subroutine bounadd(node,is,ie,val,nodeboun,ndirboun,xboun,
     &  nboun,nboun_,iamboun,iamplitude,nam,ipompc,nodempc,
     &  coefmpc,nmpc,nmpc_,mpcfree,inotr,trab,
     &  ntrans,ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,type,
     &  typeboun,nmethod,iperturb,fixed,vold,nodetrue,mi,label)
!
!     adds a boundary condition to the data base
!
      implicit none
!
      logical fixed,rottracoupling
!
      character*1 type,typeboun(*)
      character*20 labmpc(*),label
!
      integer nodeboun(*),ndirboun(*),node,is,ie,nboun,nboun_,i,j,
     &  iamboun(*),iamplitude,nam,ipompc(*),nodempc(3,*),nmpc,nmpc_,
     &  mpcfree,inotr(2,*),ntrans,ikboun(*),ilboun(*),ikmpc(*),
     &  ilmpc(*),itr,idof,newnode,number,id,idofnew,idnew,nk,nk_,
     &  mpcfreenew,nmethod,iperturb(*),ii,nodetrue,mi(*),three,kflag,
     &  iy(3),inumber,irotnode(11),irotdof(11)
!
      real*8 xboun(*),val,coefmpc(*),trab(7,*),a(3,3),co(3,*),
     &  vold(0:mi(2),*),dx(3)
!
      if(ntrans.le.0) then
         itr=0
      elseif(inotr(1,node).eq.0) then
         itr=0
      else
         itr=inotr(1,node)
      endif
!
!     checking for boundary conditions on rotational dofs of
!     distributing couplings 
!
      rottracoupling=.false.
      if((ie.ge.4).and.(ie.le.6)) then
!
!        rotational dof
!
         do ii=is,ie
            irotnode(ii)=node
            if(ii.gt.3) then
c               idof=8*(node-1)+ii+1
               idof=8*(node-1)+ii
               call nident(ikmpc,idof,nmpc,id)
               if(id.gt.0) then
                  if(ikmpc(id).eq.idof) then
                     if(labmpc(ilmpc(id))(1:14).eq.'ROTTRACOUPLING')then
                        rottracoupling=.true.
                        irotnode(ii)=
     &                     nodempc(1,nodempc(3,ipompc(ilmpc(id))))
                        irotdof(ii)=
     &                     nodempc(2,nodempc(3,ipompc(ilmpc(id))))
                        itr=0
                     endif
                  endif
               endif
            endif
         enddo
      endif
!
      loop: do ii=is,ie
!
!     change: transformations on rotations are taken into account
!     by the normal of the mean rotation MPC, not by expanding the
!     MPC in Carthesian coordinates
!
         if((itr.eq.0).or.(ii.eq.0).or.(ii.gt.3)) then
!
!        no transformation applies: simple SPC
!
            if(rottracoupling) then
               node=irotnode(ii)
               if(ii.gt.3) then
                  i=irotdof(ii)
               else
                  i=ii
               endif
            else
c               if(ii.le.3) then
               if(ii.le.6) then
                  i=ii
c               elseif(ii.eq.4) then
c                  i=5
c               elseif(ii.eq.5) then
c                  i=6
c               elseif(ii.eq.6) then
c                  i=7
               elseif(ii.eq.8) then
                  i=4
               elseif(ii.eq.11) then
                  i=0
               else
                  write(*,*) '*ERROR in bounadd: unknown DOF: ',
     &                 ii
                  call exit(201)
               endif
            endif
!
            if((fixed).and.(i.lt.5)) then
               val=vold(i,nodetrue)
            elseif(fixed) then
               write(*,*) '*ERROR in bounadd: parameter FIXED cannot'
               write(*,*) '       be used for rotations'
               call exit(201)
            endif
            idof=8*(node-1)+i
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
            if(nboun.gt.nboun_) then
               write(*,*) '*ERROR in bounadd: increase nboun_'
               call exit(201)
            endif
            if((nmethod.eq.4).and.(iperturb(1).le.1)) then
               write(*,*) '*ERROR in bounadd: in a modal dynamic step'
               write(*,*) '       new SPCs are not allowed'
               call exit(201)
            endif
            nodeboun(nboun)=node
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
            call transformatrix(trab(1,itr),co(1,node),a)
c            if(ii.le.3) then
            if(ii.le.6) then
               i=ii
c            elseif(ii.eq.4) then
c               i=5
c            elseif(ii.eq.5) then
c               i=6
c            elseif(ii.eq.6) then
c               i=7
            elseif(ii.eq.8) then
               i=4
            elseif(ii.eq.11) then
               i=0
            else
               write(*,*) '*ERROR in bounadd: unknown DOF: ',
     &              ii
               call exit(201)
            endif
            if((fixed).and.(i.lt.5)) then
               val=vold(i,nodetrue)
            elseif(fixed) then
               write(*,*) '*ERROR in bounadd: parameter FIXED cannot'
               write(*,*) '       be used for rotations'
               call exit(201)
            endif
            if(inotr(2,node).ne.0) then
               newnode=inotr(2,node)
               idofnew=8*(newnode-1)+i
               call nident(ikboun,idofnew,nboun,idnew)
               if(idnew.gt.0) then
                  if(ikboun(idnew).eq.idofnew) then
                     j=ilboun(idnew)
c
                     if(typeboun(j).ne.type) cycle
c
                     xboun(j)=val
c                     typeboun(j)=type
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
               inotr(2,node)=newnode
               idofnew=8*(newnode-1)+i
               idnew=nboun
!
!              copying the initial conditions from node into newnode
!
               do j=0,mi(2)
                  vold(j,newnode)=vold(j,node)
               enddo
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
               idof=8*(node-1)+number
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
               labmpc(nmpc)=label
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
!           check whether a dependent term was found; if none was
!           found this can be due to the fact that:
!           - all dofs were used by other MPC's
!           - the MPC coefficients were too small
!           - or a combination of both
!
            if(inumber.gt.3) then
               write(*,*) '*ERROR in bounadd'
               write(*,*) '       SPC in node',node
               write(*,*) '       and local direction',ii
               write(*,*) '       cannot be applied: all'
               write(*,*) '       degrees of freedom have'
               write(*,*) '       been used by other MPCs'
               write(*,*) '       or the coefficient is'
               write(*,*) '       too small'
               call exit(201)
            endif
!
            inumber=inumber-1
            do j=1,3
               inumber=inumber+1
               if(inumber.gt.3) inumber=1
               number=iy(inumber)
c               if(dabs(a(number,i)).lt.1.d-5) cycle
               if(dabs(a(number,i)).lt.1.d-30) cycle
               nodempc(1,mpcfree)=node
               nodempc(2,mpcfree)=number
               coefmpc(mpcfree)=a(number,i)
               mpcfree=nodempc(3,mpcfree)
               if(mpcfree.eq.0) then
                  write(*,*) '*ERROR in bounadd: increase memmpc_'
                  call exit(201)
               endif
            enddo
            nodempc(1,mpcfree)=newnode
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
c         enddo
         endif
      enddo loop
!
      return
      end

