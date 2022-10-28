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
      subroutine usermpc(ipompc,nodempc,coefmpc,
     &  labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,nk,nk_,nodeboun,ndirboun,
     &  ikboun,ilboun,nboun,nboun_,nnodes,node,co,label,typeboun,
     &  iperturb,noderef,idirref,xboun)
!
!     initializes mpc fields for a user MPC
!
      implicit none
!
      character*1 typeboun(*)
      character*20 labmpc(*),label
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,nk,nk_,
     &  ikmpc(*),idir,indexmax,indexold(3),indexnew(3),nodeold,nodemax,
     &  ilmpc(*),node,id,mpcfreeold,idof,l,nodeboun(*),iperturb(*),
     &  ndirboun(*),ikboun(*),ilboun(*),nboun,nboun_,nnodes,nodevector,
     &  index,index1,node1,i,j,nkn,idirold,idirmax,noderef,idirref,
     &  nendnode
!
      real*8 coefmpc(*),co(3,*),aa(3),dd,cgx(3),pi(3),c1,c4,c9,
     &  c10,amax,xcoef,transcoef(3),xboun(*),stdev
!
      save nodevector
!
      if(node.ne.0) then
         if(nnodes.eq.1) then
            if((label(1:7).ne.'MEANROT').and.
     &           (label(1:1).ne.'1')) then
!
!              define a new MPC
!              default for the dependent DOF direction is 1
!     
               idof=8*(node-1)+1
!     
               call nident(ikmpc,idof,nmpc,id)
               if(id.gt.0) then
                  if(ikmpc(id).eq.idof) then
                     write(*,*)'*WARNING in usermpc: DOF for node ',node
                     write(*,*) '         in direction 1 has been used'
                     write(*,*) 
     &                   '         on the dependent side of another'
                     write(*,*) '         MPC. ',label
                     write(*,*) '         constraint cannot be applied'
                     write(*,*)
                     return
                  endif
               endif
            endif
!
            nmpc=nmpc+1
            if(nmpc.gt.nmpc_) then
               write(*,*) '*ERROR in usermpc: increase nmpc_'
               call exit(201)
            endif
!
            ipompc(nmpc)=mpcfree
            labmpc(nmpc)=label
!
!           for a mean rotation MPC the dependent degree of freedom
!           is determined at the end of this routine
!
            if((label(1:7).ne.'MEANROT').and.
     &           (label(1:1).ne.'1')) then
               do l=nmpc,id+2,-1
                  ikmpc(l)=ikmpc(l-1)
                  ilmpc(l)=ilmpc(l-1)
               enddo
               ikmpc(id+1)=idof
               ilmpc(id+1)=nmpc
            endif
         endif
!
!        general case: add a term to the MPC
!
         nodempc(1,mpcfree)=node
!
!        nodevector: additional node such that:
!         - the coordinates of this node are the axis direction
!         - the 1st DOF is reserved for the mean rotation value
!
         if((labmpc(nmpc)(1:7).eq.'MEANROT').or.
     &        (labmpc(nmpc)(1:1).eq.'1')) then
            nodevector=node
            labmpc(nmpc)(1:7)='MEANROT'
         endif
!
         if(nnodes.eq.1) then
            nodempc(2,mpcfree)=1
         else
            nodempc(2,mpcfree)=0
         endif
         coefmpc(mpcfree)=0.d0
         mpcfree=nodempc(3,mpcfree)
      else
!
!        MPC definition finished: add a nonhomogeneous term
!        (the 2nd degree of freedom of an extra node is taken; 
!         the 1st degree of freedom can be taken by the user
!         for the angle dof)
!
         nk=nk+1
         if(nk.gt.nk_) then
            write(*,*) '*ERROR in usermpc: increase nk_'
            call exit(201)
         endif
!
         nodempc(1,mpcfree)=nk
         nodempc(2,mpcfree)=2
!
         coefmpc(mpcfree)=1.d0
!
         mpcfreeold=mpcfree
         mpcfree=nodempc(3,mpcfree)
         nodempc(3,mpcfreeold)=0
         idof=8*(nk-1)+2
         call nident(ikboun,idof,nboun,id)
         nboun=nboun+1
         if(nboun.gt.nboun_) then
            write(*,*) '*ERROR in usermpc: increase nboun_'
            call exit(201)
         endif
         nodeboun(nboun)=nk
         ndirboun(nboun)=2
         xboun(nboun)=0.d0
         typeboun(nboun)='U'
         do l=nboun,id+2,-1
            ikboun(l)=ikboun(l-1)
            ilboun(l)=ilboun(l-1)
         enddo
         ikboun(id+1)=idof
         ilboun(id+1)=nboun
!
!        calculating the MPC coefficients for linear applications
!
         if((labmpc(nmpc)(1:7).eq.'MEANROT').or.
     &        (labmpc(nmpc)(1:1).eq.'1')) then
            nkn=(nnodes-1)/3
            if(3*nkn.ne.nnodes-1) then
               write(*,*)
     &              '*ERROR in usermpc: MPC has wrong number of terms'
               call exit(201)
            endif
!
!           normal along the rotation axis
!
            dd=0.d0
            do i=1,3
               aa(i)=co(i,nodevector)
               dd=dd+aa(i)**2
            enddo
            dd=dsqrt(dd)
            if(dd.lt.1.d-10) then
               write(*,*) 
     &           '*ERROR in usermpc: rotation vector has zero length'
               call exit(201)
            endif
            do i=1,3
               aa(i)=aa(i)/dd
            enddo
!     
!     finding the center of gravity of the position and the
!     displacements of the nodes involved in the MPC
!
            do i=1,3
               cgx(i)=0.d0
            enddo
!
            index=ipompc(nmpc)
            do
               node=nodempc(1,index)
               if(node.eq.nodevector) exit
               do j=1,3
                  cgx(j)=cgx(j)+co(j,node)
               enddo
               index=nodempc(3,nodempc(3,nodempc(3,index)))
            enddo
!
            do i=1,3
               cgx(i)=cgx(i)/nkn
            enddo
!
!           calculating a standard deviation; this quantity will
!           serve as a limit for checking the closeness of individual
!           nodes to the center of gravity
!     
            stdev=0.d0
!
            index=ipompc(nmpc)
            do
               node=nodempc(1,index)
               if(node.eq.nodevector) exit
               do j=1,3
                  stdev=stdev+(co(j,node)-cgx(j))**2
               enddo
               index=nodempc(3,nodempc(3,nodempc(3,index)))
            enddo
            stdev=stdev/nkn
!
!           calculating the derivatives
!
!           loop over all nodes belonging to the mean rotation MPC
!
            index=ipompc(nmpc)
            do
               node=nodempc(1,index)
               if(node.eq.nodevector) exit
!
!              relative positions
!               
               do j=1,3
                  pi(j)=co(j,node)-cgx(j)
               enddo
!
!              projection on a plane orthogonal to the rotation vector
!
               c1=pi(1)*aa(1)+pi(2)*aa(2)+pi(3)*aa(3)
               do j=1,3
                  pi(j)=pi(j)-c1*aa(j)
               enddo
!
               c1=pi(1)*pi(1)+pi(2)*pi(2)+pi(3)*pi(3)
c               if(c1.lt.1.d-20) then
               if(c1.lt.stdev*1.d-10) then
                  if(label(8:9).ne.'BS') then
                     write(*,*) '*WARNING in usermpc: node ',node
                     write(*,*) '         is very close to the '
                     write(*,*) '         rotation axis through the'
                     write(*,*) '         center of gravity of'
                     write(*,*) '         the nodal cloud in a'
                     write(*,*) '         mean rotation MPC.'
                     write(*,*) '         This node is not taken'
                     write(*,*) '         into account in the MPC'
                     write(*,*)
                  endif
                  index=nodempc(3,nodempc(3,nodempc(3,index)))
                  cycle
               endif
!
               do j=1,3
                  if(j.eq.1) then
                     c4=aa(2)*pi(3)-aa(3)*pi(2)
                  elseif(j.eq.2) then
                     c4=aa(3)*pi(1)-aa(1)*pi(3)
                  else
                     c4=aa(1)*pi(2)-aa(2)*pi(1)
                  endif
                  c9=c4/c1
!
                  index1=ipompc(nmpc)
                  do
                     node1=nodempc(1,index1)
                     if(node1.eq.nodevector) exit
                     if(node1.eq.node) then
                        c10=c9*(1.d0-1.d0/real(nkn))
                     else
                        c10=-c9/real(nkn)
                     endif
                     if(j.eq.1) then
                        coefmpc(index1)=coefmpc(index1)+c10
                     elseif(j.eq.2) then
                        coefmpc(nodempc(3,index1))=
     &                  coefmpc(nodempc(3,index1))+c10
                     else
                        coefmpc(nodempc(3,nodempc(3,index1)))=
     &                  coefmpc(nodempc(3,nodempc(3,index1)))+c10
                     endif
                     index1=nodempc(3,nodempc(3,nodempc(3,index1)))
                  enddo
               enddo
               index=nodempc(3,nodempc(3,nodempc(3,index)))
            enddo
            coefmpc(index)=-nkn
!
!     assigning the degrees of freedom
!
            j=0
            index=ipompc(nmpc)
            do
               j=j+1
               if(j.gt.3) j=1
               nodempc(2,index)=j
               index=nodempc(3,index)
               if(nodempc(3,index).eq.0) exit
            enddo
!
!           look for biggest coefficient of all but the last
!           regular node. The general form of the MPC is:
!           regular nodes, angle dof and inhomogeneous dof
!           The last regular node is exempted so that it can
!           be used for translational dofs in the application of
!           moments/rotations to beams/shells.
!
!           check for the coefficients of the last regular node
!           this node is used for translational dofs in shell
!           and beam applications
!
            if(label(8:9).eq.'BS') then
               index=ipompc(nmpc)
               do i=1,nnodes-4
                  index=nodempc(3,index)
               enddo
               do i=1,3
                  transcoef(i)=coefmpc(index)
                  index=nodempc(3,index)
               enddo
            endif
!
            indexmax=0
            index=ipompc(nmpc)
            amax=1.d-5
!
!           loop over all nodes - angle node - last regular node
!
            if(label(8:9).eq.'BS') then
               nendnode=nnodes-4
            else
               nendnode=nnodes-1
            endif
            do i=1,nendnode
               if(dabs(coefmpc(index)).gt.amax) then
                  idir=nodempc(2,index)
!
!                 dependent node of MPC should not have the same
!                 sign as the corresponding dof of the translational
!                 degrees of freedom: avoids the occurence of a
!                 zero coefficient of the dependent term if both
!                 rotational and translational dofs are suppressed
!
                  if(label(8:9).eq.'BS') then
                     if(coefmpc(index)*transcoef(idir).gt.0) then
                        index=nodempc(3,index)
                        cycle
                     endif
                  endif
!
                  node=nodempc(1,index)
                  idof=8*(node-1)+idir
!     
!                 check whether the node was used in another MPC
!
                  call nident(ikmpc,idof,nmpc-1,id)
                  if(id.gt.0) then
                     if(ikmpc(id).eq.idof) then
                        index=nodempc(3,index)
                        cycle
                     endif
                  endif
!
!                 check whether the node was used in a SPC
!
                  call nident(ikboun,idof,nboun,id)
                  if(id.gt.0) then
                     if(ikboun(id).eq.idof) then
                        index=nodempc(3,index)
                        cycle
                     endif
                  endif
!
                  amax=dabs(coefmpc(index))
                  indexmax=index
               endif
!
!              after each node has been treated, a check is performed
!              whether indexmax is already nonzero. Taking too many
!              nodes into account may lead to dependent equations if
!              all rotational dofs are constrained.
!
               if((i/3)*3.eq.i) then
                  if(indexmax.ne.0) exit
               endif
!
               index=nodempc(3,index)
            enddo
!
            if(indexmax.eq.0) then
               if(idirref.ne.0) then
                  write(*,*) '*WARNING in usermpc: no MPC is '
                  write(*,*) '         generated for the mean'
                  write(*,*) '         rotation in node ',noderef
                  write(*,*) '         and direction ',idirref+3
                  write(*,*)
               else
                  write(*,*) '*WARNING in usermpc: no MPC is '
                  write(*,*) '         generated for the mean'
                  write(*,*) '         rotation definition starting '
                  write(*,*) '         with node ',noderef
                  write(*,*)
               endif
!
!              removing the MPC
!
               mpcfreeold=mpcfree
               index=ipompc(nmpc)
               ipompc(nmpc)=0
               mpcfree=index
!     
!              removing the MPC from fields nodempc and coefmpc
!     
               do
                  nodempc(1,index)=0
                  nodempc(2,index)=0
                  coefmpc(index)=0
                  if(nodempc(3,index).ne.0) then
                     index=nodempc(3,index)
                  else
                     nodempc(3,index)=mpcfreeold
                     exit
                  endif
               enddo
!     
!              decrementing nmpc
!     
               nmpc=nmpc-1
!
!              token to signify that no MPC was generated
!              (needed in gen3dboun.f)
!
               node=-1
               return
            endif
!
            nodemax=nodempc(1,indexmax)
            idirmax=nodempc(2,indexmax)
            index=ipompc(nmpc)
!
            nodeold=nodempc(1,index)
            idirold=nodempc(2,index)
!
!           exchange the node information in the MPC
!
            if(nodemax.ne.nodeold) then
!
               indexold(1)=index
               indexold(2)=nodempc(3,index)
               indexold(3)=nodempc(3,indexold(2))
!
               index=nodempc(3,indexold(3))
               do
                  if(nodempc(1,index).eq.nodevector) exit
                  if(nodempc(1,index).eq.nodemax) then
                     if(nodempc(2,index).eq.1) then
                        indexnew(1)=index
                     elseif(nodempc(2,index).eq.2) then
                        indexnew(2)=index
                     else
                        indexnew(3)=index
                     endif
                  endif
                  index=nodempc(3,index)
               enddo
!
               do j=1,3
                  node=nodempc(1,indexold(j))
                  idir=nodempc(2,indexold(j))
                  xcoef=coefmpc(indexold(j))
                  nodempc(1,indexold(j))=nodempc(1,indexnew(j))
                  nodempc(2,indexold(j))=nodempc(2,indexnew(j))
                  coefmpc(indexold(j))=coefmpc(indexnew(j))
                  nodempc(1,indexnew(j))=node
                  nodempc(2,indexnew(j))=idir
                  coefmpc(indexnew(j))=xcoef
               enddo
            endif
!
!           exchange the direction information in the MPC
!
            index=ipompc(nmpc)
            if(idirmax.ne.1) then
               indexold(1)=index
               if(idirmax.eq.2) then
                  indexnew(1)=nodempc(3,index)
               else
                  indexnew(1)=nodempc(3,nodempc(3,index))
               endif
!     
               do j=1,1
                  idir=nodempc(2,indexold(j))
                  xcoef=coefmpc(indexold(j))
                  nodempc(2,indexold(j))=nodempc(2,indexnew(j))
                  coefmpc(indexold(j))=coefmpc(indexnew(j))
                  nodempc(2,indexnew(j))=idir
                  coefmpc(indexnew(j))=xcoef
               enddo
            endif
!     
!           determining the dependent dof of the MPC
!     
            idof=8*(nodemax-1)+idirmax
            call nident(ikmpc,idof,nmpc-1,id)
            do l=nmpc,id+2,-1
               ikmpc(l)=ikmpc(l-1)
               ilmpc(l)=ilmpc(l-1)
            enddo
            ikmpc(id+1)=idof
            ilmpc(id+1)=nmpc
!     
         elseif(labmpc(nmpc)(1:4).eq.'DIST') then
            iperturb(2)=1
            write(*,*) '*INFO in usermpc: nonlinear geometric'
            write(*,*) '      effects are turned on'
            write(*,*)
            if(iperturb(1).eq.0) iperturb(1)=2
c         elseif(labmpc(nmpc)(1:3).eq.'GAP') then
c            iperturb(2)=1
c            write(*,*) '*INFO in usermpc: nonlinear geometric'
c            write(*,*) '      effects are turned on'
c            write(*,*)
c            if(iperturb(1).eq.0) iperturb(1)=2
         elseif(labmpc(nmpc)(1:4).eq.'USER') then
            iperturb(2)=1
            write(*,*) '*INFO in usermpc: nonlinear geometric'
            write(*,*) '      effects are turned on'
            write(*,*)
            if(iperturb(1).eq.0) iperturb(1)=2
         else
            write(*,*) '*ERROR in usermpc: mpc of type',labmpc(nmpc)
            write(*,*) '       is unknown'
            call exit(201)
         endif
      endif
!
      return
      end


