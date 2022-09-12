!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2017 Guido Dhondt
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
      subroutine gen3dboun(ikboun,ilboun,nboun,nboun_,nodeboun,ndirboun,
     &     xboun,iamboun,typeboun,iponoel,inoel,iponoelmax,kon,ipkon,
     &     lakon,ne,iponor,xnor,knor,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &     mpcfree,ikmpc,ilmpc,labmpc,rig,ntrans,inotr,trab,nam,nk,nk_,
     &     co,nmethod,iperturb,istep,vold,mi,ne2boun)
!     
!     connects nodes of 1-D and 2-D elements, for which SPC's were
!     defined, to the nodes of their expanded counterparts
!     
      implicit none
!     
      logical fixed
!     
      character*1 type,typeboun(*)
      character*8 lakon(*)
      character*20 labmpc(*),label
!     
      integer ikboun(*),ilboun(*),nboun,nboun_,nodeboun(*),nodeact,
     &     ndirboun(*),idim,ier,matz,ncgnodes,lstart,lend,linc,nnodes,
     &     iamboun(*),iponoel(*),inoel(3,*),iponoelmax,kon(*),ipkon(*),
     &     iponor(2,*),knor(*),ipompc(*),nodempc(3,*),nmpc,nmpc_,
     &     ikmpc(*),ilmpc(*),rig(*),ntrans,inotr(2,*),nbounold,i,node,
     &     index,ielem,j,indexe,indexk,idir,iamplitude,irotnode,nk,nk_,
     &     newnode,idof,id,mpcfreenew,k,nam,nmethod,iperturb(*),
     &     idepnodes(80),l,iexpnode,indexx,irefnode,imax,isol,
     &     nod,impc,istep,nrhs,ipiv(3),info,m,mi(*),itr,idirref,nnode,
     &     ne2boun(2,*),ispc1,ispc2,node1,node2,ne,mpcfree,ndepnodes,
     &     midfix(3,4),mpcfreeold
!     
      real*8 xboun(*),xnor(*),coefmpc(*),trab(7,*),val,co(3,*),
     &     xnoref(3),dmax,d(3,3),e(3,3,3),alpha,q(3),w(3),xn(3),
     &     a1(3),a2(3),dd,c1,c2,c3,ww,c(3,3),vold(0:mi(2),*),a(3,3),
     &     e1(3),e2(3),t1(3),b(3,3),x(3),y(3),fv1(3),dot,
     &     fv2(3),z(3,3),xi1,xi2,xi3,u(3,3),r(3,3),xnode,dot1,dot2
!     
      data d /1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      data e /0.d0,0.d0,0.d0,0.d0,0.d0,-1.d0,0.d0,1.d0,0.d0,
     &     0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,-1.d0,0.d0,0.d0,
     &     0.d0,-1.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0/
      data midfix /5,1,2,6,2,3,7,3,4,8,4,1/
!     
      label='                    '
      fixed=.false.
!     
!     remove any rotational value applied to shells
!     
      do i=1,iponoelmax
        if(ne2boun(1,i).ne.0) then
          xboun(ne2boun(1,i))=0.d0
          if(ne2boun(2,i).ne.0) then
            xboun(ne2boun(2,i))=0.d0
          endif
        endif
      enddo
!     
!     global loop over all boundary conditions
!     
      nbounold=nboun
      loop: do i=1,nbounold
        node=nodeboun(i)
        if(node.gt.iponoelmax) then
          cycle
        endif
        index=iponoel(node)
!
!     next loop is needed, since an element may have
!     been deactivated
!
        do
          if(index.eq.0) cycle loop
          ielem=inoel(1,index)
          indexe=ipkon(ielem)
          if(indexe.ge.0) exit
          index=inoel(3,index)
        enddo
!          
c        if(index.eq.0) then
c          cycle
c        endif
c        ielem=inoel(1,index)
        j=inoel(2,index)
c        indexe=ipkon(ielem)
        indexk=iponor(2,indexe+j)
        idir=ndirboun(i)
        val=xboun(i)
        if(nam.gt.0) iamplitude=iamboun(i)
!     
        if(rig(node).ne.0) then
!     
!     existing knot
!     
          if(idir.gt.3) then
            if(rig(node).lt.0) then
              write(*,*) '*ERROR in gen3dboun: in node ',node
              write(*,*) '       a rotational DOF is constrained'
              write(*,*) '       by a SPC; however, the elements'
              write(*,*) '       to which this node belongs do not'
              write(*,*) '       have rotational DOFs'
              call exit(201)
            endif
            j=idir-3
            irotnode=rig(node)
            type='B'
            call bounadd(irotnode,j,j,val,nodeboun,
     &           ndirboun,xboun,nboun,nboun_,iamboun,
     &           iamplitude,nam,ipompc,nodempc,coefmpc,
     &           nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &           ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &           type,typeboun,nmethod,iperturb,fixed,vold,
     &           irotnode,mi,label)
          endif
        else
!     
!     1. no existing knot
!     rotational dof
!     dynamics
!     => knot is created
!     
!     check for rotational DOFs defined in any but the first step
!     
!     nonlinear dynamic case and/or beams: creation of knots
!     
!     knots (expandable rigid bodies) can take rotational
!     values arbitrarily exceeding 90 degrees
!     
          if((idir.gt.3).and.((nmethod.eq.4).and.(iperturb(1).gt.1)))
     &         then
!     
!     create a knot: determine the knot
!     
            ndepnodes=0
            if(lakon(ielem)(7:7).eq.'L') then
              do k=1,3
                ndepnodes=ndepnodes+1
                idepnodes(ndepnodes)=knor(indexk+k)
              enddo
              idim=1
            elseif(lakon(ielem)(7:7).eq.'B') then
              do k=1,8
                ndepnodes=ndepnodes+1
                idepnodes(ndepnodes)=knor(indexk+k)
              enddo
              idim=3
            else
              write(*,*) 
     &             '*ERROR in gen3dboun: a rotational DOF was applied'
              write(*,*) 
     &             '*      to node',node,' without rotational DOFs'
              call exit(201)
            endif
!     
!     remove all MPC's in which the knot nodes are
!     dependent nodes
!     
            do k=1,ndepnodes
              nod=idepnodes(k)
              do l=1,3
                idof=8*(nod-1)+l
                call nident(ikmpc,idof,nmpc,id)
                if(id.gt.0) then
                  if(ikmpc(id).eq.idof) then
                    impc=ilmpc(id)
                    call mpcrem(impc,mpcfree,nodempc,nmpc,
     &                   ikmpc,ilmpc,labmpc,coefmpc,ipompc)
                  endif
                endif
              enddo
            enddo
!     
!     generate a rigid body knot
!     
            irefnode=node
            nk=nk+1
            if(nk.gt.nk_) then
              write(*,*) '*ERROR in rigidbodies: increase nk_'
              call exit(201)
            endif
            irotnode=nk
            rig(node)=irotnode
            write(27,*) 'a KNOT was generated in node ',node
            write(27,*)
            nk=nk+1
            if(nk.gt.nk_) then
              write(*,*) '*ERROR in rigidbodies: increase nk_'
              call exit(201)
            endif
            iexpnode=nk
            do k=1,ndepnodes
              call knotmpc(ipompc,nodempc,coefmpc,irefnode,
     &             irotnode,iexpnode,
     &             labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,nk,nk_,
     &             nodeboun,ndirboun,ikboun,ilboun,nboun,nboun_,
     &             idepnodes,typeboun,co,xboun,istep,k,ndepnodes,
     &             idim,e1,e2,t1)
            enddo
!     
!     determine the location of the center of gravity of
!     the section and its displacements
!     
            do l=1,3
              q(l)=0.d0
              w(l)=0.d0
            enddo
            if(ndepnodes.eq.3) then
              do k=1,ndepnodes,2
                nod=idepnodes(k)
                do l=1,3
                  q(l)=q(l)+co(l,nod)
                  w(l)=w(l)+vold(l,nod)
                enddo
              enddo
              do l=1,3
                q(l)=q(l)/2.d0
                w(l)=w(l)/2.d0
              enddo
            else
              do k=1,ndepnodes
                nod=idepnodes(k)
                do l=1,3
                  q(l)=q(l)+co(l,nod)
                  w(l)=w(l)+vold(l,nod)
                enddo
              enddo
              do l=1,3
                q(l)=q(l)/ndepnodes
                w(l)=w(l)/ndepnodes
              enddo
            endif
!     
!     check whether the displacements are zero
!     
            dd=dsqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3))
            if(dd.lt.1.d-20) then
              do l=1,3
                vold(l,irefnode)=0.d0
                vold(l,irotnode)=0.d0
                vold(l,iexpnode)=0.d0
              enddo
            else
!     
!     determine the first displacements of iexpnode
!     
              alpha=0.d0
              ncgnodes=0
              do k=1,ndepnodes
                nod=idepnodes(k)
                dd=(co(1,nod)-q(1))**2
     &               +(co(2,nod)-q(2))**2
     &               +(co(3,nod)-q(3))**2
                if(dd.lt.1.d-20) then
                  ncgnodes=ncgnodes+1
                  cycle
                endif
                alpha=alpha+dsqrt(
     &               ((co(1,nod)+vold(1,nod)-q(1)-w(1))**2
     &               +(co(2,nod)+vold(2,nod)-q(2)-w(2))**2
     &               +(co(3,nod)+vold(3,nod)-q(3)-w(3))**2)/dd)
              enddo
              if(ndepnodes-ncgnodes.gt.0) then
                alpha=alpha/(ndepnodes-ncgnodes)
              endif
!     
!     determine the displacements of irotnodes
!     
              do l=1,3
                xn(l)=0.d0
              enddo
!     
              ncgnodes=0
              do k=1,ndepnodes
                nod=idepnodes(k)
                dd=0.d0
                do l=1,3
                  a1(l)=co(l,nod)-q(l)
                  a2(l)=vold(l,nod)-w(l)
                  dd=dd+a1(l)*a1(l)
                enddo
                dd=dsqrt(dd)
                if(dd.lt.1.d-10) then
                  ncgnodes=ncgnodes+1
                  cycle
                endif
                do l=1,3
                  a1(l)=a1(l)/dd
                  a2(l)=a2(l)/dd
                enddo
                xn(1)=xn(1)+(a1(2)*a2(3)-a1(3)*a2(2))
                xn(2)=xn(2)+(a1(3)*a2(1)-a1(1)*a2(3))
                xn(3)=xn(3)+(a1(1)*a2(2)-a1(2)*a2(1))
              enddo
!     
              if(ndepnodes-ncgnodes.gt.0) then
                do l=1,3
                  xn(l)=xn(l)/(ndepnodes-ncgnodes)
                enddo
              endif
!     
              dd=0.d0
              do l=1,3
                dd=dd+xn(l)*xn(l)
              enddo
              dd=dsqrt(dd)
              do l=1,3
                xn(l)=dasin(dd/alpha)*xn(l)/dd
              enddo
!     
!     determine the displacements of irefnode
!     
              ww=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
!     
              c1=dcos(ww)
              if(ww.gt.1.d-10) then
                c2=dsin(ww)/ww
              else
                c2=1.d0
              endif
              if(ww.gt.1.d-5) then
                c3=(1.d0-c1)/ww**2
              else
                c3=0.5d0
              endif
!     
!     rotation matrix c
!     
              do k=1,3
                do l=1,3
                  r(k,l)=c1*d(k,l)+
     &                 c2*(e(k,1,l)*xn(1)+e(k,2,l)*xn(2)+
     &                 e(k,3,l)*xn(3))+c3*xn(k)*xn(l)
                enddo
              enddo
!     
!     copying the displacements
!     
              do l=1,3
                vold(l,irefnode)=w(l)
                vold(l,irotnode)=xn(l)
              enddo
              vold(1,iexpnode)=alpha-1.d0
!     
!     correction of the expansion values for beam sections
!     
              if(idim.eq.2) then
!     
!     initializing matrices b and c
!     
                do l=1,3
                  do m=1,3
                    b(l,m)=0.d0
                    c(l,m)=0.d0
                  enddo
                enddo
!     
!     the transpose of the deformation gradient:
!     c.F^T=b
!     
                do k=1,ndepnodes
                  nod=idepnodes(k)
                  do l=1,3
                    x(l)=co(l,nod)-q(l)
                    y(l)=x(l)+vold(l,nod)-w(l)
                  enddo
                  do l=1,3
                    do m=1,3
                      c(l,m)=c(l,m)+x(l)*x(m)
                      b(l,m)=b(l,m)+x(l)*y(m)
                    enddo
                  enddo
                enddo
!     
!     solving the linear equation system
!     
                m=3
                nrhs=3
                call dgesv(m,nrhs,c,m,ipiv,b,m,info)
                if(info.ne.0) then
                  write(*,*) '*ERROR in gen3dforc:'
                  write(*,*) '       singular system of equations'
                  call exit(201)
                endif
!     
!     now b=F^T
!     
!     constructing the right stretch tensor
!     U=F^T.R
!     
                do l=1,3
                  do m=l,3
                    u(l,m)=b(l,1)*r(1,m)+b(l,2)*r(2,m)+
     &                   b(l,3)*r(3,m)
                  enddo
                enddo
                u(2,1)=u(1,2)
                u(3,1)=u(1,3)
                u(3,2)=u(2,3)
!     
!     determining the eigenvalues and eigenvectors of U
!     
                m=3
                matz=1
                ier=0
                call rs(m,m,u,w,matz,z,fv1,fv2,ier)
                if(ier.ne.0) then
                  write(*,*) 
     &                 '*ERROR in knotmpc while calculating the'
                  write(*,*) '       eigenvalues/eigenvectors'
                  call exit(201)
                endif
!     
                if((dabs(w(1)-1.d0).lt.dabs(w(2)-1.d0)).and.
     &               (dabs(w(1)-1.d0).lt.dabs(w(3)-1.d0))) then
                  l=2
                  m=3
                elseif((dabs(w(2)-1.d0).lt.dabs(w(1)-1.d0)).and.
     &                 (dabs(w(2)-1.d0).lt.dabs(w(3)-1.d0))) then
                  l=1
                  m=3
                else
                  l=1
                  m=2
                endif
                xi1=datan2
     &               ((z(1,l)*e2(1)+z(2,l)*e2(2)+z(3,l)*e2(2)),
     &               (z(1,l)*e1(1)+z(2,l)*e1(2)+z(3,l)*e1(2)))
                xi2=w(l)-1.d0
                xi3=w(m)-1.d0
!     
                vold(1,iexpnode)=xi1
                vold(2,iexpnode)=xi2
                vold(3,iexpnode)=xi3
              endif
            endif
!     
!     apply the boundary condition
!     
            idir=idir-3
            type='B'
            call bounadd(irotnode,idir,idir,val,nodeboun,
     &           ndirboun,xboun,nboun,nboun_,iamboun,
     &           iamplitude,nam,ipompc,nodempc,coefmpc,
     &           nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &           ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &           type,typeboun,nmethod,iperturb,fixed,vold,
     &           irotnode,mi,label)
!     
!     check for shells whether the rotation about the normal
!     on the shell has been eliminated
!     
            if(lakon(ielem)(7:7).eq.'L') then
              indexx=iponor(1,indexe+j)
              do j=1,3
                xnoref(j)=xnor(indexx+j)
              enddo
              dmax=0.d0
              imax=0
              do j=1,3
                if(dabs(xnoref(j)).gt.dmax) then
                  dmax=dabs(xnoref(j))
                  imax=j
                endif
              enddo
!     
!     check whether a SPC suffices
!     
              if(dabs(1.d0-dmax).lt.1.d-3) then
                val=0.d0
                if(nam.gt.0) iamplitude=0
                type='R'
                call bounadd(irotnode,imax,imax,val,nodeboun,
     &               ndirboun,xboun,nboun,nboun_,iamboun,
     &               iamplitude,nam,ipompc,nodempc,coefmpc,
     &               nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &               ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &               type,typeboun,nmethod,iperturb,fixed,vold,
     &               irotnode,mi,label)
              else
!     
!     check for an unused rotational DOF
!     
                isol=0
                do l=1,3
                  idof=8*(node-1)+3+imax
                  call nident(ikboun,idof,nboun,id)
                  if((id.gt.0).and.(ikboun(id).eq.idof)) then
                    imax=imax+1
                    if(imax.gt.3) imax=imax-3
                    cycle
                  endif
                  isol=1
                  exit
                enddo
!     
!     if one of the rotational dofs was not used so far,
!     it can be taken as dependent side for fixing the
!     rotation about the normal. If all dofs were used,
!     no additional equation is needed.
!     
                if(isol.eq.1) then
                  idof=8*(irotnode-1)+imax
                  call nident(ikmpc,idof,nmpc,id)
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                    write(*,*) 
     &                   '*ERROR in gen3dboun: increase nmpc_'
                    call exit(201)
                  endif
!     
                  ipompc(nmpc)=mpcfree
                  labmpc(nmpc)='                    '
!     
                  do l=nmpc,id+2,-1
                    ikmpc(l)=ikmpc(l-1)
                    ilmpc(l)=ilmpc(l-1)
                  enddo
                  ikmpc(id+1)=idof
                  ilmpc(id+1)=nmpc
!     
                  nodempc(1,mpcfree)=irotnode
                  nodempc(2,mpcfree)=imax
                  coefmpc(mpcfree)=xnoref(imax)
                  mpcfree=nodempc(3,mpcfree)
                  imax=imax+1
                  if(imax.gt.3) imax=imax-3
                  nodempc(1,mpcfree)=irotnode
                  nodempc(2,mpcfree)=imax
                  coefmpc(mpcfree)=xnoref(imax)
                  mpcfree=nodempc(3,mpcfree)
                  imax=imax+1
                  if(imax.gt.3) imax=imax-3
                  nodempc(1,mpcfree)=irotnode
                  nodempc(2,mpcfree)=imax
                  coefmpc(mpcfree)=xnoref(imax)
                  mpcfreeold=mpcfree
                  mpcfree=nodempc(3,mpcfree)
                  nodempc(3,mpcfreeold)=0
                endif
              endif
            endif
            cycle
          endif
!     
!     2. no existing knot
!     rotational dof
!     no dynamics
!     => mean rotation MPC is created
!     
!     all cases except nonlinear dynamic case: creation
!     of meanrotation MPC's
!     
          if((idir.gt.3).and.((nmethod.ne.4).or.(iperturb(1).le.1)))
     &         then
!     
!     create a mean rotation MPC
!     advantage: more accurate since less constraining
!     disadvantage: cannot exceed 90 degrees rotation
!     
!     if a mean rotation MPC has already been created 
!     for idof, ilboun(id) contains the index of the
!     SPC created for the angle value
!     
            idof=8*(node-1)+idir
            call nident(ikboun,idof,nboun,id)
            if(ilboun(id).ne.i) cycle

            idirref=idir-3
!     
            if(lakon(ielem)(7:7).eq.'L') then
              lstart=3
              lend=1
              linc=-2
!     
!     calculating the normal vector =
!     vector along the drilling direction
!     
              do j=1,3
                xn(j)=co(j,knor(indexk+3))-co(j,knor(indexk+1))
              enddo
              dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
              do j=1,3
                xn(j)=xn(j)/dd
              enddo
!     
            elseif(lakon(ielem)(7:7).eq.'B') then
              lstart=4
              lend=1
              linc=-1
            endif
!     
!     check for transformations
!     
            if(ntrans.le.0) then
              itr=0
            elseif(inotr(1,node).eq.0) then
              itr=0
            else
              itr=inotr(1,node)
            endif
!     
!     determine a unit vector on the rotation axis
!     
            if(itr.eq.0) then
              do j=1,3
                do k=1,3
                  a(j,k)=0.d0
                enddo
                a(j,j)=1.d0
              enddo
            else
              call transformatrix(trab(1,itr),co(1,node),a)
            endif
!     
!     check whether the rotation vector does not
!     have a component along the drilling direction
!     
            if(lakon(ielem)(7:7).eq.'L') then
              dot=a(1,idirref)*xn(1)+a(2,idirref)*xn(2)+
     &             a(3,idirref)*xn(3)
!     
!     rotation vectors along the drilling direction 
!     are not taken into account
!     
!     in general, rotation vectors are projected onto
!     the shell surface
!     
              if(dabs(dot).gt.0.999d0) cycle
!     
!     check for all applied rotations in this node
!     
!     idir1 is the dof the projection of which
!     corresponds to the first unit vector in the
!     shell plane (corresponds to ispc1)
!     
!     if one rotation is applied in this node, its
!     values is stored in SPC ispc1
!     if more than one rotation is applied in this node,
!     the value of the second rotation (only two are
!     possible in the shell plane) is stored in SPC ispc2
!     
              ispc1=ne2boun(1,node)
              ispc2=ne2boun(2,node)
!     
!     idr(k),k=1,3 are the SPC's corresponding to
!     dofs 4..6 (0 if not applied)
!     
              if(ispc1.eq.0) then
!     
!     rotation not applied yet: 
!     - project rotation on shell plane
!     - create MPC
!     
!     projecting the rotation vector on the tangent plane
!     
                do k=1,3
                  a(k,idirref)=a(k,idirref)-dot*xn(k)
                enddo
                dd=0.d0
                do k=1,3
                  dd=dd+a(k,idirref)**2
                enddo
                dd=dsqrt(dd)
                do k=1,3
                  a(k,idirref)=a(k,idirref)/dd
                enddo
                val=val*dsqrt(1.d0-dot*dot)
!     
!     specific label for mean rotations for beams and
!     shells
!     
                label='MEANROTBS           '
                write(27,*) 
     &               'a MEAN ROTATION MPC was generated in node ',
     &               node
                write(27,*) '                about the axis (',
     &               a(1,idirref),',',a(2,idirref),',',a(3,idirref),')'
                write(27,*)
                nnodes=0
                do j=lstart,lend,linc
                  nodeact=knor(indexk+j)
                  do k=1,3
                    nnodes=nnodes+1
                    call usermpc(ipompc,nodempc,coefmpc,
     &                   labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &                   nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &                   nboun,nboun_,nnodes,nodeact,co,label,
     &                   typeboun,iperturb,node,idirref,xboun)
                  enddo
                enddo
!     
!     rotation value term
!     
                nodeact=nk+1
                do k=1,3
                  co(k,nodeact)=a(k,idirref)
                enddo
                nnodes=nnodes+1
                call usermpc(ipompc,nodempc,coefmpc,
     &               labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &               nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &               nboun,nboun_,nnodes,nodeact,co,label,
     &               typeboun,iperturb,node,idirref,xboun)
!     
!     inhomogeneous term
!     
                nodeact=0
                call usermpc(ipompc,nodempc,coefmpc,
     &               labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &               nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &               nboun,nboun_,nnodes,nodeact,co,label,
     &               typeboun,iperturb,node,idirref,xboun)
!     
!     end meanrotationmpc
!     
!     SPC angle term
!     
                if(nodeact.ne.-1) then
                  idir=1
                  type='B'
                  call bounadd(nk,idir,idir,val,nodeboun,
     &                 ndirboun,xboun,nboun,nboun_,iamboun,
     &                 iamplitude,nam,ipompc,nodempc,coefmpc,
     &                 nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &                 ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &                 type,typeboun,nmethod,iperturb,fixed,vold,
     &                 nk,mi,label)
!     
!     storing the index of the SPC with the angle
!     value in ne2boun
!     
                  ne2boun(1,node)=nboun
                endif
              elseif(a(1,idirref)*co(1,nodeboun(ispc1))+
     &               a(2,idirref)*co(2,nodeboun(ispc1))+
     &               a(3,idirref)*co(3,nodeboun(ispc1)).gt.
     &               0.9999d0) then
!     
!     the dof of the applied rotation corresponds
!     to the first dof in the shell plane
!     
                val=val*dsqrt(1.d0-dot*dot)
                xboun(ispc1)=xboun(ispc1)+val
              else
!     
!     applied rotation does not correspond to the
!     first coordinate direction in the shell plane
!     
!     determine first direction
!     
                node1=nodeboun(ispc1)
                do k=1,3
                  e1(k)=co(k,node1)
                enddo
!     
!     determine second direction
!     
                if(ispc2.ne.0) then
                  node2=nodeboun(ispc2)
                  do k=1,3
                    e2(k)=co(k,node2)
                  enddo
                else
!     
!     e2=n x e1
!     
                  e2(1)=xn(2)*e1(3)-xn(3)*e1(2)
                  e2(2)=xn(3)*e1(1)-xn(1)*e1(3)
                  e2(3)=xn(1)*e1(2)-xn(2)*e1(1)
!     
                  dd=dsqrt(e2(1)*e2(1)+e2(2)*e2(2)+e2(3)*e2(3))
!     
                  do k=1,3
                    e2(k)=e2(k)/dd
                  enddo
                endif
!     
!     determine cosines between rotation direction and
!     unit vectors in shell plane  
!     
                dot1=a(1,idirref)*e1(1)
     &               +a(2,idirref)*e1(2)
     &               +a(3,idirref)*e1(3)
                dot2=a(1,idirref)*e2(1)
     &               +a(2,idirref)*e2(2)
     &               +a(3,idirref)*e2(3)
!     
!     contribution to first direction
!     
                xboun(ispc1)=xboun(ispc1)+val*dot1
!     
                if(ispc2.eq.0) then
!     
!     SPC corresponding to second direction has to
!     be created
!     
                  val=val*dot2
                  do k=1,3
                    a(k,idirref)=e2(k)
                  enddo
!     
                  label='MEANROTBS           '
                  write(27,*) 
     &                 'a MEAN ROTATION MPC was generated in node ',
     &                 node
                  write(27,*) '                about the axis (',
     &                 a(1,idirref),',',a(2,idirref),',',
     &                 a(3,idirref),')'
                  write(27,*)
                  nnodes=0
                  do j=lstart,lend,linc
                    nodeact=knor(indexk+j)
                    do k=1,3
                      nnodes=nnodes+1
                      call usermpc(ipompc,nodempc,coefmpc,
     &                     labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &                     nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &                     nboun,nboun_,nnodes,nodeact,co,label,
     &                     typeboun,iperturb,node,idirref,xboun)
                    enddo
                  enddo
!     
!     rotation value term
!     
                  nodeact=nk+1
                  do k=1,3
                    co(k,nodeact)=a(k,idirref)
                  enddo
                  nnodes=nnodes+1
                  call usermpc(ipompc,nodempc,coefmpc,
     &                 labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &                 nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &                 nboun,nboun_,nnodes,nodeact,co,label,
     &                 typeboun,iperturb,node,idirref,xboun)
!     
!     inhomogeneous term
!     
                  nodeact=0
                  call usermpc(ipompc,nodempc,coefmpc,
     &                 labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &                 nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &                 nboun,nboun_,nnodes,nodeact,co,label,
     &                 typeboun,iperturb,node,idirref,xboun)
!     
!     end meanrotationmpc
!     
!     SPC angle term
!     
                  if(nodeact.ne.-1) then
                    idir=1
                    type='B'
                    call bounadd(nk,idir,idir,val,nodeboun,
     &                   ndirboun,xboun,nboun,nboun_,iamboun,
     &                   iamplitude,nam,ipompc,nodempc,coefmpc,
     &                   nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &                   ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &                   type,typeboun,nmethod,iperturb,fixed,vold,
     &                   nk,mi,label)
!     
!     storing the index of the SPC with the angle
!     value in ne2boun
!     
c     ne2boun(2,ispcref)=nboun
                    ne2boun(2,node)=nboun
                  endif
                else
!     
!     SPC corresponding to second direction in shell
!     plane has already been created
!     
                  xboun(ispc2)=xboun(ispc2)+val*dot2
                endif
              endif
              cycle
            endif
!     
!     beams             
!     
!     specific label for mean rotations for beams and
!     shells
!     
            label='MEANROTBS           '
            write(27,*) 'a MEAN ROTATION MPC was generated in node ',
     &           node
            write(27,*) '                about the axis (',
     &           a(1,idirref),',',a(2,idirref),',',a(3,idirref),')'
            write(27,*)
            nnodes=0
            do j=lstart,lend,linc
              nodeact=knor(indexk+j)
              do k=1,3
                nnodes=nnodes+1
                call usermpc(ipompc,nodempc,coefmpc,
     &               labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &               nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &               nboun,nboun_,nnodes,nodeact,co,label,
     &               typeboun,iperturb,node,idirref,xboun)
              enddo
            enddo
!     
!     rotation value term
!     
            nodeact=nk+1
            do k=1,3
              co(k,nodeact)=a(k,idirref)
            enddo
            nnodes=nnodes+1
            call usermpc(ipompc,nodempc,coefmpc,
     &           labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &           nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &           nboun,nboun_,nnodes,nodeact,co,label,
     &           typeboun,iperturb,node,idirref,xboun)
!     
!     inhomogeneous term
!     
            nodeact=0
            call usermpc(ipompc,nodempc,coefmpc,
     &           labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &           nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &           nboun,nboun_,nnodes,nodeact,co,label,
     &           typeboun,iperturb,node,idirref,xboun)
!     
!     end meanrotationmpc
!     
!     SPC angle term
!     
            if(nodeact.ne.-1) then
              idir=1
              type='B'
              call bounadd(nk,idir,idir,val,nodeboun,
     &             ndirboun,xboun,nboun,nboun_,iamboun,
     &             iamplitude,nam,ipompc,nodempc,coefmpc,
     &             nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &             ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &             type,typeboun,nmethod,iperturb,fixed,vold,
     &             nk,mi,label)
!     
!     storing the index of the SPC with the angle
!     value in ilboun(id)
!     
              ilboun(id)=nboun
            endif
!     
!     for quadratic beams: fixing the middle nodes in between the
!     end nodes
!     
            if((lakon(ielem)(7:7).eq.'B').and.
     &           (lakon(ielem)(4:4).eq.'2')) then
              do k=1,4
                newnode=knor(indexk+midfix(1,k))
                do idir=1,3
                  idof=8*(newnode-1)+idir
                  call nident(ikmpc,idof,nmpc,id)
                  if((id.le.0).or.(ikmpc(id).ne.idof)) then
                    nmpc=nmpc+1
                    if(nmpc.gt.nmpc_) then
                      write(*,*) 
     &                     '*ERROR in gen3dboun: increase nmpc_'
                      call exit(201)
                    endif
                    labmpc(nmpc)='                    '
                    ipompc(nmpc)=mpcfree
                    do j=nmpc,id+2,-1
                      ikmpc(j)=ikmpc(j-1)
                      ilmpc(j)=ilmpc(j-1)
                    enddo
                    ikmpc(id+1)=idof
                    ilmpc(id+1)=nmpc
                    nodempc(1,mpcfree)=newnode
                    nodempc(2,mpcfree)=idir
                    coefmpc(mpcfree)=2.d0
                    mpcfree=nodempc(3,mpcfree)
                    if(mpcfree.eq.0) then
                      write(*,*) 
     &                     '*ERROR in gen3dboun: increase memmpc_'
                      call exit(201)
                    endif
                    nodempc(1,mpcfree)=knor(indexk+midfix(2,k))
                    nodempc(2,mpcfree)=idir
                    coefmpc(mpcfree)=-1.d0
                    mpcfree=nodempc(3,mpcfree)
                    if(mpcfree.eq.0) then
                      write(*,*) 
     &                     '*ERROR in gen3dboun: increase memmpc_'
                      call exit(201)
                    endif
                    nodempc(1,mpcfree)=knor(indexk+midfix(3,k))
                    nodempc(2,mpcfree)=idir
                    coefmpc(mpcfree)=-1.d0
                    mpcfreenew=nodempc(3,mpcfree)
                    if(mpcfreenew.eq.0) then
                      write(*,*) 
     &                     '*ERROR in gen3dboun: increase memmpc_'
                      call exit(201)
                    endif
                    nodempc(3,mpcfree)=0
                    mpcfree=mpcfreenew
                  endif
                enddo
              enddo
            endif
!     
            cycle
          endif
!     
!     3. no existing knot
!     translational dof
!     
!     a) 2d shell element: generate MPC's
!     
!     u(n_1)+u(n_3)=2*u(n)
!     
          if(lakon(ielem)(7:7).eq.'L') then
            newnode=knor(indexk+1)
            idof=8*(newnode-1)+idir
            call nident(ikmpc,idof,nmpc,id)
            if((id.le.0).or.(ikmpc(id).ne.idof)) then
              nmpc=nmpc+1
              if(nmpc.gt.nmpc_) then
                write(*,*) 
     &               '*ERROR in gen3dboun: increase nmpc_'
                call exit(201)
              endif
              labmpc(nmpc)='                    '
              ipompc(nmpc)=mpcfree
              do j=nmpc,id+2,-1
                ikmpc(j)=ikmpc(j-1)
                ilmpc(j)=ilmpc(j-1)
              enddo
              ikmpc(id+1)=idof
              ilmpc(id+1)=nmpc
              nodempc(1,mpcfree)=newnode
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=1.d0
              mpcfree=nodempc(3,mpcfree)
              if(mpcfree.eq.0) then
                write(*,*) 
     &               '*ERROR in gen3dboun: increase memmpc_'
                call exit(201)
              endif
              nodempc(1,mpcfree)=knor(indexk+3)
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=1.d0
              mpcfree=nodempc(3,mpcfree)
              if(mpcfree.eq.0) then
                write(*,*) 
     &               '*ERROR in gen3dboun: increase memmpc_'
                call exit(201)
              endif
              nodempc(1,mpcfree)=node
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=-2.d0
              mpcfreenew=nodempc(3,mpcfree)
              if(mpcfreenew.eq.0) then
                write(*,*) 
     &               '*ERROR in gen3dboun: increase memmpc_'
                call exit(201)
              endif
              nodempc(3,mpcfree)=0
              mpcfree=mpcfreenew
            endif
!     
!     u(n_2)=u(n)
!     
            newnode=knor(indexk+2)
            idof=8*(newnode-1)+idir
            call nident(ikmpc,idof,nmpc,id)
            if((id.le.0).or.(ikmpc(id).ne.idof)) then
              nmpc=nmpc+1
              if(nmpc.gt.nmpc_) then
                write(*,*) 
     &               '*ERROR in gen3dboun: increase nmpc_'
                call exit(201)
              endif
              labmpc(nmpc)='                    '
              ipompc(nmpc)=mpcfree
              do j=nmpc,id+2,-1
                ikmpc(j)=ikmpc(j-1)
                ilmpc(j)=ilmpc(j-1)
              enddo
              ikmpc(id+1)=idof
              ilmpc(id+1)=nmpc
              nodempc(1,mpcfree)=newnode
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=1.d0
              mpcfree=nodempc(3,mpcfree)
              if(mpcfree.eq.0) then
                write(*,*) 
     &               '*ERROR in gen3dboun: increase memmpc_'
                call exit(201)
              endif
              nodempc(1,mpcfree)=node
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=-1.d0
              mpcfreenew=nodempc(3,mpcfree)
              if(mpcfreenew.eq.0) then
                write(*,*) 
     &               '*ERROR in gen3dboun: increase memmpc_'
                call exit(201)
              endif
              nodempc(3,mpcfree)=0
              mpcfree=mpcfreenew
            endif
!     
!     fixing the temperature degrees of freedom
!     
            if(idir.eq.0) then
!     
!     t(n_3)=t(n)
!     
              newnode=knor(indexk+3)
              idof=8*(newnode-1)+idir
              call nident(ikmpc,idof,nmpc,id)
              if((id.le.0).or.(ikmpc(id).ne.idof)) then
                nmpc=nmpc+1
                if(nmpc.gt.nmpc_) then
                  write(*,*) 
     &                 '*ERROR in gen3dboun: increase nmpc_'
                  call exit(201)
                endif
                labmpc(nmpc)='                    '
                ipompc(nmpc)=mpcfree
                do j=nmpc,id+2,-1
                  ikmpc(j)=ikmpc(j-1)
                  ilmpc(j)=ilmpc(j-1)
                enddo
                ikmpc(id+1)=idof
                ilmpc(id+1)=nmpc
                nodempc(1,mpcfree)=newnode
                nodempc(2,mpcfree)=idir
                coefmpc(mpcfree)=1.d0
                mpcfree=nodempc(3,mpcfree)
                if(mpcfree.eq.0) then
                  write(*,*) 
     &                 '*ERROR in gen3dboun: increase memmpc_'
                  call exit(201)
                endif
                nodempc(1,mpcfree)=node
                nodempc(2,mpcfree)=idir
                coefmpc(mpcfree)=-1.d0
                mpcfreenew=nodempc(3,mpcfree)
                if(mpcfreenew.eq.0) then
                  write(*,*) 
     &                 '*ERROR in gen3dboun: increase memmpc_'
                  call exit(201)
                endif
                nodempc(3,mpcfree)=0
                mpcfree=mpcfreenew
              endif
            endif
          elseif(lakon(ielem)(7:7).eq.'B') then
!     
!     b) 1d beam element: generate MPC's
!     
!     u(n_1)+u(n_2)+u(n_3)+u(n_4)=4*u(n)
!     
            nnode=4
            xnode=4.d0
            newnode=knor(indexk+1)
            idof=8*(newnode-1)+idir
            call nident(ikmpc,idof,nmpc,id)
            if((id.le.0).or.(ikmpc(id).ne.idof)) then
              nmpc=nmpc+1
              if(nmpc.gt.nmpc_) then
                write(*,*) 
     &               '*ERROR in gen3dboun: increase nmpc_'
                call exit(201)
              endif
              labmpc(nmpc)='                    '
              ipompc(nmpc)=mpcfree
              do j=nmpc,id+2,-1
                ikmpc(j)=ikmpc(j-1)
                ilmpc(j)=ilmpc(j-1)
              enddo
              ikmpc(id+1)=idof
              ilmpc(id+1)=nmpc
              nodempc(1,mpcfree)=newnode
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=1.d0
              mpcfree=nodempc(3,mpcfree)
              if(mpcfree.eq.0) then
                write(*,*) 
     &               '*ERROR in gen3dboun: increase memmpc_'
                call exit(201)
              endif
              do k=2,nnode
                nodempc(1,mpcfree)=knor(indexk+k)
                nodempc(2,mpcfree)=idir
                coefmpc(mpcfree)=1.d0
                mpcfree=nodempc(3,mpcfree)
                if(mpcfree.eq.0) then
                  write(*,*) 
     &                 '*ERROR in gen3dboun: increase memmpc_'
                  call exit(201)
                endif
              enddo
              nodempc(1,mpcfree)=node
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=-xnode
              mpcfreenew=nodempc(3,mpcfree)
              if(mpcfreenew.eq.0) then
                write(*,*) 
     &               '*ERROR in gen3dboun: increase memmpc_'
                call exit(201)
              endif
              nodempc(3,mpcfree)=0
              mpcfree=mpcfreenew
            endif
!     
!     fixing the temperature degrees of freedom
!     
            if(idir.eq.0) then
              do k=2,4
!     
!     t(n_k)=t(n), k=2,4
!     
                newnode=knor(indexk+k)
                idof=8*(newnode-1)+idir
                call nident(ikmpc,idof,nmpc,id)
                if((id.le.0).or.(ikmpc(id).ne.idof)) then
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                    write(*,*) 
     &                   '*ERROR in gen3dboun: increase nmpc_'
                    call exit(201)
                  endif
                  labmpc(nmpc)='                    '
                  ipompc(nmpc)=mpcfree
                  do j=nmpc,id+2,-1
                    ikmpc(j)=ikmpc(j-1)
                    ilmpc(j)=ilmpc(j-1)
                  enddo
                  ikmpc(id+1)=idof
                  ilmpc(id+1)=nmpc
                  nodempc(1,mpcfree)=newnode
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                    write(*,*) 
     &                   '*ERROR in gen3dboun: increase memmpc_'
                    call exit(201)
                  endif
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-1.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                    write(*,*) 
     &                   '*ERROR in gen3dboun: increase memmpc_'
                    call exit(201)
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
                endif
              enddo
            endif
          else
!     
!     c)       2d plane stress, plane strain or axisymmetric
!     element: MPC in all but z-direction
!     
            newnode=knor(indexk+2)
            idof=8*(newnode-1)+idir
            call nident(ikmpc,idof,nmpc,id)
            if(((id.le.0).or.(ikmpc(id).ne.idof)).and.
     &           (idir.ne.3)) then
              nmpc=nmpc+1
              if(nmpc.gt.nmpc_) then
                write(*,*) 
     &               '*ERROR in gen3dmpc: increase nmpc_'
                call exit(201)
              endif
              labmpc(nmpc)='                    '
              ipompc(nmpc)=mpcfree
              do j=nmpc,id+2,-1
                ikmpc(j)=ikmpc(j-1)
                ilmpc(j)=ilmpc(j-1)
              enddo
              ikmpc(id+1)=idof
              ilmpc(id+1)=nmpc
              nodempc(1,mpcfree)=newnode
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=1.d0
              mpcfree=nodempc(3,mpcfree)
              if(mpcfree.eq.0) then
                write(*,*) 
     &               '*ERROR in gen3dmpc: increase memmpc_'
                call exit(201)
              endif
              nodempc(1,mpcfree)=node
              nodempc(2,mpcfree)=idir
              coefmpc(mpcfree)=-1.d0
              mpcfreenew=nodempc(3,mpcfree)
              if(mpcfreenew.eq.0) then
                write(*,*) 
     &               '*ERROR in gen3dmpc: increase memmpc_'
                call exit(201)
              endif
              nodempc(3,mpcfree)=0
              mpcfree=mpcfreenew
            endif
          endif
        endif
      enddo loop
!     
      return
      end


