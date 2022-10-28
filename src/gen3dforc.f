!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine gen3dforc(ikforc,ilforc,nforc,nforc_,nodeforc,
     &     ndirforc,xforc,iamforc,ntrans,inotr,trab,rig,ipompc,nodempc,
     &     coefmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,labmpc,iponoel,inoel,
     &     iponoelmax,kon,ipkon,lakon,ne,iponor,xnor,knor,nam,nk,nk_,
     &     co,thicke,nodeboun,ndirboun,ikboun,ilboun,nboun,nboun_,
     &     iamboun,typeboun,xboun,nmethod,iperturb,istep,vold,mi,
     &     idefforc)
!     
!     connects nodes of 1-D and 2-D elements, for which 
!     concentrated forces were
!     defined, to the nodes of their expanded counterparts
!     
      implicit none
!     
      logical add,fixed,user,quadratic
!     
      character*1 type,typeboun(*)
      character*8 lakon(*)
      character*20 labmpc(*),label
!     
      integer ikforc(*),ilforc(*),nodeforc(2,*),ndirforc(*),itr,idirref,
     &     iamforc(*),idim,ier,matz,nodeact,linc,lend,lstart,nnodes,
     &     nforc,nforc_,ntrans,inotr(2,*),rig(*),ipompc(*),nodempc(3,*),
     &     nmpc,nmpc_,mpcfree,ikmpc(*),ilmpc(*),iponoel(*),inoel(3,*),
     &     iponoelmax,kon(*),ipkon(*),ne,iponor(2,*),knor(*),nforcold,
     &     i,node,index,ielem,j,indexe,indexk,nam,iamplitude,idir,
     &     irotnode,nk,nk_,newnode,idof,id,mpcfreenew,k,isector,
     &     idepnodes(80),l,iexpnode,indexx,irefnode,imax,isol,
     &     nod,impc,istep,nodeboun(*),ndirboun(*),ikboun(*),ilboun(*),
     &     nboun,nboun_,iamboun(*),nmethod,iperturb(*),nrhs,ipiv(3),
     &     mi(*),idefforc(*),nedge,idirstart,idirend,idirl,ndepnodes,
     &     mpcfreeold,info,m
!     
      real*8 xforc(*),trab(7,*),coefmpc(*),xnor(*),val,co(3,*),dot,
     &     thicke(mi(3),*),pi,xboun(*),xnoref(3),dmax,d(3,3),e(3,3,3),
     &     alpha,q(3),w(3),xn(3),a(3,3),a1(3),a2(3),dd,c1,c2,c3,ww,
     &     vold(0:mi(2),*),e1(3),e2(3),t1(3),b(3,3),x(3),y(3),fv1(3),
     &     fv2(3),z(3,3),xi1,xi2,xi3,u(3,3),r(3,3),c(3,3)
!     
      data d /1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1./
      data e /0.d0,0.d0,0.d0,0.d0,0.d0,-1.d0,0.d0,1.d0,0.d0,
     &     0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,-1.d0,0.d0,0.d0,
     &     0.d0,-1.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0/
!     
      label='                    '
      fixed=.false.
!     
      add=.false.
      user=.false.
      pi=4.d0*datan(1.d0)
      isector=0
!     
      nforcold=nforc
      loop: do i=1,nforcold
        node=nodeforc(1,i)
        if(node.gt.iponoelmax) then
          if(ndirforc(i).gt.3) then
            write(*,*) '*WARNING: in gen3dforc: node ',i,
     &           ' does not'
            write(*,*) '       belong to a beam nor shell'
            write(*,*) '       element and consequently has no'
            write(*,*) '       rotational degrees of freedom'
          endif
          cycle
        endif
        index=iponoel(node)
!     
!     next loop is needed, since an element may have
!     been deactivated
!     
        do
          if(index.eq.0) then
            if(ndirforc(i).gt.3) then
              write(*,*) '*WARNING: in gen3dforc: node ',i,
     &             ' does not'
              write(*,*) '       belong to a beam nor shell'
              write(*,*) '       element and consequently has no'
              write(*,*) '       rotational degrees of freedom'
            endif
            cycle loop
          endif
          ielem=inoel(1,index)
          indexe=ipkon(ielem)
          if(indexe.ge.0) exit
          index=inoel(3,index)
        enddo
!     
!     checking whether element is linear or quadratic
!     
        if((lakon(ielem)(4:4).eq.'6').or.
     &       (lakon(ielem)(4:4).eq.'8')) then
          quadratic=.false.
        else
          quadratic=.true.
        endif
!     
!     checking whether element is 3-sided or 4-sided
!     
        if((lakon(ielem)(4:4).eq.'6').or.
     &       (lakon(ielem)(4:5).eq.'15')) then
          nedge=3
        else
          nedge=4
        endif
!     
        j=inoel(2,index)
c     indexe=ipkon(ielem)
        indexk=iponor(2,indexe+j)
        if(nam.gt.0) iamplitude=iamforc(i)
        idir=ndirforc(i)
        val=xforc(i)
!     
        if(rig(node).ne.0) then
          if(idir.gt.3) then
            if(rig(node).lt.0) then
              write(*,*) '*ERROR in gen3dforc: in node ',node
              write(*,*) '       a rotational DOF is loaded;'
              write(*,*) '       however, the elements to which'
              write(*,*) '       this node belongs do not have'
              write(*,*) '       rotational DOFs'
              call exit(201)
            endif
c     val=xforc(i)
            k=idir-3
            irotnode=rig(node)
            call forcadd(irotnode,k,val,nodeforc,
     &           ndirforc,xforc,nforc,nforc_,iamforc,
     &           iamplitude,nam,ntrans,trab,inotr,co,
     &           ikforc,ilforc,isector,add,user,idefforc,
     &           ipompc,nodempc,nmpc,ikmpc,ilmpc,labmpc)
          endif
        else
!     
!     check for moments defined in any but the first step
!     nonlinear dynamic case: creation of knots
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
!     determine the uniform expansion
!     
            alpha=0.d0
            do k=1,ndepnodes
              nod=idepnodes(k)
              dd=(co(1,nod)-q(1))**2
     &             +(co(2,nod)-q(2))**2
     &             +(co(3,nod)-q(3))**2
              if(dd.lt.1.d-20) cycle
              alpha=alpha+dsqrt(
     &             ((co(1,nod)+vold(1,nod)-q(1)-w(1))**2
     &             +(co(2,nod)+vold(2,nod)-q(2)-w(2))**2
     &             +(co(3,nod)+vold(3,nod)-q(3)-w(3))**2)/dd)
            enddo
            alpha=alpha/ndepnodes
!     
!     determine the displacements of irotnodes
!     
            do l=1,3
              do m=1,3
                a(l,m)=0.d0
              enddo
              xn(l)=0.d0
            enddo
            do k=1,ndepnodes
              nod=idepnodes(k)
              dd=0.d0
              do l=1,3
                a1(l)=co(l,nod)-q(l)
                a2(l)=vold(l,nod)-w(l)
                dd=dd+a1(l)*a1(l)
              enddo
              dd=dsqrt(dd)
              if(dd.lt.1.d-10) cycle
              do l=1,3
                a1(l)=a1(l)/dd
                a2(l)=a2(l)/dd
              enddo
              xn(1)=xn(1)+(a1(2)*a2(3)-a1(3)*a2(2))
              xn(2)=xn(2)+(a1(3)*a2(1)-a1(1)*a2(3))
              xn(3)=xn(3)+(a1(1)*a2(2)-a1(2)*a2(1))
              do l=1,3
                do m=1,3
                  a(l,m)=a(l,m)+a1(l)*a1(m)
                enddo
              enddo
            enddo
!     
            do l=1,3
              do m=1,3
                a(l,m)=a(l,m)/ndepnodes
              enddo
              xn(l)=xn(l)/ndepnodes
              a(l,l)=1.d0-a(l,l)
            enddo
!     
            m=3
            nrhs=1
            call dgesv(m,nrhs,a,m,ipiv,xn,m,info)
            if(info.ne.0) then
              write(*,*) '*ERROR in gen3dforc:'
              write(*,*) '       singular system of equations'
              call exit(201)
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
!     rotation matrix r
!     
            do k=1,3
              do l=1,3
                r(k,l)=c1*d(k,l)+
     &               c2*(e(k,1,l)*xn(1)+e(k,2,l)*xn(2)+
     &               e(k,3,l)*xn(3))+c3*xn(k)*xn(l)
              enddo
            enddo
!     
!     copying the displacements
!     
            do l=1,3
              vold(l,irefnode)=w(l)
              vold(l,irotnode)=xn(l)
            enddo
c     vold(1,iexpnode)=alpha
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
!     solving a least squares problem to determine 
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
     &                 b(l,3)*r(3,m)
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
     &               '*ERROR in knotmpc while calculating the'
                write(*,*) '       eigenvalues/eigenvectors'
                call exit(201)
              endif
!     
              if((dabs(w(1)-1.d0).lt.dabs(w(2)-1.d0)).and.
     &             (dabs(w(1)-1.d0).lt.dabs(w(3)-1.d0))) then
                l=2
                m=3
              elseif((dabs(w(2)-1.d0).lt.dabs(w(1)-1.d0)).and.
     &               (dabs(w(2)-1.d0).lt.dabs(w(3)-1.d0))) then
                l=1
                m=3
              else
                l=1
                m=2
              endif
              xi1=datan2((z(1,l)*e2(1)+z(2,l)*e2(2)+z(3,l)*e2(2)),
     &             (z(1,l)*e1(1)+z(2,l)*e1(2)+z(3,l)*e1(2)))
              xi2=w(l)-1.d0
              xi3=w(m)-1.d0
!     
              vold(1,iexpnode)=xi1
              vold(2,iexpnode)=xi2
              vold(3,iexpnode)=xi3
            endif
!     
!     apply the moment
!     
            idir=idir-3
            val=xforc(i)
            call forcadd(irotnode,idir,val,nodeforc,
     &           ndirforc,xforc,nforc,nforc_,iamforc,
     &           iamplitude,nam,ntrans,trab,inotr,co,
     &           ikforc,ilforc,isector,add,user,idefforc,
     &           ipompc,nodempc,nmpc,ikmpc,ilmpc,labmpc)
!     
!     check for shells whether the rotation about the normal
!     on the shell has been eliminated
!     
            if(lakon(ielem)(7:7).eq.'L') then
              indexx=iponor(1,indexe+j)
              do k=1,3
                xnoref(k)=xnor(indexx+k)
              enddo
              dmax=0.d0
              imax=0
              do k=1,3
                if(dabs(xnoref(k)).gt.dmax) then
                  dmax=dabs(xnoref(k))
                  imax=k
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
     &                   '*ERROR in gen3dnor: increase nmpc_'
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
!     all cases except nonlinear dynamic case: creation
!     of mean rotation MPC's
!     
          if((idir.gt.3).and.((nmethod.ne.4).or.(iperturb(1).le.1)))then
!     
!     create a mean rotation MPC
!     advantage: more accurate since less constraining
!     disadvantage: cannot exceed 90 degrees rotation
!     
!     if a mean rotation MPC has already been created 
!     for idof, ilforc(id) contains the index of the
!     SPC created for the angle value
!     
            idof=8*(node-1)+idir
            call nident(ikforc,idof,nforc,id)
            if(ilforc(id).ne.i) cycle

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
              if(dot.gt.999d0) then
                write(*,*) '*WARNING in gen3dforc: applied'
                write(*,*) '         moment in node ',node
                write(*,*) '         and direction ',idir-1
                write(*,*) '         is nearly'
                write(*,*) 
     &               '         along the drilling'
                write(*,*) '         direction; this load is'
                write(*,*) '         taken into account'
                write(*,*)
                cycle
              endif
            endif
!     
!     specific label for mean rotations for beams and
!     shells
!     
            label='MEANROTBS           '
            write(27,*) 'a MEAN ROTATION MPC was generated in node ',
     &           node
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
              call forcadd(nk,idir,val,nodeforc,
     &             ndirforc,xforc,nforc,nforc_,iamforc,
     &             iamplitude,nam,ntrans,trab,inotr,co,
     &             ikforc,ilforc,isector,add,user,idefforc,
     &             ipompc,nodempc,nmpc,ikmpc,ilmpc,labmpc)
!     
!     storing the index of the SPC with the angle
!     value in ilboun(id)
!     
              ilforc(id)=nforc
            endif
!     
            cycle
          endif
!     
!     2d shell element: generate MPC's
!     
          if(lakon(ielem)(7:7).eq.'L') then
            newnode=knor(indexk+1)
!     
!     if a transformation applies to the node all
!     dofs have to be connected
!     
            if(ntrans.le.0) then
              idirstart=idir
              idirend=idir
            elseif(inotr(1,node).eq.0) then
              idirstart=idir
              idirend=idir
            else
              idirstart=1
              idirend=3
            endif
            do idirl=idirstart,idirend
              idof=8*(newnode-1)+idirl
              call nident(ikmpc,idof,nmpc,id)
              if((id.le.0).or.(ikmpc(id).ne.idof)) then
                nmpc=nmpc+1
                if(nmpc.gt.nmpc_) then
                  write(*,*) 
     &                 '*ERROR in gen3dforc: increase nmpc_'
                  call exit(201)
                endif
                labmpc(nmpc)='                    '
                ipompc(nmpc)=mpcfree
                do k=nmpc,id+2,-1
                  ikmpc(k)=ikmpc(k-1)
                  ilmpc(k)=ilmpc(k-1)
                enddo
                ikmpc(id+1)=idof
                ilmpc(id+1)=nmpc
!     
!     for middle nodes: u_1+u_3-2*u_node=0
!     for end nodes: -u_1+4*u_2-u_3-2*u_node=0
!     
!     u_1 corresponds to knor(indexk+1)....
!     
                nodempc(1,mpcfree)=newnode
                nodempc(2,mpcfree)=idirl
                if((j.gt.nedge).or.(.not.quadratic)) then
                  coefmpc(mpcfree)=1.d0
                else
                  coefmpc(mpcfree)=-1.d0
                endif
                mpcfree=nodempc(3,mpcfree)
                if(mpcfree.eq.0) then
                  write(*,*) 
     &                 '*ERROR in gen3dforc: increase memmpc_'
                  call exit(201)
                endif
!     
                if((j.le.nedge).and.(quadratic)) then
                  nodempc(1,mpcfree)=knor(indexk+2)
                  nodempc(2,mpcfree)=idirl
                  coefmpc(mpcfree)=4.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                    write(*,*) 
     &                   '*ERROR in gen3dforc: increase memmpc_'
                    call exit(201)
                  endif
                endif
!     
                nodempc(1,mpcfree)=knor(indexk+3)
                nodempc(2,mpcfree)=idirl
                if((j.gt.nedge).or.(.not.quadratic)) then
                  coefmpc(mpcfree)=1.d0
                else
                  coefmpc(mpcfree)=-1.d0
                endif
                mpcfree=nodempc(3,mpcfree)
                if(mpcfree.eq.0) then
                  write(*,*) 
     &                 '*ERROR in gen3dforc: increase memmpc_'
                  call exit(201)
                endif
!     
                nodempc(1,mpcfree)=node
                nodempc(2,mpcfree)=idirl
                coefmpc(mpcfree)=-2.d0
                mpcfreenew=nodempc(3,mpcfree)
                if(mpcfreenew.eq.0) then
                  write(*,*) 
     &                 '*ERROR in gen3dforc: increase memmpc_'
                  call exit(201)
                endif
!     
                nodempc(3,mpcfree)=0
                mpcfree=mpcfreenew
              endif
            enddo
          elseif(lakon(ielem)(7:7).eq.'B') then
!     
!     1d beam element: generate MPC's
!     
            newnode=knor(indexk+1)
!     
!     if a transformation applies to the node all
!     dofs have to be connected
!     
            if(ntrans.le.0) then
              idirstart=idir
              idirend=idir
            elseif(inotr(1,node).eq.0) then
              idirstart=idir
              idirend=idir
            else
              idirstart=1
              idirend=3
            endif
            do idirl=idirstart,idirend
              idof=8*(newnode-1)+idirl
              call nident(ikmpc,idof,nmpc,id)
              if((id.le.0).or.(ikmpc(id).ne.idof)) then
                nmpc=nmpc+1
                if(nmpc.gt.nmpc_) then
                  write(*,*) 
     &                 '*ERROR in gen3dforc: increase nmpc_'
                  call exit(201)
                endif
                labmpc(nmpc)='                    '
                ipompc(nmpc)=mpcfree
                do k=nmpc,id+2,-1
                  ikmpc(k)=ikmpc(k-1)
                  ilmpc(k)=ilmpc(k-1)
                enddo
                ikmpc(id+1)=idof
                ilmpc(id+1)=nmpc
                do k=1,4
                  nodempc(1,mpcfree)=knor(indexk+k)
                  nodempc(2,mpcfree)=idirl
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                    write(*,*) 
     &                   '*ERROR in gen3dforc: increase memmpc_'
                    call exit(201)
                  endif
                enddo
                nodempc(1,mpcfree)=node
                nodempc(2,mpcfree)=idirl
                coefmpc(mpcfree)=-4.d0
                mpcfreenew=nodempc(3,mpcfree)
                if(mpcfreenew.eq.0) then
                  write(*,*) 
     &                 '*ERROR in gen3dforc: increase memmpc_'
                  call exit(201)
                endif
                nodempc(3,mpcfree)=0
                mpcfree=mpcfreenew
              endif
            enddo
          else
!     
!     2d plane stress, plane strain or axisymmetric
!     element: MPC in all but z-direction
!     
            newnode=knor(indexk+2)
!     
            if(idir.eq.3) cycle
!     
!     if a transformation applies to the node all
!     dofs have to be connected
!     
            if(ntrans.le.0) then
              idirstart=idir
              idirend=idir
            elseif(inotr(1,node).eq.0) then
              idirstart=idir
              idirend=idir
            else
              idirstart=1
              idirend=2
            endif
            do idirl=idirstart,idirend
              idof=8*(newnode-1)+idirl
              call nident(ikmpc,idof,nmpc,id)
              if((id.le.0).or.(ikmpc(id).ne.idof)) then
                nmpc=nmpc+1
                if(nmpc.gt.nmpc_) then
                  write(*,*) 
     &                 '*ERROR in gen3dmpc: increase nmpc_'
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
                nodempc(2,mpcfree)=idirl
                coefmpc(mpcfree)=1.d0
                mpcfree=nodempc(3,mpcfree)
                if(mpcfree.eq.0) then
                  write(*,*) 
     &                 '*ERROR in gen3dmpc: increase memmpc_'
                  call exit(201)
                endif
                nodempc(1,mpcfree)=node
                nodempc(2,mpcfree)=idirl
                coefmpc(mpcfree)=-1.d0
                mpcfreenew=nodempc(3,mpcfree)
                if(mpcfreenew.eq.0) then
                  write(*,*) 
     &                 '*ERROR in gen3dmpc: increase memmpc_'
                  call exit(201)
                endif
                nodempc(3,mpcfree)=0
                mpcfree=mpcfreenew
              endif
            enddo
          endif
        endif
      enddo loop
!     
      return
      end
