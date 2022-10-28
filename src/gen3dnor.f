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
      subroutine gen3dnor(nk,nk_,co,iponoel,inoel,iponoelmax,kon,ipkon,
     &  lakon,ne,thicke,offset,iponor,xnor,knor,rig,iperturb,tinc,
     &  tper,tmin,tmax,ctrl,ipompc,nodempc,coefmpc,nmpc,nmpc_,mpcfree,
     &  ikmpc,ilmpc,labmpc,ikboun,ilboun,nboun,nboun_,nodeboun,ndirboun,
     &  xboun,iamboun,typeboun,nam,ntrans,inotr,trab,ikfree,ixfree,
     &  nmethod,ithermal,istep,mi,icomposite,ielmat,vold,iflagpl)
!     
!     calculates normals on 1-D and 2-D elements
!     
!     degrees of freedom:
!     0: Temperature
!     1: tra_in_x
!     2: tra_in_y
!     3: tra_in_z
!     4: pressure (only CFD) or rot_about_x (else)
!     5: rot_about_y
!     6: rot_about_z
!     
      implicit none
!     
      logical fixed,composite,beam
!     
      character*1 type,typeboun(*)
      character*8 lakon(*)
      character*20 labmpc(*),label
!     
      integer nk,nk_,iponoel(*),inoel(3,*),iponoelmax,kon(*),ipkon(*),
     &     ne,iponor(2,*),knor(*),rig(*),iperturb(*),ipompc(*),
     &     nmpc,nmpc_,mpcfree,ikmpc(*),ilmpc(*),ikboun(*),ilboun(*),
     &     nboun_,nodeboun(*),ndirboun(*),iamboun(*),nam,ntrans,
     &     isol,istep,idummy,mi(*),icomposite,ielmat(mi(3),*),nkold,
     &     i,ndepnodes,index,nexp,nnor,nel,ielem,indexe,j,iel(999),
     &     jl(999),ial(999),ifi(999),idepnodes(800),indexx,k,l,ifix,
     &     jact,ixfree,ikfree,node,nelshell,irefnode,idof,id,mpcfreeold,
     &     irotnode,imax,iamplitude,nmethod,ithermal(*),iexpnode,idim,
     &     iflagpl,nodempc(3,*),nboun,inotr(2,*),idmpc,nemin
!     
      real*8 co(3,*),thicke(mi(3),*),offset(2,*),xnor(*),tinc,tper,tmin,
     &     tmax,ctrl(*),coefmpc(*),xboun(*),trab(7,*),vold(0:mi(2),*),
     &     xno(3,999),xta(3,999),xn1(3,999),thl1(999),thl2(999),
     &     off1(999),off2(999),xi,et,coloc6(2,6),coloc8(2,8),xl(3,8),
     &     dd,xnoref(3),dot,coloc3(3),dot1,dot2,dmax,val,coloc2(2),
     &     e1(3),e2(3),t1(3)
!     
      data coloc2 /-1.d0,1.d0/
      data coloc3 /-1.d0,0.d0,1.d0/
      data coloc6 /0.d0,0.d0,1.d0,0.d0,0.d0,1.d0,0.5d0,0.d0,
     &     0.5d0,0.5d0,0.d0,0.5d0/
      data coloc8 /-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0,1.d0,
     &     0.d0,-1.d0,1.d0,0.d0,0.d0,1.d0,-1.d0,0.d0/
!     
      label='                    '
      fixed=.false.
!     
!     calculating the normals in nodes belonging to shells/beams
!     
      nkold=nk
!     
      do i=1,nkold
        ndepnodes=0
        idim=0
        index=iponoel(i)
        if(index.eq.0) cycle
!     
!     nexp indicates how many times the node was expanded
!     
        nexp=0
!     
!     nnor indicates whether the expanded nodes lie on a point
!     (nnor=0, only for plane stress, plane strain or axisymmetric
!     elements), on a line (nnor=1) or in a plane (nnor=2)
!     
        nnor=0
!     
!     locating the shell elements to which node i belongs
!     
        nel=0
        do
          if(index.eq.0) exit
          ielem=inoel(1,index)
          if((lakon(ielem)(1:1).ne.'B').and.
     &         (lakon(ielem)(1:1).ne.'T')) then
            if((lakon(ielem)(1:1).eq.'S').or.
     &           (lakon(ielem)(1:1).eq.'M')) nnor=1
            indexe=ipkon(ielem)
            nel=nel+1
            if(nel.gt.999) then
              write(*,*) '*ERROR in gen3dnor: more than 999'
              write(*,*) '       shell elements share the'
              write(*,*) '       same node'
              call exit(201)
            endif
            j=inoel(2,index)
            jl(nel)=j
            iel(nel)=ielem
            thl1(nel)=thicke(1,indexe+j)
!     
!     correcting the thickness for more than one layer (composite)
!     
            do k=2,mi(3)
              thl1(nel)=thl1(nel)+thicke(k,indexe+j)
            enddo
!     
            off1(nel)=offset(1,ielem)
          endif
          index=inoel(3,index)
        enddo
!     
        if(nel.gt.0) then
          do j=1,nel
            ial(j)=0
          enddo
!     
!     estimate the normal
!     
          do j=1,nel
            indexe=ipkon(iel(j))
            indexx=iponor(1,indexe+jl(j))
            if(indexx.ge.0) then
              do k=1,3
                xno(k,j)=xnor(indexx+k)
              enddo
              ifi(j)=1
              cycle
            else
              ifi(j)=0
            endif
!     
!     local normal on the element (Jacobian)
!     
            if((lakon(iel(j))(1:2).eq.'S3').or.
     &           (lakon(iel(j))(4:4).eq.'3')) then
              xi=coloc6(1,jl(j))
              et=coloc6(2,jl(j))
              do k=1,3
                node=kon(indexe+k)
                do l=1,3
                  xl(l,k)=co(l,node)
                enddo
              enddo
              call norshell3(xi,et,xl,xno(1,j))
            elseif((lakon(iel(j))(2:2).eq.'4').or.
     &             (lakon(iel(j))(4:4).eq.'4')) then
              xi=coloc8(1,jl(j))
              et=coloc8(2,jl(j))
              do k=1,4
                node=kon(indexe+k)
                do l=1,3
                  xl(l,k)=co(l,node)
                enddo
              enddo
              call norshell4(xi,et,xl,xno(1,j))
            elseif((lakon(iel(j))(2:2).eq.'6').or.
     &             (lakon(iel(j))(4:4).eq.'6')) then
              xi=coloc6(1,jl(j))
              et=coloc6(2,jl(j))
              do k=1,6
                node=kon(indexe+k)
                do l=1,3
                  xl(l,k)=co(l,node)
                enddo
              enddo
              call norshell6(xi,et,xl,xno(1,j))
            elseif((lakon(iel(j))(2:2).eq.'8').or.
     &             (lakon(iel(j))(4:4).eq.'8')) then
              xi=coloc8(1,jl(j))
              et=coloc8(2,jl(j))
              do k=1,8
                node=kon(indexe+k)
                do l=1,3
                  xl(l,k)=co(l,node)
                enddo
              enddo
              call norshell8(xi,et,xl,xno(1,j))
            endif
!     
            dd=dsqrt(xno(1,j)**2+xno(2,j)**2+xno(3,j)**2)
            if(dd.lt.1.d-10) then
              write(*,*) '*ERROR in gen3dnor: size of estimated'
              write(*,*) '       shell normal in node ',i,
     &             ' element ',iel(j)
              write(*,*) '       is smaller than 1.e-10'
              call exit(201)
            endif
            do k=1,3
              xno(k,j)=xno(k,j)/dd
            enddo
          enddo
!     
          do
!     
!     determining a fixed normal which was not treated yet,
!     or, if none is left, the minimum element number of all
!     elements containing node i and for which no normal was
!     determined yet
!     
!     if ial(j)=0: the normal on this element has not been
!     treated yet
!     if ial(j)=2: normal has been treated
!     
            ifix=0
            nemin=ne+1
            do j=1,nel
              if(ial(j).ne.0) cycle
              if(ifi(j).eq.1) then
                jact=j
                ifix=1
                exit
              endif
            enddo
            if(ifix.eq.0) then
              do j=1,nel
                if(ial(j).eq.0) then
                  if(iel(j).lt.nemin) then
                    nemin=iel(j)
                    jact=j
                  endif
                endif
              enddo
              if(nemin.eq.ne+1) exit
            endif
!     
            do j=1,3
              xnoref(j)=xno(j,jact)
            enddo
!     
!     determining all elements whose normal in node i makes an
!     angle smaller than 0.5 or 20 degrees with the reference normal,
!     depending whether the reference normal was given by the
!     user or is being calculated; the thickness and offset must
!     also fit.
!     
!     if ial(j)=1: normal on element is being treated now
!     
            do j=1,nel
              if(ial(j).eq.2) cycle
              if(j.eq.jact) then
                ial(jact)=1
              else
                dot=xno(1,j)*xnoref(1)+xno(2,j)*xnoref(2)+
     &               xno(3,j)*xnoref(3)
                if(ifix.eq.0) then
                  if(dot.gt.0.939693d0)then
                    if((dabs(thl1(j)-thl1(jact)).lt.1.d-10)
     &                   .and.
     &                   (dabs(off1(j)-off1(jact)).lt.1.d-10)
     &                   .and.
     &                   (lakon(iel(j))(1:1).ne.'M').and.
     &                   (lakon(iel(jact))(1:1).ne.'M').and.
     &                   ((lakon(iel(j))(1:3).eq.lakon(iel(jact))(1:3))
     &                   .or.
     &                   ((lakon(iel(j))(1:1).eq.'S').and.
     &                   (lakon(iel(jact))(1:1).eq.'S'))))
     &                   ial(j)=1
c     
                    if(dot.lt.0.999962d0) nnor=2
c     
                  else
                    if((lakon(iel(j))(1:1).eq.'S').and.
     &                   (lakon(iel(jact))(1:1).eq.'S')) then
!     
!     if the normals have the opposite
!     direction, the expanded nodes are on a
!     straight line
!     
                      if(dot.gt.-0.999962d0) then
                        nnor=2
                      else
                        write(*,*) '*INFO in gen3dnor: in node ',i, 
     &                       ' opposite normals are defined'
                      endif
                    endif
                  endif
                else
                  if(dot.gt.0.999962d0) then
                    if((dabs(thl1(j)-thl1(jact)).lt.1.d-10)
     &                   .and.
     &                   (dabs(off1(j)-off1(jact)).lt.1.d-10)
     &                   .and.
     &                   (lakon(iel(j))(1:1).ne.'M').and.
     &                   (lakon(iel(jact))(1:1).ne.'M').and.
     &                   ((lakon(iel(j))(1:3).eq.lakon(iel(jact))(1:3))
     &                   .or.
     &                   ((lakon(iel(j))(1:1).eq.'S').and.
     &                   (lakon(iel(jact))(1:1).eq.'S'))))
     &                   ial(j)=1
c     
c     if(dot.lt.0.999962) nnor=2
c     
                  else
                    if((lakon(iel(j))(1:1).eq.'S').and.
     &                   (lakon(iel(jact))(1:1).eq.'S')) then
!     
!     if the normals have the opposite
!     direction, the expanded nodes are on a
!     straight line
!     
                      if(dot.gt.-0.999962d0) then
                        nnor=2
                      else
                        write(*,*) '*INFO in gen3dnor: in node ',i,
     &                       ' opposite normals are defined'
                      endif
                    endif
                  endif
                endif
              endif
            enddo
!     
!     determining the mean normal for the selected elements
!     
            if(ifix.eq.0) then
              do j=1,3
                xnoref(j)=0.d0
              enddo
              do j=1,nel
                if(ial(j).eq.1) then
                  do k=1,3
                    xnoref(k)=xnoref(k)+xno(k,j)
                  enddo
                endif
              enddo
              dd=dsqrt(xnoref(1)**2+xnoref(2)**2+xnoref(3)**2)
              if(dd.lt.1.d-10) then
                write(*,*) '*ERROR in gen3dnor: size of'
                write(*,*) '        estimated shell normal is'
                write(*,*) '        smaller than 1.e-10'
                call exit(201)
              endif
              do j=1,3
                xnoref(j)=xnoref(j)/dd
              enddo
            endif
!     
!     check whether any elements are composite elements
!     
            composite=.false.
            if(icomposite.eq.1) then
              do j=1,nel
                if(ial(j).eq.1) then
                  ielem=iel(j)
                  if(ielmat(2,ielem).ne.0) then
                    composite=.true.
                    exit
                  endif
                endif
              enddo
            endif
!     
!     updating the pointers iponor
!     
            nexp=nexp+1
            do j=1,nel
              if(ial(j).eq.1) then
                ial(j)=2
                if(ifix.eq.0) then
                  iponor(1,ipkon(iel(j))+jl(j))=ixfree
                elseif(j.ne.jact) then
                  iponor(1,ipkon(iel(j))+jl(j))=
     &                 iponor(1,ipkon(iel(jact))+jl(jact)) 
                endif
                iponor(2,ipkon(iel(j))+jl(j))=ikfree
              endif
            enddo
!     
!     storing the normal in xnor and generating 3 nodes
!     for knor
!     
            if(ifix.eq.0) then
              do j=1,3
                xnor(ixfree+j)=xnoref(j)
              enddo
              ixfree=ixfree+3
            endif
!     
            do k=1,3
              nk=nk+1
              if(nk.gt.nk_) then
                write(*,*) '*ERROR in gen3dnor: increase nk_'
                call exit(201)
              endif
              knor(ikfree+k)=nk
!     
!     for plane stress, plane strain and axisymmetric
!     elements only the middle node is included in the
!     rigid body definition
!     
              if((lakon(iel(jact))(2:2).ne.'P').and.
     &             (lakon(iel(jact))(2:2).ne.'A').and.
     &             (lakon(iel(jact))(1:1).ne.'M')) then
                idepnodes(ndepnodes+1)=nk
                ndepnodes=ndepnodes+1
                idim=max(idim,1)
              elseif(k.eq.2) then
                if (lakon(iel(jact))(1:1).ne.'M') then
!     
!     plane stress/strain axisymmetric
!     
                  idepnodes(ndepnodes+1)=nk
                  ndepnodes=ndepnodes+1
                  idim=max(idim,1)
                else
!     
!     membrane: insert hinge
!     
                  call gen3dmembrane(ipompc,nodempc,coefmpc,nmpc,
     &                 nmpc_,mpcfree,ikmpc,ilmpc,labmpc,nk,
     &                 ithermal,i)
                endif
              endif
            enddo
            ikfree=ikfree+3
!     
!     extra nodes for composite shells
!     
            if(composite) then
              do l=1,mi(3)
                do k=1,3
                  nk=nk+1
                  if(nk.gt.nk_) then
                    write(*,*) '*ERROR in gen3dnor: increase nk_'
                    call exit(201)
                  endif
                  knor(ikfree+k)=nk
                enddo
                ikfree=ikfree+3
              enddo
            endif
!     
          enddo
        endif
!     
        nelshell=nel+1
!     
!     locating the beam elements to which node i belongs
!     
        index=iponoel(i)
        do
          if(index.eq.0) exit
          ielem=inoel(1,index)
          if((lakon(ielem)(1:1).eq.'B').or.
     &         (lakon(ielem)(1:1).eq.'T')) then
            if(lakon(ielem)(1:1).eq.'B') then
              beam=.true.
            else
              beam=.false.
            endif
            indexe=ipkon(ielem)
            nel=nel+1
            if(nel.gt.999) then
              write(*,*) '*ERROR in gen3dnor: more than 999'
              write(*,*) '        beam/shell elements share'
              write(*,*) '        the same node'
              call exit(201)
            endif
            j=inoel(2,index)
            jl(nel)=j
            iel(nel)=ielem
            thl1(nel)=thicke(1,indexe+j)
            thl2(nel)=thicke(2,indexe+j)
            off1(nel)=offset(1,ielem)
            off2(nel)=offset(2,ielem)
          endif
          index=inoel(3,index)
        enddo
!     
        if(nel.ge.nelshell) then
          if(beam) nnor=2
          do j=nelshell,nel
            ial(j)=0
          enddo
!     
!     estimate the normal
!     
          do j=nelshell,nel
            if((lakon(iel(j))(3:3).eq.'1').or.
     &           (lakon(iel(j))(4:4).eq.'2')) then
!     
!     linear beam or truss element
!     
              xi=coloc2(jl(j))
              indexe=ipkon(iel(j))
              do k=1,2
                node=kon(indexe+k)
                do l=1,3
                  xl(l,k)=co(l,node)
                enddo
              enddo
!     
!     determining the tangent vector xta
!     
              do k=1,3
                xta(k,j)=(xl(k,2)-xl(k,1))/2.d0
              enddo
            else
!     
!     quadratic beam element
!     
              xi=coloc3(jl(j))
              indexe=ipkon(iel(j))
              do k=1,3
                node=kon(indexe+k)
                do l=1,3
                  xl(l,k)=co(l,node)
                enddo
              enddo
!     
!     determining the tangent vector xta
!     
              do k=1,3
                xta(k,j)=(xi-0.5d0)*xl(k,1)-2.d0*xi*xl(k,2)+
     &               (xi+0.5d0)*xl(k,3)
              enddo
            endif
!     
            dd=dsqrt(xta(1,j)**2+xta(2,j)**2+xta(3,j)**2)
            if(dd.lt.1.d-10) then
              write(*,*) '*ERROR in gen3dnor: size of estimated'
              write(*,*)'       beam tangent in node ',i,' element '
     &             ,iel(j)
              write(*,*) '       is smaller than 1.e-10'
              call exit(201)
            endif
            do k=1,3
              xta(k,j)=xta(k,j)/dd
            enddo
!     
!     check whether normal was defined with *NORMAL and
!     determine vector n1
!     
            if(iponor(1,indexe+jl(j)).ge.0) then
              indexx=iponor(1,indexe+jl(j))
              if(dabs(xnor(indexx+4)**2+xnor(indexx+5)**2+
     &             xnor(indexx+6)**2-1.d0).lt.1.d-5) then
                do k=1,3
                  xno(k,j)=xnor(indexx+3+k)
                enddo
                ifi(j)=1
                cycle
              endif
              ifi(j)=0
              do k=1,3
                xn1(k,j)=xnor(indexx+k)
              enddo
            else
              ifi(j)=0
              xn1(1,j)=0.d0
              xn1(2,j)=0.d0
              xn1(3,j)=-1.d0
            endif
!     
!     normal (=n2) = xta x xn1
!     
            xno(1,j)=xta(2,j)*xn1(3,j)-xta(3,j)*xn1(2,j)
            xno(2,j)=xta(3,j)*xn1(1,j)-xta(1,j)*xn1(3,j)
            xno(3,j)=xta(1,j)*xn1(2,j)-xta(2,j)*xn1(1,j)
            dd=dsqrt(xno(1,j)**2+xno(2,j)**2+xno(3,j)**2)
            if(dd.lt.1.d-10) then
              write(*,*) '*ERROR in gen3dnor: size of estimated'
              write(*,*)'       beam normal in 2-direction in node '
     &             ,i,' element ',iel(j)
              write(*,*) '       is smaller than 1.e-10'
              call exit(201)
            endif
            do k=1,3
              xno(k,j)=xno(k,j)/dd
            enddo
          enddo
!     
          do
!     
!     determining a fixed normal which was not treated yet,
!     or, if none is left, the minimum element number of all
!     elements containing node i and for which no normal was
!     determined yet
!     
            ifix=0
            nemin=ne+1
            do j=nelshell,nel
              if(ial(j).ne.0) cycle
              if(ifi(j).eq.1) then
                jact=j
                ifix=1
                exit
              endif
            enddo
            if(ifix.eq.0) then
              do j=nelshell,nel
                if(ial(j).eq.0) then
                  if(iel(j).lt.nemin) then
                    nemin=iel(j)
                    jact=j
                  endif
                endif
              enddo
              if(nemin.eq.ne+1) exit
            endif
!     
!     the reference normal is the one on the element with the
!     smallest element number
!     
            do j=1,3
              xnoref(j)=xno(j,jact)
            enddo
!     
!     determining all elements whose normal in node i makes an
!     angle smaller than 0.5 or 20 degrees with the reference normal,
!     depending whether the reference normal was given by the
!     user or is being calculated; the thickness and offset must
!     also fit.
!     
            do j=nelshell,nel
              if(j.eq.jact) then
                ial(jact)=1
              else
                dot1=xno(1,j)*xnoref(1)+xno(2,j)*xnoref(2)+
     &               xno(3,j)*xnoref(3)
                dot2=xta(1,j)*xta(1,jact)+xta(2,j)*xta(2,jact)+
     &               xta(3,j)*xta(3,jact)
                if(ifix.eq.0) then
                  if((dot1.gt.0.939693d0).and.
     &                 (dot2.gt.0.939693d0)) then
                    if((dabs(thl1(j)-thl1(jact)).lt.1.d-10)
     &                   .and.
     &                   (dabs(thl2(j)-thl2(jact)).lt.1.d-10)
     &                   .and.
     &                   (dabs(off1(j)-off1(jact)).lt.1.d-10)
     &                   .and.
     &                   (dabs(off2(j)-off2(jact)).lt.1.d-10)
     &                   .and.
     &                   (lakon(iel(j))(1:1).ne.'T').and.
     &                   (lakon(iel(jact))(1:1).ne.'T').and.
     &                   (lakon(iel(j))(8:8).eq.lakon(iel(jact))(8:8)))
     &                   ial(j)=1
                  endif
                else
                  if((dot1.gt.0.999962d0).and.
     &                 (dot2.gt.0.999962d0)) then
                    if((dabs(thl1(j)-thl1(jact)).lt.1.d-10)
     &                   .and.
     &                   (dabs(thl2(j)-thl2(jact)).lt.1.d-10)
     &                   .and.
     &                   (dabs(off1(j)-off1(jact)).lt.1.d-10)
     &                   .and.
     &                   (dabs(off2(j)-off2(jact)).lt.1.d-10)
     &                   .and.
     &                   (lakon(iel(j))(1:1).ne.'T').and.
     &                   (lakon(iel(jact))(1:1).ne.'T').and.
     &                   (lakon(iel(j))(8:8).eq.lakon(iel(jact))(8:8)))
     &                   ial(j)=1
                  endif
                endif
              endif
            enddo
!     
!     determining the mean normal for the selected elements
!     
            if(ifix.eq.0) then
              do j=1,3
                xnoref(j)=0.d0
              enddo
              do j=nelshell,nel
                if(ial(j).eq.1) then
                  do k=1,3
                    xnoref(k)=xnoref(k)+xno(k,j)
                  enddo
                endif
              enddo
            endif
!     
!     calculating the mean tangent
!     
            do j=nelshell,nel
              if((ial(j).eq.1).and.(j.ne.jact)) then
                do k=1,3
                  xta(k,jact)=xta(k,jact)+xta(k,j)
                enddo
              endif
            enddo
            dd=dsqrt(xta(1,jact)**2+xta(2,jact)**2+xta(3,jact)**2)
            if(dd.lt.1.d-10) then
              write(*,*) '*ERROR in gen3dnor: size of mean'
              write(*,*)'    beam tangent is smaller than 1.e-10'
              call exit(201)
            endif
            do k=1,3
              xta(k,jact)=xta(k,jact)/dd
            enddo
!     
!     taking care that the mean normal is orthogonal towards
!     the mean tangent
!     
            dd=xnoref(1)*xta(1,jact)+xnoref(2)*xta(2,jact)+
     &           xnoref(3)*xta(3,jact)
            do j=1,3
              xnoref(j)=xnoref(j)-dd*xta(j,jact)
            enddo
            dd=dsqrt(xnoref(1)**2+xnoref(2)**2+xnoref(3)**2)
            if(dd.lt.1.d-10) then
              write(*,*) '*ERROR in gen3dnor: size of'
              write(*,*) '        estimated beam normal is'
              write(*,*) '        smaller than 1.e-10'
              call exit(201)
            endif
            do j=1,3
              xnoref(j)=xnoref(j)/dd
            enddo
!     
!     calculating xn1 = xn2  x tangent              
!     
            xn1(1,jact)=xnoref(2)*xta(3,jact)-xnoref(3)*xta(2,jact)
            xn1(2,jact)=xnoref(3)*xta(1,jact)-xnoref(1)*xta(3,jact)
            xn1(3,jact)=xnoref(1)*xta(2,jact)-xnoref(2)*xta(1,jact)
!     
!     storing the normals in xnor and generating 8 nodes
!     for knor
!     
            nexp=nexp+1
            do j=nelshell,nel
              if(ial(j).eq.1) then
                ial(j)=2
                if(ifix.eq.0) then
                  iponor(1,ipkon(iel(j))+jl(j))=ixfree
                else
                  iponor(1,ipkon(iel(j))+jl(j))=
     &                 iponor(1,ipkon(iel(jact))+jl(jact))
                endif
                iponor(2,ipkon(iel(j))+jl(j))=ikfree
              endif
            enddo
!     
            do j=1,3
              xnor(ixfree+j)=xn1(j,jact)
            enddo
            do j=1,3
              xnor(ixfree+3+j)=xnoref(j)
            enddo
            ixfree=ixfree+6
            if(lakon(iel(jact))(1:1).ne.'T') then
!     
!     beam
!     
              do k=1,8
                nk=nk+1
                if(nk.gt.nk_) then
                  write(*,*) '*ERROR in gen3dnor: increase nk_'
                  call exit(201)
                endif
                knor(ikfree+k)=nk
                idepnodes(ndepnodes+k)=nk
              enddo
              ndepnodes=ndepnodes+8
              idim=max(idim,3)
            else
!     
!     truss
!     
              do k=1,8
                nk=nk+1
                if(nk.gt.nk_) then
                  write(*,*) '*ERROR in gen3dnor: increase nk_'
                  call exit(201)
                endif
                knor(ikfree+k)=nk
!     
!     assigning coordinates (necessary for the 
!     no-rotation-about-the-truss-axis-MPC created
!     in gen3dtruss
!     
                if(k.eq.1) then
                  do j=1,3
                    co(j,nk)=co(j,i)
     &                   -thl1(jact)*xn1(j,jact)*(.5d0+off1(jact))
     &                   +thl2(jact)*xnoref(j)*(.5d0-off2(jact))
                  enddo
                elseif(k.eq.2) then
                  do j=1,3
                    co(j,nk)=co(j,i)
     &                   -thl1(jact)*xn1(j,jact)*(.5d0+off1(jact))
     &                   -thl2(jact)*xnoref(j)*(.5d0+off2(jact))
                  enddo
                elseif(k.eq.3) then
                  do j=1,3
                    co(j,nk)=co(j,i)
     &                   +thl1(jact)*xn1(j,jact)*(.5d0-off1(jact))
     &                   -thl2(jact)*xnoref(j)*(.5d0+off2(jact))
                  enddo
                elseif(k.eq.4) then
                  do j=1,3
                    co(j,nk)=co(j,i)
     &                   +thl1(jact)*xn1(j,jact)*(.5d0-off1(jact))
     &                   +thl2(jact)*xnoref(j)*(.5d0-off2(jact))
                  enddo
                endif
              enddo
!     
              call gen3dtruss(ipompc,nodempc,coefmpc,nmpc,
     &             nmpc_,mpcfree,ikmpc,ilmpc,labmpc,nk,
     &             ithermal,i,nodeboun,ndirboun,ikboun,ilboun,
     &             nboun,nboun_,typeboun,xboun,xta,jact,co,
     &             knor,ntrans,inotr,trab,vold,mi,nmethod,nk_,
     &             nam,iperturb,ikfree,iamboun,iflagpl)
            endif
            ikfree=ikfree+8
          enddo
        endif
!     
!     check whether the user has specified rotational degrees
!     of freedom (in that case rig(i)=-1 was assigned in 
!     subroutine gen3delem); if so, a rigid MPC must be defined
!     
c     if(rig(i).ne.0) then
c     rig(i)=0
c     if(nexp.le.1) then
c     nexp=2
c     endif
c     endif
!     
!     storing the expanded nodes
!     
!     
!     generate rigid MPC's if necessary
!     
        if(nexp.gt.1) then
          irefnode=i
!     
          rig(i)=-1
!     
          if(ithermal(2).ne.2) then
            if(nnor.eq.0) then
!     
!     the node belongs to plane stress, plane strain
!     or axisymmetric elements only. These are only linked
!     through the node in the midplane: the nodes
!     coincide; only DOF1 and DOF2 are linked.
!     rig(i)=-1 to indicate that a knot has
!     been generated without rotational node
!     
              write(27,*) 'a KNOT without rotation was generated'
              write(27,*) '       in node ',i
              write(27,*)
!     
              do k=1,ndepnodes
                node=idepnodes(k)
                do j=1,2
                  idof=8*(node-1)+j
                  call nident(ikmpc,idof,nmpc,id)
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                    write(*,*) 
     &                   '*ERROR in rigidmpc: increase nmpc_'
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
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=j
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  nodempc(1,mpcfree)=irefnode
                  nodempc(2,mpcfree)=j
                  coefmpc(mpcfree)=-1.d0
                  mpcfreeold=mpcfree
                  mpcfree=nodempc(3,mpcfree)
                  nodempc(3,mpcfreeold)=0
                enddo
              enddo
            else
!     
!     generate a rigid body knot; rig(i) contains the
!     rotational node of the knot
!     
              nk=nk+1
              if(nk.gt.nk_) then
                write(*,*) '*ERROR in rigidbodies: increase nk_'
                call exit(201)
              endif
              irotnode=nk
              rig(i)=irotnode
              write(27,*) 'a KNOT was generated in node ',i
              write(27,*)
              nk=nk+1
              if(nk.gt.nk_) then
                write(*,*) '*ERROR in rigidbodies: increase nk_'
                call exit(201)
              endif
              iexpnode=nk
              do k=1,ndepnodes
                call knotmpc(ipompc,nodempc,coefmpc,irefnode,
     &               irotnode,iexpnode,
     &               labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,nk,nk_,
     &               nodeboun,ndirboun,ikboun,ilboun,nboun,nboun_,
     &               idepnodes,typeboun,co,xboun,istep,k,
     &               ndepnodes,idim,e1,e2,t1)
              enddo
            endif
          endif
!     
!     MPC's for the temperature DOFs
!     
          if(ithermal(2).ge.2) then
            do k=1,ndepnodes
              node=idepnodes(k)
              idof=8*(node-1)
              call nident(ikmpc,idof,nmpc,id)
              nmpc=nmpc+1
              if(nmpc.gt.nmpc_) then
                write(*,*) 
     &               '*ERROR in gen3dnor: increase nmpc_'
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
              nodempc(1,mpcfree)=node
              nodempc(2,mpcfree)=0
              coefmpc(mpcfree)=1.d0
              mpcfree=nodempc(3,mpcfree)
              nodempc(1,mpcfree)=irefnode
              nodempc(2,mpcfree)=0
              coefmpc(mpcfree)=-1.d0
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
              nodempc(3,mpcfreeold)=0
            enddo
          endif
!     
          if((nnor.eq.1).and.(ithermal(2).ne.2)) then
!     
!     generate an additional SPC or MPC for rigid body nodes
!     lying on a line to prevent rotation about the
!     line
!     
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
!     
!     check whether the dof has been used as dependent
!     term of a MPC
!     
              idof=8*(i-1)+3+imax
              call nident(ikmpc,idof,nmpc,idmpc)
!     
              if((idmpc.le.0).or.(ikmpc(idmpc).ne.idof)) then
                call bounadd(irotnode,imax,imax,val,nodeboun,
     &               ndirboun,xboun,nboun,nboun_,iamboun,
     &               iamplitude,nam,ipompc,nodempc,coefmpc,
     &               nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &               ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &               type,typeboun,nmethod,iperturb,fixed,vold,
     &               idummy,mi,label)
              endif
            else
!     
!     check whether the rotational degree of freedom
!     imax is fixed through a SPC
!     
              isol=0
              do l=1,3
                idof=8*(i-1)+3+imax
                call nident(ikboun,idof,nboun,id)
                call nident(ikmpc,idof,nmpc,idmpc)
                if(((idmpc.gt.0).and.(ikmpc(idmpc).eq.idof)).or.
     &               ((id.gt.0).and.(ikboun(id).eq.idof)).or.
     &               (dabs(xnoref(imax)).lt.1.d-2)) then
c     changed 24.11.2014
c     &                   (dabs(xnoref(imax)).lt.1.d-10)) then
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
     &                 '*ERROR in gen3dnor: increase nmpc_'
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
        endif
      enddo
!     
      return
      end
