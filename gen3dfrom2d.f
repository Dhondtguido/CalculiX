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
      subroutine gen3dfrom2d(i,kon,ipkon,lakon,ne,iponor,xnor,knor,
     &  thicke,offset,ntrans,inotr,trab,ikboun,ilboun,nboun,nboun_,
     &  nodeboun,ndirboun,xboun,iamboun,typeboun,ipompc,nodempc,coefmpc,
     &  nmpc,nmpc_,mpcfree,ikmpc,ilmpc,labmpc,nk,nk_,co,rig,nmethod,
     &  iperturb,ithermal,mi,nam,icomposite,ielmat,vold,orname,orab,
     &  norien,norien_,ielorien)
!
!     expands 2d element i into a 3d element
!
!     generates additional MPC's for plane stress, plane strain and
!     axisymmetric elements
!
      implicit none
!
      logical axial,fixed,quadratic
!
      character*1 type,typeboun(*)
      character*8 lakon(*)
      character*20 labmpc(*),label
      character*80 orname(*)
!
      integer kon(*),ipkon(*),ne,iponor(2,*),knor(*),ntrans,
     &  inotr(2,*),nodefixz,norien,norien_,norieninput,iflag,
     &  ikboun(*),ilboun(*),nboun,nboun_,nodeboun(*),ndirboun(*),
     &  iamboun(*),nam,ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,
     &  ikmpc(*),ilmpc(*),nk,nk_,i,rig(*),nmethod,iperturb(*),ishift,
     &  indexe,j,nodel(8),indexx,indexk,k,nedge,nodes(3,8),nodec(3,8),
     &  iamplitude,l,newnode,idir,idof,id,m,mpcfreenew,node,ithermal(*),
     &  jmin,jmax,idummy,mi(*),indexc,indexl,icomposite,ielmat(mi(3),*),
     &  nope,neworien(0:norien),iorien,ipos,ielorien(mi(3),*)
! 
      real*8 xnor(*),thicke(mi(3),*),offset(2,*),trab(7,*),xboun(*),
     &  coefmpc(*),co(3,*),thicks(8),xnors(3,8),dc,ds,val,orab(7,*),
     &  x,y,thickness(8),vold(0:mi(2),*),xi,et,xl(3,8),xno(3),xs(3,7),
     &  shp(7,8),p(3),a(3,3),dd,e1(3),e2(3)
!
      norieninput=norien
      iflag=2
!
      label='                    '
      fixed=.false.
!
      if(lakon(i)(1:2).eq.'CA') then
         axial=.true.
      else
         axial=.false.
      endif
!
      indexe=ipkon(i)
!
!           localizing the nodes, thicknesses and normals for the
!           2-D element
!            
      if((lakon(i)(1:2).eq.'S3').or.
     &     (lakon(i)(4:4).eq.'3')) then
         quadratic=.false.
         nedge=3
         nope=3
         ishift=6
      elseif((lakon(i)(2:2).eq.'4').or.
     &     (lakon(i)(4:4).eq.'4')) then
         quadratic=.false.
         nedge=4
         nope=4
!
!        CPE4R, CPS4R, CAX4R and S4R elements are expanded into C3D8R
!
         if((lakon(i)(3:3).eq.'R').or.
     &      (lakon(i)(5:5).eq.'R')) then
            ishift=8
!
!           S4 elements are expanded into C3D8I
!
         elseif(lakon(i)(1:1).eq.'S') then
            ishift=11
!
!           CPE4, CPS4, CAX4 and M3D4 elements are expanded into C3D8
!
         else
            ishift=8
         endif
      elseif((lakon(i)(2:2).eq.'6').or.
     &     (lakon(i)(4:4).eq.'6')) then
         quadratic=.true.
         nedge=3
         nope=6
         ishift=15
      elseif((lakon(i)(2:2).eq.'8').or.
     &     (lakon(i)(4:4).eq.'8')) then
         quadratic=.true.
         nedge=4
         nope=8
         ishift=20
      endif
!
c      do j=1,2*nedge
      do j=1,nope
         nodel(j)=kon(indexe+j)
         kon(indexe+ishift+j)=nodel(j)
         indexk=iponor(2,indexe+j)
         thicks(j)=0.d0
         do k=1,mi(3)
            thicks(j)=thicks(j)+thicke(k,indexe+j)
         enddo
         do k=1,3
            nodes(k,j)=knor(indexk+k)
         enddo
      enddo
!
!           generating the 3-D element topology for shell and plane
!           stress/strain elements
!
      if(lakon(i)(1:2).ne.'CA') then
c         do j=1,2*nedge
         do j=1,nope
            indexx=iponor(1,indexe+j)
            do k=1,3
               xnors(k,j)=xnor(indexx+k)
            enddo
            if(ntrans.gt.0) then
               do k=1,3
                  inotr(1,nodes(k,j))=inotr(1,nodel(j))
               enddo
            endif
         enddo
!
!        vertex nodes
!
         do k=1,nedge
            kon(indexe+k)=nodes(1,k)
!
            do j=1,3
               co(j,nodes(1,k))=co(j,nodel(k))
     &              -thicks(k)*xnors(j,k)*(.5d0+offset(1,i))
            enddo
         enddo
         do k=1,nedge
            kon(indexe+nedge+k)=nodes(3,k)
            do j=1,3
               co(j,nodes(3,k))=co(j,nodel(k))
     &              +thicks(k)*xnors(j,k)*(.5d0-offset(1,i))
            enddo
         enddo
!
!        middle nodes for quadratic expansion
!
         if(quadratic) then
            do k=nedge+1,2*nedge
               kon(indexe+nedge+k)=nodes(1,k)
               do j=1,3
                  co(j,nodes(1,k))=co(j,nodel(k))
     &                 -thicks(k)*xnors(j,k)*(.5d0+offset(1,i))
               enddo
            enddo
            do k=nedge+1,2*nedge
               kon(indexe+2*nedge+k)=nodes(3,k)
               do j=1,3
                  co(j,nodes(3,k))=co(j,nodel(k))
     &                 +thicks(k)*xnors(j,k)*(.5d0-offset(1,i))
               enddo
            enddo
            do k=1,nedge
               kon(indexe+4*nedge+k)=nodes(2,k)
               do j=1,3
                  co(j,nodes(2,k))=co(j,nodel(k))
     &                 -thicks(k)*xnors(j,k)*offset(1,i)
               enddo
            enddo
!
!           generating coordinates for the expanded nodes which
!           are not used by the C3D20(R) element (needed for the
!           determination of the knot dimension)
!
            do k=nedge+1,2*nedge
               do j=1,3
                  co(j,nodes(2,k))=co(j,nodel(k))
     &                 -thicks(k)*xnors(j,k)*offset(1,i)
               enddo
            enddo
         else
!
!           for linear elements the coordinates of the nodes
!           in the middle may be needed in case of knots (to
!           determine the coefficients of the equations)
!
            do k=1,nedge
               do j=1,3
                  co(j,nodes(2,k))=co(j,nodel(k))
     &                 -thicks(k)*xnors(j,k)*offset(1,i)
               enddo
            enddo
         endif
!
!        generating the layer geometry for composite shells
!        only for S8R elements
!
         if(lakon(i)(8:8).eq.'C') then
            do k=1,2*nedge
               thickness(k)=0.d0
            enddo
            indexc=indexe+2*nedge
            indexl=0
            do l=1,mi(3)
               if(ielmat(l,i).eq.0) exit
               indexc=indexc+ishift
               indexl=indexl+3
!
               do j=1,2*nedge
                  do k=1,3
                     nodec(k,j)=knor(iponor(2,indexe+j)+indexl+k)
                  enddo
               enddo
!
!              nodes 1 to 4 (vertex nodes face 1)
!
               do k=1,nedge
                  kon(indexc+k)=nodec(1,k)
                  do j=1,3
                     co(j,nodec(1,k))=co(j,nodel(k))+
     &                    (thickness(k)-thicks(k)*(.5d0+offset(1,i)))
     &                    *xnors(j,k)
                  enddo
               enddo
!
!              nodes 5 to 8 (vertex nodes face 2)
!
               do k=1,nedge
                  kon(indexc+nedge+k)=nodec(3,k)
                  do j=1,3
                     co(j,nodec(3,k))=co(j,nodel(k))+
     &                   (thickness(k)+thicke(l,indexe+k)
     &                    -thicks(k)*(.5d0+offset(1,i)))
     &                   *xnors(j,k)
                  enddo
               enddo
!
!              nodes 9 to 12 (middle nodes face 1)
!
               do k=nedge+1,2*nedge
                  kon(indexc+nedge+k)=nodec(1,k)
                  do j=1,3
                     co(j,nodec(1,k))=co(j,nodel(k))+
     &                    (thickness(k)-thicks(k)*(.5d0+offset(1,i)))
     &                     *xnors(j,k)
                  enddo
               enddo
!
!              nodes 13 to 16 (middle nodes face 2)
!
               do k=nedge+1,2*nedge
                  kon(indexc+2*nedge+k)=
     &              nodec(3,k)
                  do j=1,3
                     co(j,nodec(3,k))=co(j,nodel(k))+
     &                   (thickness(k)+thicke(l,indexe+k)
     &                    -thicks(k)*(.5d0+offset(1,i)))
     &                   *xnors(j,k)
                  enddo
               enddo
!
!              nodes 17 to 20 (middle nodes in between face 1 and 2)
!
               do k=1,nedge
                  kon(indexc+4*nedge+k)=
     &                  nodec(2,k)
                  do j=1,3
                     co(j,nodec(2,k))=co(j,nodel(k))+
     &                    (thickness(k)+thicke(l,indexe+k)/2.d0
     &                    -thicks(k)*(.5d0+offset(1,i)))
     &                    *xnors(j,k)
                  enddo
               enddo
!
               do k=1,2*nedge
                  thickness(k)=thickness(k)+thicke(l,indexe+k)
               enddo
            enddo
         endif
      else
!
!           generating the 3-D element topology for axisymmetric elements
!
         dc=dcos(thicks(1)/2.d0)
         ds=dsin(thicks(1)/2.d0)
         do j=1,nedge
            indexk=iponor(2,indexe+j)
            x=co(1,nodel(j))
            y=co(2,nodel(j))
!
            node=knor(indexk+1)
            co(1,node)=x*dc
            co(2,node)=y
            co(3,node)=-x*ds
            kon(indexe+j)=node
!
            if(quadratic) then
               node=knor(indexk+2)
               co(1,node)=x
               co(2,node)=y
               co(3,node)=0.d0
               kon(indexe+4*nedge+j)=node
            else
!
!              for linear elements the coordinates of the nodes
!              in the middle may be needed in case of knots (to
!              determine the coefficients of the equations)
!
               node=knor(indexk+2)
               co(1,node)=x
               co(2,node)=y
               co(3,node)=0.d0
            endif
!
            node=knor(indexk+3)
            co(1,node)=x*dc
            co(2,node)=y
            co(3,node)=x*ds
            kon(indexe+nedge+j)=node
         enddo
!
         if(quadratic) then
            do j=nedge+1,2*nedge
               indexk=iponor(2,indexe+j)
               x=co(1,nodel(j))
               y=co(2,nodel(j))
!     
               node=knor(indexk+1)
               co(1,node)=x*dc
               co(2,node)=y
               co(3,node)=-x*ds
               kon(indexe+nedge+j)=node
!     
               node=knor(indexk+3)
               co(1,node)=x*dc
               co(2,node)=y
               co(3,node)=x*ds
               kon(indexe+2*nedge+j)=node
            enddo
         endif
      endif
!
!           additional SPC's due to plane strain/plane stress/axisymmetric
!           conditions
!
      do j=1,nope
         if((lakon(i)(1:1).ne.'S').and.(lakon(i)(1:1).ne.'M')) then
!
!                    fixing the middle plane
!
c            if(rig(nodel(j)).gt.0) cycle
            if(rig(nodel(j)).gt.0) then
               nodefixz=nodel(j)
            else
               nodefixz=nodes(2,j)
            endif
!
            if(ithermal(2).ne.2) then
               val=0.d0
               k=3
               if(nam.gt.0) iamplitude=0
               type='M'
c               call bounadd(nodes(2,j),k,k,val,nodeboun,
               call bounadd(nodefixz,k,k,val,nodeboun,
     &              ndirboun,xboun,nboun,nboun_,iamboun,
     &              iamplitude,nam,ipompc,nodempc,coefmpc,
     &              nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &              ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,
     &              labmpc,type,typeboun,nmethod,iperturb,
     &              fixed,vold,idummy,mi,label)
            endif
!     
!           specifying that the side planes do the same
!           as the middle plane (in all directions for
!           plane strain and axisymmetric elements, in the
!           plane for plane stress elements)
!
            if(ithermal(2).le.1) then
               jmin=1
               jmax=3
            elseif(ithermal(2).eq.2) then
               jmin=0
               jmax=0
            else
               jmin=0
               jmax=3
            endif
!
            do l=1,3,2
               newnode=nodes(l,j)
               do idir=jmin,jmax
                  if((idir.eq.3).and.(lakon(i)(1:3).eq.'CPS')) then
!
!                    for linear elements and plane stress it is not
!                    sufficient to set the z-displacement of the 
!                    midplane nodes to zero, since these are not used.
!                    Instead, the sum of the z-displacement of the
!                    boundary plane nodes are set to zero:
!                    w_1+w_3=0
!
                     if(((lakon(i)(4:4).eq.'3').or.
     &                   (lakon(i)(4:4).eq.'4')).and.
     &                  (l.eq.1)) then
                        idof=8*(newnode-1)+3
                        call nident(ikmpc,idof,nmpc,id)
                        if((id.le.0).or.(ikmpc(id).ne.idof)) then
                           nmpc=nmpc+1
                           if(nmpc.gt.nmpc_) then
                              write(*,*) 
     &                          '*ERROR in gen3dfrom2d: increase nmpc_'
                              call exit(201)
                           endif
                           labmpc(nmpc)='                    '
                           ipompc(nmpc)=mpcfree
                           do m=nmpc,id+2,-1
                              ikmpc(m)=ikmpc(m-1)
                              ilmpc(m)=ilmpc(m-1)
                           enddo
                           ikmpc(id+1)=idof
                           ilmpc(id+1)=nmpc
                           nodempc(1,mpcfree)=newnode
                           nodempc(2,mpcfree)=3
                           coefmpc(mpcfree)=1.d0
                           mpcfree=nodempc(3,mpcfree)
                           if(mpcfree.eq.0) then
                              write(*,*) 
     &                        '*ERROR in gen3dfrom2d: increase memmpc_'
                              call exit(201)
                           endif
                           nodempc(1,mpcfree)=nodes(3,j)
                           nodempc(2,mpcfree)=3
                           coefmpc(mpcfree)=1.d0
                           mpcfreenew=nodempc(3,mpcfree)
                           if(mpcfreenew.eq.0) then
                              write(*,*) 
     &                        '*ERROR in gen3dfrom2d: increase memmpc_'
                              call exit(201)
                           endif
                           nodempc(3,mpcfree)=0
                           mpcfree=mpcfreenew
                        endif
                     else
                        cycle
                     endif
                  endif
                  idof=8*(newnode-1)+idir
                  call nident(ikmpc,idof,nmpc,id)
                  if((id.le.0).or.(ikmpc(id).ne.idof)) then
                     nmpc=nmpc+1
                     if(nmpc.gt.nmpc_) then
                        write(*,*) 
     &                       '*ERROR in gen3dfrom2d: increase nmpc_'
                        call exit(201)
                     endif
                     labmpc(nmpc)='                    '
                     ipompc(nmpc)=mpcfree
                     do m=nmpc,id+2,-1
                        ikmpc(m)=ikmpc(m-1)
                        ilmpc(m)=ilmpc(m-1)
                     enddo
                     ikmpc(id+1)=idof
                     ilmpc(id+1)=nmpc
                     nodempc(1,mpcfree)=newnode
                     nodempc(2,mpcfree)=idir
                     coefmpc(mpcfree)=1.d0
                     mpcfree=nodempc(3,mpcfree)
                     if(mpcfree.eq.0) then
                        write(*,*) 
     &                       '*ERROR in gen3dfrom2d: increase memmpc_'
                        call exit(201)
                     endif
                     nodempc(1,mpcfree)=nodes(2,j)
                     if((lakon(i)(2:2).eq.'A').and.(idir.eq.3))
     &                    then
                        nodempc(2,mpcfree)=1
                     else
                        nodempc(2,mpcfree)=idir
                     endif
                     if(lakon(i)(2:2).eq.'A') then
                        if(idir.eq.1) then
                           coefmpc(mpcfree)=-dc
                        elseif(idir.eq.3) then
                           if(l.eq.1) then
                              coefmpc(mpcfree)=ds
                           else
                              coefmpc(mpcfree)=-ds
                           endif
                        else
                           coefmpc(mpcfree)=-1.d0
                        endif
                     else
                        coefmpc(mpcfree)=-1.d0
                     endif
                     mpcfreenew=nodempc(3,mpcfree)
                     if(mpcfreenew.eq.0) then
                        write(*,*) 
     &                       '*ERROR in gen3dfrom2d: increase memmpc_'
                        call exit(201)
                     endif
                     nodempc(3,mpcfree)=0
                     mpcfree=mpcfreenew
                  endif
               enddo
            enddo
         endif
      enddo
!
!     create a local coordinate system for each layer in a shell
!     element
!
      if(lakon(i)(1:1).eq.'S') then
!
!        initialize neworien
!
         do j=0,norieninput
            neworien(j)=0
         enddo
!
!        loop over the layers
!
         do j=1,mi(3)
!
!           check whether layer exists (i.e. whether a material
!           was assigned to the layer)
!            
            if(ielmat(j,i).le.0) exit
!
            iorien=max(0,ielorien(j,i))
!
!           check whether this orientation already occurred within
!           this element (i.e. for a previous layer)
!
            if(neworien(iorien).gt.0) then
               ielorien(j,i)=neworien(iorien)
               cycle
            endif
!
!           calculate the coordinates of the point at
!           (xi,et)=(0,0) and the normal at that point
!
            xi=0.d0
            et=0.d0
!
            if(lakon(i)(2:2).eq.'3') then
               nope=3
            elseif(lakon(i)(2:2).eq.'4') then
               nope=4
            elseif(lakon(i)(2:2).eq.'6') then
               nope=6
            elseif(lakon(i)(2:2).eq.'8') then
               nope=8
            else
               exit
            endif
!
            do k=1,nope
               node=kon(indexe+ishift+k)
               do l=1,3
                  xl(l,k)=co(l,node)
               enddo
            enddo
!
            if(nope.eq.3) then
               call shape3tri(xi,et,xl,xno,xs,shp,iflag)
            elseif(nope.eq.4) then
               call shape4q(xi,et,xl,xno,xs,shp,iflag)
            elseif(nope.eq.6) then
               call shape6tri(xi,et,xl,xno,xs,shp,iflag)
            else
               call shape8q(xi,et,xl,xno,xs,shp,iflag)
            endif
!
            dd=dsqrt(xno(1)*xno(1)+xno(2)*xno(2)+xno(3)*xno(3))
            do l=1,3
               xno(l)=xno(l)/dd
            enddo
c            write(*,*) 'gen3dfrom2d xnor',(xnor(l),l=1,3)
!
!           coordinates of the point at (xi,et)=(0.,0.)
!
            do l=1,3
               p(l)=0.d0
               do k=1,nope
                  p(l)=p(l)+shp(4,k)*xl(l,k)
               enddo
            enddo
!
            if(iorien.eq.0) then
!
!              unit matrix
!
               do k=1,3
                  do l=1,3
                     a(k,l)=0.d0
                  enddo
                  a(k,k)=1.d0
               enddo
            else
               call transformatrix(orab(1,iorien),p,a)
            endif
!
!           dd=e1.n
!            
            dd=a(1,1)*xno(1)+a(2,1)*xno(2)+a(3,1)*xno(3)
c            write(*,*) 'gen3dfrom2d ',i,dd
!
!           check whether e1 is within 0.1° from n
!           e1 is the local 1-direction in the shell plane
!
            if(dabs(dd).gt.0.999999999536d0) then
!
!              project the z-axis
!
c               write(*,*) 'gen3dfrom2d project z'
               dd=a(1,3)*xno(1)+a(2,3)*xno(2)+a(3,3)*xno(3)
               do l=1,3
                  e1(l)=a(l,3)-dd*xno(l)
               enddo
            else
!
!              project the x-axis
!
c               write(*,*) 'gen3dfrom2d project x'
               do l=1,3
                  e1(l)=a(l,1)-dd*xno(l)
               enddo
            endif
!
            dd=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
            do l=1,3
               e1(l)=e1(l)/dd
            enddo
!
!           e2=nxe1
!
            e2(1)=xno(2)*e1(3)-e1(2)*xno(3)
            e2(2)=xno(3)*e1(1)-e1(3)*xno(1)
            e2(3)=xno(1)*e1(2)-e1(1)*xno(2)
!
            norien=norien+1
            if(norien.gt.norien_) then
               write(*,*) '*ERROR in gen3dfrom2d: increase norien_'
               call exit(201)
            endif
!
!           creating the name of the new orientation =
!           old name + _shell_ + element number (10 digits;
!           this corresponds to the maximum number of digits
!           for element numbers in the frd-format)
!
            if(iorien.gt.0) then
               ipos=index(orname(iorien),' ')
               orname(norien)(1:ipos-1)=orname(iorien)(1:ipos-1)
            else
               ipos=1
            endif
            orname(norien)(ipos:ipos+6)='_shell_'
            write(orname(norien)(ipos+7:ipos+16),'(i10.10)') i
            do l=ipos+17,80
               orname(norien)(l:l)=' '
            enddo
!
            do l=1,3
               orab(l,norien)=e1(l)
               orab(l+3,norien)=e2(l)
            enddo
            orab(7,norien)=1.d0
c            write(*,*) 'gen3dfrom2d o',(orab(l,iorien),l=1,7)
c            write(*,*) 'gen3dfrom2d n',(orab(l,norien),l=1,7)
!
            neworien(iorien)=norien
            ielorien(j,i)=norien
!
         enddo
!
      endif
!            
      return
      end


