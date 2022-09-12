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
      subroutine generatecycmpcs(tolloc,co,nk,ipompc,nodempc,
     &     coefmpc,nmpc,ikmpc,ilmpc,mpcfree,rcs,zcs,ics,nr,nz,
     &     rcs0,zcs0,labmpc,mcs,triangulation,csab,xn,yn,zn,phi,
     &     noded,ncsnodes,rcscg,rcs0cg,zcscg,zcs0cg,nrcg,nzcg,jcs,lcs,
     &     kontri,straight,ne,ipkon,kon,lakon,ifacetet,inodface,
     &     vold,nef,mi,indepset,ithermal,icount)
!     
!     generate cyclic mpc's
!     
      implicit none
!     
      logical triangulation,interpolation,multistage
!     
      character*1 c
      character*3 m1,m2,m3
      character*5 p0,p1,p2,p3,p7,p9999
      character*8 lakon(*)
      character*20 labmpc(*),label
      character*81 indepset
      character*132 fntria
!     
      integer ipompc(*),nodempc(3,*),nneigh,ne,ipkon(*),kon(*),
     &     j,k,nk,nmpc,mpcfree,ics(*),nterms,ncyclicsymmetrymodel,
     &     nr(*),nz(*),noded,nodei,ikmpc(*),ilmpc(*),kontri(3,*),
     &     number,idof,ndir,node,ncsnodes,id,mpcfreeold,icount,
     &     mcs,nrcg(*),nzcg(*),jcs(*),lcs(*),nodef(8),
     &     netri,ifacetet(*),inodface(*),lathyp(3,6),inum,one,i,
     &     noden(10),ier,ipos,nef,mi(*),ilen,ithermal(2)
!     
      real*8 tolloc,co(3,*),coefmpc(*),rcs(*),zcs(*),rcs0(*),zcs0(*),
     &     csab(7),xn,yn,zn,xap,yap,zap,rp,zp,al(3,3),ar(3,3),phi,
     &     x2,y2,z2,x3,y3,z3,rcscg(*),rcs0cg(*),zcscg(*),zcs0cg(*),
     &     straight(9,*),ratio(8),vold(0:mi(2),*)
!     
      save netri,ncyclicsymmetrymodel
!     
      data ncyclicsymmetrymodel /0/
!     
!     latin hypercube positions in a 3 x 3 matrix
!     
      data lathyp /1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1/
!     
      multistage=.false.
      nneigh=10
!     
      xap=co(1,noded)-csab(1)
      yap=co(2,noded)-csab(2)
      zap=co(3,noded)-csab(3)
!     
      zp=xap*xn+yap*yn+zap*zn
      rp=dsqrt((xap-zp*xn)**2+(yap-zp*yn)**2+(zap-zp*zn)**2)
!     
      call near2d(rcs0,zcs0,rcs,zcs,nr,nz,rp,zp,ncsnodes,noden,nneigh)
      node=noden(1)
      nodei=abs(ics(noden(1)))
!     
!     check whether node is on axis
!     
      if(nodei.eq.noded) then
        return
      endif
!     
      interpolation=.false.
!     
      if(rp.gt.1.d-10) then
        x2=(xap-zp*xn)/rp
        y2=(yap-zp*yn)/rp
        z2=(zap-zp*zn)/rp
        x3=yn*z2-y2*zn
        y3=x2*zn-xn*z2
        z3=xn*y2-x2*yn
      endif
!     
      if((tolloc.ge.0.d0).and.
     &     (tolloc.le.dsqrt((rp-rcs0(node))**2+(zp-zcs0(node))**2)))
     &     then
!     
!     the nodal positions on the dependent and independent
!     sides of the mpc's do no agree: interpolation is
!     necessary. 
!     
!     
        interpolation=.true.
!     
        if(.not.triangulation) then
          call triangulate(ics,rcs0,zcs0,ncsnodes,
     &         rcscg,rcs0cg,zcscg,zcs0cg,nrcg,nzcg,jcs,kontri,
     &         straight,ne,ipkon,kon,lakon,lcs,netri,ifacetet,
     &         inodface)
          triangulation=.true.
!     
c     fntria(1:28)='TriMasterCyclicSymmetryModel'
c     ncyclicsymmetrymodel=ncyclicsymmetrymodel+1
c     if(ncyclicsymmetrymodel.lt.10) then
c     write(fntria(29:29),'(i1)')ncyclicsymmetrymodel
c     fntria(30:33)='.frd'
c     ipos=34
c     elseif(ncyclicsymmetrymodel.lt.100) then
c     write(fntria(29:30),'(i2)')ncyclicsymmetrymodel
c     fntria(31:34)='.frd'
c     ipos=35
c     else
c     write(*,*) '*ERROR in generatecycmpcs: no more than'
c     write(*,*) '       99 cyclic symmetry model cards'
c     write(*,*) '       allowed'
c     call exit(201)
c     endif
c     do i=ipos,132
c     fntria(i:i)=' '
c     enddo
c     
          ilen=index(indepset,' ')
          fntria(1:3)='Tri'
          do j=4,ilen+2
            fntria(j:j)=indepset(j-3:j-3)
          enddo
          fntria(ilen+3:ilen+6)='.frd'
          do j=ilen+7,132
            fntria(j:j)=' '
          enddo
!     
c     open(70,file=fntria,status='unknown')
c     c='C'
c     m1=' -1'
c     m2=' -2'
c     m3=' -3'
c     p0='    0'
c     p1='    1'
c     p2='    2'
c     p3='    3'
c     p7='    7'
c     p9999=' 9999'
c     one=1
c     write(70,'(a5,a1)') p1,c
c     write(70,'(a5,a1,67x,i1)') p2,c,one
c     do i=1,nk
c     write(70,'(a3,i10,1p,3e12.5)') m1,i,(co(j,i),j=1,3)
c     enddo
c     write(70,'(a3)') m3
c     write(70,'(a5,a1,67x,i1)') p3,c,one
c     do i=1,netri
c     write(70,'(a3,i10,2a5)')m1,i,p7,p0
c     write(70,'(a3,3i10)') m2,(kontri(j,i),j=1,3)
c     enddo
c     write(70,'(a3)') m3
c     write(70,'(a5)') p9999
c     close(70)
!     
        endif
!     
        label='CYCLIC              '
        if(mcs.lt.10) then
          write(label(7:7),'(i1)') mcs
        elseif(mcs.lt.100) then
          write(label(7:8),'(i2)') mcs
        else
          write(*,*)'*ERROR in generatecycmpcs: no more than 99'
          write(*,*)'       cyclic symmetry definitions allowed'
          call exit(201)
        endif
!     
        nodei=nk+1
!     
!     copying the initial conditions for the new node
!     
        do i=0,mi(2)
          vold(i,nodei)=vold(i,noded)
        enddo
!     
        co(1,nodei)=csab(1)+zp*xn+rp*(x2*dcos(phi)+x3*dsin(phi))
        co(2,nodei)=csab(2)+zp*yn+rp*(y2*dcos(phi)+y3*dsin(phi))
        co(3,nodei)=csab(3)+zp*zn+rp*(z2*dcos(phi)+z3*dsin(phi))
!     
        ier=0
!     
        call linkdissimilar(co,csab,
     &       rcscg,rcs0cg,zcscg,zcs0cg,nrcg,nzcg,straight,
     &       nodef,ratio,nterms,rp,zp,netri,
     &       nodei,ifacetet,inodface,noded,xn,yn,
     &       zn,ier,multistage,icount)
!     
        if(ier.ne.0) then
          return
        endif
!     
      else
        if(ics(node).lt.0) return
      endif
!     
!     generating the mechanical MPC's; the generated MPC's are for
!     nodal diameter 0. For other nodal diameters the MPC's are
!     changed implicitly in mastructcs and mafillsmcs
!     
      call transformatrix(csab,co(1,noded),al)
      call transformatrix(csab,co(1,nodei),ar)
!     
!     checking for latin hypercube positions in matrix al none of
!     which are zero
!     
      do inum=1,6
        if((dabs(al(lathyp(1,inum),1)).gt.1.d-3).and.
     &       (dabs(al(lathyp(2,inum),2)).gt.1.d-3).and.
     &       (dabs(al(lathyp(3,inum),3)).gt.1.d-3)) exit
      enddo
!     
      do ndir=1,3
!     
!     determining which direction to use for the
!     dependent side: should not occur on the dependent
!     side in another MPC and should have a nonzero
!     coefficient
!     
        number=lathyp(ndir,inum)
        idof=8*(noded-1)+number
        call nident(ikmpc,idof,nmpc,id)
        if(id.gt.0) then
          if(ikmpc(id).eq.idof) then
            write(*,*) '*WARNING in generatecycmpcs: cyclic MPC in no
     &de'
            write(*,*) '         ',noded,' and direction ',ndir
            write(*,*) '         cannot be created: the'
            write(*,*) '         DOF in this node is already used'
            cycle
          endif
        endif
        number=number-1
!     
        nmpc=nmpc+1
        labmpc(nmpc)='CYCLIC              '
        if(mcs.lt.10) then
          write(labmpc(nmpc)(7:7),'(i1)') mcs
        elseif(mcs.lt.100) then
          write(labmpc(nmpc)(7:8),'(i2)') mcs
        else
          write(*,*)'*ERROR in generatecycmpcs: no more than 99'
          write(*,*)'       cyclic symmetry definitions allowed'
          call exit(201)
        endif
        ipompc(nmpc)=mpcfree 
!     
!     updating ikmpc and ilmpc
!     
        do j=nmpc,id+2,-1
          ikmpc(j)=ikmpc(j-1)
          ilmpc(j)=ilmpc(j-1)
        enddo
        ikmpc(id+1)=idof
        ilmpc(id+1)=nmpc
!     
        do j=1,3
          number=number+1
          if(number.gt.3) number=1
          if(dabs(al(number,ndir)).lt.1.d-5) cycle
          nodempc(1,mpcfree)=noded
          nodempc(2,mpcfree)=number
          coefmpc(mpcfree)=al(number,ndir)
          mpcfree=nodempc(3,mpcfree)
          if(mpcfree.eq.0) then
            write(*,*)'*ERROR in generatecycmpcs: increase memmpc_'
            call exit(201)
          endif
        enddo
        do j=1,3
          number=number+1
          if(number.gt.3) number=1
          if(dabs(ar(number,ndir)).lt.1.d-5) cycle
          if(.not.interpolation) then
            nodempc(1,mpcfree)=nodei
            nodempc(2,mpcfree)=number
            coefmpc(mpcfree)=-ar(number,ndir)
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
            if(mpcfree.eq.0) then
              write(*,*)
     &             '*ERROR in generatecycmpcs: increase memmpc_'
              call exit(201)
            endif
          else
            do k=1,nterms
              nodempc(1,mpcfree)=nodef(k)
              nodempc(2,mpcfree)=number
              coefmpc(mpcfree)=-ar(number,ndir)*ratio(k)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
              if(mpcfree.eq.0) then
                write(*,*) '*ERROR in generatecycmpcs: increase nmp
     &c_'
                call exit(201)
              endif
            enddo
          endif
        enddo
        nodempc(3,mpcfreeold)=0
      enddo
!     
      nmpc=nmpc+1
      labmpc(nmpc)='CYCLIC              '
      if(mcs.lt.10) then
        write(labmpc(nmpc)(7:7),'(i1)') mcs
      elseif(mcs.lt.100) then
        write(labmpc(nmpc)(7:8),'(i2)') mcs
      else
        write(*,*)'*ERROR in generatecycmpcs: no more than 99'
        write(*,*)'       cyclic symmetry definitions allowed'
        call exit(201)
      endif
      ipompc(nmpc)=mpcfree 
      idof=8*(noded-1)
      call nident(ikmpc,idof,nmpc-1,id)
      if(id.gt.0) then
        if(ikmpc(id).eq.idof) then
          write(*,*) '*ERROR in generatecycmpcs: temperature'
          write(*,*) '       in node',noded,'is already used'
          call exit(201)
        endif
      endif
!     
!     updating ikmpc and ilmpc
!     
      do j=nmpc,id+2,-1
        ikmpc(j)=ikmpc(j-1)
        ilmpc(j)=ilmpc(j-1)
      enddo
      ikmpc(id+1)=idof
      ilmpc(id+1)=nmpc
!     
      nodempc(1,mpcfree)=noded
      nodempc(2,mpcfree)=0
      coefmpc(mpcfree)=1.d0
      mpcfree=nodempc(3,mpcfree)
      if(mpcfree.eq.0) then
        write(*,*)'*ERROR in generatecycmpcs: increase memmpc_'
        call exit(201)
      endif
      if(.not.interpolation) then
        nodempc(1,mpcfree)=nodei
        nodempc(2,mpcfree)=0
        coefmpc(mpcfree)=-1.d0
        mpcfreeold=mpcfree
        mpcfree=nodempc(3,mpcfree)
        if(mpcfree.eq.0) then
          write(*,*)'*ERROR in generatecycmpcs: increase memmpc_'
          call exit(201)
        endif
      else
        do k=1,nterms
          nodempc(1,mpcfree)=nodef(k)
          nodempc(2,mpcfree)=0
          coefmpc(mpcfree)=-ratio(k)
          mpcfreeold=mpcfree
          mpcfree=nodempc(3,mpcfree)
          if(mpcfree.eq.0) then
            write(*,*)'*ERROR in generatecycmpcs: increase memmpc_'
            call exit(201)
          endif
        enddo
      endif
      nodempc(3,mpcfreeold)=0
!     
!     generating the static pressure and possibly turbulence MPC's
!     for 3D fluid calculations; 
!     the generated MPC's are for nodal diameter 0. 
!     
      if(nef.gt.0) then
        do i=4,mi(2)
          nmpc=nmpc+1
          labmpc(nmpc)='CYCLIC              '
          if(mcs.lt.10) then
            write(labmpc(nmpc)(7:7),'(i1)') mcs
          elseif(mcs.lt.100) then
            write(labmpc(nmpc)(7:8),'(i2)') mcs
          else
            write(*,*)'*ERROR in generatecycmpcs: no more than 99'
            write(*,*)'       cyclic symmetry definitions allowed'
            call exit(201)
          endif
          ipompc(nmpc)=mpcfree 
          idof=8*(noded-1)+i
          call nident(ikmpc,idof,nmpc-1,id)
          if(id.gt.0) then
            if(ikmpc(id).eq.idof) then
              if(i.eq.4) then
                write(*,*) '*ERROR in generatecycmpcs: pressure'
                write(*,*) '       in node',noded,'is already used'
                call exit(201)
              else
                write(*,*)
     &               '*ERROR in generatecycmpcs: turbulence parameter',
     &               i-4
                write(*,*) '       in node',noded,'is already used'
                call exit(201)
              endif
            endif
          endif
!     
!     updating ikmpc and ilmpc
!     
          do j=nmpc,id+2,-1
            ikmpc(j)=ikmpc(j-1)
            ilmpc(j)=ilmpc(j-1)
          enddo
          ikmpc(id+1)=idof
          ilmpc(id+1)=nmpc
!     
          nodempc(1,mpcfree)=noded
          nodempc(2,mpcfree)=i
          coefmpc(mpcfree)=1.d0
          mpcfree=nodempc(3,mpcfree)
          if(mpcfree.eq.0) then
            write(*,*)'*ERROR in generatecycmpcs: increase memmpc_'
            call exit(201)
          endif
          if(.not.interpolation) then
            nodempc(1,mpcfree)=nodei
            nodempc(2,mpcfree)=i
            coefmpc(mpcfree)=-1.d0
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
            if(mpcfree.eq.0) then
              write(*,*)'*ERROR in generatecycmpcs: increase memmpc_'
              call exit(201)
            endif
          else
            do k=1,nterms
              nodempc(1,mpcfree)=nodef(k)
              nodempc(2,mpcfree)=i
              coefmpc(mpcfree)=-ratio(k)
              mpcfreeold=mpcfree
              mpcfree=nodempc(3,mpcfree)
              if(mpcfree.eq.0) then
                write(*,*)
     &               '*ERROR in generatecycmpcs: increase memmpc_'
                call exit(201)
              endif
            enddo
          endif
          nodempc(3,mpcfreeold)=0
        enddo
      endif
!     
      return
      end

