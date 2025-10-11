!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2024 Guido Dhondt
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
      subroutine multistages(set,istartset,iendset,
     &     ialset,nset,tieset,tietol,co,nk,ipompc,nodempc,
     &     coefmpc,nmpc,ikmpc,ilmpc,mpcfree,xind,yind,ics,nx,ny,
     &     xind0,yind0,cs,labmpc,ntie,mcs,rcscg,rcs0cg,zcscg,
     &     zcs0cg,nrcg,nzcg,jcs,kontri,straight,ne,ipkon,kon,
     &     lakon,lcs,ifacetet,inodface,jobnamec,nmethod)
!     
      implicit none
!     
!     creates multistage MPC's: connection of dissimilar cyclic
!     symmetric segments 
!     
!     author: Conrad Mottl
!     
      logical cylindrical,left,right,multistage
!     
      character*8 lakon(*)
      character*20 labmpc(*)
      character*81 set(*),leftset,rightset,tieset(3,*),temp,indepties,
     &     indeptiet
      character*132 jobnamec(*)
      character*256 fn
!     
      integer istartset(*),iendset(*),ialset(*),ipompc(*),nodempc(3,*),
     &     nset,i,j,k,nk,nmpc,mpcfree,ics(*),l,ikmpc(*),ilmpc(*),n,
     &     lcs(*),kflag,ncsnodes,mcs,ntie,nrcg(*),nzcg(*),jcs(*),
     &     kontri(3,*),ne,ipkon(*),kon(*),ifacetet(*),inodface(*),
     &     nodele(4),noderi(4),indexe,nope,ipos,nelem,ilen,nprint,
     &     indcs,node_cycle,itemp(4),nx(*),ny(*),netri,noderi0,
     &     nodef(8),nterms,kseg,k2,ndir,idof,number,id,mpcfreeold,
     &     lathyp(3,6),inum,ier,icount,imcs,nmethod,ksegmin,ksegmax
!     
      real*8 tolloc,co(3,* ),coefmpc(*),xind(*),yind(*),xind0(*),
     &     yind0(*),dd,xap,yap,zap,tietol(4,*),cs(18,*),xdep,ydep,
     &     phi,rcscg(*),rcs0cg(*),zcscg(*),zcs0cg(*),zp,rp,scale,
     &     straight(9,*),T(3,3),csab(7),ratio(8),Tinv(3,3),
     &     coord(3),pnod(3),T2D(3,3),phi0,al(3,3),ar(3,3),
     &     rind,dxmax,dxmin,drmax,drmin,pi,phi_min,xsegle(2),
     &     xsegri(2),xtemp(2),coglob(3),coloc(3),corotloc(3),
     &     corotglob(3),radius,angle,phi_max
!
      integer,dimension(:),allocatable::ksegmcs
!     
      data lathyp /1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1/
!     
      pi=4.d0*atan(1.d0)
      multistage=.true.
!     
!     opening a file to store the nodes which are not connected
!     
      ilen=index(jobnamec(1),char(0))-1
      fn=jobnamec(1)(1:ilen)//'_WarnNodeMissMultiStage.nam'
      open(40,file=fn,status='unknown')
      write(40,*) '*NSET,NSET=WarnNodeMultiStage'
      icount=0
!     
!     Find the TIE numbers which describe multistage connections
!     
      do i=1,ntie
        if (tieset(1,i)(81:81).eq.'M') then
          tieset(1,i)(81:81)=' '
!     
!     Creating of the fields nodele and noderi with specific 
!     information for each field
!     
!     l = left     
!     r = right
!     
!     node(le,ri)(1)=number of a node belonging to the left/right side
!     node(le,ri)(2)=setnumber
!     node(le,ri)(3)=number of an element belonging to the left/right side
!     node(le,ri)(4)=cyclic symmetry parts number of the left/right side
!
!     xseg(le,ri)(1)=value of N on the left/right side
!     xseg(le,ri)(2)=value of NPHYS on the left/right side
!          
!     Defining the left and right node sets      
!     
          leftset=tieset(2,i)
          rightset=tieset(3,i)
!     
          if (tietol(1,i).eq.-1.d0) then
            tolloc=0.1d0
            write(*,*) 
     &           '*INFO in multistages: no tolerance was defined'
            write(*,*) '      in the *TIE option; a tolerance of ',
     &           tolloc
            write(*,*) '      will be used'
            write(*,*)
          else
            tolloc=tietol(1,i)
          endif
!     
          call cident81(set,leftset,nset,j)
          if(j.gt.0) then
            if(leftset.eq.set(j)) then
              nodele(1)=ialset(istartset(j))
              nodele(2)=j         
            endif
          endif
          call cident81(set,rightset,nset,j)
          if(j.gt.0) then
            if(rightset.eq.set(j)) then
              noderi(1)=ialset(istartset(j))
              noderi(2)=j         
            endif
          endif
!     
!     identifying an element on either side of the
!     multistage connection
!     
          left=.false.
          right=.false.
          loop: do k=1,ne
            indexe=ipkon(k)
            if(indexe.lt.0) cycle
!     
!     number of nodes belonging to the element
!     
            if(lakon(k)(1:5).eq.'C3D8I') then
              nope=11
            elseif(lakon(k)(4:4).eq.'2') then
              nope=20
            elseif(lakon(k)(4:4).eq.'8') then
              nope=8
            elseif(lakon(k)(4:5).eq.'10') then
              nope=10
            elseif(lakon(k)(4:4).eq.'4') then
              nope=4
            elseif(lakon(k)(4:5).eq.'13') then
              nope=13
            elseif(lakon(k)(4:4).eq.'5') then
              nope=5
            elseif(lakon(k)(4:5).eq.'15') then
              nope=15
            elseif(lakon(k)(4:4).eq.'6') then
              nope=6
            elseif(lakon(k)(1:2).eq.'ES') then
              read(lakon(k)(8:8),'(i1)') nope
              nope=nope+1
            endif
!     
            do l=indexe+1,indexe+nope
              if(.not.left) then
                if(nodele(1).eq.kon(l))then
                  nodele(3)=k
                  left=.true.
                endif
              endif
              if(.not.right) then
                if(noderi(1).eq.kon(l))then
                  noderi(3)=k
                  right=.true.
                endif
              endif
              if(left.and.right) exit loop
            enddo
          enddo loop
!     
!         determining the number of sectors and the cyclic symmetry     
!         number for the left/right side
!     
!         assumption: number of parts is a natural number     
!     
          do j=1,mcs
            do l=istartset(int(cs(13,j))),iendset(int(cs(13,j)))
              if (ialset(l).eq.nodele(3)) then
                xsegle(1)=cs(1,j)
                xsegle(2)=cs(18,j)
                nodele(4)=j
              elseif (ialset(l).eq.noderi(3)) then
                xsegri(1)=cs(1,j)
                xsegri(2)=cs(18,j)
                noderi(4)=j
              endif   
            enddo
            csab(1)=cs(6,j)
            csab(2)=cs(7,j)
            csab(3)=cs(8,j)
            csab(4)=cs(9,j)
            csab(5)=cs(10,j)
            csab(6)=cs(11,j)
            csab(7)=-1.d0
          enddo
!     
!     Sorting such that rightset has the smaller angle;
!     it is taken is independent side for the multistage MPC's     
!     
          if (xsegle(2).ge.xsegri(2)) then
            indcs=nodele(4)
            phi0=(2.d0*pi)/xsegle(2)
            temp=rightset;
            rightset=leftset
            leftset=temp
            do j=1,4
              itemp(j)=noderi(j)
              noderi(j)=nodele(j)
              nodele(j)=itemp(j)
            enddo
            do j=1,2
              xtemp(j)=xsegri(j)
              xsegri(j)=xsegle(j)
              xsegle(j)=xtemp(j)
            enddo
          else
            indcs=noderi(4)
            phi0=(2.d0*pi)/xsegri(2)
          endif
!     
!         scaling factor for the connecting MPC's in static calculations    
!     
          scale=xsegle(2)/xsegri(2)
!     
!     Looking for a node on the independent cyclic symmetry side of
!     the right side of the multistage tie    
!     
          indepties=tieset(3,int(cs(17,indcs)))
          indeptiet=indepties
          ipos=index(indepties,' ')
          indepties(ipos:ipos)='S'
          indeptiet(ipos:ipos)='T'
          do j=1,nset
            if(indepties.eq.set(j)) then
!     
!     nodal independent surface
!     
              node_cycle=ialset(istartset(j))
              exit
            elseif(indeptiet.eq.set(j)) then
!     
!     facial independent surface
!     
              nelem=int(ialset(istartset(j))/10)
              node_cycle=kon(ipkon(nelem)+1)
              exit
            endif
          enddo
!     
!     Defining a transformation matrix from the global system into a local
!     system with the x-direction along the rotation axis 
!     
          T(1,1)=csab(4)-csab(1)
          T(1,2)=csab(5)-csab(2)
          T(1,3)=csab(6)-csab(3)
          dd=dsqrt(T(1,1)*T(1,1)+T(1,2)*T(1,2)+T(1,3)*T(1,3))
          T(1,1)=T(1,1)/dd
          T(1,2)=T(1,2)/dd
          T(1,3)=T(1,3)/dd
!     
!     second local direction is orthogonal to the rotation axis
!     
          xap=co(1,node_cycle)-csab(1)
          yap=co(2,node_cycle)-csab(2)
          zap=co(3,node_cycle)-csab(3)
!     
          zp=xap*T(1,1)+yap*T(1,2)+zap*T(1,3)
          rp=((xap-T(1,1)*zp)**2+(yap-zp*T(1,2))**2+
     &         (zap-zp*T(1,3))**2)
          rp=dsqrt(rp)
!     
!     Performing the vector product for the third local direction
!     
          if(rp.gt.1.d-10) then
            T(2,1)=(xap-zp*T(1,1))/rp
            T(2,2)=(yap-zp*T(1,2))/rp
            T(2,3)=(zap-zp*T(1,3))/rp
            T(3,1)=T(1,2)*T(2,3)-T(2,2)*T(1,3)
            T(3,2)=T(2,1)*T(1,3)-T(1,1)*T(2,3)
            T(3,3)=T(1,1)*T(2,2)-T(2,1)*T(1,2)
          endif
!     
!     (global coordinates)=T^T.(local coordinates)    
!     (local coordinates)=T.(global coordinates)
!     
!     Inverting the rotary matrix transform from local to global
!     
          call invert3D(T,Tinv,3)
!     
!         mapping the independent nodes into the local system;
!         based on the number of sectors a circumferential base angle     
!         is determined. Since the master surface can be slanted
!         (changing phi-value for increasing radius for an axial
!          multistage tie or changing phi-value for increasing 
!          axial coordinate values for a radial multistage tie)
!         it may cover a multiple (n) of the base angle. For n=0
!         the surface does not exceed tha base angle, otherwise it
!         exceeds the base angle n times
!
          l=0
          do j=istartset(noderi(2)),iendset(noderi(2))
            l=l+1
            pnod(1)=co(1,ialset(j))-csab(1)
            pnod(2)=co(2,ialset(j))-csab(2)
            pnod(3)=co(3,ialset(j))-csab(3)
            call Mprod(T,pnod,coord,3)
            xind(l)=coord(2)
            yind(l)=coord(3) 
            nx(l)=l
            ny(l)=l
            ics(l)=ialset(j)
            xind0(l)=xind(l)
            yind0(l)=yind(l)
            rind=dsqrt(coord(2)**2+coord(3)**2)
            phi=datan2(-coord(3),coord(2))
            if (l.gt.1) then
              dxmax=max(dxmax,dabs(coord(1)))
              drmax=max(drmax,dabs(rind))
              
              dxmin=min(dxmin,dabs(coord(1)))
              drmin=min(drmin,dabs(rind))
              phi_min=min(phi_min,phi)
              phi_max=max(phi_max,phi)
            else
              dxmax=dabs(coord(1))
              phi_min=phi
              phi_max=phi
              drmax=rind
              dxmin=dabs(coord(1))
              drmin=rind
            endif
          enddo
!
!         allowing for 0.1 % sector angle error
!
          n=floor(xsegri(2)*(phi_max-phi_min)/(2.d0*pi)-0.001d0)
!     
          cylindrical=.false.
          if ((dxmax-dxmin).ge.(drmax-drmin)) then
            l=0
            n=0
            do j=istartset(noderi(2)),iendset(noderi(2))
              l=l+1
              pnod(1)=co(1,ialset(j))-csab(1)
              pnod(2)=co(2,ialset(j))-csab(2)
              pnod(3)=co(3,ialset(j))-csab(3)
              call Mprod(T,pnod,coord,3)
              xind(l)=coord(1)
              if(l.eq.1) radius=dsqrt(coord(2)**2+coord(3)**2)
              angle=datan2(-coord(3),coord(2))
              yind(l)=radius*angle
              n=max(n,floor(xsegri(2)*(angle-phi_min)/
     &             (2.d0*pi)-0.001d0))
              nx(l)=l
              ny(l)=l
              ics(l)=ialset(j)
              xind0(l)=xind(l)
              yind0(l)=yind(l)
            enddo
            cylindrical=.true.
          endif
!     
!     Sorting the coordinates for the further consideration
!     
          ncsnodes=l
          kflag=2
          call dsort(xind,nx,ncsnodes,kflag)
          call dsort(yind,ny,ncsnodes,kflag)
!     
          call triangulate(ics,xind0,yind0,ncsnodes,
     &         rcscg,rcs0cg,zcscg,zcs0cg,nrcg,nzcg,jcs,kontri,
     &         straight,ne,ipkon,kon,lakon,lcs,netri,ifacetet,
     &         inodface)
!     
!     field containing the cyclic symmetry part number corresponding
!     to each value of kseg
!     
          allocate(ksegmcs(nint(xsegri(2))+1))
          do j=1,nint(xsegri(2))+1
            ksegmcs(j)=0
          enddo
!     
!         loop over all dependent nodes   
!     
          do j=istartset(nodele(2)),iendset(nodele(2))
            coglob(1)=co(1,ialset(j))-csab(1)
            coglob(2)=co(2,ialset(j))-csab(2)
            coglob(3)=co(3,ialset(j))-csab(3)
!     
!     change to local coordinates "coord"
!     
            call Mprod(T,coglob,coloc,3)
!     
!     Determining the phase shift for the tie constraints
!     
            phi=datan2(-coloc(3),coloc(2))
!     
!     kseg is the number of times the base sector has to be shifted
!     in order to contain the node at stake    
!     
            ksegmax=floor(xsegri(2)*(phi-phi_min)/(2.d0*pi))
            ksegmin=ksegmax-n
!
            do kseg=ksegmin,ksegmax
              T2D(1,1)=1.d0
              T2D(1,2)=0.d0
              T2D(1,3)=0.d0
              T2D(2,1)=0.d0
              T2D(2,2)=dcos(-kseg*phi0)
              T2D(2,3)=dsin(-kseg*phi0)
              T2D(3,1)=0.d0
              T2D(3,2)=-dsin(-kseg*phi0)
              T2D(3,3)=dcos(-kseg*phi0)
!     
!     Rotating the dependent nodes by the number of 
!     segments to "pnod" in local coordinates
!     
              call Mprod(T2D,coloc,corotloc,3)
!     
!     copying the local coordinates to the local variables
!     
              if (cylindrical) then
!     
!               radial connection    
!     
                xdep=corotloc(1)
                ydep=radius*datan2(-corotloc(3),corotloc(2))
              else
!     
!               axial connection     
!     
                xdep=corotloc(2)
                ydep=corotloc(3)
              endif
              noderi0=nk+1
!     
!     transforming back to global coordinates
!     
              call Mprod(Tinv,corotloc,corotglob,3)
              co(1,noderi0)=corotglob(1)
              co(2,noderi0)=corotglob(2)
              co(3,noderi0)=corotglob(3)
!     
!     determining the coefficients of the multistage MPC's
!     
              ier=0
              nprint=-1
              call linkdissimilar(co,csab,
     &             rcscg,rcs0cg,zcscg,zcs0cg,nrcg,nzcg,
     &             straight,nodef,ratio,nterms,xdep,ydep,netri,
     &             noderi0,ifacetet,inodface,ialset(j),
     &             T(1,1),T(1,2),T(1,3),ier,multistage,nprint)
!     
!     if a corresponding face was found: go to next slave node
!     
              if(ier.eq.0) exit
            enddo
!
            if(ier.lt.0) then
              write(*,*)
     &             '*WARNING in multistages: no suitable partner'
              write(*,*) '         face found for node',ialset(j),'.'
              write(*,*) 
     &             '         Nodes belonging to the best partner face:'
              write(*,*) (nodef(k2),k2=1,nterms)
              write(*,*) 
              write(40,*) ialset(j)
              icount=icount+1
            endif
            if(kseg.lt.0) kseg=kseg+nint(xsegri(2))
!     
!     scaling the independent degrees of freedom of the
!     multistage equations for non-frequency calculations
!     
            if(nmethod.ne.2) then
              do k2=1,nterms
                ratio(k2)=ratio(k2)*scale
              enddo
            endif
!     
!     all sector shifts for multistage MPC's are stored
!     
            if(kseg.eq.1) then
              imcs=noderi(4)
            elseif(kseg.gt.1) then
              if(ksegmcs(kseg).gt.0) then
                imcs=ksegmcs(kseg)
              else
                mcs=mcs+1
                ksegmcs(kseg)=mcs
                imcs=mcs
                cs(1,mcs)=xsegri(1)/kseg
                do k=2,17
                  cs(k,mcs)=cs(k,noderi(4))
                enddo
                cs(13,mcs)=0.5d0
                cs(18,mcs)=xsegri(2)/kseg
              endif
            endif
!     
!     the tie MPC relationships apply in cylindrical coordinates
!     
            call transformatrix(csab,co(1,ialset(j)),al)
            call transformatrix(csab,co(1,noderi0),ar)   
!     
!     checking for latin hypercube positions in matrix al none of
!     which are zero
!     
            do inum=1,6
              if((dabs(al(lathyp(1,inum),1)).gt.1.d-3).and.
     &             (dabs(al(lathyp(2,inum),2)).gt.1.d-3).and.
     &             (dabs(al(lathyp(3,inum),3)).gt.1.d-3)) exit
            enddo
!     
            do ndir=1,3
              nmpc=nmpc+1
              ipompc(nmpc)=mpcfree 
              if(kseg.gt.0) then
                labmpc(nmpc)='CYCLIC              '
                if (imcs.lt.10) then
                  write(labmpc(nmpc)(7:7),'(i1)') imcs
                elseif(imcs.lt.100) then
                  write(labmpc(nmpc)(7:8),'(i2)') imcs
                elseif(imcs.lt.1000) then
                  write(labmpc(nmpc)(7:9),'(i3)') imcs
                else
                  write(*,*)
     &                 '*ERROR in multistages: no more than 999'
                  write(*,*)
     &                 '       cyclic symmetry definitions allowed'
                  call exit(201)
                endif
!     
                number=lathyp(ndir,inum)
!     
!     determining which direction to use for the
!     dependent side: should not occur on the dependent
!     side in another MPC and should have a nonzero
!     coefficient
!     
                idof=8*(ialset(j)-1)+number
                call nident(ikmpc,idof,nmpc-1,id)
                if(id.gt.0) then
                  if(ikmpc(id).eq.idof) then
                    write(*,*) 
     &                   '*WARNING in multistages: cyclic MPC in node'
                    write(*,*) '         ',ialset(j),
     &                   ' and direction',ndir
                    write(*,*) '         cannot be created: the'
                    write(*,*) 
     &                   '         DOF in this node is already used'
                    nmpc=nmpc-1
                    cycle
                  endif
                endif
!     
                number=number-1
!     
!     updating ikmpc and ilmpc
!     
                do k=nmpc,id+2,-1
                  ikmpc(k)=ikmpc(k-1)
                  ilmpc(k)=ilmpc(k-1)
                enddo
                ikmpc(id+1)=idof
                ilmpc(id+1)=nmpc
                do k=1,3
                  number=number+1
                  if(number.gt.3) number=1
                  if(dabs(al(number,ndir)).lt.1.d-5) cycle
                  nodempc(1,mpcfree)=ialset(j)
                  nodempc(2,mpcfree)=number
                  coefmpc(mpcfree)=al(number,ndir)
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                    write(*,*)
     &                   '*ERROR in multistages: increase memmpc_'
                    call exit(201)
                  endif
                enddo
                do k=1,3
                  number=number+1
                  if(number.gt.3) number=1
                  if(dabs(ar(number,ndir)).lt.1.d-5) cycle
                  do k2=1,nterms
                    if (dabs(ratio(k2)).gt.1.d-10) then
                      nodempc(1,mpcfree)=nodef(k2)
                      nodempc(2,mpcfree)=number
                      coefmpc(mpcfree)=
     %                     -ar(number,ndir)*ratio(k2)
                      mpcfreeold=mpcfree
                      mpcfree=nodempc(3,mpcfree)
                    endif
                    if(mpcfree.eq.0) then
                      write(*,*) 
     &                     '*ERROR in multistages: increase memmpc_'
                      call exit(201)
                    endif
                  enddo
                enddo
                nodempc(3,mpcfreeold)=0
                phi=0.d0
              else
!     
!     Starting the creation of the standard MPC
!     
!     determining which direction to use for the
!     dependent side: should not occur on the dependent
!     side in another MPC and should have a nonzero
!     coefficient
!     
                labmpc(nmpc)='                    '
                idof=8*(ialset(j)-1)+ndir
                call nident(ikmpc,idof,nmpc-1,id)
                if(id.gt.0) then
                  if(ikmpc(id).eq.idof) then
                    write(*,*) 
     &                   '*WARNING in multistages: cyclic MPC in node'
                    write(*,*) '         ',ialset(j),
     &                   ' and direction',ndir
                    write(*,*) '         cannot be created: the'
                    write(*,*) 
     &                   '         DOF in this node is already used'
                    nmpc=nmpc-1
                    cycle
                  endif
                endif
!     
!     updating ikmpc and ilmpc
!     
                do k=nmpc,id+2,-1
                  ikmpc(k)=ikmpc(k-1)
                  ilmpc(k)=ilmpc(k-1)
                enddo
                ikmpc(id+1)=idof
                ilmpc(id+1)=nmpc
                
                nodempc(1,mpcfree)=ialset(j)
                nodempc(2,mpcfree)=ndir
                coefmpc(mpcfree)=-1          
                mpcfree=nodempc(3,mpcfree)
                if(mpcfree.eq.0) then
                  write(*,*)
     &                 '*ERROR in multistages: increase memmpc_'
                  call exit(201)
                endif
                
                do k2=1,nterms
                  if (dabs(ratio(k2)).gt.1.d-6) then
                    nodempc(1,mpcfree)=nodef(k2)
                    nodempc(2,mpcfree)=ndir
                    coefmpc(mpcfree)=ratio(k2)
                    mpcfreeold=mpcfree
                    mpcfree=nodempc(3,mpcfree)
                  endif
                  if(mpcfree.eq.0) then
                    write(*,*) 
     &                   '*ERROR in multistages: increase memmpc_'
                    call exit(201)
                  endif
                enddo  
              endif
              nodempc(3,mpcfreeold)=0
            enddo
!     
          enddo                 !Loop over nodes on dependent side
!     
          deallocate(ksegmcs)
!     
        endif
      enddo                     !Loop over ties
!
      if(icount.gt.0) then
        write(*,*) '*INFO multistages:'
        write(*,*) '      failed nodes are stored in file'
        write(*,*) '      ',fn(1:ilen+27)
        write(*,*) '      This file can be loaded into'
        write(*,*) '      an active cgx-session by typing'
        write(*,*) 
     &       '      read ',fn(1:ilen+27),' inp'
        write(*,*)
        close(40)
      else
        close(40,status='delete')
      endif
!     
      return
      end
!     
      subroutine Mprod(M,v_in,v_out,size)
!     
      implicit none
!     
      integer size,i,j
      real*8 M(size,size),v_in(size),v_out(size),line
!     
      do i=1,size
        line=0
        do j=1,size
          line=M(i,j)*v_in(j)+line
        enddo
        v_out(i)=line
      enddo
!     
      return
      end
!     
      subroutine invert3D(M,Minv,size)
!     
      implicit none
!     
      integer i,j,size
      real*8 detA,M(3,3),Minv(3,3)
!     
      detA=M(1,1)*M(2,2)*M(3,3)+M(2,1)*M(3,2)*M(1,3)+
     &     M(3,1)*M(1,2)*M(2,3)-
     &     M(3,1)*M(2,2)*M(1,3)-M(1,1)*M(3,2)*M(2,3)-
     &     M(2,1)*M(1,2)*M(3,3)
      
      Minv(1,1)=M(2,2)*M(3,3)-M(3,2)*M(2,3)
      Minv(2,1)=M(3,1)*M(2,3)-M(2,1)*M(3,3)
      Minv(3,1)=M(2,1)*M(3,2)-M(3,1)*M(2,2)
      
      Minv(1,2)=M(3,2)*M(1,3)-M(1,2)*M(3,3)
      Minv(2,2)=M(1,1)*M(3,3)-M(3,1)*M(1,3)
      Minv(3,2)=M(3,1)*M(1,2)-M(1,1)*M(3,2)
      
      Minv(1,3)=M(1,2)*M(2,3)-M(2,2)*M(1,3)
      Minv(2,3)=M(2,1)*M(1,3)-M(1,1)*M(2,3)
      Minv(3,3)=M(1,1)*M(2,2)-M(2,1)*M(1,2)
      do i=1,size
        do j=1,size
          Minv(i,j)=(1/detA)*Minv(i,j)
        enddo
      enddo
!     
      return
      end
      
