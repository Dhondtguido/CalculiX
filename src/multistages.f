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
      subroutine multistages(nkon,set,istartset,iendset,
     &     ialset,nset,tieset,tietol,co,nk,ipompc,nodempc,
     &     coefmpc,nmpc,nmpc_,ikmpc,ilmpc,mpcfree,xind,yind,ics,nx,ny,
     &     xind0,yind0,ncs_,cs,labmpc,ntie,mcs,rcscg,rcs0cg,zcscg,
     &     zcs0cg,nrcg,nzcg,jcs,kontri,straight,ne,ipkon,kon,
     &     lakon,lcs,ifacetet,inodface,jobnamec)
!     
      implicit none
!     
!     creates multistage MPC's: connection of dissimilar cyclic
!     symmetric segments 
!     
!     author: Konrad Mottl
!     
      logical nodesonaxis,cylindrical,replace,left,right,multistage
!     
      character*8 lakon(*)
      character*20 labmpc(*)
      character*81 set(*),leftset,rightset,tieset(3,*),temp,indepties,
     &     indeptiet
      character*132 jobnamec(*)
      character*256 fn
!     
      integer istartset(*),iendset(*),ialset(*),ipompc(*),nodempc(3,*),
     &     nset,i,j,k,nk,nmpc,nmpc_,mpcfree,ics(*),l,ikmpc(*),ilmpc(*),
     &     lcs(*),kflag,ncsnodes,ncs_,mcs,ntie,nrcg(*),nzcg(*),jcs(*),
     &     kontri(3,*),ne,ipkon(*),kon(*),ifacetet(*),inodface(*),
     &     nodel(5),noder(5),nkon,indexe,nope,ipos,nelem,ilen,
     &     indcs,node_cycle,itemp(5),nx(*),ny(*),netri,noder0,
     &     nodef(8),nterms,kseg,k2,ndir,idof,number,id,mpcfreeold,
     &     lathyp(3,6),inum,ier,icount
!     
      real*8 tolloc,co(3,* ),coefmpc(*),xind(*),yind(*),xind0(*),
     &     yind0(*),dd,xap,yap,zap,tietol(4,*),cs(17,*),xp,yp,
     &     phi,rcscg(*),rcs0cg(*),zcscg(*),zcs0cg(*),zp,rp,
     &     straight(9,*),T(3,3),csab(7),ratio(8),Tinv(3,3),
     &     coord(3),node(3),T2D(3,3),phi0,al(3,3),ar(3,3),
     &     rind,dxmax,dxmin,drmax,drmin,pi,phi_min
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
!     **********Starting the main loop over all multistage ties*********
!     Creating of the fields nodel and noder with specific 
!     information for each node
!     node(l,r)(1)=nodenumber
!     node(l,r)(2)=setnumber
!     node(l,r)(3)=elementnumber
!     node(l,r)(4)=segments of cyclic symmetry part
!     node(l,r)(5)=cyclic symmetry parts number
!     
          replace=.true.
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
!     Copying the nodes of each Set of the Tie
!     
c     do j=1,nset
c     if(leftset.eq.set(j)) then
c     nodel(1)=ialset(istartset(j))
c     nodel(2)=j         
c     elseif(rightset.eq.set(j)) then
c     noder(1)=ialset(istartset(j))
c     noder(2)=j         
c     endif    
c     enddo
          call cident81(set,leftset,nset,j)
          if(j.gt.0) then
            if(leftset.eq.set(j)) then
              nodel(1)=ialset(istartset(j))
              nodel(2)=j         
            endif
          endif
          call cident81(set,rightset,nset,j)
          if(j.gt.0) then
            if(rightset.eq.set(j)) then
              noder(1)=ialset(istartset(j))
              noder(2)=j         
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
              if(nodel(1).eq.kon(l))then
                nodel(3)=k
                left=.true.
              endif
            endif
            if(.not.right) then
              if(noder(1).eq.kon(l))then
                noder(3)=k
                right=.true.
              endif
            endif
            if(left.and.right) exit loop
          enddo
        enddo loop
c     !     
c     !     Looking for the cyclic symmetry the nodes belong to
c     !     
c     !     
c     do k=1,ne-1
c     write(*,*) 'multistages ',k,ipkon(k),ipkon(k+1)
c     do l=ipkon(k)+1,ipkon(k+1)
c     if (nodel(1).eq.kon(l)) then
c     nodel(3)=k
c     elseif (noder(1).eq.kon(l)) then
c     noder(3)=k
c     endif
c     enddo   
c     enddo
c     !     
c     !     For the last element a different loop is needed 
c     !     
c     do l=ipkon(ne),nkon
c     if (nodel(1).eq.kon(l)) then
c     nodel(3)=k
c     !     
c     elseif (noder(1).eq.kon(l)) then
c     noder(3)=k
c     endif
c     enddo   
!     
        do j=1,mcs
          do l=istartset(int(cs(13,j))),iendset(int(cs(13,j)))
            if (ialset(l).eq.nodel(3)) then
              nodel(4)=cs(1,j)
              nodel(5)=j
            elseif (ialset(l).eq.noder(3)) then
              noder(4)=cs(1,j)
              noder(5)=j
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
!     Sorting such that rightset is independent with lesser angle
!     and noder is in the independent set
!     
        if (nodel(4).ge.noder(4)) then
          indcs=nodel(5)
          phi0=(2.d0*pi)/nodel(4)
          temp=rightset;
          rightset=leftset
          leftset=temp
          do j=1,5
            itemp(j)=noder(j)
            noder(j)=nodel(j)
            nodel(j)=itemp(j)
          enddo
        else
          indcs=noder(5)
          phi0=(2.d0*pi)/noder(4)
        endif
!     
!     Looking for a node on the independent cyclic symmetry side
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
!     Defining the rotary matrix for the tie level
!     
        T(1,1)=csab(4)-csab(1)
        T(1,2)=csab(5)-csab(2)
        T(1,3)=csab(6)-csab(3)
        dd=dsqrt(T(1,1)*T(1,1)+T(1,2)*T(1,2)+T(1,3)*T(1,3))
        T(1,1)=T(1,1)/dd
        T(1,2)=T(1,2)/dd
        T(1,3)=T(1,3)/dd
!     
!     Defining the Position of the Leftnode, which contains the angle boundary
!     for the second parameter for the tie level
!     
        xap=co(1,node_cycle)-csab(1)
        yap=co(2,node_cycle)-csab(2)
        zap=co(3,node_cycle)-csab(3)
!     
        zp=xap*T(1,1)+yap*T(1,2)+zap*T(1,3)
        rp=((xap-T(1,1)*zp)**2+(yap-zp*T(1,2))**2+
     &       (zap-zp*T(1,3))**2)
        rp=dsqrt(rp)
!     
!     Performing the vector product for the third axis of the rotary matrix
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
!     Inverting the rotary matrix to rotate the node back afterwards
!     
        call invert3D(T,Tinv,3)
!     
!     Writing a secondary rotary matrix which rotates the nodes on
!     the tie level with the angle phi
!     
        l=0
        do j=istartset(noder(2)),iendset(noder(2))
          l=l+1
          node(1)=co(1,ialset(j))-csab(1)
          node(2)=co(2,ialset(j))-csab(2)
          node(3)=co(3,ialset(j))-csab(3)
          call Mprod(T,node,coord,3)
          xind(l)=coord(2)
          yind(l)=coord(3) 
          nx(l)=l
          ny(l)=l
          ics(l)=ialset(j)
          xind0(l)=xind(l)
          yind0(l)=yind(l)
          rind=dsqrt(coord(2)**2+coord(3)**2)
          phi=datan2(-coord(3),coord(2))
          if (l.gt.1.d0) then
            dxmax=max(dxmax,dabs(coord(1)))
            drmax=max(drmax,dabs(rind))
            
            dxmin=min(dxmin,dabs(coord(1)))
            drmin=min(drmin,dabs(rind))
            phi_min=min(phi_min,phi)
          else
            dxmax=dabs(coord(1))
            phi_min=phi
            drmax=rind
            dxmin=dabs(coord(1))
            drmin=rind
          endif
        enddo
!     
        cylindrical=.false.
        if ((dxmax-dxmin).ge.(drmax-drmin)) then
          l=0.d0
          do j=istartset(noder(2)),iendset(noder(2))
            l=l+1
            node(1)=co(1,ialset(j))-csab(1)
            node(2)=co(2,ialset(j))-csab(2)
            node(3)=co(3,ialset(j))-csab(3)
            call Mprod(T,node,coord,3)
            xind(l)=coord(2)
            yind(l)=datan2(-coord(3),coord(2))
            nx(l)=l
            ny(l)=l
            ics(l)=ialset(j)
            xind0(l)=xind(l)
            yind0(l)=yind(l)
          enddo
          cylindrical=.true.
c     write(*,*) 'Multistage Tie',tieset(1,i)(1:80),
c     &              'is horizontally triangulated'
        else
c     write(*,*) 'Multistage Tie',tieset(1,i)(1:80),
c     &               'is vertically triangulated'
        endif
c     
c     Sorting the coordinates for the further consideration
c     
        ncsnodes=l
        kflag=2
        call dsort(xind,nx,ncsnodes,kflag)
        call dsort(yind,ny,ncsnodes,kflag)
!     
        call triangulate(ics,xind0,yind0,ncsnodes,
     &       rcscg,rcs0cg,zcscg,zcs0cg,nrcg,nzcg,jcs,kontri,
     &       straight,ne,ipkon,kon,lakon,lcs,netri,ifacetet,
     &       inodface)
!     
        do j=istartset(nodel(2)),iendset(nodel(2))
          node(1)=co(1,ialset(j))-csab(1)
          node(2)=co(2,ialset(j))-csab(2)
          node(3)=co(3,ialset(j))-csab(3)
!     
!     Rotating into local coordinates "coord"
!     
          call Mprod(T,node,coord,3)
!     
!     Determining the number of segments for the rotation
!     
          phi=datan2(-coord(3),coord(2))
          if (phi.gt.(-1.d-5+phi_min)) then
            kseg=int(noder(4)*0.5d0*(phi-phi_min)/pi)
          else 
            kseg=int(noder(4)*0.5d0*(2.d0*pi+(phi-phi_min))/pi)
          endif
!     
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
!     segments to "node" in local coordinates
!     
          call Mprod(T2D,coord,node,3)
!     
!     Copying the local coordinates to the local variables
!     
          if (cylindrical) then
            yp=node(1)
            zp=phi
          else
            xp=node(1)
            yp=node(2)
            zp=node(3)
          endif
          noder0=nk+1
!     
!     Rotating back to global coordinates to "coord"
!     
          call Mprod(Tinv,node,coord,3)
          co(1,noder0)=coord(1)
          co(2,noder0)=coord(2)
          co(3,noder0)=coord(3)
!     
          ier=0
          call linkdissimilar(co,csab,
     &         rcscg,rcs0cg,zcscg,zcs0cg,nrcg,nzcg,
     &         straight,nodef,ratio,nterms,yp,zp,netri,
     &         noder0,ifacetet,inodface,ialset(j),
     &         T(1,1),T(1,2),T(1,3),ier,multistage,icount)
!     
          if ((kseg.eq.(noder(4))-1.d0).and.(replace)) then
            cs(1,mcs+1)=-noder(4)
            do k=2,17
              cs(k,mcs+1)=cs(k,noder(5))
            enddo
            cs(13,mcs+1)=0.5d0
            mcs=mcs+1
            replace=.false.
          endif
!     
!     
          call transformatrix(csab,co(1,ialset(j)),al)
          call transformatrix(csab,co(1,noder0),ar)   
!     
!     checking for latin hypercube positions in matrix al none of
!     which are zero
!     
          do inum=1,6
            if((dabs(al(lathyp(1,inum),1)).gt.1.d-3).and.
     &           (dabs(al(lathyp(2,inum),2)).gt.1.d-3).and.
     &           (dabs(al(lathyp(3,inum),3)).gt.1.d-3)) exit
          enddo
!     
          do ndir=1,3
            nmpc=nmpc+1
            ipompc(nmpc)=mpcfree 
            if((kseg.ne.0).or.(kseg.eq.(noder(4))-1.d0)) then
              labmpc(nmpc)='CYCLIC              '
              if (kseg.eq.1) then
                if (noder(5).lt.10) then
                  write(labmpc(nmpc)(7:7),'(i1)') noder(5)
                else
                  write(labmpc(nmpc)(7:8),'(i2)') noder(5)
                endif
              endif
              if (kseg.eq.(noder(4))-1.d0) then
                if (mcs.lt.10) then
                  write(labmpc(nmpc)(7:7),'(i1)') mcs
                elseif (mcs.lt.100) then
                  write(labmpc(nmpc)(7:8),'(i2)') mcs
                else
                  write(*,*)
     &                 '*ERROR in multistages: no more than 49'
                  write(*,*)
     &                 '       cyclic symmetry definitions allowed'
                  call exit(201)
                endif
              endif
              number=ndir-1 
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
     &                 '*WARNING in multistages: cyclic MPC in node'
                  write(*,*) '         ',ialset(j),
     &                 ' and direction',ndir
                  write(*,*) '         cannot be created: the'
                  write(*,*) 
     &                 '         DOF in this node is already used'
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
     &                 '*ERROR in multistages: increase memmpc_'
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
     %                   -ar(number,ndir)*ratio(k2)
                    mpcfreeold=mpcfree
                    mpcfree=nodempc(3,mpcfree)
                  endif
                  if(mpcfree.eq.0) then
                    write(*,*) 
     &                   '*ERROR in multistages: increase memmpc_'
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
     &                 '*WARNING in multistages: cyclic MPC in node'
                  write(*,*) '         ',ialset(j),
     &                 ' and direction',ndir
                  write(*,*) '         cannot be created: the'
                  write(*,*) 
     &                 '         DOF in this node is already used'
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
     &               '*ERROR in multistages: increase memmpc_'
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
     &                 '*ERROR in multistages: increase memmpc_'
                  call exit(201)
                endif
              enddo  
            endif
            nodempc(3,mpcfreeold)=0
          enddo
!     
        enddo                   !Loop over nodes on dependent side
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
!     *********Ending the main loop over all multistage ties***********
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
      
