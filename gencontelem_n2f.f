!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine gencontelem_n2f(tieset,ntie,itietri,ne,ipkon,kon,
     &     lakon,cg,straight,ifree,
     &     koncont,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,ielmat,
     &     elcon,istep,iinc,iit,ncmat_,ntmat_,
     &     nmethod,mi,imastop,nslavnode,islavnode,islavsurf,
     &     itiefac,areaslav,iponoels,inoels,springarea,set,nset,
     &     istartset,
     &     iendset,ialset,tietol,reltime,filab,nasym,xnoels,icutb,ne0,
     &     jobnamef,mortar,auw,jqw,iroww,nzsw,imastnode,nmastnode)
!     
!     generate contact elements for the slave contact nodes
!     
      implicit none
!     
      character*8 lakon(*)
      character*33 cfile
      character*81 tieset(3,*),slavset,set(*),noset,setname
      character*87 filab(*)
      character*132 jobnamef(*)
!     
      integer ntie,ifree,nasym,icutb,ne0,mortar,jqw(*),iroww(*),
     &     itietri(2,ntie),ipkon(*),kon(*),koncont(4,*),ne,node,
     &     neigh(1),iflag,kneigh,i,j,k,l,isol,iset,idummy(1),
     &     itri,ll,kflag,n,nx(*),ny(*),istep,iinc,mi(*),nzsw,
     &     nz(*),nstart,ielmat(mi(3),*),imat,ifaceq(8,6),ifacet(6,4),
     &     ifacew1(4,5),ifacew2(8,5),nelem,jface,indexe,iit,
     &     nface,nope,nodef(9),ncmat_,ntmat_,index1,indexel,
     &     nmethod,ifaces,jfaces,irelslavface,number(4),lenset,
     &     imastop(3,*), itriangle(100),ntriangle,ntriangle_,itriold,
     &     itrinew,id,nslavnode(*),islavnode(*),islavsurf(2,*),
     &     itiefac(2,*),iponoels(*),inoels(2,*),konl(26),nelems,m,
     &     mint2d,nopes,ipos,nset,istartset(*),iendset(*),m1,m2,
     &     ialset(*),length,iorder(4),imastnode(*),nmastnode(*),
     &     jbasis,kbasis,filedebug
!     
      real*8 cg(3,*),straight(16,*),co(3,*),vold(0:mi(2),*),p(3),
     &     dist,xo(*),yo(*),zo(*),x(*),y(*),z(*),c0coef,auw(*),
     &     beta,c0,elcon(0:ncmat_,ntmat_,*),weight,t1(3),t2(3),dt1,
     &     areaslav(*),springarea(2,*),xl2(3,9),area,xi,et,shp2(7,9),
     &     xs2(3,2),xsj2(3),adjust,tietol(4,*),reltime,xns(3,3),
     &     clear,ratio(9),pl(3,9),
     &     pproj(3),al(3),xn(3),xm(3),dm,xnoels(*)
!     
!     nodes per face for hex elements
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
!     
!     nodes per face for tet elements
!     
      data ifacet /1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/
!     
!     nodes per face for linear wedge elements
!     
      data ifacew1 /1,3,2,0,
     &     4,5,6,0,
     &     1,2,5,4,
     &     2,3,6,5,
     &     3,1,4,6/
!     
!     nodes per face for quadratic wedge elements
!     
      data ifacew2 /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     3,1,4,6,9,13,12,15/
!     
!     flag for shape functions
!     
      data iflag /2/
      data indexel /0/
!     
      save indexel
!     
      include "gauss.f"
!     
!     debugging file
!     
      filedebug=0
      if(filedebug.gt.0) then
        do i=1,132
          if(jobnamef(1)(i:i).eq.' ') exit
        enddo
        i=i-1
        cfile=jobnamef(1)(1:i)//'.wb'
        open(88,file=cfile,status='unknown',position='append')
        setname(1:15)='contactelements'
        lenset=15
        number(1)=istep
        number(2)=iinc
        number(3)=icutb+1
        number(4)=iit
        if(number(4).le.0) number(4)=1
        do i=1,4
          setname(lenset+1:lenset+1)='_'
          if(i.eq.1) then
            setname(lenset+2:lenset+3)='st'
          elseif(i.eq.2) then
            setname(lenset+2:lenset+3)='in'
          elseif(i.eq.3) then
            setname(lenset+2:lenset+3)='at'
          else
            setname(lenset+2:lenset+3)='it'
          endif
          if(number(i).lt.10) then
            write(setname(lenset+4:lenset+4),'(i1)') number(i)
            lenset=lenset+4
          elseif(number(i).lt.100) then
            write(setname(lenset+4:lenset+5),'(i2)') number(i)
            lenset=lenset+5
          elseif(number(i).lt.1000) then
            write(setname(lenset+4:lenset+6),'(i3)') number(i)
            lenset=lenset+6
          else
            write(*,*) '*ERROR in gencontelem_n2f: no more than 1000'
            write(*,*) '       steps/increments/cutbacks/iterations'
            write(*,*) '       allowed (for output in '
            write(*,*) '       [JOBNAME].wb)'
            stop
          endif
        enddo
      endif
!     
!     opening a file to store the contact spring elements
!     
      if(filab(1)(3:3).eq.'C') then
c     iteller=iteller+1
!     
        do i=1,132
          if(jobnamef(1)(i:i).eq.' ') exit
        enddo
        i=i-1
        cfile=jobnamef(1)(1:i)//'.cel'
        open(27,file=cfile,status='unknown',position='append')
!     
        setname(1:15)='contactelements'
        lenset=15
        number(1)=istep
        number(2)=iinc
        number(3)=icutb+1
        number(4)=iit
        if(number(4).le.0) number(4)=1
        do i=1,4
          setname(lenset+1:lenset+1)='_'
          if(i.eq.1) then
            setname(lenset+2:lenset+3)='st'
          elseif(i.eq.2) then
            setname(lenset+2:lenset+3)='in'
          elseif(i.eq.3) then
            setname(lenset+2:lenset+3)='at'
          else
            setname(lenset+2:lenset+3)='it'
          endif
          if(number(i).lt.10) then
            write(setname(lenset+4:lenset+4),'(i1)') number(i)
            lenset=lenset+4
          elseif(number(i).lt.100) then
            write(setname(lenset+4:lenset+5),'(i2)') number(i)
            lenset=lenset+5
          elseif(number(i).lt.1000) then
            write(setname(lenset+4:lenset+6),'(i3)') number(i)
            lenset=lenset+6
          else
            write(*,*) '*ERROR in gencontelem_n2f: no more than 1000'
            write(*,*) '       steps/increments/cutbacks/iterations'
            write(*,*) '       allowed (for output in '
            write(*,*) '       contactelements.inp)'
            stop
          endif
        enddo
      endif
!     
!     massless contact: initialization of nzsw (contact matrix W_b)
!     
      if(mortar.eq.-1) nzsw=0
!     
      do i=1,ntie
        if(tieset(1,i)(81:81).ne.'C') cycle
        kneigh=1
        slavset=tieset(2,i)
        imat=int(tietol(2,i))
        c0coef=elcon(4,1,imat)
!     
!     check whether an adjust node set has been defined
!     only checked at the start of the first step
!     
c     if((istep.eq.1).and.(iinc.eq.1).and.(iit.le.0)) then
        if((istep.eq.1).and.(iit.lt.0)) then
          iset=0
          if(tieset(1,i)(1:1).ne.' ') then
            noset(1:80)=tieset(1,i)(1:80)
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
c     do iset=1,nset
c     if(set(iset).eq.noset) exit
c     enddo
            call cident81(set,noset,nset,id)
            iset=nset+1
            if(id.gt.0) then
              if(noset.eq.set(id)) then
                iset=id
              endif
            endif
            kflag=1
            call isortii(ialset(istartset(iset)),idummy,
     &           iendset(iset)-istartset(iset)+1,kflag)
          endif
        endif
!     
!     determine the area of the slave surfaces
!     
        do l = itiefac(1,i), itiefac(2,i)
          ifaces = islavsurf(1,l)
          nelems = int(ifaces/10)
          jfaces = ifaces - nelems*10
!     
!     Decide on the max integration points number, just consider 2D situation 
!     
          if(lakon(nelems)(4:5).eq.'8R') then
            mint2d=1
            nopes=4
            nope=8
          elseif(lakon(nelems)(4:4).eq.'8') then
            mint2d=4
            nopes=4
            nope=8
          elseif(lakon(nelems)(4:6).eq.'20R') then
            mint2d=4
            nopes=8
            nope=20
          elseif(lakon(nelems)(4:5).eq.'20') then
            mint2d=9
            nopes=8
            nope=20
          elseif(lakon(nelems)(4:5).eq.'10') then
            mint2d=3
            nopes=6
            nope=10
          elseif(lakon(nelems)(4:4).eq.'4') then
            mint2d=1
            nopes=3
            nope=4
!     
!     treatment of wedge faces
!     
          elseif(lakon(nelems)(4:4).eq.'6') then
            mint2d=1
            nope=6
            if(jfaces.le.2) then
              nopes=3
            else
              nopes=4
            endif
          elseif(lakon(nelems)(4:5).eq.'15') then
            nope=15
            if(jfaces.le.2) then
              mint2d=3
              nopes=6
            else
              mint2d=4
              nopes=8
            endif
          endif
!     
!     actual position of the nodes belonging to the
!     slave surface
!     
          do j=1,nope
            konl(j)=kon(ipkon(nelems)+j)
          enddo
!     
          if((nope.eq.20).or.(nope.eq.8)) then
            do m=1,nopes
              do j=1,3
                xl2(j,m)=co(j,konl(ifaceq(m,jfaces)))+
     &               vold(j,konl(ifaceq(m,jfaces)))
              enddo
            enddo
          elseif((nope.eq.10).or.(nope.eq.4)) then
            do m=1,nopes
              do j=1,3
                xl2(j,m)=co(j,konl(ifacet(m,jfaces)))+
     &               vold(j,konl(ifacet(m,jfaces)))
              enddo
            enddo
          elseif(nope.eq.15) then
            do m=1,nopes
              do j=1,3
                xl2(j,m)=co(j,konl(ifacew2(m,jfaces)))+
     &               vold(j,konl(ifacew2(m,jfaces)))
              enddo
            enddo
          else
            do m=1,nopes
              do j=1,3
                xl2(j,m)=co(j,konl(ifacew1(m,jfaces)))+
     &               vold(j,konl(ifacew1(m,jfaces)))
              enddo
            enddo
          endif
!     
!     calculating the area of the slave face
!     
          area=0.d0
          do m=1,mint2d
            if((lakon(nelems)(4:5).eq.'8R').or.
     &           ((lakon(nelems)(4:4).eq.'6').and.(nopes.eq.4))) then
              xi=gauss2d1(1,m)
              et=gauss2d1(2,m)
              weight=weight2d1(m)
            elseif((lakon(nelems)(4:4).eq.'8').or.
     &             (lakon(nelems)(4:6).eq.'20R').or.
     &             ((lakon(nelems)(4:5).eq.'15').and.
     &             (nopes.eq.8))) then
              xi=gauss2d2(1,m)
              et=gauss2d2(2,m)
              weight=weight2d2(m)
            elseif(lakon(nelems)(4:4).eq.'2') then
              xi=gauss2d3(1,m)
              et=gauss2d3(2,m)
              weight=weight2d3(m)
            elseif((lakon(nelems)(4:5).eq.'10').or.
     &             ((lakon(nelems)(4:5).eq.'15').and.
     &             (nopes.eq.6))) then
              xi=gauss2d5(1,m)
              et=gauss2d5(2,m)
              weight=weight2d5(m)
            elseif((lakon(nelems)(4:4).eq.'4').or.
     &             ((lakon(nelems)(4:4).eq.'6').and.
     &             (nopes.eq.3))) then
              xi=gauss2d4(1,m)
              et=gauss2d4(2,m)
              weight=weight2d4(m)
            endif
!     
            if(nopes.eq.8) then
              call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.4) then
              call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.6) then
              call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            else
              call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            endif
            area=area+weight*dsqrt(xsj2(1)**2+xsj2(2)**2+
     &           xsj2(3)**2)
          enddo
          areaslav(l)=area
        enddo
!     
!     search a master face for each slave node and generate a contact
!     spring element if successful
!     
        nstart=itietri(1,i)-1
        n=itietri(2,i)-nstart
        if(n.lt.kneigh) kneigh=n
        do j=1,n
          xo(j)=cg(1,nstart+j)
          x(j)=xo(j)
          nx(j)=j
          yo(j)=cg(2,nstart+j)
          y(j)=yo(j)
          ny(j)=j
          zo(j)=cg(3,nstart+j)
          z(j)=zo(j)
          nz(j)=j
        enddo
        kflag=2
        call dsort(x,nx,n,kflag)
        call dsort(y,ny,n,kflag)
        call dsort(z,nz,n,kflag)
!     
        do j=nslavnode(i)+1,nslavnode(i+1)
          if(mortar.eq.-1)then
            if((istep.eq.1).and.(iit.lt.0)) springarea(2,j)=0.d0
          else
            if(iit.le.0) springarea(2,j)=0.d0
          endif
          node=islavnode(j)
!     
!     calculating the area corresponding to the
!     slave node; is made up of the area
!     of the neighboring slave faces
!     
          area=0.d0
          index1=iponoels(node)
          do
            if(index1.eq.0) exit
            irelslavface=inoels(1,index1)
            if((itiefac(1,i).le.irelslavface).and.
     &           (irelslavface.le.itiefac(2,i))) then
              area=area+areaslav(irelslavface)*
     &             xnoels(index1)
            endif
            index1=inoels(2,index1)
          enddo
!     
          do k=1,3
            p(k)=co(k,node)+vold(k,node)
          enddo
!     
!     determining the kneigh neighboring master contact
!     triangle centers of gravity
!     
          call near3d(xo,yo,zo,x,y,z,nx,ny,nz,p(1),p(2),p(3),
     &         n,neigh,kneigh)
!     
          isol=0
!     
          itriold=0
          itri=neigh(1)+itietri(1,i)-1
          ntriangle=0
          ntriangle_=100
!     
          loop1: do
          do l=1,3
            ll=4*l-3
            dist=straight(ll,itri)*p(1)+
     &           straight(ll+1,itri)*p(2)+
     &           straight(ll+2,itri)*p(3)+
     &           straight(ll+3,itri)
!     
!     1.d-6 was increased to 1.d-3 on 19/04/2012
!     this is important for 2d-calculations or
!     calculations for which structures fit exactly
!     at their boundaries
!     
            if(dist.gt.1.d-3*dsqrt(area)) then
              itrinew=imastop(l,itri)
              if(itrinew.eq.0) then
c     write(*,*) '**border reached'
                exit loop1
              elseif((itrinew.lt.itietri(1,i)).or.
     &               (itrinew.gt.itietri(2,i))) then
c     write(*,*) '**border reached'
                exit loop1
              elseif(itrinew.eq.itriold) then
c     write(*,*) '**solution in between triangles'
                isol=itri
                exit loop1
              else
                call nident(itriangle,itrinew,ntriangle,id)
                if(id.gt.0) then
                  if(itriangle(id).eq.itrinew) then
c     write(*,*) '**circular path;no solution'
                    exit loop1
                  endif
                endif
                ntriangle=ntriangle+1
                if(ntriangle.gt.ntriangle_) then
c     write(*,*) '**too many iterations'
                  exit loop1
                endif
                do k=ntriangle,id+2,-1
                  itriangle(k)=itriangle(k-1)
                enddo
                itriangle(id+1)=itrinew
                itriold=itri
                itri=itrinew
                cycle loop1
              endif
            elseif(l.eq.3) then
c     write(*,*) '**regular solution'
              isol=itri
              exit loop1
            endif
          enddo
        enddo loop1
!     
!     check whether distance is larger than c0:
!     no element is generated
!     
        if(isol.ne.0) then
!     
!     determining the clearance
!     
!     identifying the element face to which the
!     triangle belongs
!     
          nelem=int(koncont(4,itri)/10.d0)
          jface=koncont(4,itri)-10*nelem
!     
          indexe=ipkon(nelem)
          if(lakon(nelem)(4:5).eq.'20') then
            nopes=8
            nface=6
          elseif(lakon(nelem)(4:4).eq.'8') then
            nopes=4
            nface=6
          elseif(lakon(nelem)(4:5).eq.'10') then
            nopes=6
            nface=4
          elseif(lakon(nelem)(4:4).eq.'4') then
            nopes=3
            nface=4
          elseif(lakon(nelem)(4:5).eq.'15') then
            if(jface.le.2) then
              nopes=6
            else
              nopes=8
            endif
            nface=5
            nope=15
          elseif(lakon(nelem)(4:4).eq.'6') then
            if(jface.le.2) then
              nopes=3
            else
              nopes=4
            endif
            nface=5
            nope=6
          else
            cycle
          endif
!     
!     determining the nodes of the face
!     
          if(nface.eq.4) then
            do k=1,nopes
              nodef(k)=kon(indexe+ifacet(k,jface))
            enddo
          elseif(nface.eq.5) then
            if(nope.eq.6) then
              do k=1,nopes
                nodef(k)=kon(indexe+ifacew1(k,jface))
              enddo
            elseif(nope.eq.15) then
              do k=1,nopes
                nodef(k)=kon(indexe+ifacew2(k,jface))
              enddo
            endif
          elseif(nface.eq.6) then
            do k=1,nopes
              nodef(k)=kon(indexe+ifaceq(k,jface))
            enddo
          endif
!     
!     orthogonal projection on the element face
!     
          do k=1,nopes
            do l=1,3
              pl(l,k)=co(l,nodef(k))+vold(l,nodef(k))
            enddo
          enddo
          do l=1,3
            pproj(l)=p(l)
          enddo
!     
!     project node pproj to the face 
!     
          call attach_2d(pl,pproj,nopes,ratio,dist,xi,et)
!     
          do l=1,3
            al(l)=p(l)-pproj(l)
          enddo
!     
!     determining the jacobian vector on the surface 
!     
          if(nopes.eq.8) then
            call shape8q(xi,et,pl,xm,xs2,shp2,iflag)
          elseif(nopes.eq.4) then
            call shape4q(xi,et,pl,xm,xs2,shp2,iflag)
          elseif(nopes.eq.6) then
            call shape6tri(xi,et,pl,xm,xs2,shp2,iflag)
          else
            call shape3tri(xi,et,pl,xm,xs2,shp2,iflag)
          endif
!     
!     normal on the surface
!     
          dm=dsqrt(xm(1)*xm(1)+xm(2)*xm(2)+xm(3)*xm(3))
          do l=1,3
            xn(l)=xm(l)/dm
          enddo
!     
!     massless contact: filling matrix Wb
!     including normal vectors on master side
!     
          if(mortar.eq.-1) then
!     
!     determining a local system (cf. springforc_n2f.f)
!     
            if(1.d0-dabs(xn(1)).lt.1.5231d-6) then       
!     
!     calculating the local directions on master surface
!     
              t1(1)=-xn(3)*xn(1)
              t1(2)=-xn(3)*xn(2)
              t1(3)=1.d0-xn(3)*xn(3)
            else
              t1(1)=1.d0-xn(1)*xn(1)
              t1(2)=-xn(1)*xn(2)
              t1(3)=-xn(1)*xn(3)
            endif
!     
!     normalize
!     
            dt1=dsqrt(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))
            do m1=1,3
              t1(m1)=t1(m1)/dt1
            enddo
            t2(1)=xn(2)*t1(3)-xn(3)*t1(2)
            t2(2)=xn(3)*t1(1)-xn(1)*t1(3)
            t2(3)=xn(1)*t1(2)-xn(2)*t1(1)           
!     
!     storing in a matrix xns=[xn t1 t2]
!     
            do m1=1,3
              xns(m1,1)=xn(m1)
              xns(m1,2)=t1(m1)
              xns(m1,3)=t2(m1)
            enddo
!     
!     sorting the nodes of the master face in
!     ascending order
!     
            length=nmastnode(i+1)-nmastnode(i)
            do k=1,nopes
              iorder(k)=k
            enddo
            kflag=2
            call isortii(nodef,iorder,nopes,kflag)
!     
            if(filedebug.gt.0)then
! Writing Wb matrix entries
              write(88,102) setname(1:lenset)
 102          format('*DEBUGDATA,TYPE=WB_COL,ELSET=',A)
              write(88,*) nopes ! tot master nodes associated to slave node
              write(88,*) node  ! slave node
              write(88,'(I12)')  nodef(1:nopes) ! list master nodes
              write(88,'(E12.6)')  -ratio(iorder) ! list slave nodes
!     write(*,'(E12.6)')  -ratio(iorder)! list slave nodes
            endif
!     
            jbasis=3*(j-1)
c     write(*,*) 'jbasis', jbasis
            do m2=1,3
c     write(*,*)'ix@jqw',(jbasis+m2),'jqw(ix)',(nzsw+1)
!     Setting new column
              jqw(jbasis+m2)=nzsw+1
!     Filling for DEPENDENT
              do m1=1,3
                nzsw=nzsw+1
                auw(nzsw)=xns(m1,m2)
                iroww(nzsw)=jbasis+m1
              enddo
!     
!     Filling for INDEPENDENT
              do k=1,nopes
                call nident(imastnode(nmastnode(i)+1),
     &               nodef(k),length,id)
                kbasis=3*(nslavnode(ntie+1)+nmastnode(i)+id-1)
                do m1=1,3
                  nzsw=nzsw+1
                  auw(nzsw)=-ratio(iorder(k))*xns(m1,m2)
                  iroww(nzsw)=kbasis+m1
                enddo
              enddo
!     
            enddo
            if(mortar.ne.-1) cycle ! ommitted to capture clearances???
          endif
!     
!     distance from surface along normal (= clearance)
!     
          clear=al(1)*xn(1)+al(2)*xn(2)+al(3)*xn(3)
c     if((istep.eq.1).and.(iit.lt.0.d0)) then
c     if(clear.lt.0.d0) then
c     springarea(2,j)=clear
c     endif
c     endif
c     if(nmethod.eq.1) then
c     clear=clear-springarea(2,j)*(1.d0-reltime)
c     endif
!     
!     check for an adjust parameter (only at the start
!     of the first step)
!     
          if((istep.eq.1).and.(iit.lt.0)) then
            if(iset.ne.0) then
!     
!     check whether node belongs to the adjust node
!     set
!     
              call nident(ialset(istartset(iset)),node,
     &             iendset(iset)-istartset(iset)+1,id)
              if(id.gt.0) then
                if(ialset(istartset(iset)+id-1).eq.node) then
                  do k=1,3
                    co(k,node)=co(k,node)-
     &                   clear*straight(12+k,itri)
                  enddo
                  clear=0.d0
                endif
              endif
            elseif(dabs(tietol(1,i)).ge.2.d0) then
!     
!     adjust parameter
!     
              adjust=dabs(tietol(1,i))-2.d0
              if(clear.le.adjust) then
                do k=1,3
                  co(k,node)=co(k,node)-
     &                 clear*straight(12+k,itri)
                enddo
                clear=0.d0
              endif
            endif
          endif
!     
!     storing  the clearance if penetration takes place,
!     i.e. clear <0 at the start of the first step
!     
          if((istep.eq.1).and.(iit.lt.0.d0)) then
            if(clear.lt.0.d0) then
              springarea(2,j)=clear
            endif
          endif
!     
          if(nmethod.eq.1) then
            clear=clear-springarea(2,j)*(1.d0-reltime)
          endif
!     
!     Storing clearance into springareas for contact gap evaluations
          if(mortar.eq.-1) then
            if((istep.eq.1).and.(iit.lt.0)) then
              springarea(2,j)=clear ! working only for MASSLESS
c     
c     ! TODO CMT: Taking it into account for N2F contact:
!     if(tieset(3,i).eq.0.d0) then
!     springarea(2,j)=clear
!     else
!     !TODO CMT: user defined clearance as char, convert to real
!     springarea(2,j)=tieset(3,i)
!     endif
            endif
            cycle
          endif
!     
          if(int(elcon(3,1,imat)).eq.1) then
!     
!     exponential overclosure
!     
            beta=elcon(1,1,imat)
            c0=dlog(100.d0)/beta
          else
!     
!     linear or tabular overclosure
!     
            if(dabs(area).gt.0.d0) then
              c0=c0coef*dsqrt(area)
            else
              c0=1.d-10
            endif
          endif
          if(clear.gt.c0) then
            isol=0
!     
          endif
        elseif(mortar.eq.-1) then
!
!         taking care of a zero column in the W-matrix for massless contact
!
          jbasis=3*(j-1)
          do m2=1,3
            jqw(jbasis+m2)=nzsw+1
          enddo
          cycle
        endif
!     
        if(isol.ne.0) then
!     
!     plane spring
!     
          ne=ne+1
          ipkon(ne)=ifree
          lakon(ne)='ESPRNGC '
          ielmat(1,ne)=imat
!     
!     nasym indicates whether at least one contact
!     spring elements exhibits friction in the present
!     step. If so, nasym=1, else nasym=0; nasym=1
!     triggers the asymmetric equation solver
!     
          if(ncmat_.ge.7) then
            if(elcon(6,1,imat).gt.0) then
              nasym=1
            endif
          endif
          if(ncmat_.ge.8) then
            if(elcon(8,1,imat).gt.0) then
              nasym=1
            endif
          endif
!     
!     storing the area corresponding to the slave node
!     
          springarea(1,j)=area
!     
          do k=1,nopes
            kon(ifree+k)=nodef(k)
          enddo
          ifree=ifree+nopes+1
          kon(ifree)=node
          ifree=ifree+1
          kon(ifree)=j
!     
          write(lakon(ne)(8:8),'(i1)') nopes
!     
          indexel=indexel+1
!     
          if((nopes.eq.3).or.(nopes.eq.6)) then
            if(filab(1)(3:3).eq.'C') then
              write(27,100) setname(1:lenset)
 100          format('*ELEMENT,TYPE=C3D4,ELSET=',A)
              write(27,*) ne0+indexel,',',nodef(1),',',
     &             nodef(2),',',nodef(3),',',node
            endif
          else
            if(filab(1)(3:3).eq.'C') then
              write(27,101) setname(1:lenset)
 101          format('*ELEMENT,TYPE=C3D6,ELSET=',A)
              write(27,*) ne0+indexel,',',nodef(2),',',node,
     &             ',',
     &             nodef(3),',',nodef(1),',',node,',',nodef(4)
            endif
          endif
        endif
!     
      enddo
      enddo
!     
!     closing the file containing the contact elements
!     
      if(filab(1)(3:3).eq.'C') then
        close(27)
      endif
      if(filedebug.gt.0) close(88)
!     
!     complete the definition of jqw for massless contact
!     
      if(mortar.eq.-1) then
        jqw(3*nslavnode(ntie+1)+1)=nzsw+1
      endif
!     
!     DEBUG PRINT WB MATRIX
      if((filab(1)(3:3).eq.'C').and.(mortar.eq.-1))then
!     nmastnode(ntie+1) ! QUESTION: total master nodes?
        open(unit=31416,file='WbPre_colsort.csv')
        do i=1,nslavnode(ntie+1)
          write(31416,*) islavnode(i) ,',N'
          write(31416,*) islavnode(i) ,',T1'
          write(31416,*) islavnode(i) ,',T2'
        enddo
        close(31416)

        open(unit=31416,file='WbPre_size.csv')
        write(31416,*) (nmastnode(ntie+1)+nslavnode(ntie+1))*3
        write(31416,*) nslavnode(ntie+1)*3
        close(31416)

        open(unit=31416,file='WbPre_auw.csv')
        do i=1,nzsw
          write(31416,*) auw(i)
        enddo
        close(31416)

        open(unit=31416,file='WbPre_iroww.csv')
        do i=1,nzsw
          write(31416,*) iroww(i)
        enddo
        close(31416)

        open(unit=31416,file='WbPre_jqw.csv')
        do i=1,(nslavnode(ntie+1)*3)+1
          write(31416,*) jqw(i)
        enddo
        close(31416)
      endif
!     
      return
      end
