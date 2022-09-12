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
      subroutine extrapolatecontact(yi,yn,ipkon,inum,kon,lakon,nfield,
     &  nk,ne,mi,ndim,co,cflag,vold,force,pslavsurf,islavact,islavnode,
     &  nslavnode,ntie,islavsurf,ielprop,prop,ielmat,ne0)
!
!     extrapolates contact values at the integration points to the 
!     nodes (for face-to-face penalty contact)
!
      implicit none
!
      logical force
!
      character*1 cflag
      character*8 lakon(*),lakonl
!
      integer ipkon(*),inum(*),kon(*),mi(*),ne,indexe,nope,
     &  nonei20(3,12),nfield,nonei10(3,6),nk,i,j,k,l,ndim,
     &  nonei15(3,9),m,iflag,jj,indexc,islavsurf(2,*),ll,
     &  indexcj,ifacew1(4,5),islavnode(*),nslavnode(*),ntie,
     &  nlayer,nopeexp,ifacew2(8,5),ifaceq(8,6),ifacet(6,4),
     &  nopespring,ifaces,nopespringj,ifacej,jfaces,n,nelems,
     &  idgn,idgnr,idglda,idgip(4),idgldb,idginfo,igauss,islavact(*),
     &  konl(26),node,nopes,ielprop(*),ielmat(mi(3),*),ne0,
     &  nodes(8)
!
      real*8 yi(ndim,mi(1),*),yn(nfield,*),field(nfield,20*mi(3)),
     &  co(3,*),xi,et,vold(0:mi(2),*),xs2(3,7),xsj2(3),shp2(7,8),
     &  aa(4,4),bb(4,6),pl(3,8),pslavsurf(3,*),xslavnor(3,nk),
     &  cc(3,4),dd(3,6),c_limit(2,nfield),nodepos(4,2),xn(3),
     &  t1(3),t2(3),trac(3),xquad(2,9),xtri(2,7),xl2s(3,9),
     &  stn(6,nk),dt1,dl,a2(6,2),a4(4,4),a27(20,27),a9(6,9),
     &  a8(8,8),prop(*)
!
      data nonei10 /5,1,2,6,2,3,7,3,1,8,1,4,9,2,4,10,3,4/
!
      data nonei15 /7,1,2,8,2,3,9,3,1,10,4,5,11,5,6,12,6,4,
     &  13,1,4,14,2,5,15,3,6/
!
      data nonei20 /9,1,2,10,2,3,11,3,4,12,4,1,
     &  13,5,6,14,6,7,15,7,8,16,8,5,
     &  17,1,5,18,2,6,19,3,7,20,4,8/
!
!     nodes per face for hex elements
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
!
!     nodes per face for tet elements
!
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
!
!     nodes per face for linear wedge elements
!
      data ifacew1 /1,3,2,0,
     &             4,5,6,0,
     &             1,2,5,4,
     &             2,3,6,5,
     &             3,1,4,6/
!
!     nodes per face for quadratic wedge elements
!
      data ifacew2 /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             3,1,4,6,9,13,12,15/
!
      data iflag /2/
!
      if(nfield.eq.0) return
!
      do i=1,nk
         inum(i)=0
      enddo
!
      do i=1,nk
         do j=1,nfield
            yn(j,i)=0.d0
         enddo
      enddo
!
      i=0
      do
         i=i+1
         if(i.gt.ne) exit
!     
         if(ipkon(i).lt.0) cycle
         indexc=ipkon(i)
!     
         if((lakon(i)(1:1).ne.'E').or.(lakon(i)(7:7).ne.'C')) cycle
         nopespring=kon(indexc)
         ifaces=islavsurf(1,kon(indexc+nopespring+2))
!     
         nelems=int(ifaces/10.d0)
         lakonl=lakon(nelems)
!     
!        determining all contact elements belonging to slave face
!        "ifaces"
!
         n=i
!     
         do
            n=n+1
            if(n.gt.ne) exit
            indexcj=ipkon(n)
            nopespringj=kon(indexcj)
            ifacej=islavsurf(1,kon(indexcj+nopespringj+2))
            if(ifaces.ne.ifacej) exit
         enddo
         n=n-1
         jfaces=ifaces-10*int(ifaces/10.d0)
         indexe=ipkon(nelems)
!     
         if(lakonl(4:4).eq.'2') then
            nope=20
         elseif(lakonl(4:4).eq.'8') then
            nope=8
         elseif(lakonl(4:5).eq.'10') then
            nope=10
         elseif(lakonl(4:4).eq.'4') then
            nope=4
         elseif(lakonl(4:5).eq.'15') then
            nope=15
         elseif(lakonl(4:4).eq.'6') then
            nope=6
         elseif((lakon(i)(1:1).eq.'E').and.
     &          ((lakon(i)(7:7).eq.'A').or.
     &           (lakon(i)(7:7).eq.'2'))) then
            inum(kon(indexe+1))=inum(kon(indexe+1))+1
            inum(kon(indexe+2))=inum(kon(indexe+2))+1
            cycle
         else
            cycle
         endif
!     
         if((lakonl(4:5).eq.'20').or.(lakonl(4:4).eq.'8').or.
     &        (((lakonl(4:5).eq.'15').or.(lakonl(4:4).eq.'6')).and.
     &        (jfaces.gt.2))) then
            if(lakonl(7:8).ne.'LC') then
               field(1:nfield,1:20)=0.d0
               do k=1,4
                  do l=1,4
                     aa(k,l)=0.d0
                  enddo
                  do l=1,nfield
                     bb(k,l)=0.d0
                  enddo
               enddo
               nodepos(1,1)=-1.d0
               nodepos(1,2)=-1.d0
               nodepos(2,1)=1.d0
               nodepos(2,2)=-1.d0
               nodepos(3,1)=1.d0
               nodepos(3,2)=1.d0
               nodepos(4,1)=-1.d0
               nodepos(4,2)=1.d0
               do j=i,n
                  do k=1,nfield
                     if(j.eq.i) then
                        c_limit(1,k)=yi(k,1,j)
                        c_limit(2,k)=yi(k,1,j)
                     endif
                     if(c_limit(1,k).lt.yi(k,1,j)) then
                        c_limit(1,k)=yi(k,1,j)
                     endif
                     if(c_limit(2,k).gt.yi(k,1,j)) then
                        c_limit(2,k)=yi(k,1,j)
                     endif
                  enddo
                  indexcj=ipkon(j)
                  nopespringj=kon(indexcj)
                  igauss=kon(indexcj+nopespringj+1)
                  xi=pslavsurf(1,igauss)
                  et=pslavsurf(2,igauss)
                  if((n-i).lt.2) exit
                  if((n-i).eq.2) then
                     aa(j+1-i,1)=xi
                     aa(j+1-i,2)=et
                     do k=1,nfield
                        dd(j+1-i,k)=yi(k,1,j)
                     enddo
                  elseif((n-i).eq.3) then
                     call shape4q(xi,et,pl,xsj2,xs2,shp2,iflag)
                     do k=1,4
                        aa(j+1-i,k)=shp2(4,k)
                     enddo
                     do k=1,nfield
                        bb(j+1-i,k)=yi(k,1,j)
                     enddo
                  else
                     call shape4q(xi,et,pl,xsj2,xs2,shp2,iflag)
                     do k=1,4
                        do l=k,4
                           aa(k,l)=aa(k,l)+shp2(4,k)*shp2(4,l)
                        enddo
                     enddo
                     do k=1,4
                        do l=1,nfield
                           bb(k,l)=bb(k,l)+shp2(4,k)*yi(l,1,j)
                        enddo
                     enddo
                  endif
               enddo
               if(n.eq.i) then
                  do k=1,4
                     do l=1,nfield
                        bb(k,l)=bb(k,l)+yi(l,1,n)
                     enddo
                  enddo   
               elseif((n-i).eq.1) then
                  do j=i,n
                     do k=1,4
                        do l=1,nfield
                           bb(k,l)=bb(k,l)+yi(l,1,j)*0.5d0
                        enddo
                     enddo
                  enddo
               elseif((n-i).eq.2) then
                  do k=1,4
                     do l=1,nfield
                        call plane_eq(aa(1,1),aa(1,2),dd(1,l),
     &                       aa(2,1),aa(2,2),dd(2,l),aa(3,1),
     &                       aa(3,2),dd(3,l),nodepos(k,1),
     &                       nodepos(k,2),bb(k,l))
                     enddo
                  enddo
               else
                  if((n-i).ne.3) then
                     do k=1,4
                        do l=1,k-1
                           aa(k,l)=aa(l,k)
                        enddo
                     enddo
                  endif
                  idgn=4
                  idgnr=nfield
                  idglda=4
                  idgldb=4
                  call dgesv(idgn,idgnr,aa,idglda,idgip,bb,idgldb,
     &                         idginfo)
               endif
               if((lakonl(4:4).eq.'6').or.(lakonl(4:5).eq.'15')) then
                  do j=1,4
                     do k=1,nfield
                        if((c_limit(1,k).gt.bb(j,k)).and.
     &                     (c_limit(2,k).lt.bb(j,k))) then
                           field(k,ifacew1(j,jfaces))=bb(j,k)
                        endif
                        if(c_limit(1,k).lt.bb(j,k)) then
                           field(k,ifacew1(j,jfaces))=c_limit(1,k)
                        endif
                        if(c_limit(2,k).gt.bb(j,k)) then
                           field(k,ifacew1(j,jfaces))=c_limit(2,k)
                        endif   
                     enddo
                  enddo
               else
                  do j=1,4
                     do k=1,nfield
                        if((c_limit(1,k).gt.bb(j,k)).and.
     &                     (c_limit(2,k).lt.bb(j,k))) then
                           field(k,ifaceq(j,jfaces))=bb(j,k)
                        endif
                        if(c_limit(1,k).lt.bb(j,k)) then
                           field(k,ifaceq(j,jfaces))=c_limit(1,k)
                        endif
                        if(c_limit(2,k).gt.bb(j,k)) then
                           field(k,ifaceq(j,jfaces))=c_limit(2,k)
                        endif   
                     enddo
                  enddo
               endif
            endif
         elseif((lakonl(4:5).eq.'10').or.(lakonl(4:4).eq.'4').or.
     &           (((lakonl(4:5).eq.'15').or.(lakonl(4:4).eq.'6')).and.
     &           (jfaces.le.2))) then
            field(1:nfield,1:15)=0.d0
            if(lakonl(7:8).ne.'LC') then
               do k=1,3
                  do l=1,3
                     cc(k,l)=0.d0
                  enddo
                  do l=1,nfield
                     dd(k,l)=0.d0
                  enddo
               enddo
               do j=i,n
                  do k=1,nfield
                     if(j.eq.i) then
                        c_limit(1,k)=yi(k,1,j)
                        c_limit(2,k)=yi(k,1,j)
                     endif
                     if(c_limit(1,k).lt.yi(k,1,j)) then
                        c_limit(1,k)=yi(k,1,j)
                     endif
                     if(c_limit(2,k).gt.yi(k,1,j)) then
                        c_limit(2,k)=yi(k,1,j)
                     endif
                  enddo
                  indexcj=ipkon(j)
                  nopespringj=kon(indexcj)
                  igauss=kon(indexcj+nopespringj+1)
                  xi=pslavsurf(1,igauss)
                  et=pslavsurf(2,igauss)
                  if((n-i).lt.2) exit
                  if((n-i).eq.2) then
                     call shape3tri(xi,et,pl,xsj2,xs2,shp2,iflag)
                     do k=1,3
                        cc(j+1-i,k)=shp2(4,k)
                     enddo
                     do k=1,nfield
                        dd(j+1-i,k)=yi(k,1,j)
                     enddo
                  else
                     call shape3tri(xi,et,pl,xsj2,xs2,shp2,iflag)
                     do k=1,3
                        do l=1,3
                           cc(k,l)=cc(k,l)+shp2(4,k)*shp2(4,l)
                        enddo
                     enddo
                     do k=1,3
                        do l=1,nfield
                           dd(k,l)=dd(k,l)+shp2(4,k)*yi(l,1,j)
                        enddo
                     enddo
                  endif
               enddo
               if(n.eq.i) then
                  do k=1,3
                     do l=1,nfield
                        dd(k,l)=dd(k,l)+yi(l,1,j)
                     enddo
                  enddo   
               elseif((n-i).eq.1) then
                  do j=i,n
                     do k=1,3
                        do l=1,nfield
                           dd(k,l)=dd(k,l)+yi(l,1,j)*0.5d0
                        enddo
                     enddo
                  enddo
               else
                  idgn=3
                  idgnr=nfield
                  idglda=3
                  idgldb=3
                  call dgesv(idgn,idgnr,cc,idglda,idgip,dd,idgldb,
     &                         idginfo)
               endif
               if((lakonl(4:4).eq.'6').or.(lakonl(4:5).eq.'15')) then
                  do j=1,3
                     do k=1,nfield
                        if((c_limit(1,k).gt.dd(j,k)).and.
     &                     (c_limit(2,k).lt.dd(j,k))) then
                           field(k,ifacew1(j,jfaces))=dd(j,k)
                        endif
                        if(c_limit(1,k).lt.dd(j,k)) then
                           field(k,ifacew1(j,jfaces))=c_limit(1,k)
                        endif
                        if(c_limit(2,k).gt.dd(j,k)) then
                           field(k,ifacew1(j,jfaces))=c_limit(2,k)
                        endif   
                     enddo
                  enddo
               else
                  do j=1,3
                     do k=1,nfield
                        if((c_limit(1,k).gt.dd(j,k)).and.
     &                       (c_limit(2,k).lt.dd(j,k))) then
                           field(k,ifacet(j,jfaces))=dd(j,k)
                        endif
                        if(c_limit(1,k).lt.dd(j,k)) then
                           field(k,ifacet(j,jfaces))=c_limit(1,k)
                        endif
                        if(c_limit(2,k).gt.dd(j,k)) then
                           field(k,ifacet(j,jfaces))=c_limit(2,k)
                        endif   
                     enddo
                  enddo
               endif
            endif  
         endif
!     
!     determining the field values in the midside nodes of the
!     slave face
!     
         if(lakonl(4:5).eq.'20') then
            if(lakonl(7:8).ne.'LC') then
               do j=5,8
                  do k=1,nfield
                     field(k,ifaceq(j,jfaces))=
     &                   (field(k,nonei20(2,ifaceq(j,jfaces)-8))+
     &                    field(k,nonei20(3,ifaceq(j,jfaces)-8)))/2.d0
                  enddo
               enddo
            else
               nlayer=0
               do j=1,mi(3)
                  if(ielmat(j,i).gt.0) then
                     nlayer=nlayer+1
                  else
                     exit
                  endif
               enddo
               do m=1,nlayer
                  jj=20*(m-1)
                  do j=9,20
                     do k=1,nfield
                        field(k,jj+j)=(field(k,jj+nonei20(2,j-8))
     &                       +field(k,jj+nonei20(3,j-8)))/2.d0
                     enddo
                  enddo
               enddo
            endif
         elseif(lakonl(4:5).eq.'10') then
            do j=4,6
               do k=1,nfield
                  field(k,ifacet(j,jfaces))=
     &                (field(k,nonei10(2,ifacet(j,jfaces)-4))+
     &                 field(k,nonei10(3,ifacet(j,jfaces)-4)))/2.d0
               enddo
            enddo
         elseif(lakonl(4:5).eq.'15') then
            if(lakonl(7:8).ne.'LC') then
               if(jfaces.le.2) then
                  do j=4,6
                     do k=1,nfield
                        field(k,ifacew2(j,jfaces))=
     &                      (field(k,nonei15(2,ifacew2(j,jfaces)-6))+
     &                      field(k,nonei15(3,ifacew2(j,jfaces)-6)))
     &                      /2.d0
                     enddo
                  enddo
               else
                  do j=5,8
                     do k=1,nfield
                        field(k,ifacew2(j,jfaces))=
     &                      (field(k,nonei15(2,ifacew2(j,jfaces)-6))+
     &                      field(k,nonei15(3,ifacew2(j,jfaces)-6)))
     &                      /2.d0
                     enddo
                  enddo
               endif
            else
               nlayer=0
               do j=1,mi(3)
                  if(ielmat(j,i).gt.0) then
                     nlayer=nlayer+1
                  else
                     exit
                  endif
               enddo
               do m=1,nlayer
                  jj=15*(m-1)
                  do j=7,15
                     do k=1,nfield
                        field(k,jj+j)=(field(k,jj+nonei15(2,j-6))
     &                     +field(k,jj+nonei15(3,j-6)))/2.d0
                     enddo
                  enddo
               enddo
            endif
         endif
!     
!     transferring the field values into yn
!     
         if(lakonl(7:8).ne.'LC') then
            do j=1,nope
               do k=1,nfield-2
                  yn(k,kon(indexe+j))=yn(k,kon(indexe+j))+
     &                 field(k,j)
               enddo
!     
!     interchanging positions of two last rows of yn  
!     
               yn(nfield-1,kon(indexe+j))=yn(nfield-1,kon(indexe+j))+
     &              field(nfield,j)
               yn(nfield,kon(indexe+j))=yn(nfield,kon(indexe+j))+
     &              field(nfield-1,j)
               inum(kon(indexe+j))=inum(kon(indexe+j))+1
            enddo
         else
            do j=1,nope*nlayer
               do k=1,nfield
                  yn(k,kon(indexe+nopeexp+j))=
     &                 yn(k,kon(indexe+nopeexp+j))+field(k,j)
               enddo
               inum(kon(indexe+nopeexp+j))=inum(kon(indexe+nopeexp+j))+1
            enddo
         endif
!     
c     Bernhardi start
c     incompatible modes elements
         if(lakonl(1:5).eq.'C3D8I') then
            do j=1,3
               do k=1,nfield
                  yn(k,kon(indexe+nope+j))=0.0d0
               enddo
c               inum(kon(indexe+nope+j))=inum(kon(indexe+nope+j))+1
            enddo
         endif
c     Bernhardi end
!     
         i=n
!     
      enddo
!     
!     taking the mean of nodal contributions coming from different
!     elements having the node in common
!     
      do i=1,nk
         if(inum(i).gt.0) then
            do j=1,nfield
               yn(j,i)=yn(j,i)/inum(i)
            enddo
         endif
      enddo
!     
!     zeroing nonactive nodes
!  
c      do i=1,nslavnode(ntie+1)
c         if(islavact(i).ne.1) then
c            do j=1,nfield
c               yn(j,islavnode(i))=0.d0
c            enddo
c         endif
c      enddo
!     
!     for 1d and 2d elements only:
!     finding the solution in the original nodes
!     
      if((cflag.ne.' ').and.(cflag.ne.'E')) then
         call map3dto1d2d(yn,ipkon,inum,kon,lakon,nfield,nk,ne,cflag,co,
     &        vold,force,mi,ielprop,prop)
      endif
!     
      return
      end
