!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine rearrangecfd(ne,ipkon,lakon,ielmat,ielorien,norien,
     &     nef,ipkonf,lakonf,ielmatf,ielorienf,mi,nelold,nelnew,nkold,
     &     nknew,nk,nkf,konf,nkonf,nmpc,ipompc,nodempc,coefmpc,memmpc_,
     &     nmpcf,ipompcf,nodempcf,coefmpcf,memmpcf,nboun,nodeboun,
     &     ndirboun,xboun,nbounf,nodebounf,ndirbounf,xbounf,nload,
     &     nelemload,sideload,xload,nloadf,nelemloadf,sideloadf,
     &     xloadf,ipobody,ipobodyf,kon,nkftot,co,cof,vold,voldf,
     &     ikbounf,ilbounf,ikmpcf,ilmpcf,iambounf,iamloadf,iamboun,
     &     iamload,xbounold,xbounoldf,xbounact,xbounactf,xloadold,
     &     xloadoldf,xloadact,xloadactf,inotr,inotrf,nam,ntrans,
     &     nbody)
!
!     renumbering nodes and elements for fluids such that no
!     gaps result
!     
      implicit none
!
      logical fluid
!     
      character*8 lakon(*),lakonf(*)
      character*20 sideload(*),sideloadf(*)
!     
      integer ne,ipkon(*),mi(*),ielmat(mi(3),*),ielorien(mi(3),*),
     &     norien,nef,ipkonf(*),ielmatf(mi(3),*),ielorienf(mi(3),*),
     &     nelold(*),nelnew(*),nkold(*),nknew(*),nk,nkf,i,j,k,indexe,
     &     nope,node,konf(*),nkonf,nmpc,ipompc(*),nodempc(3,*),index,
     &     memmpc_,nmpcf,ipompcf(*),nodempcf(3,*),memmpcf,nboun,nam,
     &     nodeboun(*),ndirboun(*),nbounf,nodebounf(*),ndirbounf(*),
     &     nload,nelemload(2,*),nloadf,nelemloadf(2,*),ipobody(2,*),
     &     ipobodyf(2,*),kon(*),nelem,nkftot,kflag,ikbounf(*),
     &     ikmpcf(*),ilmpcf(*),ndir,iambounf(*),iamloadf(2,*),ntrans,
     &     iamboun(*),iamload(2,*),inotr(2,*),inotrf(2,*),nbody,
     &     ilbounf(*)
!     
      real*8 coefmpc(*),coefmpcf(*),xboun(*),xbounf(*),xload(2,*),
     &     xloadf(2,*),co(3,*),cof(3,*),vold(0:mi(2),*),
     &     voldf(0:mi(2),*),xbounold(*),xbounoldf(*),xbounact(*),
     &     xbounactf(*),xloadold(2,*),xloadoldf(2,*),xloadact(2,*),
     &     xloadactf(2,*)
!
      kflag=2
!
!     rearranging the elements
!
      nef=0
      do i=1,ne
        if(ipkon(i).lt.0) cycle
        if(lakon(i)(1:1).ne.'F') cycle
        nef=nef+1
        nelold(nef)=i
        nelnew(i)=nef
        ipkonf(nef)=ipkon(i)
        lakonf(nef)=lakon(i)
        do j=1,mi(3)
          ielmatf(j,nef)=ielmat(j,i)
        enddo
        if(norien.gt.0) then
          do j=1,mi(3)
            ielorienf(j,nef)=ielorien(j,i)
          enddo
        endif
      enddo
!
!     rearranging the nodes belonging to fluid elements
!     setting nknew to 1 for all used fluid nodes
!
      do i=1,nef
        nope=ichar(lakonf(i)(4:4))-48
        indexe=ipkonf(i)
        do j=1,nope
          node=kon(indexe+j)
          nknew(node)=1
        enddo
      enddo
!
      nkf=0
      do i=1,nk
        if(nknew(i).eq.1) then
          nkf=nkf+1
          nknew(i)=nkf
          nkold(nkf)=i
          do j=1,3
            cof(j,nkf)=co(j,i)
          enddo
          do j=0,mi(2)
            voldf(j,nkf)=vold(j,i)
          enddo
        endif
      enddo
      nkftot=nkf
!
!     setting up konf
!
      nkonf=0
      do i=1,nef
        nope=ichar(lakonf(i)(4:4))-48
        indexe=ipkonf(i)
        ipkonf(i)=nkonf
        do j=1,nope
          nkonf=nkonf+1
          node=kon(indexe+j)
          konf(nkonf)=nknew(node)
        enddo
      enddo
!
!     identifying the fluid MPC's
!     adapting the node numbers
!
      nmpcf=0
      memmpcf=0
      do i=1,nmpc
        index=ipompc(i)
        fluid=.false.
        do
          if(index.eq.0) exit
          node=nodempc(1,index)
          if(nknew(node).gt.0) then
            fluid=.true.
            exit
          endif
          index=nodempc(3,index)
        enddo
!
        if(fluid) then
          index=ipompc(i)
          nmpcf=nmpcf+1
          ipompcf(nmpcf)=memmpcf+1
          do
            memmpcf=memmpcf+1
            node=nodempc(1,index)
            if(nknew(node).eq.0) then
              nkftot=nkftot+1
              nknew(node)=nkftot
              nkold(nkftot)=node
              do j=1,3
                cof(j,nkftot)=co(j,node)
              enddo
              do j=0,mi(2)
                voldf(j,nkftot)=vold(j,node)
              enddo
            endif
            nodempcf(1,memmpcf)=nknew(node)
            nodempcf(2,memmpcf)=nodempc(2,index)
            coefmpcf(memmpcf)=coefmpc(index)
            index=nodempc(3,index)
            if(index.eq.0) then
              nodempcf(3,memmpcf)=0
              exit
            else
              nodempcf(3,memmpcf)=memmpcf+1
            endif
          enddo
          index=ipompcf(nmpcf)
          node=nodempcf(1,index)
          ndir=nodempcf(2,index)
          ikmpcf(nmpcf)=8*(node-1)+ndir
          ilmpcf(nmpcf)=nmpcf
        endif
      enddo
      call isortii(ikmpcf,ilmpcf,nmpcf,kflag)
!
!     identifying the fluid SPC's
!     adapting the node numbers
!
      nbounf=0
      do i=1,nboun
        node=nodeboun(i)
        if(nknew(node).ne.0) then
          nbounf=nbounf+1
          nodebounf(nbounf)=nknew(node)
          ndirbounf(nbounf)=ndirboun(i)
          if(nam.gt.0) iambounf(nbounf)=iamboun(i)
          ikbounf(nbounf)=8*(nodebounf(nbounf)-1)+ndirbounf(nbounf)
          ilbounf(nbounf)=nbounf
          xbounf(nbounf)=xboun(i)
          xbounoldf(nbounf)=xbounold(i)
          xbounactf(nbounf)=xbounact(i)
        endif
      enddo
      call isortii(ikbounf,ilbounf,nbounf,kflag)
!
!     rearranging distributed load
!
      nloadf=0
      do i=1,nload
        nelem=nelemload(1,i)
        if(nelnew(nelem).eq.0) cycle
        nloadf=nloadf+1
        nelemloadf(1,nloadf)=nelnew(nelem)
        node=nelemload(2,i)
        if(node.gt.0) then
          if(nknew(node).eq.0) then
            nkftot=nkftot+1
            nknew(node)=nkftot
            nkold(nkftot)=node
            do j=1,3
              cof(j,nkftot)=co(j,node)
            enddo
            do j=0,mi(2)
              voldf(j,nkftot)=vold(j,node)
            enddo
          endif
          nelemloadf(2,nloadf)=nknew(node)
        endif
        sideloadf(nloadf)=sideload(i)
        if(nam.gt.0) then
          iamloadf(1,nloadf)=iamload(1,i)
          iamloadf(2,nloadf)=iamload(2,i)
        endif
        xloadf(1,nloadf)=xload(1,i)
        xloadf(2,nloadf)=xload(2,i)
        xloadoldf(1,nloadf)=xloadold(1,i)
        xloadoldf(2,nloadf)=xloadold(2,i)
        xloadactf(1,nloadf)=xloadact(1,i)
        xloadactf(2,nloadf)=xloadact(2,i)
      enddo
!
!     transformations
!
      if(ntrans.gt.0) then
        do i=1,nkftot
          inotrf(1,i)=inotr(1,nkold(i))
          if(inotrf(2,nkold(i)).eq.0) then
            inotrf(2,i)=inotr(2,nkold(i))
          else
            inotrf(2,i)=nknew(inotr(2,nkold(i)))
          endif
        enddo
      endif
!
!     rearranging body loads
!
      if(nbody.gt.0) then
        do i=1,nef
          do j=1,2
            ipobodyf(j,i)=ipobody(j,nelold(i))
          enddo
        enddo
      endif
!     
      return
      end
