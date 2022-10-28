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
      subroutine basis(x,y,z,xo,yo,zo,nx,ny,nz,planfa,ifatet,nktet,
     &     netet,field,nfield,cotet,kontyp,ipkon,kon,iparent,xp,yp,zp,
     &     value,ratio,iselect,nselect,istartset,iendset,ialset,
     &     imastset,ielemnr,nterms,konl,nelem,loopa,dist)
!
!     interpolates fields from a master mesh onto a new node at an
!     arbitrary position p. The master mesh must consist of volumetric
!     elements. These are called the parent elements. They are supposed
!     to be subdivided into 1 or more linear (= 4-noded) tetrahedral
!     elements.
!     The search to which parent element position p belongs is done
!     using the tetrangulation. It permits to find the appropriate
!     tetrahedron to which p belongs and through the parent relationship
!     to which parent element      
!
!     INPUT
!      
!     x(*)            x-coordinates of the center of gravity of the
!                     tetrahedra in the tetrangulation of the master mesh,
!                     = xo, sorted in ascending order
!     y(*)            y-coordinates of the center of gravity of the
!                     tetrahedra in the tetrangulation of the master mesh,
!                     = yo, sorted in ascending order
!     z(*)            z-coordinates of the center of gravity of the
!                     tetrahedra in the tetrangulation of the master mesh,
!                     = zo, sorted in ascending order
!     xo(*)           x-coordinates of the center of gravity of the
!                     tetrahedra in the tetrangulation of the master mesh,
!     yo(*)           y-coordinates of the center of gravity of the
!                     tetrahedra in the tetrangulation of the master mesh,
!     zo(*)           z-coordinates of the center of gravity of the
!                     tetrahedra in the tetrangulation of the master mesh,
!     nx(i)           position in xo of entry i of x
!     ny(i)           position in yo of entry i of y
!     nz(i)           position in zo of entry i of z
!     planfa(4,i)     coefficients of the plane of tetrahedral face i:    
!                     planfa(1)*x+planfa(2)*y+planfa(3)*z+planfa(4)=0
!     ifatet(1..4,i)  faces of tetrahedron i
!     nktet           highest node number in the tetrahedral mesh
!     netet           number of tets in the tetrahedral mesh = highest
!                     element number
!     field(i,nktet)  nodal values of field i in the tet nodes
!     nfield          maximum number of available field in the master mesh
!     cotet(3,*)      coordinates of the tet nodes
!     kontyp(i)       type of parent element i (convention for the types as
!                     in CalculiX GraphiX)
!     ipkon(i)        entry before the start of the topology description
!                     for parent element i (= list of nodes belonging to the
!                     element) in field kon
!     kon(*)          field with the topology of all parent elements
!     iparent(i)      parent element of tetrahedral element i;      
!     xp              x-coordinate for the location in which the field are
!                     to be interpolated
!     yp              y-coordinate for the location in which the field are
!                     to be interpolated
!     zp              z-coordinate for the location in which the field are
!                     to be interpolated
!     iselect(*)      selected fields for interpolation. Must be numbers
!                     in the interval [1,nfield]
!     nselect         number of selected fields for interpolation. Must be
!                     smaller or equal to nfield
!     istartset(i)    starting position in field ialset for set i
!     iendset(i)      ending position in field ialset for set i
!     ialset(*)       field containing all set information
!     imastset        element set consisting of parent element within
!                     which the interpolation is to be done. This
!                     allows to limit the interpolation to a subset
!                     of parent elements. This may be important in
!                     case of a discontinuous master field
!     ielemnr(i)      number of parent element i as used in set imastset
!     loopa           number of loops used in attached_3d.f. The accuracy
!                     of the projection onto the parent element in local
!                     coordinates is 10**(-loopa)
!
!     OUTPUT      
!
!     value(nfield)   interpolated values (max. nfield values)
!     ratio(nterms)   coefficients for the interpolation within the parent
!                     element
!     nterms          number of coeffients in the interpolation. This is
!                     also the number of nodes belonging to the parent
!                     element
!     konl(nterms)    node numbers of the parent element ielemnr
!     nelem           number of the parent element to which p belongs.
!     dist            distance from p to the parent element to which
!                     it belongs. Ideally this is a very small number,
!                     however, if p is external to the master mesh this      
!                     number may be finite
!     xp              x-coordinate of the closest location within the
!                     master mesh
!     yp              y-coordinate of the closest location within the
!                     master mesh
!     zp              z-coordinate of the closest location within the
!                     master mesh
!      
      implicit none
!     
!     100 nearest nodes: node(100),idummy1(100),idummy2(100),iparentel(100)
!     
      integer ifatet(4,*),id,node(netet),near,nx(*),ny(*),nz(*),
     &     ifs,iface,i,konl(20),nfield,nktet,ielement,ielmax,netet,k,
     &     j,iselect(nselect),nselect,kontyp(*),ipkon(*),kon(*),
     &     nterms,indexe,nelem,konl_opt(20),idummy1(netet),
     &     iparentel(netet),nparentel,kflag,ii,inside,nterms_opt,
     &     istartset(*),iendset(*),ialset(*),imastset,ielemnr(*),
     &     ielementnr,nlength,nearset,nelem_opt,loopa,ielement2,
     &     iparent(*),idummy2(netet)
!     
      real*8 cotet(3,*),planfa(4,*),dface,dfacemax,tolerance,
     &     dist,field(nfield,nktet),ratio(20),pneigh(3,20),pnode(3),xi,
     &     et,ze,dist_opt,ratio_opt(20),xp,yp,zp,x(*),y(*),z(*),xo(*),
     &     yo(*),zo(*),value(nfield)
!     
      tolerance=1.d-6
!     
      if(imastset.ne.0) then
        nlength=iendset(imastset)-istartset(imastset)+1
      endif
!     
!     look for the global element encompassing point (xp,yp,zp) 
!     
      do ii=1,4
!     
!     three/four-level algorithm
!     - start with the search for the nearest neighboring element
!     - if this element is not the right one search for the 10
!     nearest neighbors
!     - if no fitting element was found, search for the 100 nearest
!     neighbors
!     - fourth level (only if a master element set is given)
!     if none of the 100 nearest neighbors belongs to the given
!     master element set, loop over all elements
!     
        if(ii.eq.1) then
          near=1
        elseif(ii.eq.2) then
          near=min(10,netet)
        elseif(ii.eq.3) then
          near=min(100,netet)
        else
          write(*,*) '*WARNING in basis: more than 100'
          write(*,*) '         neighbors are needed in order'
          write(*,*) '         to find a global element to'
          write(*,*) '         which the local node belongs;'
          write(*,*) '         the global element set may not'
          write(*,*) '         be correctly selected by the'
          write(*,*) '         user;'
          write(*,*) '         global coordinates of the local'
          write(*,*) '         node: ',xp,yp,zp
          near=netet
        endif
c     write(*,*) 'basis ',ii
!     
        if(ii.lt.4) then
!     
!     looking for the nearest linear tetrahedral element
!     
          call near3d(xo,yo,zo,x,y,z,nx,ny,nz,xp,yp,zp,netet,
     &         node,near)
        endif
!     
        inside=0
        nearset=0
        ielmax=0
!     
        do i=1,near
          if(ii.lt.4) then
            ielement=node(i)
          else
            ielement=i
          endif
!     
!     check whether the linear tetrahedral element belongs to the 
!     right set, if a set is specified
!     
!     ielement is a tetrahedral element number
!     iparent(ielement) is a running number from 1 to the
!     total number of elements
!     ielemnr(iparent(ielement)) is the real element number
!     (there can be gaps in the numbering)
!     
          if(imastset.ne.0) then
            ielementnr=ielemnr(iparent(ielement))
            call nident(ialset(istartset(imastset)),ielementnr,
     &           nlength,id)
            if(id.le.0) cycle
            if(ialset(istartset(imastset)+id-1).ne.ielementnr) cycle
            nearset=nearset+1
            node(nearset)=node(i)
          else
            nearset=near
          endif
!     
          dface=0.d0
          do j=1,4
!     
!     the face number in ifatet is negative, if the equation
!     of the face plane is such, that its value in the 
!     remaining node of the tetrahedron is negative
!     
            ifs=ifatet(j,ielement)
            iface=abs(ifs)
            dist=planfa(1,iface)*xp+planfa(2,iface)*yp
     &           +planfa(3,iface)*zp+planfa(4,iface)
            if(dist*ifs.lt.-1.d-10*iface) then
              dface=dface+dist*ifs/iface
            endif
          enddo
c          write(*,*) 'basis dface 1',dface,ielemnr(iparent(ielement))
          if(dface.gt.-1.d-10) then
            inside=1
!     
!     check whether the other parent elements are in the
!     master set (the fact that the node is inside a linear
!     tetrahedron does not mean that the node is within the
!     parent element; if not, all parent elements are checked;
!     therefore the other parent elements have
!     to be checked on whether they belong to the master set;
!     
            if(imastset.ne.0) then
              do j=i+1,near
                if(ii.lt.4) then
                  ielement2=node(j)
                else
                  ielement2=j
                endif
                ielementnr=ielemnr(iparent(ielement2))
                call nident(ialset(istartset(imastset)),ielementnr,
     &               nlength,id)
                if(id.le.0) cycle
                if(ialset(istartset(imastset)+id-1).ne.ielementnr)
     &               cycle
                nearset=nearset+1
                node(nearset)=node(j)
              enddo
            else
              nearset=near
            endif
!     
            exit
          endif
c     if(i.eq.1) then
          if(ielmax.eq.0) then
            dfacemax=dface
            ielmax=ielement
          else
            if(dface.gt.dfacemax) then
              ielmax=ielement
              dfacemax=dface
            endif
          endif
        enddo
!     
!     if a global element set was defined, check whether an
!     appropriate element was found
!     
        if(imastset.ne.0) then
          if(nearset.eq.0) then
            if(ii.lt.4) then
              cycle
            else
              write(*,*) '*ERROR: no suitable global element found'
              write(*,*) '        for location (',xp,yp,zp,')'
              call exit(201)
            endif
          endif
        endif
!     
!     if no linear tetrahedral element was found, the linear tetrahedral
!     element with the smallest dfacemax (in absolute value; summed distance)
!     is taken 
!     
        if(inside.eq.0) then
          ielement=ielmax
        endif
!     
        nelem=iparent(ielement)
c     write(*,*) 'basis element ',nelem,kontyp(nelem)
        if(kontyp(nelem).eq.1) then
          nterms=8
        elseif(kontyp(nelem).eq.2) then
          nterms=6
        elseif(kontyp(nelem).eq.3) then
          nterms=4
        elseif(kontyp(nelem).eq.4) then
          nterms=20
        elseif(kontyp(nelem).eq.5) then
          nterms=15
        elseif(kontyp(nelem).eq.6) then
          nterms=10
        else
          nterms=0
        endif
        indexe=ipkon(nelem)
!     
!     modify order of connectivity ccx <-> cgx
!     
        if(kontyp(nelem).eq.4) then
          do i=1,12
            konl(i)=kon(indexe+i)
          enddo
          do i=13,16
            konl(i+4)=kon(indexe+i)
          enddo
          do i=17,20
            konl(i-4)=kon(indexe+i)
          enddo
        elseif(kontyp(nelem).eq.5) then
          do i=1,9
            konl(i)=kon(indexe+i)
          enddo
          do i=10,12
            konl(i+3)=kon(indexe+i)
          enddo
          do i=13,15
            konl(i-3)=kon(indexe+i)
          enddo
        else
          do i=1,nterms
            konl(i)=kon(indexe+i)
          enddo
        endif
!     
!     nodes of master element
!     
        do i=1,nterms
          do j=1,3
            pneigh(j,i)=cotet(j,konl(i))
          enddo
        enddo
!     
!     slave node
!     
        pnode(1)=xp
        pnode(2)=yp
        pnode(3)=zp
!     
!     attaching slave node to master element
!     
        if(nterms.ne.0) then
          call attach_3d(pneigh,pnode,nterms,ratio,dist,xi,et,ze,
     &         loopa)
        endif
c        write(*,*) 'basis dist 2',dist
c     write(*,123) xi,et,ze,dist
c     123     format('basis ',4(1x,e11.4))
!     
!     checking the parent elements of the "best" tetrahedra
!     in case the distance between slave node and location of
!     interpolation exceeds "tolerance"
!     
        if(dist.gt.tolerance) then
          do i=1,nterms
            konl_opt(i)=konl(i)
            ratio_opt(i)=ratio(i)
          enddo
          nelem_opt=nelem
          nterms_opt=nterms
          dist_opt=dist
!     
!     sorting the parent elements (different linear tetrahedrals
!     may have the same parent element)
!     
          do i=1,nearset
            idummy1(i)=iparent(node(i))
          enddo
          kflag=1
          call isortii(idummy1,idummy2,nearset,kflag)
          nparentel=0
          do i=1,nearset
            if(idummy1(i).eq.nelem) cycle
            if(i.eq.1) then
              iparentel(1)=idummy1(1)
              nparentel=1
            else
              if(idummy1(i).ne.idummy1(i-1)) then
                nparentel=nparentel+1
                iparentel(nparentel)=idummy1(i)
              endif
            endif
          enddo
!     
          do k=1,nparentel
!     
!     slave node
!     
            pnode(1)=xp
            pnode(2)=yp
            pnode(3)=zp
            nelem=iparentel(k)
            if(kontyp(nelem).eq.1) then
              nterms=8
            elseif(kontyp(nelem).eq.2) then
              nterms=6
            elseif(kontyp(nelem).eq.3) then
              nterms=4
            elseif(kontyp(nelem).eq.4) then
              nterms=20
            elseif(kontyp(nelem).eq.5) then
              nterms=15
            elseif(kontyp(nelem).eq.6) then
              nterms=10
            else
              nterms=0
            endif
            indexe=ipkon(nelem)
!     
!     modify order of connectivity ccx <-> cgx
!     
            if(kontyp(nelem).eq.4) then
              do i=1,12
                konl(i)=kon(indexe+i)
              enddo
              do i=13,16
                konl(i+4)=kon(indexe+i)
              enddo
              do i=17,20
                konl(i-4)=kon(indexe+i)
              enddo
            elseif(kontyp(nelem).eq.5) then
              do i=1,9
                konl(i)=kon(indexe+i)
              enddo
              do i=10,12
                konl(i+3)=kon(indexe+i)
              enddo
              do i=13,15
                konl(i-3)=kon(indexe+i)
              enddo
            else
              do i=1,nterms
                konl(i)=kon(indexe+i)
              enddo
            endif
!     
!     nodes of master element
!     
            do i=1,nterms
              do j=1,3
                pneigh(j,i)=cotet(j,konl(i))
              enddo
            enddo
!     
!     attaching slave node to master element
!     
            if(nterms.ne.0) then
              call attach_3d(pneigh,pnode,nterms,ratio,dist,
     &             xi,et,ze,loopa)
            endif
c      write(*,*) 'basis ii 3',ii,nelem,dist
!     
!     check whether the present element yields better results
!     
            if(dist.lt.dist_opt) then
c     write(*,123) xi,et,ze,dist
              do i=1,nterms
                konl_opt(i)=konl(i)
                ratio_opt(i)=ratio(i)
              enddo
              nelem_opt=nelem
              nterms_opt=nterms
              dist_opt=dist
            endif
            if(dist.lt.tolerance) exit
          enddo
!     
!     storing the optimal configuration
!     
          nelem=nelem_opt
          nterms=nterms_opt
c     
          dist=dist_opt
c     
          do i=1,nterms
            konl(i)=konl_opt(i)
            ratio(i)=ratio_opt(i)
          enddo
        endif
!     
        if((ii.eq.3).or.(dist.lt.tolerance)) exit
!     
      enddo
!
!     calculate the coordinates of the closest location
!
      xp=0.d0
      yp=0.d0
      zp=0.d0
      do j=1,nterms
        xp=xp+ratio(j)*cotet(1,konl(j))
        yp=yp+ratio(j)*cotet(2,konl(j))
        zp=zp+ratio(j)*cotet(3,konl(j))
      enddo
!     
!     interpolating the fields
!     
c      write(*,*) 'basis end 4',dist,tolerance,near
      do k=1,nselect
        i=iselect(k)
        value(k)=0.d0
        do j=1,nterms
!     
!     due to round-off the values can be slightly outside the
!     (0,1) range           
!     
c          ratio(j)=max(0.d0,ratio(j))
c          ratio(j)=min(1.d0,ratio(j))
!          
          value(k)=value(k)+ratio(j)*field(i,konl(j))
c     write(*,*) 'basis j',j,konl(j),field(i,konl(j))
        enddo
c     write(*,*) 'basis select',k,i,value(k)
      enddo
!     
      return
      end
      
