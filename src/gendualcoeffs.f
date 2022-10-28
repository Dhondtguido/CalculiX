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
!     
!     Determining the coefficients of the dual shape functions on
!     the slave surface
!     see phd thesis Sitzmann Chapter 3.3. (p.34) and Chapter 4.3.
!     
!     
!     [in]     islavact  (i) indicates, if slave node i is active (=-3 no-slave-node, =-2 no-LM-node, =-1 no-gap-node, =0 inactive node, =1 sticky node, =2 slipping/active node)
!     [in]     islavsurf  islavsurf(1,i) slaveface i islavsurf(2,i) pointer into imastsurf and pmastsurf
!     [in]     itiefac    pointer into field islavsurf: (1,i) beginning slave_i (2,i) end of slave_i
!     [in]     islavnode  fields containing nodes of slace surfaces
!     [in]     nslavnode  (i) for contraint i pointer into field islavnode
!     [in]     pslavsurf  integration points and weights on slave side 
!     [in] pslavdual (:,i)coefficients \f$ \alpha_{ij}\f$, \f$ 1,j=1,..8\f$ for dual shape functions for face i
!     [in] pslavdualpg (:,i)coefficients \f$ \alpha_{ij}\f$, \f$ 1,j=1,..8\f$ for Petrov-Galerkin shape functions for face i
!     [in]  iflagdualquad   flag indicating what mortar contact is used (=1 quad-lin, =2 quad-quad, =3 PG quad-lin, =4 PG quad-quad)
!     
      subroutine gendualcoeffs(tieset,ntie,ipkon,kon,lakon,co,vold,
     &     islavact,islavsurf,itiefac,islavnode,nslavnode,
     &     mi,pslavsurf,pslavdual,pslavdualpg,iflagdualquad)
!     
!     Determining the coefficients of the dual shape functions on
!     the slave surface
!     
!     Author: Sitzmann,Saskia ;
!     
      implicit none
!     
      logical checkbiorthogonality, checknorm
!     
      character*8 lakon(*)
      character*81 tieset(3,*),slavset
!     
      integer ntie,ipkon(*),kon(*),iflag,kneigh,i,ii,j,k,l,
     &     ipos,islavact(*),nope,islavsurf(2,*),islavnode(*),
     &     nslavnode(ntie+1),itiefac(2,*),ifaces,nelems,jfaces,mi(*),
     &     mint2d,m,nopes,konl(20),id,nopes2,ipiv(8),info,ipnt,ifac,
     &     getlocno,lnode(2,8),nnogap,n1,n2,iscontr(8*8),
     &     imcontr(8*8),jj,locm,locs,nodesf,nodem,ifs,
     &     ifm,ns,indexf,idummy,modf,iflagdualquad
!     
      real*8 co(3,*),vold(0:mi(2),*),ets,xis,weight,xl2s(3,8),xsj2(3),
     &     shp2(7,8),shp2s(7,8),dx,help,xs2(3,7),
     &     pslavdual(64,*),diag_els(8),m_els(36),contribution,work(8),
     &     contr(8*8),xs2m(3,7),xsj2m(3),shp2m(7,8),etm,xim,
     &     pslavsurf(3,*),pslavdualpg(64,*)
!     
!     
!     
      data iflag /2/
!     
      checkbiorthogonality=.false.
      checknorm=.false.
!     
!     loop over all ties
!     
      do i=1,ntie
        if(tieset(1,i)(81:81).ne.'C') cycle
        kneigh=1     
        slavset=tieset(2,i)
        ipos=index(slavset,' ')
        if(slavset(ipos-1:ipos-1).eq.'S') cycle
!     
!     ntri: number of triangles in the triangulation of the master
!     surface corresponding to tie i
!     
        do l=itiefac(1,i),itiefac(2,i)
          if(checkbiorthogonality)then
            write(*,*)'face',l
          endif
          ifaces=islavsurf(1,l)
          nelems=int(ifaces/10)
          jfaces=ifaces-nelems*10
          call getnumberofnodes(nelems,jfaces,lakon,nope,
     &         nopes,idummy) 
!     
!     initialization for Dualshape Coefficient matrix
!     
          do j=1,8
            diag_els(j)=0.0
            do k=1,j
              m_els(k+((j-1)*j/2))=0.0
            enddo
          enddo
!     
!     actual position of the nodes belonging to the
!     slave surface
!     
          do j=1,nope
            konl(j)=kon(ipkon(nelems)+j)
          enddo
          nnogap=0
!     
!     check for nogap or nolm Nodes
!     here redistribution of Lagrange Multiplier contribution
!     has to be applied
!     see Phd thesis Sitzmann Chapter 4.3.
!     
          do m=1,nopes
            ifac=getlocno(m,jfaces,nope)
            lnode(1,m)=konl(ifac)
            call nident(islavnode(nslavnode(i)+1),
     &           konl(ifac),(nslavnode(i+1)-nslavnode(i)),id)
            if(islavnode(nslavnode(i)+id).eq.konl(ifac)) then
              lnode(2,m)=islavact(nslavnode(i)+id)
!     
!     for quad-lin mortar method only corner nodes carry
!     Lagrange Multiplier contribution
!     
              if((nopes.eq.8).and. 
     &             ((iflagdualquad.eq.1).or.(iflagdualquad.eq.3)))
     &             then
                if((lnode(2,m).lt.0).and.(m.le.4))then
                  nnogap=nnogap+1
                endif
                if(m.gt.4) then
                  islavact(nslavnode(i)+id)=-2
                  lnode(2,m)=islavact(nslavnode(i)+id)
                endif
              else if((nopes.eq.6).and. 
     &               ((iflagdualquad.eq.1).or.(iflagdualquad.eq.3)))then
                if((lnode(2,m).lt.0).and.(m.le.3)) then
                  nnogap=nnogap+1
                endif
                if(m.gt.3) then
                  islavact(nslavnode(i)+id)=-2
                  lnode(2,m)=islavact(nslavnode(i)+id)
                endif
              else
                if(lnode(2,m).lt.0) nnogap=nnogap+1
              endif             
              if(checknorm)then
                write(*,*) 'node',lnode(1,m),lnode(2,m)
              endif
            else
              write(*,*)'createbd: node',konl(ifac)
              write(*,*)'was not catalogued properly in islavnode'
              call exit(201)
            endif  
            
            do j=1,3
              xl2s(j,m)=co(j,konl(ifac))+
     &             vold(j,konl(ifac))
            enddo
          enddo
          if(((iflagdualquad.eq.1).or.(iflagdualquad.eq.3)))then
            if(nopes.eq.6)then
              nopes2=3
            else if(nopes.eq.8)then
              nopes2=4
            else
              nopes2=nopes
            endif
          else
            nopes2=nopes
          endif
!     
          mint2d=islavsurf(2,l+1)-islavsurf(2,l)
          if(mint2d.eq.0) cycle
          indexf=islavsurf(2,l)
!     
!     Calculate the Mass matrix for compilation of the dualshapefunction 
!     pslavdual(16,*)
!     create Z, W (see Phd thesis Sitzmann equation (3.34))
!     
!     loop over integration points
!     
          do m=1,mint2d
            xis=pslavsurf(1,indexf+m)
            ets=pslavsurf(2,indexf+m)
            weight=pslavsurf(3,indexf+m)
            ns=l
            iflag=2
!     
!     construct dual basis functions for transformed standart
!     basis functions
!     
            if(nopes.eq.8) then
              if((iflagdualquad.eq.2).or.(iflagdualquad.eq.4))then
                call shape8qtilde(xis,ets,xl2s,xsj2,xs2,shp2,iflag)
              else
                call shape8qtilde_lin(xis,ets,xl2s,xsj2,xs2,
     &               shp2,iflag)
              endif
            elseif(nopes.eq.4) then
              call shape4q(xis,ets,xl2s,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.6) then
              if((iflagdualquad.eq.2).or.(iflagdualquad.eq.4))then
                call shape6tritilde(xis,ets,xl2s,xsj2,xs2,shp2,
     &               iflag)
              else
                call shape6tritilde_lin(xis,ets,xl2s,xsj2,xs2,
     &               shp2,iflag)
              endif
            else
              call shape3tri(xis,ets,xl2s,xsj2,xs2,shp2,iflag)
            endif 
            dx=dsqrt(xsj2(1)**2+xsj2(2)**2+
     &           xsj2(3)**2)
            ipnt=0
!     
!     column
!     
            do j=1,nopes2
              diag_els(j)=diag_els(j)+shp2(4,j)*weight
     &             *dx
              do k=1,j
                m_els(k+((j-1)*j/2))=m_els(k+((j-1)*j/2))
     &               +shp2(4,k)*shp2(4,j)
     &               *weight
     &               *dx  
              enddo
            enddo
          enddo
!     
!     compute inverse of me_ls
!     factorisation
!     
          if(checknorm .or.checkbiorthogonality)then
            write(*,*)'diag_els',l
            write(*,105)(diag_els(j),j=1,nopes2)
            write(*,*)'m_els'
            do j=1,nopes2
              write(*,105) (m_els(k+((j-1)*j/2)),k=1,j)
            enddo               
          endif
          call dsptrf('U',nopes2,m_els,ipiv,info)
!     
!     inverse
!     
          call dsptri('U',nopes2,m_els,ipiv,work,info)
!     
!     stack of pslavdual multiplication with diag_els
!     A=W*Z^-1
!     
          do k=1,64
            pslavdual(k,l)=0.0
            if(iflagdualquad.gt.2)then
              pslavdualpg(k,l)=0.0
            endif
          enddo
          do k=1,nopes2
            if(iflagdualquad.gt.2)then
              pslavdualpg((k-1)*8+k,l)=1.0
            endif
            do j=1,nopes2
              if(k.le.j)then
                pslavdual((k-1)*8+j,l)=diag_els(k)
     &               *m_els(k+((j-1)*j/2))
              else
                pslavdual((k-1)*8+j,l)=diag_els(k)
     &               *m_els(j+((k-1)*k/2))
              endif
            enddo
          enddo 
!     
          if(checknorm)then
            write(*,*)'pslavdual1,face',l
            do ii=1,nopes2
              help=0.0
              do j=1,nopes2
                help=help+pslavdual((ii-1)*8+j,l)
              enddo
              write(*,105) (pslavdual((ii-1)*8+j,l),j=1,nopes2)
            enddo
          endif
!     
!     remove Lagrange Multiplier contribution for nogap and noLM nodes
!     phd thesis Sitzmann Chapter 4.3.
!     
          if((nnogap.gt.0).and.(nopes2.eq.4)) then
            if(nnogap.eq.1) then
              do ii=1,nopes2
                if(lnode(2,ii).lt.0) exit
              enddo
              n1=modf(nopes2,ii-1)
              n2=modf(nopes2,ii+1)
              do jj=1,nopes2
                pslavdual((n1-1)*8+jj,l)=pslavdual((n1-1)*8+jj,l)
     &               +0.5*pslavdual((ii-1)*8+jj,l)
                pslavdual((n2-1)*8+jj,l)=pslavdual((n2-1)*8+jj,l)
     &               +0.5*pslavdual((ii-1)*8+jj,l)
                if(iflagdualquad.gt.2)then
                  pslavdualpg((n1-1)*8+jj,l)=
     &                 pslavdualpg((n1-1)*8+jj,l)
     &                 +0.5*pslavdualpg((ii-1)*8+jj,l)
                  pslavdualpg((n2-1)*8+jj,l)=
     &                 pslavdualpg((n2-1)*8+jj,l)
     &                 +0.5*pslavdualpg((ii-1)*8+jj,l)
                endif
              enddo
              do jj=1,nopes2
                pslavdual((ii-1)*8+jj,l)=0.0
                if(iflagdualquad.gt.2)then
                  pslavdualpg((ii-1)*8+jj,l)=0.0
                endif
              enddo
            elseif(nnogap.eq.2)then
              do ii=1,nopes2
                if(lnode(2,ii).lt.0) then
                  n1=modf(nopes2,ii-1)
                  n2=modf(nopes2,ii+1)
                  if(lnode(2,n1).lt.0) n1=n2
                  if(lnode(2,n2).lt.0) n2=n1
                  do jj=1,nopes2
                    pslavdual((n1-1)*8+jj,l)=
     &                   pslavdual((n1-1)*8+jj,l)
     &                   +0.5*pslavdual((ii-1)*8+jj,l)
                    pslavdual((n2-1)*8+jj,l)=
     &                   pslavdual((n2-1)*8+jj,l)
     &                   +0.5*pslavdual((ii-1)*8+jj,l)
                    if(iflagdualquad.gt.2)then
                      pslavdualpg((n1-1)*8+jj,l)=
     &                     pslavdualpg((n1-1)*8+jj,l)
     &                     +0.5*pslavdualpg((ii-1)*8+jj,l)
                      pslavdualpg((n2-1)*8+jj,l)=
     &                     pslavdualpg((n2-1)*8+jj,l)
     &                     +0.5*pslavdualpg((ii-1)*8+jj,l)
                    endif
                  enddo
                  do jj=1,nopes2
                    pslavdual((ii-1)*8+jj,l)=0.0
                    if(iflagdualquad.gt.2)then
                      pslavdualpg((ii-1)*8+jj,l)=0.0
                    endif
                  enddo
                endif
              enddo
            elseif(nnogap.eq.3)then
              do ii=1,nopes2
                if(lnode(2,ii).gt.-1) n1=ii
              enddo
              do ii=1,nopes2
                if(lnode(2,ii).lt.0) then
                  do jj=1,nopes2
                    pslavdual((n1-1)*8+jj,l)=
     &                   pslavdual((n1-1)*8+jj,l)
     &                   +pslavdual((ii-1)*8+jj,l)
                    if(iflagdualquad.gt.2)then
                      pslavdualpg((n1-1)*8+jj,l)=
     &                     pslavdualpg((n1-1)*8+jj,l)
     &                     +pslavdualpg((ii-1)*8+jj,l)
                    endif 
                  enddo
                  do jj=1,nopes2
                    pslavdual((ii-1)*8+jj,l)=0.0
                    if(iflagdualquad.gt.2)then
                      pslavdualpg((ii-1)*8+jj,l)=0.0
                    endif
                  enddo
                endif
              enddo
            endif
          endif
          if((nnogap.gt.0).and.(nopes2.eq.3)) then
            if(nnogap.eq.1) then
              do ii=1,nopes2
                if(lnode(2,ii).lt.0) exit
              enddo
              n1=modf(nopes2,ii-1)
              n2=modf(nopes2,ii+1)
              do jj=1,nopes2
                pslavdual((n1-1)*8+jj,l)=pslavdual((n1-1)*8+jj,l)
     &               +0.5*pslavdual((ii-1)*8+jj,l)
                pslavdual((n2-1)*8+jj,l)=pslavdual((n2-1)*8+jj,l)
     &               +0.5*pslavdual((ii-1)*8+jj,l)
                if(iflagdualquad.gt.2)then
                  pslavdualpg((n1-1)*8+jj,l)=
     &                 pslavdualpg((n1-1)*8+jj,l)
     &                 +0.5*pslavdualpg((ii-1)*8+jj,l)
                  pslavdualpg((n2-1)*8+jj,l)=
     &                 pslavdualpg((n2-1)*8+jj,l)
     &                 +0.5*pslavdualpg((ii-1)*8+jj,l)
                endif
              enddo
              do jj=1,nopes2
                pslavdual((ii-1)*8+jj,l)=0.0
                if(iflagdualquad.gt.2)then
                  pslavdualpg((ii-1)*8+jj,l)=0.0
                endif
              enddo
            elseif(nnogap.eq.2)then
              do ii=1,nopes2
                if(lnode(2,ii).gt.-1) n1=ii
              enddo
              do ii=1,nopes2
                if(lnode(2,ii).lt.0) then
                  do jj=1,nopes2
                    pslavdual((n1-1)*8+jj,l)=
     &                   pslavdual((n1-1)*8+jj,l)
     &                   +pslavdual((ii-1)*8+jj,l)
                    if(iflagdualquad.gt.2)then
                      pslavdualpg((n1-1)*8+jj,l)=
     &                     pslavdualpg((n1-1)*8+jj,l)
     &                     +pslavdualpg((ii-1)*8+jj,l)
                    endif
                  enddo
                  do jj=1,nopes2
                    pslavdual((ii-1)*8+jj,l)=0.0
                    if(iflagdualquad.gt.2)then
                      pslavdualpg((ii-1)*8+jj,l)=0.0
                    endif
                  enddo
                endif
              enddo                  
            endif
          elseif((nnogap.gt.0).and.(nopes2.eq.8))then
            if(nnogap.eq.1)then
              do ii=1,nopes2
                if(lnode(2,ii).lt.0) exit
              enddo
              if(ii.le.4)then
                n1=modf(4,ii-1)+4
                n2=ii+4
              else
                n1=ii-4
                n2=modf(4,ii-4+1)
              endif
              do jj=1,nopes2
                pslavdual((n1-1)*8+jj,l)=pslavdual((n1-1)*8+jj,l)
     &               +0.5*pslavdual((ii-1)*8+jj,l)
                pslavdual((n2-1)*8+jj,l)=pslavdual((n2-1)*8+jj,l)
     &               +0.5*pslavdual((ii-1)*8+jj,l)
                if(iflagdualquad.gt.2)then
                  pslavdualpg((n1-1)*8+jj,l)=
     &                 pslavdualpg((n1-1)*8+jj,l)
     &                 +0.5*pslavdualpg((ii-1)*8+jj,l)
                  pslavdualpg((n2-1)*8+jj,l)=
     &                 pslavdualpg((n2-1)*8+jj,l)
     &                 +0.5*pslavdualpg((ii-1)*8+jj,l)
                endif
              enddo
              do jj=1,nopes2
                pslavdual((ii-1)*8+jj,l)=0.0
                if(iflagdualquad.gt.2)then
                  pslavdualpg((ii-1)*8+jj,l)=0.0
                endif
              enddo                  
            elseif((nnogap.gt.1).and.(nnogap.lt.7))then
              do ii=1,nopes2
                if(lnode(2,ii).lt.0) then                        
                  if(ii.le.4)then
                    n1=modf(4,ii-1)+4
                    n2=ii+4
                  else
                    n1=ii-4
                    n2=modf(4,ii-4+1)
                  endif
 10               if((lnode(2,n1).lt.0).and.(lnode(2,n2).lt.0))then
                    if(n1.le.4)then
                      n1=modf(4,n1-1)+4
                    else
                      n1=n1-4
                    endif
                    if(n2.le.4)then
                      n2=n2+4
                    else
                      n2=modf(4,n2-4+1)
                    endif
                    go to 10
                  else
                    if(lnode(2,n1).lt.0) n1=n2
                    if(lnode(2,n2).lt.0) n2=n1
                    do jj=1,nopes2
                      pslavdual((n1-1)*8+jj,l)=
     &                     pslavdual((n1-1)*8+jj,l)
     &                     +0.5*pslavdual((ii-1)*8+jj,l)
                      pslavdual((n2-1)*8+jj,l)=
     &                     pslavdual((n2-1)*8+jj,l)
     &                     +0.5*pslavdual((ii-1)*8+jj,l)
                      if(iflagdualquad.gt.2)then
                        pslavdualpg((n1-1)*8+jj,l)=
     &                       pslavdualpg((n1-1)*8+jj,l)
     &                       +0.5*pslavdualpg((ii-1)*8+jj,l)
                        pslavdualpg((n2-1)*8+jj,l)=
     &                       pslavdualpg((n2-1)*8+jj,l)
     &                       +0.5*pslavdualpg((ii-1)*8+jj,l)
                      endif
                    enddo
                  endif
                  do jj=1,nopes2
                    pslavdual((ii-1)*8+jj,l)=0.0
                    if(iflagdualquad.gt.2)then
                      pslavdualpg((ii-1)*8+jj,l)=0.0
                    endif
                  enddo
                endif
              enddo
            elseif(nnogap.eq.7)then
              do ii=1,nopes2
                if(lnode(2,ii).gt.-1) n1=ii
              enddo
              do ii=1,nopes2
                if(lnode(2,ii).lt.0) then
                  do jj=1,nopes2
                    pslavdual((n1-1)*8+jj,l)=
     &                   pslavdual((n1-1)*8+jj,l)
     &                   +pslavdual((ii-1)*8+jj,l)
                    if(iflagdualquad.gt.2)then
                      pslavdualpg((n1-1)*8+jj,l)=
     &                     pslavdualpg((n1-1)*8+jj,l)
     &                     +pslavdualpg((ii-1)*8+jj,l)
                    endif
                  enddo
                  do jj=1,nopes2
                    pslavdual((ii-1)*8+jj,l)=0.0
                    if(iflagdualquad.gt.2)then
                      pslavdualpg((ii-1)*8+jj,l)=0.0
                    endif
                  enddo
                endif
              enddo                  
            endif
          elseif((nnogap.gt.0).and.(nopes2.eq.6))then
            if(nnogap.eq.1)then
              do ii=1,nopes2
                if(lnode(2,ii).lt.0) exit
              enddo
              if(ii.le.3)then
                n1=modf(3,ii-1)+3
                n2=ii+3
              else
                n1=ii-3
                n2=modf(3,ii-3+1)
              endif
              do jj=1,nopes2
                pslavdual((n1-1)*8+jj,l)=pslavdual((n1-1)*8+jj,l)
     &               +0.5*pslavdual((ii-1)*8+jj,l)
                pslavdual((n2-1)*8+jj,l)=pslavdual((n2-1)*8+jj,l)
     &               +0.5*pslavdual((ii-1)*8+jj,l)
                if(iflagdualquad.gt.2)then
                  pslavdualpg((n1-1)*8+jj,l)=
     &                 pslavdualpg((n1-1)*8+jj,l)
     &                 +0.5*pslavdualpg((ii-1)*8+jj,l)
                  pslavdualpg((n2-1)*8+jj,l)=
     &                 pslavdualpg((n2-1)*8+jj,l)
     &                 +0.5*pslavdualpg((ii-1)*8+jj,l)
                endif
              enddo
              do jj=1,nopes2
                pslavdual((ii-1)*8+jj,l)=0.0
                if(iflagdualquad.gt.2)then
                  pslavdualpg((ii-1)*8+jj,l)=0.0
                endif
              enddo                   
            elseif((nnogap.gt.1).and.(nnogap.lt.5))then
              do ii=1,nopes2
                if(lnode(2,ii).lt.0) then                        
                  if(ii.le.3)then
                    n1=modf(3,ii-1)+3
                    n2=ii+3
                  else
                    n1=ii-3
                    n2=modf(3,ii-3+1)
                  endif
 20               if((lnode(2,n1).lt.0).and.(lnode(2,n2).lt.0))then
                    if(n1.le.3)then
                      n1=modf(3,n1-1)+3
                    else
                      n1=n1-3
                    endif
                    if(n2.le.3)then
                      n2=n2+3
                    else
                      n2=modf(3,n2-3+1)
                    endif
                    go to 20
                  else
                    if(lnode(2,n1).lt.0) n1=n2
                    if(lnode(2,n2).lt.0) n2=n1
                    do jj=1,nopes2
                      pslavdual((n1-1)*8+jj,l)=
     &                     pslavdual((n1-1)*8+jj,l)
     &                     +0.5*pslavdual((ii-1)*8+jj,l)
                      pslavdual((n2-1)*8+jj,l)=
     &                     pslavdual((n2-1)*8+jj,l)
     &                     +0.5*pslavdual((ii-1)*8+jj,l)
                      if(iflagdualquad.gt.2)then
                        pslavdualpg((n1-1)*8+jj,l)=
     &                       pslavdualpg((n1-1)*8+jj,l)
     &                       +0.5*pslavdualpg((ii-1)*8+jj,l)
                        pslavdualpg((n2-1)*8+jj,l)=
     &                       pslavdualpg((n2-1)*8+jj,l)
     &                       +0.5*pslavdualpg((ii-1)*8+jj,l)
                      endif
                    enddo
                  endif
                  do jj=1,nopes2
                    pslavdual((ii-1)*8+jj,l)=0.0
                    if(iflagdualquad.gt.2)then
                      pslavdualpg((ii-1)*8+jj,l)=0.0
                    endif
                  enddo
                endif
              enddo
            elseif(nnogap.eq.5)then
              do ii=1,nopes2
                if(lnode(2,ii).gt.-1) n1=ii
              enddo
              do ii=1,nopes2
                if(lnode(2,ii).lt.0) then
                  do jj=1,nopes2
                    pslavdual((n1-1)*8+jj,l)=
     &                   pslavdual((n1-1)*8+jj,l)
     &                   +pslavdual((ii-1)*8+jj,l)
                    if(iflagdualquad.gt.2)then
                      pslavdualpg((n1-1)*8+jj,l)=
     &                     pslavdualpg((n1-1)*8+jj,l)
     &                     +pslavdualpg((ii-1)*8+jj,l)
                    endif
                  enddo
                  do jj=1,nopes2
                    pslavdual((ii-1)*8+jj,l)=0.0
                    if(iflagdualquad.gt.2)then
                      pslavdualpg((ii-1)*8+jj,l)=0.0
                    endif
                  enddo
                endif
              enddo                  
            endif
          endif
          if(checknorm)then
            write(*,*)'pslavdual2,face',l
            do ii=1,8
              help=0.0
              do j=1,8
                help=help+pslavdual((ii-1)*8+j,l)
              enddo
              write(*,105) (pslavdual((ii-1)*8+j,l),j=1,8)
            enddo
          endif
 105      format(8(1x,e15.4))
        enddo
!     
!     FIRST SLAVE SURFACE LOOP DONE
!     
        if(checkbiorthogonality)then
          do l=itiefac(1,i),itiefac(2,i)
            ifaces=islavsurf(1,l)
            nelems=int(ifaces/10)
            jfaces=ifaces-nelems*10
            call getnumberofnodes(nelems,jfaces,lakon,nope,
     &           nopes,idummy) 
!     
!     actual position of the nodes belonging to the
!     slave surface
!     
            do j=1,nope
              konl(j)=kon(ipkon(nelems)+j)
            enddo
            do m=1,nopes
              ifac=getlocno(m,jfaces,nope)               
              do j=1,3
                xl2s(j,m)=co(j,konl(ifac))+
     &               vold(j,konl(ifac))
              enddo
            enddo
            mint2d=islavsurf(2,l+1)-islavsurf(2,l)
            if(mint2d.eq.0) cycle
            indexf=islavsurf(2,l)
            do m=1,mint2d
              xis=pslavsurf(1,indexf+m)
              ets=pslavsurf(2,indexf+m)
              weight=pslavsurf(3,indexf+m)
              ns=l
              iflag=2
              if(nopes.eq.8) then
                if((iflagdualquad.eq.2).or.(iflagdualquad.eq.4))then
                  call dualshape8qtilde(xis,ets,xl2s,xsj2,xs2,
     &                 shp2s,ns,pslavdual,iflag)
                else
                  call dualshape8qtilde_lin(xis,ets,xl2s,xsj2,
     &                 xs2,shp2s,ns,pslavdual,iflag)
                endif
              elseif(nopes.eq.4) then
                call dualshape4q(xis,ets,xl2s,xsj2,xs2,shp2s,ns,
     &               pslavdual,iflag)
              elseif(nopes.eq.6) then
                if((iflagdualquad.eq.2).or.(iflagdualquad.eq.4))then
                  call dualshape6tritilde
     &                 (xis,ets,xl2s,xsj2,xs2,shp2s,ns,
     &                 pslavdual,iflag)
                else
                  call dualshape6tritilde_lin
     &                 (xis,ets,xl2s,xsj2,xs2,shp2s,ns,
     &                 pslavdual,iflag)
                endif
              else
                call dualshape3tri(xis,ets,xl2s,xsj2,xs2,shp2s,ns,
     &               pslavdual,iflag)
              endif
              xim=pslavsurf(1,indexf+m)
              etm=pslavsurf(2,indexf+m)
!     
              if(nopes.eq.8) then
                if((iflagdualquad.eq.2).or.(iflagdualquad.eq.4))then
                  call shape8qtilde(xim,etm,xl2s,xsj2m,xs2m,
     &                 shp2m,iflag)
                else
                  call shape4q(xim,etm,xl2s,xsj2m,xs2m,shp2m,
     &                 iflag)
                  dx=dsqrt(xsj2m(1)**2+xsj2m(2)**2+xsj2m(3)**2)
                  call shape8qtilde_lin(xim,etm,xl2s,xsj2m,
     &                 xs2m,shp2m,iflag) 
                  dx=dsqrt(xsj2m(1)**2+xsj2m(2)**2+xsj2m(3)**2)
                endif
              elseif(nopes.eq.4) then
                call shape4q(xim,etm,xl2s,xsj2m,xs2m,shp2m,iflag)
              elseif(nopes.eq.6) then
                if((iflagdualquad.eq.2).or.(iflagdualquad.eq.4))then
                  call shape6tritilde(xim,etm,xl2s,xsj2m,
     &                 xs2m,shp2m,iflag)
                else
                  call shape6tritilde_lin(xim,etm,xl2s,xsj2m,
     &                 xs2m,shp2m,iflag)
                endif
              else
                call shape3tri(xim,etm,xl2s,xsj2m,xs2m,shp2m,iflag)
              endif
              if(m.eq.1)then
                do j=1,8
                  do jj=1,8
                    contr(8*(j-1)+jj)=0.0
                  enddo
                enddo
              endif 
              dx=dsqrt(xsj2m(1)**2+xsj2m(2)**2+xsj2m(3)**2)
              do j=1,nopes
                ifs=getlocno(j,jfaces,nope)
                nodesf=kon(ipkon(nelems)+ifs)
                locs=j 
                do jj=1,nopes
                  ifm=getlocno(jj,jfaces,nope)
                  nodem=kon(ipkon(nelems)+ifm)
                  locm=jj
                  contribution=shp2s(4,locs)*shp2m(4,locm)
     &                 *weight
     &                 *dx
                  contr(8*(j-1)+jj)=contr(8*(j-1)+jj)
     &                 +contribution
                  iscontr(8*(j-1)+jj)=nodesf
                  imcontr(8*(j-1)+jj)=nodem
                  contribution=0.d0
                enddo
              enddo
            enddo 
            if(checkbiorthogonality) then
              write(*,*) 'checkbiorth.: contri,iscontr,imcontr',l
              do j=1,nopes 
                write(*,105)(contr(8*(j-1)+jj),jj=1,nopes)
              enddo 
            endif
          enddo
        endif
!     
        checknorm=.false.
        if(checknorm)then
          do l=itiefac(1,i),itiefac(2,i)
            ifaces=islavsurf(1,l)
            nelems=int(ifaces/10)
            jfaces=ifaces-nelems*10
            call getnumberofnodes(nelems,jfaces,lakon,nope,
     &           nopes,idummy)
            mint2d=islavsurf(2,l+1)-islavsurf(2,l)
            if(mint2d.eq.0) cycle
            indexf=islavsurf(2,l)
!     
            do j=1,nope
              konl(j)=kon(ipkon(nelems)+j)
            enddo
            do m=1,nopes
              do j=1,3
                ifac=getlocno(m,jfaces,nope)
                xl2s(j,m)=co(j,konl(ifac))+
     &               vold(j,konl(ifac))
              enddo
            enddo 
            if((iflagdualquad.eq.1).or.(iflagdualquad.eq.3))then
              if(nopes.eq.6)then
                nopes=3
              else if(nopes.eq.8)then
                nopes=4
              endif
            endif
            do m=1,mint2d
              xis=pslavsurf(1,indexf+m)
              ets=pslavsurf(2,indexf+m)
              weight=pslavsurf(3,indexf+m)
              ns=l
              iflag=2
              if(nopes.eq.8) then
                call dualshape8qtilde(xis,ets,xl2s,xsj2,xs2,
     &               shp2s,ns,pslavdual,iflag)
              elseif(nopes.eq.4) then
                call dualshape4q(xis,ets,xl2s,xsj2,xs2,shp2s,ns,
     &               pslavdual,iflag)
              elseif(nopes.eq.6) then
                call dualshape6tritilde(xis,ets,xl2s,xsj2,xs2,
     &               shp2s,ns,pslavdual,iflag)
              else
                call dualshape3tri(xis,ets,xl2s,xsj2,xs2,shp2s,ns,
     &               pslavdual,iflag)
              endif
              if(m.eq.1)then
                do j=1,nopes
                  contr(j)=0.0
                enddo
              endif 
              dx= dsqrt(xsj2(1)**2+xsj2(2)**2+xsj2(3)**2)
!     
              do j=1,nopes
                ifs=getlocno(j,jfaces,nope)    
                nodesf=kon(ipkon(nelems)+ifs)
                locs=j
                contribution=shp2s(4,locs)*
     &               weight  
     &               *dx
                contr(j)=contr(j)
     &               +contribution
                iscontr(j)=nodesf
                contribution=0.d0
              enddo   
            enddo  
            write(*,*) 'checknorm: contri,iscontr,imcontr',l
            do j=1,nopes
              write(*,*)contr(j),iscontr(j)
            enddo 
!     
          enddo
        endif
      enddo
c      write(*,*) 'gendualcoeffs.f'
c      do j=1,64
c        write(*,*) j,pslavdual(j,1)
c      enddo
!     
      return
      end
