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
!     subroutine to calculate all coupling matrix combinations and gap contributions of p,q for current slave face l and store it into contri.
!     calculate 
!     \f$ \tilde{B}_d[p,q]=-<\tilde{\Phi}_q,\Psi_p> Id_d\f$ and \f$ \tilde{D}_d \f$.
!     see Phd-thesis Sitzmann Chapter 4
!     
!     Author:Saskia Sitzmann
!     
!     [in] ict               current tie
!     [in] l                 current slave face
!     [in] gapmints		(i) gap between slave surface and master surface in integration point i
!     [in] islavsurf         islavsurf(1,i) slaveface i islavsurf(2,i) pointer into imastsurf and pmastsurf
!     [in] imastsurf         index of masterface corresponding to integration point i
!     [in] pmastsurf         field storing position xil, etal and weight for integration point on master side
!     [out] contr            field containing B_d contributions for current face
!     [out] iscontr          (i) slave node  of contribution(i)
!     [out] imcontr          (i) master node of contribution(i)
!     [out] dcontr            field containing D_d contributions for current face
!     [out] idcontr1          (i) slave node of contribution(i)
!     [out] idcontr2          (i) master node of contribution(i)
!     [out] gcontr            field containing gap contributions for current face
!     [out] igcontr          (i) nodesf of contribution(i)
!     [in] pslavsurf         field storing position xil, etal and weight for integration point on slave side
!     [in] pslavdual         (:,i)coefficients \f$ \alpha_{ij}\f$, \f$ 1,j=1,..8\f$ for dual shape functions for face i
!     [in] nslavnode	(i)pointer into field isalvnode for contact tie i
!     [in] islavnode	field storing the nodes of the slave surface
!     [in] nmastnode	(i)pointer into field imastnode for contact tie i
!     [in] imastnode	field storing the nodes of the master surfaces
!     [out] icounter	counter variable for contr
!     [out] icounter2      	counter variable for dcontr
!     [in] islavact		(i) indicates, if slave node i is active (=-3 no-slave-node, =-2 no-LM-node, =-1 no-gap-node, =0 inactive node, =1 sticky node, =2 slipping/active node) 
!     [in]  iflagdualquad   flag indicating what mortar contact is used (=1 quad-lin, =2 quad-quad, =3 PG quad-lin, =4 PG quad-quad)
!     
      subroutine createbd(ict,l,ipkon,kon,lakon,co, vold, gapmints,
     &     islavsurf,imastsurf,pmastsurf,contr,iscontr,imcontr,
     &     dcontr,idcontr1,idcontr2,gcontr,igcontr,mi,
     &     pslavsurf,pslavdual,nslavnode,islavnode,nmastnode,imastnode,
     &     icounter,icounter2,islavact,iflagdualquad)
!     
      implicit none
!     
      logical debug
!     
      character*8 lakon(*)
!     
      integer ipkon(*),kon(*),konl(20),iflag,m,l,j,jj,
     &     indexe,islavsurf(2,*),
     &     imastsurf(*),ifaces,nelemens,jfaces,ifacem,
     &     mint2d,indexf,nopes1,nopes2,nodem,nodesf,
     &     locs,locm,mi(*),ns,mint2dloc1,mint2dloc2,
     &     ifs,ifm,nope1,nope2, iscontr(*),imcontr(*),getlocno,
     &     jfacem,nelemenm,icounter,idummy,ifac,idcontr1(*),idcontr2(*),
     &     igcontr(*),icounter2,nslavnode(*),islavnode(*),ict,id,
     &     nmastnode(*),imastnode(*),islavact(*),iflagdualquad
!     
      real*8 pmastsurf(2,*),co(3,*),gapmints(*),
     &     vold(0:mi(2),*),weight,dx,help,
     &     ets,xis,etm,xim,xl2s(3,8),xsj2s(3),xsj2s2(3),
     &     shp2s(7,8),xs2s(3,7),xl2m(3,8),xsj2m(3),shp2m(7,8),xs2m(3,7),
     &     contribution,pslavsurf(3,*),pslavdual(64,*),contr(*), 
     &     dcontr(*),dcontribution,gcontr(*),gcontribution,
     &     shp2s2(7,8),xs2s2(3,7)
!     
      debug=.false.
      contribution = 0.d0
      dcontribution = 0.d0
      gcontribution = 0.d0
      icounter=0
      icounter2=0
      ifaces=islavsurf(1,l)
      nelemens = int(ifaces/10)
      jfaces = ifaces - nelemens*10
      indexe = ipkon(nelemens)
      ict=ict+1
      if(debug)write(*,*) 'createbd:l=',l,'tie',ict
      call getnumberofnodes(nelemens,jfaces,lakon,nope1,nopes1,idummy)
      mint2d=islavsurf(2,l+1)-islavsurf(2,l)
      if(mint2d.eq.0) return
      indexf=islavsurf(2,l)
      if(debug)write(*,*) 'createbd:mint2d',mint2d
!     
!     loop over all nodesf of current slave face
!     
      do j=1,nope1
        konl(j)=kon(ipkon(nelemens)+j)
      enddo
      do m=1,nopes1
        ifac=getlocno(m,jfaces,nope1)
        do j=1,3
          xl2s(j,m)=co(j,konl(ifac))+
     &         vold(j,konl(ifac))       
        enddo
      enddo
!     
      do j=1,nopes1
        do jj=1,nopes1
          dcontr(icounter2+nopes1*(j-1)+jj)=0.0
        enddo        
        gcontr(icounter2+j)=0.0
      enddo
!     
      mint2dloc1=1
      mint2dloc2=1
      help=0.d0
!     
!     loop over all integration points created for current slave face
!     
      do 
        if (mint2dloc2.ge.mint2d) exit
!     
!     find current master face
!     
        ifacem=imastsurf(indexf+mint2dloc1)
        nelemenm= int (ifacem/10);
        jfacem=ifacem-nelemenm*10
        call getnumberofnodes(nelemenm,jfacem,lakon,
     &       nope2,nopes2,idummy)
!     
!     find number of integration points belonging to master face     
!     
        do
          if(ifacem.eq.imastsurf(indexf+mint2dloc2+1))then
            mint2dloc2=mint2dloc2+1
            if(mint2dloc2.eq.mint2d) exit     
          else
            exit
          endif  
        enddo
        if(debug)write(*,*) 'createbd:MF,loc1, loc2',ifacem,
     &       mint2dloc1,mint2dloc2       
        help=0.0
        do m=mint2dloc1,mint2dloc2
          xis=pslavsurf(1,indexf+m)
          ets=pslavsurf(2,indexf+m)
          weight=pslavsurf(3,indexf+m)
          ns=l
          iflag = 2
          if(nopes1.eq.8) then
            if(iflagdualquad.eq.2 .or. iflagdualquad.eq.4)then
              call dualshape8qtilde(xis,ets,xl2s,xsj2s,xs2s,shp2s,
     &             ns,pslavdual,iflag)
            else
              call dualshape8qtilde_lin(xis,ets,xl2s,xsj2s,xs2s,
     &             shp2s,ns,pslavdual,iflag)
            endif
          elseif(nopes1.eq.4) then
            call dualshape4q(xis,ets,xl2s,xsj2s,xs2s,shp2s,ns,
     &           pslavdual,iflag)
          elseif(nopes1.eq.6) then
            if(iflagdualquad.eq.2 .or. iflagdualquad.eq.4)then
              call dualshape6tritilde(xis,ets,xl2s,xsj2s,xs2s,shp2s,
     &             ns,pslavdual,iflag)
            else
              call dualshape6tritilde_lin(xis,ets,xl2s,xsj2s,xs2s,
     &             shp2s,ns,pslavdual,iflag)       
            endif
          else
            call dualshape3tri(xis,ets,xl2s,xsj2s,xs2s,shp2s,ns,
     &           pslavdual,iflag)
          endif
          if(nopes1.eq.8) then
            if(iflagdualquad.eq.2 .or. iflagdualquad.eq.4)then
              call shape8qtilde(xis,ets,xl2s,xsj2s2,xs2s2,
     &             shp2s2,iflag)
            else
              call shape8qtilde_lin(xis,ets,xl2s,xsj2s2,xs2s2,
     &             shp2s2,iflag)
            endif
          elseif(nopes1.eq.4) then
            call shape4q(xis,ets,xl2s,xsj2s2,xs2s2,
     &           shp2s2,iflag)
          elseif(nopes1.eq.6) then
            if(iflagdualquad.eq.2 .or. iflagdualquad.eq.4)then
              call shape6tritilde(xis,ets,xl2s,xsj2s2,xs2s2,
     &             shp2s2,iflag)
            else
              call shape6tritilde_lin(xis,ets,xl2s,xsj2s2,xs2s2,
     &             shp2s2,iflag)
            endif
          else
            call shape3tri(xis,ets,xl2s,xsj2s2,xs2s2,
     &           shp2s2,iflag)
          endif    
          xim = pmastsurf(1,indexf+m)
          etm = pmastsurf(2,indexf+m)
!     
          if(nopes2.eq.8) then
            call shape8q(xim,etm,xl2m,xsj2m,xs2m,shp2m,iflag)
          elseif(nopes2.eq.4) then
            call shape4q(xim,etm,xl2m,xsj2m,xs2m,shp2m,iflag)
          elseif(nopes2.eq.6) then
            call shape6tri(xim,etm,xl2m,xsj2m,xs2m,shp2m,iflag)
          else
            call shape3tri(xim,etm,xl2m,xsj2m,xs2m,shp2m,iflag)
          endif
          if(m.eq.mint2dloc1)then
            do j=1,nopes1
              do jj=1,nopes2
                contr(icounter+nopes2*(j-1)+jj)=0.0
              enddo
            enddo
          endif
          dx=dsqrt(xsj2s(1)**2+xsj2s(2)**2+xsj2s(3)**2)   
          do j=1,nopes1
            ifs=getlocno(j,jfaces,nope1) 
            nodesf=kon(ipkon(nelemens)+ifs)
            locs=j
            gcontribution=shp2s(4,locs)
     &           *(gapmints(indexf+m))
     &           *weight
     &           *dx
            gcontr(icounter2+j)=gcontr(icounter2+j)+gcontribution
            if(m.eq.1)then
              call nident(islavnode(nslavnode(ict)+1),
     &             nodesf,(nslavnode(ict+1)-nslavnode(ict)),id)
              if(islavnode(nslavnode(ict)+id).eq.nodesf) then
                igcontr(icounter2+j)=nslavnode(ict)+id
              else
                write(*,*)'createbd: node',nodesf
                write(*,*)'was not catalogued properly in', 
     &               'islavnode'
                call exit(201)
              endif
            endif
            do jj=1,nopes1
              ifs=getlocno(jj,jfaces,nope1) 
              nodem=kon(ipkon(nelemens)+ifs)
              if(m.eq.1)then
                call nident(islavnode(nslavnode(ict)+1),
     &               nodesf,(nslavnode(ict+1)-nslavnode(ict)),id)
                if(islavnode(nslavnode(ict)+id).eq.nodesf) then
                  idcontr1(icounter2+nopes1*(j-1)+jj)=
     &                 nslavnode(ict)+id
                else
                  write(*,*)'createbd: node',nodesf
                  write(*,*)'was not catalogued properly in', 
     &                 'islavnode'
                  call exit(201)
                endif
!     
                call nident(islavnode(nslavnode(ict)+1),
     &               nodem,(nslavnode(ict+1)-nslavnode(ict)),id)
                if(islavnode(nslavnode(ict)+id).eq.nodem) then
                  idcontr2(icounter2+nopes1*(j-1)+jj)=
     &                 nslavnode(ict)+id
                else
                  write(*,*)'createbd: node',nodem
                  write(*,*)'was not catalogued properly in', 
     &                 'islavnode'
                  call exit(201)
                endif
!     
              endif
              dcontribution=shp2s(4,locs)*shp2s2(4,jj)  *weight*dx
              dcontr(icounter2+nopes1*(j-1)+jj)=
     &             dcontr(icounter2+nopes1*(j-1)+jj)+dcontribution
            enddo
            do jj=1,nopes2
              ifm=getlocno(jj,jfacem,nope2)
              nodem=kon(ipkon(nelemenm)+ifm)
              locm=jj
              contribution=shp2s(4,locs)*shp2m(4,locm)
     &             *weight
     &             *dx
              contr(icounter+nopes2*(j-1)+jj)=
     &             contr(icounter+nopes2*(j-1)+jj)+contribution
              if(m.eq.mint2dloc1)then
                call nident(islavnode(nslavnode(ict)+1),
     &               nodesf,(nslavnode(ict+1)-nslavnode(ict)),id)
                if(islavnode(nslavnode(ict)+id).eq.nodesf) then
                  iscontr(icounter+nopes2*(j-1)+jj)=
     &                 nslavnode(ict)+id
                else
                  write(*,*)'createbd: node',nodesf
                  write(*,*)'was not catalogued properly in', 
     &                 ' islavnode'
                  call exit(201)
                endif                            
                call nident(imastnode(nmastnode(ict)+1),
     &               nodem,(nmastnode(ict+1)-nmastnode(ict)),id)
                if(imastnode(nmastnode(ict)+id).eq.nodem) then
                  imcontr(icounter+nopes2*(j-1)+jj)=
     &                 nmastnode(ict)+id
                else
                  write(*,*)'createbd: node',nodem
                  write(*,*)'was not catalogued properly in', 
     &                 ' imastnode',nmastnode(ict)+id,
     &                 imastnode(nmastnode(ict)+id),
     &                 nmastnode(ict)+1,
     &                 nmastnode(ict+1)
                  call exit(201)
                endif
!     
              endif
              contribution=0.d0
            enddo
          enddo
!     
        enddo
        mint2dloc1=mint2dloc2+1
        mint2dloc2=mint2dloc1
        icounter=icounter+nopes1*nopes2
      enddo
      icounter2=icounter2+nopes1*nopes1
!     
      debug=.false.
      if(debug )then
        write(*,*) 'createbd: gcontri,idcontr',l
        do j=1, nopes1
          write(*,*) gcontr(j),igcontr(j),islavact(igcontr(j))
        enddo 
      endif
      ict=ict-1
      return 
      end
