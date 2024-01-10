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
!     author: Yannick Muller
!     
      subroutine gen3dchanoutput(co,nup,ndo,nelem,g,dg,lakon,s0,sqrts0,
     &     sqrts0up,sqrts0do,nknet,conet,nenet,konnet,vnet,bup,bdo,
     &     hup,hdo,uup,udo,thup,thdo,thetaup,thetado,jumpup,sfr,hfr,
     &     jumpdo,sba,hba,ndata,dl,ha,rho,xflow,i)
!     
      implicit none
!     
      character*8 lakon(*)
!     
      integer nup,ndo,i,k,m,nelem,nknet,nenet,konnet(*),indexe,
     &     jumpup(*),jumpdo(*),index,ndata
!     
      real*8 e1(3),e2(3),e3(3),e3up(3),e3do(3),dd,co(3,*),g(3),dg,
     &     d13,aa,bb,cc,disc,al,conet(3,*),vnet(3,*),bup,bdo,hup,
     &     hdo,uup,udo,thup,thdo,thetaup,thetado,sfr(*),hfr(*),
     &     sba(*),hba(*),h1,h2,area1,area2,u1,u2,fr1,fr2,z1,z2,
     &     th1,th2,codo(3),coup(3),dl,s0,sqrts0,sqrts0up,sqrts0do,
     &     ha,rho,xflow
!     
!     generating a three-dimensional version of the network
!     
!     unit vector e1 along the bottom of the channel
!     
      do k=1,3
        e1(k)=co(k,ndo)-co(k,nup)
      enddo
      dd=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
      do k=1,3
        e1(k)=e1(k)/dd
      enddo
!     
!     unit vector e3 in the plane consisting of e1 and the gravity
!     vector and orthogonal to the bottom specified by s0 except for
!     the downstream end of a sluice gate and weir and the upstream
!     end of a reservoir
!     
      do k=1,3
        e3(k)=-g(k)/dg
      enddo
      d13=e1(1)*e3(1)+e1(2)*e3(2)+e1(3)*e3(3)
!     
      if((lakon(nelem)(6:7).eq.'SG').or.
     &     (lakon(nelem)(6:7).eq.'WE')) then
        do k=1,3
          e3up(k)=e3(k)
        enddo
      else
        aa=s0*s0
        bb=aa*d13
        cc=d13*d13-sqrts0up*sqrts0up
        disc=dsqrt(bb**2-aa*cc)
        al=(-bb+disc)/aa
        if(al.lt.0.d0) then
          al=(-bb-disc)/aa
        endif
        do k=1,3
          e3up(k)=e1(k)+al*e3(k)
        enddo
        dd=dsqrt(e3up(1)*e3up(1)+e3up(2)*e3up(2)+e3up(3)*e3up(3))
        do k=1,3
          e3up(k)=e3up(k)/dd
        enddo
      endif
!     
      if(lakon(nelem)(6:7).eq.'RE') then
        do k=1,3
          e3do(k)=e3(k)
        enddo
      else
        aa=s0*s0
        bb=aa*d13
        cc=d13*d13-sqrts0do*sqrts0do
        disc=dsqrt(bb**2-aa*cc)
        al=(-bb+disc)/aa
        if(al.lt.0.d0) then
          al=(-bb-disc)/aa
        endif
        do k=1,3
          e3do(k)=e1(k)+al*e3(k)
        enddo
        dd=dsqrt(e3do(1)*e3do(1)+e3do(2)*e3do(2)+e3do(3)*e3do(3))
        do k=1,3
          e3do(k)=e3do(k)/dd
        enddo
      endif
!     
!     unit vector e2 along the bottom in width direction
!     
      e2(1)=e3(2)*e1(3)-e3(3)*e1(2)
      e2(2)=e3(3)*e1(1)-e3(1)*e1(3)
      e2(3)=e3(1)*e1(2)-e3(2)*e1(1)
      dd=dsqrt(e2(1)*e2(1)+e2(2)*e2(2)+e2(3)*e2(3))
      do k=1,3
        e2(k)=e2(k)/dd
      enddo
!     
!     generating nodes and elements
!     
      indexe=8*nenet
      if((lakon(nelem)(6:7).ne.'  ').and.
     &     (lakon(nelem)(6:7).ne.'DS')) then
        nenet=nenet+1
!     
!     node 1
!     
        nknet=nknet+1
        if(lakon(nelem)(6:7).eq.'RE') then
          do k=1,3
            conet(k,nknet)=co(k,nup)-e2(k)*bup/2.d0
          enddo
          vnet(1,nknet)=hdo
          vnet(2,nknet)=uup/dsqrt(dg*hup*sqrts0)
          vnet(3,nknet)=thup
          konnet(indexe+1)=nknet
        else
          do k=1,3
            conet(k,nknet)=co(k,nup)-e2(k)*bup/2.d0
          enddo
          vnet(1,nknet)=hup
          vnet(2,nknet)=uup/dsqrt(dg*hup*sqrts0up)
          vnet(3,nknet)=thup
          konnet(indexe+1)=nknet
        endif
!     
!     node 2
!     
        nknet=nknet+1
        do k=1,3
          conet(k,nknet)=co(k,ndo)-e2(k)*bdo/2.d0
        enddo
        vnet(1,nknet)=hdo
        vnet(2,nknet)=udo/dsqrt(dg*hdo*sqrts0do)
        vnet(3,nknet)=thdo
        konnet(indexe+2)=nknet
!     
!     node 3
!     
        nknet=nknet+1
        do k=1,3
          conet(k,nknet)=co(k,ndo)+e2(k)*bdo/2.d0
        enddo
        vnet(1,nknet)=hdo
        vnet(2,nknet)=udo/dsqrt(dg*hdo*sqrts0do)
        vnet(3,nknet)=thdo
        konnet(indexe+3)=nknet
!     
!     node 4
!     
        nknet=nknet+1
        if(lakon(nelem)(6:7).eq.'RE') then
          do k=1,3
            conet(k,nknet)=co(k,nup)+e2(k)*bup/2.d0
          enddo
          vnet(1,nknet)=hdo
          vnet(2,nknet)=uup/dsqrt(dg*hup*sqrts0)
          vnet(3,nknet)=thup
          konnet(indexe+4)=nknet
        else
          do k=1,3
            conet(k,nknet)=co(k,nup)+e2(k)*bup/2.d0
          enddo
          vnet(1,nknet)=hup
          vnet(2,nknet)=uup/dsqrt(dg*hup*sqrts0up)
          vnet(3,nknet)=thup
          konnet(indexe+4)=nknet
        endif
!     
!     node 5
!     
        nknet=nknet+1
        if(lakon(nelem)(6:7).eq.'RE') then
          do k=1,3
            conet(k,nknet)=co(k,nup)+e3up(k)*hdo
     &           -e2(k)*(bup/2.d0+dtan(thetaup)*hdo)
          enddo
          vnet(1,nknet)=hdo
          vnet(2,nknet)=uup/dsqrt(dg*hup*sqrts0)
          vnet(3,nknet)=thup
          konnet(indexe+5)=nknet
        else
          do k=1,3
            conet(k,nknet)=co(k,nup)+e3up(k)*hup
     &           -e2(k)*(bup/2.d0+dtan(thetaup)*hup)
          enddo
          vnet(1,nknet)=hup
          vnet(2,nknet)=uup/dsqrt(dg*hup*sqrts0up)
          vnet(3,nknet)=thup
          konnet(indexe+5)=nknet
        endif
!     
!     node 6
!     
        nknet=nknet+1
        if((lakon(nelem)(6:7).eq.'SG').and.(hdo.gt.ha)) then
          do k=1,3
            conet(k,nknet)=co(k,ndo)+e3do(k)*hup/sqrts0
     &           -e2(k)*(bdo/2.d0+dtan(thetado)*hup/sqrts0)
          enddo
          vnet(1,nknet)=hup
          vnet(2,nknet)=udo/dsqrt(dg*hdo*sqrts0)
          vnet(3,nknet)=thdo
          konnet(indexe+6)=nknet
        else
          do k=1,3
            conet(k,nknet)=co(k,ndo)+e3do(k)*hdo
     &           -e2(k)*(bdo/2.d0+dtan(thetado)*hdo)
          enddo
          vnet(1,nknet)=hdo
          vnet(2,nknet)=udo/dsqrt(dg*hdo*sqrts0do)
          vnet(3,nknet)=thdo
          konnet(indexe+6)=nknet
        endif
!     
!     node 7
!     
        nknet=nknet+1
        if((lakon(nelem)(6:7).eq.'SG').and.(hdo.gt.ha)) then
          do k=1,3
            conet(k,nknet)=co(k,ndo)+e3do(k)*hup/sqrts0
     &           +e2(k)*(bdo/2.d0+dtan(thetado)*hup/sqrts0)
          enddo
          vnet(1,nknet)=hup
          vnet(2,nknet)=udo/dsqrt(dg*hdo*sqrts0)
          vnet(3,nknet)=thdo
          konnet(indexe+7)=nknet
        else
          do k=1,3
            conet(k,nknet)=co(k,ndo)+e3do(k)*hdo
     &           +e2(k)*(bdo/2.d0+dtan(thetado)*hdo)
          enddo
          vnet(1,nknet)=hdo
          vnet(2,nknet)=udo/dsqrt(dg*hdo*sqrts0do)
          vnet(3,nknet)=thdo
          konnet(indexe+7)=nknet
        endif
!     
!     node 8
!     
        nknet=nknet+1
        if(lakon(nelem)(6:7).eq.'RE') then
          do k=1,3
            conet(k,nknet)=co(k,nup)+e3up(k)*hdo
     &           +e2(k)*(bup/2.d0+dtan(thetaup)*hdo)
          enddo
          vnet(1,nknet)=hdo
          vnet(2,nknet)=uup/dsqrt(dg*hup*sqrts0)
          vnet(3,nknet)=thup
          konnet(indexe+8)=nknet
        else
          do k=1,3
            conet(k,nknet)=co(k,nup)+e3up(k)*hup
     &           +e2(k)*(bup/2.d0+dtan(thetaup)*hup)
          enddo
          vnet(1,nknet)=hup
          vnet(2,nknet)=uup/dsqrt(dg*hup*sqrts0up)
          vnet(3,nknet)=thup
          konnet(indexe+8)=nknet
        endif
      elseif(lakon(nelem)(6:7).eq.'  ') then
!     
!     straight channel
!     
        index=(i-1)*ndata
!     
!     frontwater curve
!     
        do m=1,jumpup(i)-1
          indexe=8*nenet
          nenet=nenet+1
          do k=1,3
            coup(k)=(dl-sfr(index+m))/dl*co(k,nup)
     &           +sfr(index+m)/dl*co(k,ndo)
          enddo
          h1=hfr(index+m)
          area1=(bup+h1*dtan(thetaup))*h1
          u1=xflow/(rho*area1)
          fr1=u1/dsqrt(dg*h1*sqrts0)
          z1=(-g(1)*coup(1)-g(2)*coup(2)-g(3)*coup(3))/dg
          th1=h1*sqrts0+u1**2/(2.d0*dg)+z1
!     
          do k=1,3
            codo(k)=(dl-sfr(index+m+1))/dl*co(k,nup)
     &           +sfr(index+m+1)/dl*co(k,ndo)
          enddo
          h2=hfr(index+m+1)
          area2=(bup+h2*dtan(thetaup))*h2
          u2=xflow/(rho*area2)
          fr2=u2/dsqrt(dg*h2*sqrts0)
          z2=(-g(1)*codo(1)-g(2)*codo(2)-g(3)*codo(3))/dg
          th2=h2*sqrts0+u2**2/(2.d0*dg)+z2
!     
!     node 1
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=coup(k)-e2(k)*bup/2.d0
          enddo
          vnet(1,nknet)=h1
          vnet(2,nknet)=fr1
          vnet(3,nknet)=th1
          konnet(indexe+1)=nknet
!     
!     node 2
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=codo(k)-e2(k)*bdo/2.d0
          enddo
          vnet(1,nknet)=h2
          vnet(2,nknet)=fr2
          vnet(3,nknet)=th2
          konnet(indexe+2)=nknet
!     
!     node 3
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=codo(k)+e2(k)*bdo/2.d0
          enddo
          vnet(1,nknet)=h2
          vnet(2,nknet)=fr2
          vnet(3,nknet)=th2
          konnet(indexe+3)=nknet
!     
!     node 4
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=coup(k)+e2(k)*bup/2.d0
          enddo
          vnet(1,nknet)=h1
          vnet(2,nknet)=fr1
          vnet(3,nknet)=th1
          konnet(indexe+4)=nknet
!     
!     node 5
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=coup(k)+e3up(k)*h1
     &           -e2(k)*(bup/2.d0+dtan(thetaup)*h1)
          enddo
          vnet(1,nknet)=h1
          vnet(2,nknet)=fr1
          vnet(3,nknet)=th1
          konnet(indexe+5)=nknet
!     
!     node 6
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=codo(k)+e3do(k)*h2
     &           -e2(k)*(bdo/2.d0+dtan(thetado)*h2)
          enddo
          vnet(1,nknet)=h2
          vnet(2,nknet)=fr2
          vnet(3,nknet)=th2
          konnet(indexe+6)=nknet
!     
!     node 7
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=codo(k)+e3do(k)*h2
     &           +e2(k)*(bdo/2.d0+dtan(thetado)*h2)
          enddo
          vnet(1,nknet)=h2
          vnet(2,nknet)=fr2
          vnet(3,nknet)=th2
          konnet(indexe+7)=nknet
!     
!     node 8
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=coup(k)+e3up(k)*h1
     &           +e2(k)*(bup/2.d0+dtan(thetaup)*h1)
          enddo
          vnet(1,nknet)=h1
          vnet(2,nknet)=fr1
          vnet(3,nknet)=th1
          konnet(indexe+8)=nknet
        enddo
!     
!     jump
!     
        if((jumpup(i).gt.0).and.(jumpdo(i).lt.ndata+1)) then
          indexe=8*nenet
          nenet=nenet+1
          do k=1,3
            coup(k)=(dl-sfr((i-1)*ndata+jumpup(i)))/dl*co(k,nup)
     &           +sfr((i-1)*ndata+jumpup(i))/dl*co(k,ndo)
          enddo
          h1=hfr((i-1)*ndata+jumpup(i))
          area1=(bup+h1*dtan(thetaup))*h1
          u1=xflow/(rho*area1)
          fr1=u1/dsqrt(dg*h1*sqrts0)
          z1=(-g(1)*coup(1)-g(2)*coup(2)-g(3)*coup(3))/dg
          th1=h1*sqrts0+u1**2/(2.d0*dg)+z1
!     
          do k=1,3
            codo(k)=(dl-sba((i-1)*ndata+jumpdo(i)))/dl*co(k,nup)
     &           +sba((i-1)*ndata+jumpdo(i))/dl*co(k,ndo)
          enddo
          h2=hba((i-1)*ndata+jumpdo(i))
          area2=(bup+h2*dtan(thetaup))*h2
          u2=xflow/(rho*area2)
          fr2=u2/dsqrt(dg*h2*sqrts0)
          z2=(-g(1)*codo(1)-g(2)*codo(2)-g(3)*codo(3))/dg
          th2=h2*sqrts0+u2**2/(2.d0*dg)+z2
!     
!     node 1
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=coup(k)-e2(k)*bup/2.d0
          enddo
          vnet(1,nknet)=h1
          vnet(2,nknet)=fr1
          vnet(3,nknet)=th1
          konnet(indexe+1)=nknet
!     
!     node 2
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=codo(k)-e2(k)*bdo/2.d0
          enddo
          vnet(1,nknet)=h2
          vnet(2,nknet)=fr2
          vnet(3,nknet)=th2
          konnet(indexe+2)=nknet
!     
!     node 3
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=codo(k)+e2(k)*bdo/2.d0
          enddo
          vnet(1,nknet)=h2
          vnet(2,nknet)=fr2
          vnet(3,nknet)=th2
          konnet(indexe+3)=nknet
!     
!     node 4
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=coup(k)+e2(k)*bup/2.d0
          enddo
          vnet(1,nknet)=h1
          vnet(2,nknet)=fr1
          vnet(3,nknet)=th1
          konnet(indexe+4)=nknet
!     
!     node 5
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=coup(k)+e3up(k)*h1
     &           -e2(k)*(bup/2.d0+dtan(thetaup)*h1)
          enddo
          vnet(1,nknet)=h1
          vnet(2,nknet)=fr1
          vnet(3,nknet)=th1
          konnet(indexe+5)=nknet
!     
!     node 6
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=codo(k)+e3do(k)*h2
     &           -e2(k)*(bdo/2.d0+dtan(thetado)*h2)
          enddo
          vnet(1,nknet)=h2
          vnet(2,nknet)=fr2
          vnet(3,nknet)=th2
          konnet(indexe+6)=nknet
!     
!     node 7
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=codo(k)+e3do(k)*h2
     &           +e2(k)*(bdo/2.d0+dtan(thetado)*h2)
          enddo
          vnet(1,nknet)=h2
          vnet(2,nknet)=fr2
          vnet(3,nknet)=th2
          konnet(indexe+7)=nknet
!     
!     node 8
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=coup(k)+e3up(k)*h1
     &           +e2(k)*(bup/2.d0+dtan(thetaup)*h1)
          enddo
          vnet(1,nknet)=h1
          vnet(2,nknet)=fr1
          vnet(3,nknet)=th1
          konnet(indexe+8)=nknet
        endif
!     
!     backwater curve
!     
        do m=jumpdo(i),ndata-1
          indexe=8*nenet
          nenet=nenet+1
          do k=1,3
            coup(k)=(dl-sba(index+m))/dl*co(k,nup)
     &           +sba(index+m)/dl*co(k,ndo)
          enddo
          h1=hba(index+m)
          area1=(bup+h1*dtan(thetaup))*h1
          u1=xflow/(rho*area1)
          fr1=u1/dsqrt(dg*h1*sqrts0)
          z1=(-g(1)*coup(1)-g(2)*coup(2)-g(3)*coup(3))/dg
          th1=h1*sqrts0+u1**2/(2.d0*dg)+z1
!     
          do k=1,3
            codo(k)=(dl-sba(index+m+1))/dl*co(k,nup)
     &           +sba(index+m+1)/dl*co(k,ndo)
          enddo
          h2=hba(index+m+1)
          area2=(bup+h2*dtan(thetaup))*h2
          u2=xflow/(rho*area2)
          fr2=u2/dsqrt(dg*h2*sqrts0)
          z2=(-g(1)*codo(1)-g(2)*codo(2)-g(3)*codo(3))/dg
          th2=h2*sqrts0+u2**2/(2.d0*dg)+z2
!     
!     node 1
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=coup(k)-e2(k)*bup/2.d0
          enddo
          vnet(1,nknet)=h1
          vnet(2,nknet)=fr1
          vnet(3,nknet)=th1
          konnet(indexe+1)=nknet
!     
!     node 2
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=codo(k)-e2(k)*bdo/2.d0
          enddo
          vnet(1,nknet)=h2
          vnet(2,nknet)=fr2
          vnet(3,nknet)=th2
          konnet(indexe+2)=nknet
!     
!     node 3
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=codo(k)+e2(k)*bdo/2.d0
          enddo
          vnet(1,nknet)=h2
          vnet(2,nknet)=fr2
          vnet(3,nknet)=th2
          konnet(indexe+3)=nknet
!     
!     node 4
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=coup(k)+e2(k)*bup/2.d0
          enddo
          vnet(1,nknet)=h1
          vnet(2,nknet)=fr1
          vnet(3,nknet)=th1
          konnet(indexe+4)=nknet
!     
!     node 5
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=coup(k)+e3up(k)*h1
     &           -e2(k)*(bup/2.d0+dtan(thetaup)*h1)
          enddo
          vnet(1,nknet)=h1
          vnet(2,nknet)=fr1
          vnet(3,nknet)=th1
          konnet(indexe+5)=nknet
!     
!     node 6
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=codo(k)+e3do(k)*h2
     &           -e2(k)*(bdo/2.d0+dtan(thetado)*h2)
          enddo
          vnet(1,nknet)=h2
          vnet(2,nknet)=fr2
          vnet(3,nknet)=th2
          konnet(indexe+6)=nknet
!     
!     node 7
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=codo(k)+e3do(k)*h2
     &           +e2(k)*(bdo/2.d0+dtan(thetado)*h2)
          enddo
          vnet(1,nknet)=h2
          vnet(2,nknet)=fr2
          vnet(3,nknet)=th2
          konnet(indexe+7)=nknet
!     
!     node 8
!     
          nknet=nknet+1
          do k=1,3
            conet(k,nknet)=coup(k)+e3up(k)*h1
     &           +e2(k)*(bup/2.d0+dtan(thetaup)*h1)
          enddo
          vnet(1,nknet)=h1
          vnet(2,nknet)=fr1
          vnet(3,nknet)=th1
          konnet(indexe+8)=nknet
        enddo
      endif
!     
      return
      end
