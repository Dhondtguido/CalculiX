!     
!     convection with the walls: contribution to the energy equations
!     
      do i=1,nload
        if(sideload(i)(3:4).eq.'FC') then
          nelem=nelemload(1,i)
          index=ipkon(nelem)
          if(index.lt.0) cycle
          lakonl=lakon(nelem)
          node=nelemload(2,i)
          ieq=nacteq(0,node)
          if(ieq.eq.0) then 
            cycle
          endif
!     
          call nident(itg,node,ntg,id)
!     
!     calculate the area
!     
          read(sideload(i)(2:2),'(i1)') ig
!     
!     number of nodes and integration points in the face
!     
          if(lakonl(4:4).eq.'2') then
            nope=20
            nopes=8
          elseif(lakonl(4:4).eq.'8') then
            nope=8
            nopes=4
          elseif(lakonl(4:5).eq.'10') then
            nope=10
            nopes=6
          elseif(lakonl(4:4).eq.'4') then
            nope=4
            nopes=3
          elseif(lakonl(4:5).eq.'15') then
            nope=15
          else
            nope=6
          endif
!     
          if(lakonl(4:5).eq.'8R') then
            mint2d=1
          elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R'))
     &           then
            if(lakonl(7:7).eq.'A') then
              mint2d=2
            else
              mint2d=4
            endif
          elseif(lakonl(4:4).eq.'2') then
            mint2d=9
          elseif(lakonl(4:5).eq.'10') then
            mint2d=3
          elseif(lakonl(4:4).eq.'4') then
            mint2d=1
          endif
!     
          if(lakonl(4:4).eq.'6') then
            mint2d=1
            if(ig.le.2) then
              nopes=3
            else
              nopes=4
            endif
          endif
          if(lakonl(4:5).eq.'15') then
            if(ig.le.2) then
              mint2d=3
              nopes=6
            else
              mint2d=4
              nopes=8
            endif
          endif
!     
!     connectivity of the element
!     
          do k=1,nope
            konl(k)=kon(index+k)
          enddo
!     
!     coordinates of the nodes belonging to the face
!     
          if((nope.eq.20).or.(nope.eq.8)) then
            do k=1,nopes
              tl2(k)=v(0,konl(ifaceq(k,ig)))
              do j=1,3
                xl2(j,k)=co(j,konl(ifaceq(k,ig)))+
     &               v(j,konl(ifaceq(k,ig)))
              enddo
            enddo
          elseif((nope.eq.10).or.(nope.eq.4)) then
            do k=1,nopes
              tl2(k)=v(0,konl(ifacet(k,ig)))
              do j=1,3
                xl2(j,k)=co(j,konl(ifacet(k,ig)))+
     &               v(j,konl(ifacet(k,ig)))
              enddo
            enddo
          else
            do k=1,nopes
              tl2(k)=v(0,konl(ifacew(k,ig)))
              do j=1,3
                xl2(j,k)=co(j,konl(ifacew(k,ig)))+
     &               v(j,konl(ifacew(k,ig)))
              enddo
            enddo
          endif
!     
!     integration to obtain the area and the mean
!     temperature
!     
          do m=1,mint2d
            if((lakonl(4:5).eq.'8R').or.
     &           ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
              xi=gauss2d1(1,m)
              et=gauss2d1(2,m)
              weight=weight2d1(m)
            elseif((lakonl(4:4).eq.'8').or.
     &             (lakonl(4:6).eq.'20R').or.
     &             ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
              xi=gauss2d2(1,m)
              et=gauss2d2(2,m)
              weight=weight2d2(m)
            elseif(lakonl(4:4).eq.'2') then
              xi=gauss2d3(1,m)
              et=gauss2d3(2,m)
              weight=weight2d3(m)
            elseif((lakonl(4:5).eq.'10').or.
     &             ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
              xi=gauss2d5(1,m)
              et=gauss2d5(2,m)
              weight=weight2d5(m)
            elseif((lakonl(4:4).eq.'4').or.
     &             ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
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
!     
            dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &           xsj2(3)*xsj2(3))
            areaj=dxsj2*weight
!     
            temp=0.d0
            do k=1,3
              coords(k)=0.d0
            enddo
            do j=1,nopes
              temp=temp+tl2(j)*shp2(4,j)
              do k=1,3
                coords(k)=coords(k)+xl2(k,j)*shp2(4,j)
              enddo
            enddo
!     
            sinktemp=v(0,node)
            if(sideload(i)(5:6).ne.'NU') then
              h(1)=xloadact(1,i)
            else
              read(sideload(i)(2:2),'(i1)') jltyp
              jltyp=jltyp+10
              call film(h,sinktemp,temp,istep,
     &             iinc,tvar,nelem,m,coords,jltyp,field,nfield,
     &             sideload(i),node,areaj,v,mi,ipkon,kon,lakon,
     &             iponoel,inoel,ielprop,prop,ielmat,shcon,nshcon,
     &             rhcon,nrhcon,ntmat_,cocon,ncocon,
     &             ipobody,xbodyact,ibody,heatnod,heatfac)
            endif
            if(lakonl(5:7).eq.'0RA') then
              term=2.d0*((temp-sinktemp)*h(1)+heatnod)*dxsj2*weight
              bc(ieq)=bc(ieq)+term
              qat=qat+dabs(term)
              nalt=nalt+1
            else
              term=((temp-sinktemp)*h(1)+heatnod)*dxsj2*weight
              bc(ieq)=bc(ieq)+term
              qat=qat+dabs(term)
              nalt=nalt+1
            endif
          enddo
        endif
      enddo
!     
!     prescribed heat generation: contribution to the energy equations        
!     
      do i=1,ntg
        node=itg(i)
        idof=8*(node-1)
        call nident(ikforc,idof,nforc,id)
        if(id.gt.0) then
          if(ikforc(id).eq.idof) then
            ieq=nacteq(0,node)
            if(ieq.ne.0) then
              term=xforcact(ilforc(id))
              bc(ieq)=bc(ieq)+term
              qat=qat+dabs(term)
              nalt=nalt+1
            endif
            cycle
          endif
        endif
      enddo
