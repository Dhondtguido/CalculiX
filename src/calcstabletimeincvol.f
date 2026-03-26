!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2021 Guido Dhondt
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
      subroutine calcstabletimeincvol(ne0,elcon,nelcon,
     &     rhcon,nrhcon,alcon,nalcon,orab,ntmat_,ithermal,alzero,
     &     plicon,nplicon,plkcon,nplkcon,npmat_,mi,dtime,
     &     xstiff,ncmat_,vold,ielmat,t0,t1,
     &     matname,lakon,wavespeed,nmat,ipkon,co,kon,dtvol,alpha,
     &     smscale,dtset,mscalmethod,mortar,jobnamef,iperturb)
!     
!     **************
!     ------------Wavespeed calculation------------------------CARLO MT
!     Calculates the propagation wave speed in a material, selecting
!     appropriate procedure between isotropic, single crystals, and
!     anisotropic materials. All other cases of orthotropy are treated
!     as anisotropic.
!     
!     Based on the procedure in:
!     C. Lane. The Development of a 2D Ultrasonic Array Inspection
!     for Single Crystal Turbine Blades.
!     Switzerland: Springer International Publishing, 2014.
!     
!     ------------Critical time step calculation---------------CARLO MT
!     Calculates the critical time increment (CTI) based on the Courant
!     Criterion for Explicit Dynamics calculations. Temperature is
!     assumed averaged from the centroid of the element and material
!     wave propagation speeds must be calculated before.
!     
!     ------------Mass Scaling -----------------------------CATHARINA Czech
!     
!     mscalmethod<0: no explicit dynamics
!     0: explicit dynamics, no scaling
!     1: explicit dynamics, volumetric scaling active
!     2: explicit dynamics, contact scaling active
!     3: explicit dynamics, volumetric and contact scaling
!     
      implicit none
!     
      character*8 lakon(*),lakonl
      character*80 matname(*),amat
      character*132 jobnamef(*)
      character*256 fn
!     
      integer i,j,i1,nelem,ne0,nelcon(2,*),nrhcon(*),nalcon(2,*),imat,
     &     ntmat_,ithermal(*),mattyp,ihyper,istiff,kode,mi(*),kk,
     &     nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),ncmat_,iorth,
     &     ielmat(mi(3),*),nope,iorien,ipkon(*),null,ilen,
     &     konl(26),nopered,npmat_,nmat,kon(*),indexe,iflag,nopes,
     &     nfaces,ig,ifaceq(8,6),ifacet(6,4),ifacew(8,5),
     &     mscalmethod,icount,mortar,iperturb(*)
!     
      real*8 stiff(21),wavespeed(*),rhcon(0:1,ntmat_,*),volfac,
     &     alcon(0:6,ntmat_,*),coords(3),orab(7,*),rho,alzero(*),
     &     t0l,t1l,elconloc(ncmat_),eth(6),plicon(0:2*npmat_,ntmat_,*),
     &     plkcon(0:2*npmat_,ntmat_,*),plconloc(802),dtime,
     &     xstiff(27,mi(1),*),elcon(0:ncmat_,ntmat_,*),aa,bb,cc,
     &     t0(*),t1(*),shp(4,20),vold(0:mi(2),*),tt,
     &     e,un,wavspd,xi,et,ze,weight,co(3,*),xl(3,20),xsj,
     &     xl2(3,9),xsj2(3),xs2(3,7),shp2(7,9),hmin,area,
     &     volume,dtvol,safefac,alpha,bet,gam,critom,damping,
     &     elemfac,quadfac,smscale(*),dtset   
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/
      data ifacew /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     4,6,3,1,12,15,9,13/
!     
      include "gauss.f"
!     
      iflag=2
      dtvol=1.d30
!     
!     the safety factor is smaller (i.e. smaller time step) for
!     massless contact due to the different time scheme (the
!     alpha method is not used for massless contact)
!     
      if(mortar.ne.-1) then
        safefac=0.80d0/1.3d0
c     safefac=0.80d0
      else
c     safefac=0.80d0/1.3d0
        safefac=0.5d0
      endif
!      
      quadfac=0.3d0
!     
      damping=0
      icount=0
!     
      bet=(1.d0-alpha)*(1.d0-alpha)/4.d0
      gam=0.5d0-alpha
!     
!     Calculation of Omega Critical
!     Om_cr=dt*freq_max
!     
!     formula (2.477) with +sign from Dhondt, Wiley(2004), the Finite      
!        Element Method for Three-Dimensional Thermomechanical Applications.
!     factor 0.98 in critom from miranda?
!      
      aa=1.d0-alpha*(1.d0+alpha*(2.d0+alpha))
      bb=2.d0*damping*(1.d0+alpha*(1.d0+2.d0*alpha))
      cc=-4.d0
      critom=0.98d0*(-bb+dsqrt(bb*bb-aa*cc))/aa
!     
      write(*,*)'Calculating Material Wave Speeds...'
      write(*,*)
!     
      null=0
!     
!     Initialization of wavespeed
!     
      do i=1,nmat
        wavespeed(i)=-1.d0
      enddo
!     
!     ** DO per element
      do nelem=1,ne0
        if(ipkon(nelem).lt.0) cycle
c        write(*,*) 'calcstabletimeincvol ',nelem
!     
        lakonl=lakon(nelem)
        imat=ielmat(1,nelem)
        amat=matname(imat)
        wavspd=wavespeed(imat)
        indexe=ipkon(nelem)
!     
        if(lakonl(1:2).ne.'ES') then
!     
!     orientation is not important for the wave speed
!     
          iorien=0
!     
          if(nelcon(1,imat).lt.0) then
            ihyper=1
          else
            ihyper=0
          endif
!     
!     ------------Shape Functions --------------------------------
!
          elemfac=1.d0
!
          if(lakon(nelem)(4:5).eq.'20')then
            nope=20
            nopes=8
            nfaces=6
            elemfac=quadfac
          elseif(lakon(nelem)(1:5).eq.'C3D8I')then
            nope=8
            nopes=4
            nfaces=6
            elemfac=0.125d0
          elseif(lakon(nelem)(4:5).eq.'8R') then
            nope=8
            nopes=4
            nfaces=6
            elemfac=0.625d0
          elseif(lakon(nelem)(4:4).eq.'8') then
            nope=8
            nopes=4
            nfaces=6
c            if(mortar.eq.-1) elemfac=0.9d0
            elemfac=0.9d0
          elseif(lakon(nelem)(4:5).eq.'10')then
            nope=10
            nopes=6
            nfaces=4
            elemfac=quadfac
          elseif(lakon(nelem)(4:4).eq.'4') then
            nope=4
            nopes=3
            nfaces=4
          elseif(lakon(nelem)(4:5).eq.'15')then
            nope=15
            nfaces=5
            elemfac=quadfac
          elseif(lakon(nelem)(4:4).eq.'6') then
            nope=6
            nfaces=5
          else
            cycle
          endif
!     
!     Find center of the element for avg temp value on the element  to
!     get properties later
!     if HEX
!     
          if((lakon(nelem)(4:5).eq.'20').or.
     &         (lakon(nelem)(4:4).eq.'8')) then
            xi=0.d0
            et=0.d0
            ze=0.d0
            weight=8.d0
!     if TET
          elseif((lakon(nelem)(4:5).eq.'10').or.
     &           (lakon(nelem)(4:4).eq.'4')) then
            xi=gauss3d4(1,1)    !all integration points mint3d
            et=gauss3d4(2,1)
            ze=gauss3d4(3,1)
            weight=weight3d4(1)
!     elseif WEDGES
          elseif((lakonl(4:5).eq.'15').or.
     &           (lakonl(4:4).eq.'6'))then
            xi=1.d0/3.d0
            et=1.d0/3.d0
            ze=0.d0
            weight=1.d0
          endif
!     
!     computation of the coordinates of the local nodes
          do i=1,nope
            konl(i)=kon(indexe+i)
            do j=1,3
              xl(j,i)=co(j,konl(i))
            enddo
          enddo
!
          if(nope.eq.20)then
            call shape20h(xi,et,ze,xl,xsj,shp,iflag)
          elseif(nope.eq.8) then
            call shape8h(xi,et,ze,xl,xsj,shp,iflag)
          elseif(nope.eq.10)then
            call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
          elseif(nope.eq.4) then
            call shape4tet (xi,et,ze,xl,xsj,shp,iflag)
          elseif(nope.eq.15)then
            call shape15w(xi,et,ze,xl,xsj,shp,iflag)
          else
            call shape6w(xi,et,ze,xl,xsj,shp,iflag)
          endif
!     
!     ------------Wavespeed calcualtion--------------------------------
!     calculating the temperature in the integration point
!     
          kk=1                  ! Only 1 integration point is considered, CENTER
          t0l=0.d0
          t1l=0.d0
          if(ithermal(1).eq.1) then
            if(lakonl(4:5).eq.'8 ') then
              do i1=1,nope
                t0l=t0l+t0(konl(i1))/8.d0
                t1l=t1l+t1(konl(i1))/8.d0
              enddo
            elseif(lakonl(4:6).eq.'20 ')then
              nopered=20
              call lintemp(t0,konl,nopered,kk,t0l)
              call lintemp(t1,konl,nopered,kk,t1l)
            elseif(lakonl(4:6).eq.'10T') then
              call linscal10(t0,konl,t0l,null,shp)
              call linscal10(t1,konl,t1l,null,shp)
            else
              do i1=1,nope
                t0l=t0l+shp(4,i1)*t0(konl(i1))
                t1l=t1l+shp(4,i1)*t1(konl(i1))
              enddo
            endif
          elseif(ithermal(1).ge.2) then
            if(lakonl(4:5).eq.'8 ') then
              do i1=1,nope
                t0l=t0l+t0(konl(i1))/8.d0
                t1l=t1l+vold(0,konl(i1))/8.d0
              enddo
            elseif(lakonl(4:6).eq.'20 ')then
              nopered=20
              call lintemp_th0(t0,konl,nopered,kk,t0l,mi)
              call lintemp_th1(vold,konl,nopered,kk,t1l,mi)
            elseif(lakonl(4:6).eq.'10T') then
              call linscal10(t0,konl,t0l,null,shp)
              call linscal10(vold,konl,t1l,mi(2),shp)
            else
              do i1=1,nope
                t0l=t0l+shp(4,i1)*t0(konl(i1))
                t1l=t1l+shp(4,i1)*vold(0,konl(i1))
              enddo
            endif
          endif
          tt=t1l-t0l
!     
          kode=nelcon(1,imat)
!     
!     material data
!     
          istiff=1
          call materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &         imat,amat,iorien,coords,orab,ntmat_,stiff,rho,
     &         nelem,ithermal,alzero,mattyp,t0l,t1l,
     &         ihyper,istiff,elconloc,eth,kode,plicon,
     &         nplicon,plkcon,nplkcon,npmat_,
     &         plconloc,mi(1),dtime,kk,
     &         xstiff,ncmat_,iperturb)
!     
          if(mattyp.eq.1) then
            e=stiff(1)
            un=stiff(2)
            wavspd=max(wavspd,
     &           dsqrt((e*(1-un))/((1+un)*(1-2*un)*rho)))
          elseif(mattyp.eq.2) then
!     
!           single crystal
!     
            iorth=1
            call anisomaxwavspd(stiff,rho,iorth,wavspd)
          elseif(mattyp.eq.3) then
            iorth=0
            call anisomaxwavspd(stiff,rho,iorth,wavspd)
          endif
!     
          wavespeed(imat)=wavspd
!
!         volume=area*h*volfac
!
          volume=weight*xsj
          volfac=1.d0
          if((lakon(nelem)(4:5).eq.'10').or.
     &         (lakon(nelem)(4:4).eq.'4')) then
            volfac=3.d0
          endif
!     
          hmin=1.d30
!     
!     DO over sides
          do ig=1,nfaces
            if(lakon(nelem)(4:4).eq.'6')then
              if(ig.le.2)then
                nopes=3
              else
                nopes=4
                volfac=2.d0
              endif
            elseif(lakon(nelem)(4:5).eq.'15')then
              if(ig.le.2)then
                nopes=6
              else
                nopes=8
                volfac=2.d0
              endif
            endif
!
            if((nope.eq.20).or.(nope.eq.8))then
              do i=1,nopes
                do j=1,3
                  xl2(j,i)=co(j,konl(ifaceq(i,ig)))
                enddo
              enddo
            elseif((nope.eq.10).or.(nope.eq.4))then
              do i=1,nopes
                do j=1,3
                  xl2(j,i)=co(j,konl(ifacet(i,ig)))
                enddo
              enddo
            else
              do i=1,nopes
                do j=1,3
                  xl2(j,i)=co(j,konl(ifacew(i,ig)))
                enddo
              enddo
            endif
!     
            if((nopes.eq.4).or.(nopes.eq.8))then
              xi=0.d0
              et=0.d0
              weight=4.d0
            else
              xi=1.d0/3.d0
              et=1.d0/3.d0
              weight=0.5d0
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
            area=weight*dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &           xsj2(3)*xsj2(3))

            hmin=min(hmin,(volume*volfac/area))
!     
          enddo
!     ENDDO over sides
!     
!     scaling of element: time increment required by element
!
          smscale(nelem)=critom*hmin*elemfac/(2.d0*wavspd)
!     
!     smallest dtvol
          dtvol=min(dtvol,smscale(nelem))
        endif                   !endif(lakonl(1:2).ne.'ES')
!     
      enddo
!     ** ENDDO per element
!     
      do i=1,nmat
        write(*,*) 'Wave Speed for mat. ',i,wavespeed(i)
      enddo
      write(*,*)
!     
!     ------------Mass Scaling ----------------------------------------
!     mscalmethod=1: selective mass scaling SMS
!
      if((dtvol.lt.dtset/safefac))then
        dtset=dtset/safefac
        mscalmethod=1
!     
        ilen=index(jobnamef(1),' ')-1
        fn=jobnamef(1)(1:ilen)//'_WarnElementMassScaled.nam'
        open(40,file=fn,status='unknown')
        write(40,*) '*ELSET,ELSET=MassScaled'
!     
        do nelem=1,ne0
          if(ipkon(nelem).lt.0) cycle
!     scaling matrix
          if(smscale(nelem).ge.dtset)then
            smscale(nelem)=0.d0
          else
            smscale(nelem)=((dtset/smscale(nelem))**2)-1 !beta
            icount=icount+1
            write(40,*) nelem
          endif
        enddo
!     
        if(icount.gt.0) then
          write(*,*) 'Selective Mass Scaling is active'
          write(*,*)
          write(*,*) 'Scaling factor of time increment:',dtset/dtvol
          write(*,*) 'Overall mass is not changed:'
          write(*,*) 'Manipulation of M matrix by beta (maximum) =',
     &         (((dtset/dtvol)**2)-1)
          write(*,*) 'In total ',icount,'elements were scaled'
          write(*,*) 
          write(*,*) '*INFO in calcstabletimeincvol:'
          write(*,*) '      scaled elements are stored in file'
          write(*,*) '      ',fn(1:ilen+26)
          write(*,*) '      This file can be loaded into'
          write(*,*) '      an active cgx-session by typing'
          write(*,*)
     &         '      read ',fn(1:ilen+26),' inp'
          write(*,*)
          close(40)
        else
          close(40,status='delete')
        endif
!     
        if(dtset/dtvol.gt.1.d10) then
          write(*,*) '*ERROR in calcstabletimeincvol:'
          write(*,*) '       time increment is too small;'
          write(*,*) '       divergence'
          call exit(201)
        endif
!     
        dtset=dtset*safefac
      endif                     !endif (dtvol.lt.dtset)
!     
      dtvol=dtvol*safefac
!     
      return
      end
