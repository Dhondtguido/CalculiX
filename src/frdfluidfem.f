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
      subroutine frdfluidfem(co,nk,kon,ipkon,lakon,ne,v,vold,
     &     kode,time,ielmat,matname,nnstep,
     &     physcon,filab,inomat,ntrans,inotr,trab,mi,stn,qfn,istep,
     &     sa,compressible,nactdoh,yy,jyy,ithermal,shcon,nshcon,
     &     rhcon,nrhcon,ntmat_,vcon,depth,xg,nodfreesurf,dt,shockcoef,
     &     iexplicit,nkold,nelold)
!     
!     stores the results in frd format
!     
      implicit none
!     
      character*1 c
      character*3 m1,m2,m3,m4,m5
      character*5 p0,p1,p2,p3,p4,p5,p6,p8,p10,p11,p12
      character*8 lakon(*),date,newclock,fmat
      character*10 clock
      character*20 newdate
      character*80 matname(*)
      character*87 filab(*)
      character*132 text
!     
      integer kon(*),nk,ne,kode,i,j,k,ipkon(*),indexe,inomat(*),mi(*),
     &     one,ielmat(mi(3),*),null,nnstep,inotr(2,*),ntrans,nout,
     &     istep,compressible,nactdoh(*),jyy(*),imat,ithermal(*),
     &     nshcon(*),ntmat_,nrhcon(*),nodfreesurf(*),iexplicit,nkold(*),
     &     nelold(*)
!     
      real*8 co(3,*),v(nk,0:mi(2)),time,vold(0:mi(2),*),
     &     vcon(nk,0:mi(2)),
     &     stn(6,*),qfn(3,*),pi,oner,physcon(*),depth(*),xg(3),dt(*),
     &     trab(7,*),a(3,3),sa(*),yy(*),xn(3),xt(3),ddn,ddt,temp,dvi,
     &     rho,shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*),tau,yplus,
     &     utau,uplus,shockcoef
!     
      save nout
!     
      kode=kode+1
      pi=4.d0*datan(1.d0)
!     
      c='C'
!     
      m1=' -1'
      m2=' -2'
      m3=' -3'
      m4=' -4'
      m5=' -5'
!     
      p0='    0'
      p1='    1'
      p2='    2'
      p3='    3'
      p4='    4'
      p5='    5'
      p6='    6'
      p8='    8'
      p10='   10'
      p11='   11'
      p12='   12'
!     
      if(time.le.0.d0) then
        fmat(1:8)='(e12.5) '
      elseif((dlog10(time).ge.0.d0).and.(dlog10(time).lt.11.d0)) then
        fmat(1:5)='(f12.'
        write(fmat(6:7),'(i2)') 11-int(dlog10(time)+1.d0)
        fmat(8:8)=')'
      else
        fmat(1:8)='(e12.5) '
      endif
!     
      null=0
      one=1
      oner=1.d0
!     
      if(kode.eq.1) then
!     
        write(13,'(a5,a1)') p1,c
        call date_and_time(date,clock)
        newdate(1:20)='                    '
        newdate(1:2)=date(7:8)
        newdate(3:3)='.'
        if(date(5:6).eq.'01') then
          newdate(4:11)='january.'
          newdate(12:15)=date(1:4)
        elseif(date(5:6).eq.'02') then
          newdate(4:12)='february.'
          newdate(13:16)=date(1:4)
        elseif(date(5:6).eq.'03') then
          newdate(4:9)='march.'
          newdate(10:13)=date(1:4)
        elseif(date(5:6).eq.'04') then
          newdate(4:9)='april.'
          newdate(10:13)=date(1:4)
        elseif(date(5:6).eq.'05') then
          newdate(4:7)='may.'
          newdate(8:11)=date(1:4)
        elseif(date(5:6).eq.'06') then
          newdate(4:8)='june.'
          newdate(9:12)=date(1:4)
        elseif(date(5:6).eq.'07') then
          newdate(4:8)='july.'
          newdate(9:12)=date(1:4)
        elseif(date(5:6).eq.'08') then
          newdate(4:10)='august.'
          newdate(11:14)=date(1:4)
        elseif(date(5:6).eq.'09') then
          newdate(4:13)='september.'
          newdate(14:17)=date(1:4)
        elseif(date(5:6).eq.'10') then
          newdate(4:11)='october.'
          newdate(12:15)=date(1:4)
        elseif(date(5:6).eq.'11') then
          newdate(4:12)='november.'
          newdate(13:16)=date(1:4)
        elseif(date(5:6).eq.'12') then
          newdate(4:12)='december.'
          newdate(13:16)=date(1:4)
        endif
        newclock(1:2)=clock(1:2)
        newclock(3:3)=':'
        newclock(4:5)=clock(3:4)
        newclock(6:6)=':'
        newclock(7:8)=clock(5:6)
        write(13,'(a5,''UUSER'')') p1
        write(13,'(a5,''UDATE'',14x,a20)') p1,newdate
        write(13,'(a5,''UTIME'',14x,a8)') p1,newclock
        write(13,'(a5,''UHOST'')') p1
        write(13,'(a5,''UPGM               CalculiX'')') p1
        write(13,'(a5,''UDIR'')') p1
        write(13,'(a5,''UDBN'')') p1
!     
!     storing the coordinates of the nodes
!     
        write(13,'(a5,a1,67x,i1)') p2,c,one
!     
        nout=0
        do i=1,nk
          if(inomat(i).le.0) cycle
          nout=nout+1
          write(13,100) m1,nkold(i),(co(j,i),j=1,3)
        enddo
!     
        write(13,'(a3)') m3
!     
!     storing the element topology
!     
        write(13,'(a5,a1,67x,i1)') p3,c,one
!     
        do i=1,ne
!     
          if(ipkon(i).lt.0) cycle
          indexe=ipkon(i)
          if(lakon(i)(4:4).eq.'2') then
            if((lakon(i)(7:7).eq.' ').or.
     &           (lakon(i)(7:7).eq.'H')) then
              write(13,'(a3,i10,3a5)') m1,nelold(i),p4,p0,
     &             matname(ielmat(1,i))(1:5)
              write(13,'(a3,10i10)') m2,(nkold(kon(indexe+j)),j=1,10)
              write(13,'(a3,10i10)') m2,(nkold(kon(indexe+j)),j=11,12),
     &             (nkold(kon(indexe+j)),j=17,19),nkold(kon(indexe+20)),
     &             (nkold(kon(indexe+j)),j=13,16)
            elseif(lakon(i)(7:7).eq.'B') then
              write(13,'(a3,i10,3a5)')m1,nelold(i),p12,p0,
     &             matname(ielmat(1,i))(1:5)
              write(13,'(a3,3i10)') m2,nkold(kon(indexe+21)),
     &             nkold(kon(indexe+23)),
     &             nkold(kon(indexe+22))
            else
              write(13,'(a3,i10,3a5)')m1,nelold(i),p10,p0,
     &             matname(ielmat(1,i))(1:5)
              write(13,'(a3,8i10)') m2,(nkold(kon(indexe+20+j)),j=1,8)
            endif
          elseif(lakon(i)(4:4).eq.'8') then
            write(13,'(a3,i10,3a5)') m1,nelold(i),p1,p0,
     &           matname(ielmat(1,i))(1:5)
            write(13,'(a3,8i10)') m2,(nkold(kon(indexe+j)),j=1,8)
          elseif(lakon(i)(4:5).eq.'10') then
            write(13,'(a3,i10,3a5)') m1,nelold(i),p6,p0,
     &           matname(ielmat(1,i))(1:5)
            write(13,'(a3,10i10)') m2,(nkold(kon(indexe+j)),j=1,10)
          elseif(lakon(i)(4:4).eq.'4') then
            write(13,'(a3,i10,3a5)') m1,nelold(i),p3,p0,
     &           matname(ielmat(1,i))(1:5)
            write(13,'(a3,4i10)') m2,(nkold(kon(indexe+j)),j=1,4)
          elseif(lakon(i)(4:5).eq.'15') then
            if((lakon(i)(7:7).eq.' ')) then
              write(13,'(a3,i10,3a5)') m1,nelold(i),p5,p0,
     &             matname(ielmat(1,i))(1:5)
              write(13,'(a3,10i10)') m2,(nkold(kon(indexe+j)),j=1,9),
     &             nkold(kon(indexe+13))
              write(13,'(a3,5i10)') m2,(nkold(kon(indexe+j)),j=14,15),
     &             (nkold(kon(indexe+j)),j=10,12)
            else
              write(13,'(a3,i10,3a5)') m1,nelold(i),p8,p0,
     &             matname(ielmat(1,i))(1:5)
              write(13,'(a3,6i10)') m2,(nkold(kon(indexe+15+j)),j=1,6)
            endif
          elseif(lakon(i)(4:4).eq.'6') then
            write(13,'(a3,i10,3a5)') m1,nelold(i),p2,p0,
     &           matname(ielmat(1,i))(1:5)
            write(13,'(a3,6i10)') m2,(nkold(kon(indexe+j)),j=1,6)
          elseif(lakon(i)(1:1).eq.'D') then
            if((nkold(kon(indexe+1)).eq.0).or.
     &           (nkold(kon(indexe+3)).eq.0)) cycle
            write(13,'(a3,i10,3a5)')m1,nelold(i),p12,p0,
     &           matname(ielmat(1,i))(1:5)
            write(13,'(a3,3i10)') m2,nkold(kon(indexe+1)),
     &           nkold(kon(indexe+3)),
     &           nkold(kon(indexe+2))
          elseif(lakon(i)(1:1).eq.'E') then
            write(13,'(a3,i10,3a5)')m1,nelold(i),p11,p0,
     &           matname(ielmat(1,i))(1:5)
            write(13,'(a3,2i10)') m2,(nkold(kon(indexe+j)),j=1,2)
          endif
!     
        enddo
!     
        write(13,'(a3)') m3
!     
      endif
!     
!     storing the velocities in the nodes
!     
      if((nnstep.eq.1).or.(nnstep.eq.3)) then
        do i=37,132
          text(i:i)=' '
        enddo
        if(iexplicit.eq.1) then
          text(1:16)='    1PSHOCKCOEF='
          write(text(17:36),'(f20.13)') shockcoef
          write(13,'(a132)') text
        endif
        text='    1PSTEP'
        write(text(25:36),'(i12)') kode
        write(13,'(a132)') text
!     
        text=
     & '  100CL       .00000E+00                                 3    1'
        text(75:75)='1'
        write(text(25:36),'(i12)') nout
        write(text(8:12),'(i5)') 100+kode
        write(text(13:24),fmat) time
        write(text(59:63),'(i5)') kode
        write(13,'(a132)') text
        text=' -4  DFVEL       4    1'
        write(13,'(a132)') text
        text=' -5  V1          1    2    1    0'
        write(13,'(a132)') text
        text=' -5  V2          1    2    2    0'
        write(13,'(a132)') text
        text=' -5  V3          1    2    3    0'
        write(13,'(a132)') text
        text=' -5  ALL         1    2    0    0    1ALL'
        write(13,'(a132)') text
!     
        do i=1,nk
          if(inomat(i).le.0) cycle
          write(13,100) m1,nkold(i),(v(i,j),j=1,3)
        enddo
!     
        write(13,'(a3)') m3
!     
!     
!     storing the static pressure in the nodes
!     
      elseif(nnstep.eq.2) then
        do i=37,132
          text(i:i)=' '
        enddo
        if(iexplicit.eq.1) then
          text(1:16)='    1PSHOCKCOEF='
          write(text(17:36),'(f20.13)') shockcoef
          write(13,'(a132)') text
        endif
        text='    1PSTEP'
        write(text(25:36),'(i12)') kode
        write(13,'(a132)') text
!     
        text=
     & '  100CL       .00000E+00                                 3    1'
        text(75:75)='1'
        write(text(25:36),'(i12)') nout
        write(text(8:12),'(i5)') 100+kode
        write(text(13:24),fmat) time
        write(text(59:63),'(i5)') kode
        write(13,'(a132)') text
        text=' -4  DDENSIT     1    1'
        write(13,'(a132)') text
        text=' -5  DRHO        1    1    0    0'
        write(13,'(a132)') text
!     
        do i=1,nk
          if(inomat(i).le.0) cycle
          write(13,100) m1,nkold(i),v(i,4)
        enddo
!     
        write(13,'(a3)') m3
!     
!     storing the static temperature in the nodes
!     
      elseif(nnstep.eq.4) then
        do i=37,132
          text(i:i)=' '
        enddo
        if(iexplicit.eq.1) then
          text(1:16)='    1PSHOCKCOEF='
          write(text(17:36),'(f20.13)') shockcoef
          write(13,'(a132)') text
        endif
        text='    1PSTEP'
        write(text(25:36),'(i12)') kode
        write(13,'(a132)') text
!     
        text=
     & '  100CL       .00000E+00                                 3    1'
        text(75:75)='1'
        write(text(25:36),'(i12)') nout
        write(text(8:12),'(i5)') 100+kode
        write(text(13:24),fmat) time
        write(text(59:63),'(i5)') kode
        write(13,'(a132)') text
        text=' -4  DENERGY     1    1'
        write(13,'(a132)') text
        text=' -5  DRE         1    1    0    0'
        write(13,'(a132)') text
!     
        do i=1,nk
          if(inomat(i).le.0) cycle
          write(13,100) m1,nkold(i),v(i,0)
        enddo
!     
        write(13,'(a3)') m3
!     
!     storing the turbulence parameters in the nodes
!     
      elseif(nnstep.eq.5) then
        do i=37,132
          text(i:i)=' '
        enddo
        if(iexplicit.eq.1) then
          text(1:16)='    1PSHOCKCOEF='
          write(text(17:36),'(f20.13)') shockcoef
          write(13,'(a132)') text
        endif
        text='    1PSTEP'
        write(text(25:36),'(i12)') kode
        write(13,'(a132)') text
!     
        text=
     & '  100CL       .00000E+00                                 3    1'
        text(75:75)='1'
        write(text(25:36),'(i12)') nout
        write(text(8:12),'(i5)') 100+kode
        write(text(13:24),fmat) time
        write(text(59:63),'(i5)') kode
        write(13,'(a132)') text
        text=' -4  DTURB1      1    1'
        write(13,'(a132)') text
        text=' -5  K           1    1    0    0'
        write(13,'(a132)') text
!     
        do i=1,nk
          if(inomat(i).le.0) cycle
          write(13,100) m1,nkold(i),v(i,5)
        enddo
!     
        write(13,'(a3)') m3
!     
        do i=37,132
          text(i:i)=' '
        enddo
        if(iexplicit.eq.1) then
          text(1:16)='    1PSHOCKCOEF='
          write(text(17:36),'(f20.13)') shockcoef
          write(13,'(a132)') text
        endif
        text='    1PSTEP'
        write(text(25:36),'(i12)') kode
        write(13,'(a132)') text
!     
        text=
     & '  100CL       .00000E+00                                 3    1'
        text(75:75)='1'
        write(text(25:36),'(i12)') nout
        write(text(8:12),'(i5)') 100+kode
        write(text(13:24),fmat) time
        write(text(59:63),'(i5)') kode
        write(13,'(a132)') text
        text=' -4  DTURB2      1    1'
        write(13,'(a132)') text
        text=' -5  OM          1    1    0    0'
        write(13,'(a132)') text
!     
        do i=1,nk
          if(inomat(i).le.0) cycle
          write(13,100) m1,nkold(i),v(i,6)
        enddo
!     
        write(13,'(a3)') m3
!     
      elseif(nnstep.eq.6) then
!     
        if(iexplicit.eq.1) then
          do i=37,132
            text(i:i)=' '
          enddo
          if(iexplicit.eq.1) then
            text(1:16)='    1PSHOCKCOEF='
            write(text(17:36),'(f20.13)') shockcoef
            write(13,'(a132)') text
          endif
          text='    1PSTEP'
          write(text(25:36),'(i12)') kode
          write(13,'(a132)') text
!     
          text=
     & '  100CL       .00000E+00                                 3    1'
          text(75:75)='1'
          write(text(25:36),'(i12)') nout
          write(text(8:12),'(i5)') 100+kode
          write(text(13:24),fmat) time
          write(text(59:63),'(i5)') kode
          write(13,'(a132)') text
          text=' -4  SA3DF       1    1'
          write(13,'(a132)') text
          text=' -5  SA          1    1    0    0'
          write(13,'(a132)') text
!     
          do i=1,nk
            if(nactdoh(i).le.0) then
              write(13,100) m1,nkold(i),1.d0*null
              cycle
            endif
            if(inomat(i).le.0) cycle
            write(13,100) m1,nkold(i),sa(nactdoh(i))
          enddo
!     
          write(13,'(a3)') m3
        endif
!     
        if(filab(34)(1:4).eq.'VF  ') then
          do i=37,132
            text(i:i)=' '
          enddo
          if(iexplicit.eq.1) then
            text(1:16)='    1PSHOCKCOEF='
            write(text(17:36),'(f20.13)') shockcoef
            write(13,'(a132)') text
          endif
          text='    1PSTEP'
          write(text(25:36),'(i12)') kode
          write(13,'(a132)') text
!     
          text=
     & '  100CL       .00000E+00                                 3    1'
          text(75:75)='1'
          write(text(25:36),'(i12)') nout
          write(text(8:12),'(i5)') 100+kode
          write(text(13:24),fmat) time
          write(text(59:63),'(i5)') kode
          write(13,'(a132)') text
          text=' -4  V3DF        4    1'
          write(13,'(a132)') text
          text=' -5  V1          1    2    1    0'
          write(13,'(a132)') text
          text=' -5  V2          1    2    2    0'
          write(13,'(a132)') text
          text=' -5  V3          1    2    3    0'
          write(13,'(a132)') text
          text=' -5  ALL         1    2    0    0    1ALL'
          write(13,'(a132)') text
!     
          if((ntrans.eq.0).or.(filab(34)(6:6).eq.'G')) then
            do i=1,nk
              if(inomat(i).le.0) cycle
              write(13,100) m1,nkold(i),(vold(j,i),j=1,3)
c     write(13,100) m1,nkold(i),(vcon(i,j),j=1,3)
            enddo
          else
            do i=1,nk
              if(inomat(i).le.0) cycle
              if(inotr(1,i).eq.0) then
                write(13,100) m1,nkold(i),(vold(j,i),j=1,3)
              else
                call transformatrix(trab(1,inotr(1,i)),co(1,i),a)
                write(13,100) m1,nkold(i),
     &               vold(1,i)*a(1,1)+vold(2,i)*a(2,1)+vold(3,i)*a(3,1),
     &               vold(1,i)*a(1,2)+vold(2,i)*a(2,2)+vold(3,i)*a(3,2),
     &               vold(1,i)*a(1,3)+vold(2,i)*a(2,3)+vold(3,i)*a(3,3)
              endif
            enddo
          endif
!     
          write(13,'(a3)') m3
        endif
!     
        if(filab(35)(1:4).eq.'PSF ') then
          do i=37,132
            text(i:i)=' '
          enddo
          if(iexplicit.eq.1) then
            text(1:16)='    1PSHOCKCOEF='
            write(text(17:36),'(f20.13)') shockcoef
            write(13,'(a132)') text
          endif
          text='    1PSTEP'
          write(text(25:36),'(i12)') kode
          write(13,'(a132)') text
!     
          text=
     & '  100CL       .00000E+00                                 3    1'
          text(75:75)='1'
          write(text(25:36),'(i12)') nout
          write(text(8:12),'(i5)') 100+kode
          write(text(13:24),fmat) time
          write(text(59:63),'(i5)') kode
          write(13,'(a132)') text
          text=' -4  PS3DF       1    1'
          write(13,'(a132)') text
          text=' -5  PS          1    1    0    0'
          write(13,'(a132)') text
!     
          do i=1,nk
            if(inomat(i).le.0) cycle
            write(13,100) m1,nkold(i),vold(4,i)
c     write(13,100) m1,nkold(i),vcon(i,4)
          enddo
!     
          write(13,'(a3)') m3
        endif
!     
        if(filab(49)(1:4).eq.'DEPF') then
          do i=37,132
            text(i:i)=' '
          enddo
          if(iexplicit.eq.1) then
            text(1:16)='    1PSHOCKCOEF='
            write(text(17:36),'(f20.13)') shockcoef
            write(13,'(a132)') text
          endif
          text='    1PSTEP'
          write(text(25:36),'(i12)') kode
          write(13,'(a132)') text
!     
          text=
     & '  100CL       .00000E+00                                 3    1'
          text(75:75)='1'
          write(text(25:36),'(i12)') nout/2
          write(text(8:12),'(i5)') 100+kode
          write(text(13:24),fmat) time
          write(text(59:63),'(i5)') kode
          write(13,'(a132)') text
          text=' -4  DISP        4    1'
          write(13,'(a132)') text
          text=' -5  D1          1    2    1    0'
          write(13,'(a132)') text
          text=' -5  D2          1    2    2    0'
          write(13,'(a132)') text
          text=' -5  D3          1    2    3    0'
          write(13,'(a132)') text
          text=' -5  ALL         1    2    0    0    1ALL'
          write(13,'(a132)') text
!     
          do i=1,nk
            if(inomat(i).le.0) cycle
            if(nodfreesurf(i).ne.1) cycle
            write(13,100) m1,nkold(i),
     &           (-(vcon(i,4)-depth(i))*xg(j),j=1,3)
c     write(13,100) m1,nkold(i),(depth(i),j=1,3)
          enddo
!     
          write(13,'(a3)') m3
        endif
!     
        if(filab(50)(1:4).eq.'DTF ') then
          do i=37,132
            text(i:i)=' '
          enddo
          if(iexplicit.eq.1) then
            text(1:16)='    1PSHOCKCOEF='
            write(text(17:36),'(f20.13)') shockcoef
            write(13,'(a132)') text
          endif
          text='    1PSTEP'
          write(text(25:36),'(i12)') kode
          write(13,'(a132)') text
!     
          text=
     & '  100CL       .00000E+00                                 3    1'
          text(75:75)='1'
          write(text(25:36),'(i12)') nout
          write(text(8:12),'(i5)') 100+kode
          write(text(13:24),fmat) time
          write(text(59:63),'(i5)') kode
          write(13,'(a132)') text
          text=' -4  DTIMF       1    1'
          write(13,'(a132)') text
          text=' -5  DT          1    1    0    0'
          write(13,'(a132)') text
!     
          do i=1,nk
            if(inomat(i).le.0) cycle
            write(13,100) m1,nkold(i),dt(i)
          enddo
!     
          write(13,'(a3)') m3
        endif
!     
        if(filab(36)(1:4).eq.'TSF ') then
          do i=37,132
            text(i:i)=' '
          enddo
          if(iexplicit.eq.1) then
            text(1:16)='    1PSHOCKCOEF='
            write(text(17:36),'(f20.13)') shockcoef
            write(13,'(a132)') text
          endif
          text='    1PSTEP'
          write(text(25:36),'(i12)') kode
          write(13,'(a132)') text
!     
          text=
     & '  100CL       .00000E+00                                 3    1'
          text(75:75)='1'
          write(text(25:36),'(i12)') nout
          write(text(8:12),'(i5)') 100+kode
          write(text(13:24),fmat) time
          write(text(59:63),'(i5)') kode
          write(13,'(a132)') text
          text=' -4  TS3DF       1    1'
          write(13,'(a132)') text
          text=' -5  TS          1    1    0    0'
          write(13,'(a132)') text
!     
          do i=1,nk
            if(inomat(i).le.0) cycle
            write(13,100) m1,nkold(i),vold(0,i)
c     write(13,100) m1,nkold(i),vcon(i,0)
          enddo
!     
          write(13,'(a3)') m3
        endif
!     
!     v(*,0) contains kappa (cp/cv)
!     v(*,1) contains the Mach number
!     
        if(filab(23)(1:4).eq.'MACH') then
          do i=37,132
            text(i:i)=' '
          enddo
          if(iexplicit.eq.1) then
            text(1:16)='    1PSHOCKCOEF='
            write(text(17:36),'(f20.13)') shockcoef
            write(13,'(a132)') text
          endif
          text='    1PSTEP'
          write(text(25:36),'(i12)') kode
          write(13,'(a132)') text
!     
          text=
     & '  100CL       .00000E+00                                 3    1'
          text(75:75)='1'
          write(text(25:36),'(i12)') nout
          write(text(8:12),'(i5)') 100+kode
          write(text(13:24),fmat) time
          write(text(59:63),'(i5)') kode
          write(13,'(a132)') text
          text=' -4  M3DF        1    1'
          write(13,'(a132)') text
          text=' -5  MACH        1    1    0    0'
          write(13,'(a132)') text
!     
          do i=1,nk
            if(inomat(i).le.0) cycle
            write(13,100) m1,nkold(i),v(i,1)
          enddo
!     
          write(13,'(a3)') m3
        endif
!     
        if(filab(38)(1:4).eq.'TTF ') then
          do i=37,132
            text(i:i)=' '
          enddo
          if(iexplicit.eq.1) then
            text(1:16)='    1PSHOCKCOEF='
            write(text(17:36),'(f20.13)') shockcoef
            write(13,'(a132)') text
          endif
          text='    1PSTEP'
          write(text(25:36),'(i12)') kode
          write(13,'(a132)') text
!     
          text=
     & '  100CL       .00000E+00                                 3    1'
          text(75:75)='1'
          write(text(25:36),'(i12)') nout
          write(text(8:12),'(i5)') 100+kode
          write(text(13:24),fmat) time
          write(text(59:63),'(i5)') kode
          write(13,'(a132)') text
          text=' -4  TT3DF       1    1'
          write(13,'(a132)') text
          text=' -5  TT          1    1    0    0'
          write(13,'(a132)') text
!     
          do i=1,nk
            if(inomat(i).le.0) cycle
            write(13,100) m1,nkold(i),
     &           vold(0,i)*(1.d0+(v(i,0)-1.d0)/2*v(i,1)**2)
          enddo
!     
          write(13,'(a3)') m3
        endif
!     
        if(filab(37)(1:4).eq.'PTF ') then
          do i=37,132
            text(i:i)=' '
          enddo
          if(iexplicit.eq.1) then
            text(1:16)='    1PSHOCKCOEF='
            write(text(17:36),'(f20.13)') shockcoef
            write(13,'(a132)') text
          endif
          text='    1PSTEP'
          write(text(25:36),'(i12)') kode
          write(13,'(a132)') text
!     
          text=
     & '  100CL       .00000E+00                                 3    1'
          text(75:75)='1'
          write(text(25:36),'(i12)') nout
          write(text(8:12),'(i5)') 100+kode
          write(text(13:24),fmat) time
          write(text(59:63),'(i5)') kode
          write(13,'(a132)') text
          text=' -4  PT3DF       1    1'
          write(13,'(a132)') text
          text=' -5  PT          1    1    0    0'
          write(13,'(a132)') text
!     
          do i=1,nk
            if(inomat(i).le.0) cycle
            write(13,100) m1,nkold(i),vold(4,i)*
     &         (1.d0+(v(i,0)-1.d0)/2*v(i,1)**2)**(v(i,0)/(v(i,0)-1.d0))
          enddo
!     
          write(13,'(a3)') m3
        endif
!     
!     storing the total stresses in the nodes
!     
        if(filab(39)(1:4).eq.'SF  ') then
          do i=37,132
            text(i:i)=' '
          enddo
          if(iexplicit.eq.1) then
            text(1:16)='    1PSHOCKCOEF='
            write(text(17:36),'(f20.13)') shockcoef
            write(13,'(a132)') text
          endif
          text='    1PSTEP'
          write(text(25:36),'(i12)') kode
          write(text(49:60),'(i12)') istep
          write(13,'(a132)') text
!     
          text=
     & '  100CL       .00000E+00                                 3    1'
          text(75:75)='1'
          write(text(25:36),'(i12)') nout
          write(text(8:12),'(i5)') 100+kode
          write(text(13:24),fmat) time
          write(text(59:63),'(i5)') kode
          write(13,'(a132)') text
          text=' -4  STRESS      6    1'
          write(13,'(a132)') text
          text=' -5  SXX         1    4    1    1'
          write(13,'(a132)') text
          text=' -5  SYY         1    4    2    2'
          write(13,'(a132)') text
          text=' -5  SZZ         1    4    3    3'
          write(13,'(a132)') text
          text=' -5  SXY         1    4    1    2'
          write(13,'(a132)') text
          text=' -5  SYZ         1    4    2    3'
          write(13,'(a132)') text
          text=' -5  SZX         1    4    3    1'
          write(13,'(a132)') text
          do i=1,nk
            if(inomat(i).le.0) cycle
            write(13,100) m1,nkold(i),(stn(j,i)-vold(4,i),j=1,3),
     &           stn(4,i),stn(6,i),stn(5,i)
          enddo
          write(13,'(a3)') m3
        endif
!     
!     storing the viscous stresses in the nodes
!     
        if(filab(41)(1:4).eq.'SVF ') then
          do i=37,132
            text(i:i)=' '
          enddo
          if(iexplicit.eq.1) then
            text(1:16)='    1PSHOCKCOEF='
            write(text(17:36),'(f20.13)') shockcoef
            write(13,'(a132)') text
          endif
          text='    1PSTEP'
          write(text(25:36),'(i12)') kode
          write(13,'(a132)') text
!     
          text=
     & '  100CL       .00000E+00                                 3    1'
          text(75:75)='1'
          write(text(25:36),'(i12)') nout
          write(text(8:12),'(i5)') 100+kode
          write(text(13:24),fmat) time
          write(text(59:63),'(i5)') kode
          write(13,'(a132)') text
          text=' -4  VSTRES      6    1'
          write(13,'(a132)') text
          text=' -5  SXX         1    4    1    1'
          write(13,'(a132)') text
          text=' -5  SYY         1    4    2    2'
          write(13,'(a132)') text
          text=' -5  SZZ         1    4    3    3'
          write(13,'(a132)') text
          text=' -5  SXY         1    4    1    2'
          write(13,'(a132)') text
          text=' -5  SYZ         1    4    2    3'
          write(13,'(a132)') text
          text=' -5  SZX         1    4    3    1'
          write(13,'(a132)') text
          do i=1,nk
            if(inomat(i).le.0) cycle
            write(13,100) m1,nkold(i),(stn(j,i),j=1,4),
     &           stn(6,i),stn(5,i)
          enddo
          write(13,'(a3)') m3
        endif
!     
!     storing the heat flux in the nodes
!     
        if(filab(40)(1:4).eq.'HFLF') then
          do i=37,132
            text(i:i)=' '
          enddo
          if(iexplicit.eq.1) then
            text(1:16)='    1PSHOCKCOEF='
            write(text(17:36),'(f20.13)') shockcoef
            write(13,'(a132)') text
          endif
          text='    1PSTEP'
          write(text(25:36),'(i12)') kode
          write(13,'(a132)') text
!     
          text=
     & '  100CL       .00000E+00                                 3    1'
          text(75:75)='1'
          write(text(25:36),'(i12)') nout
          write(text(8:12),'(i5)') 100+kode
          write(text(13:24),fmat) time
          write(text(59:63),'(i5)') kode
          write(13,'(a132)') text
          text=' -4  FLUX        4    1'
          write(13,'(a132)') text
          text=' -5  F1          1    2    1    0'
          write(13,'(a132)') text
          text=' -5  F2          1    2    2    0'
          write(13,'(a132)') text
          text=' -5  F3          1    2    3    0'
          write(13,'(a132)') text
          text=' -5  ALL         1    2    0    0    1ALL'
          write(13,'(a132)') text
          do i=1,nk
            if(inomat(i).le.0) cycle
            write(13,100) m1,nkold(i),(qfn(j,i),j=1,3)
          enddo
          write(13,'(a3)') m3
        endif
!     
        if(filab(24)(1:4).eq.'CP  ') then
          do i=37,132
            text(i:i)=' '
          enddo
          if(iexplicit.eq.1) then
            text(1:16)='    1PSHOCKCOEF='
            write(text(17:36),'(f20.13)') shockcoef
            write(13,'(a132)') text
          endif
          text='    1PSTEP'
          write(text(25:36),'(i12)') kode
          write(13,'(a132)') text
!     
          text=
     & '  100CL       .00000E+00                                 3    1'
          text(75:75)='1'
          write(text(25:36),'(i12)') nout
          write(text(8:12),'(i5)') 100+kode
          write(text(13:24),fmat) time
          write(text(59:63),'(i5)') kode
          write(13,'(a132)') text
          text=' -4  CP3DF       1    1'
          write(13,'(a132)') text
          text=' -5  CP          1    1    0    0'
          write(13,'(a132)') text
!     
          do i=1,nk
            if(inomat(i).le.0) cycle
            write(13,100) m1,nkold(i),(vold(4,i)-physcon(6))*2.d0/
     &           (physcon(7)*physcon(5)**2)
          enddo
!     
          write(13,'(a3)') m3
        endif
!     
        if(filab(25)(1:4).eq.'TURB') then
          do i=37,132
            text(i:i)=' '
          enddo
          if(iexplicit.eq.1) then
            text(1:16)='    1PSHOCKCOEF='
            write(text(17:36),'(f20.13)') shockcoef
            write(13,'(a132)') text
          endif
          text='    1PSTEP'
          write(text(25:36),'(i12)') kode
          write(13,'(a132)') text
!     
          text=
     & '  100CL       .00000E+00                                 3    1'
          text(75:75)='1'
          write(text(25:36),'(i12)') nout
          write(text(8:12),'(i5)') 100+kode
          write(text(13:24),fmat) time
          write(text(59:63),'(i5)') kode
          write(13,'(a132)') text
          text=' -4  TURB3DF     5    1'
          write(13,'(a132)') text
          text=' -5  K           1    1    1    0'
          write(13,'(a132)') text
          text=' -5  OM          1    2    2    0'
          write(13,'(a132)') text
          text=' -5  NUT         1    3    3    0'
          write(13,'(a132)') text
          text=' -5  Y+          1    4    4    0'
          write(13,'(a132)') text
          text=' -5  U+          1    5    5    0'
          write(13,'(a132)') text
!     
          do i=1,nk
            if(inomat(i).le.0) cycle
            imat=inomat(i)
!     
!     for the determination of y+ calculate:
!     xn: approximate normal to wall
!     xt: approximate tangent to wall in velocity direction
!     jyy(i)=0 can occur if no solid surfaces were defined in the
!     input deck
!     
            if((i.ne.jyy(i)).and.(jyy(i).ne.0)) then
              k=jyy(i)
              temp=vold(0,k)
              do j=1,3
                xn(j)=co(j,i)-co(j,k)
                xt(j)=vold(j,i)
              enddo
              ddn=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
              ddt=dsqrt(xt(1)*xt(1)+xt(2)*xt(2)+xt(3)*xt(3))
              if(ddt.gt.1.d-30) then
                do j=1,3
                  xn(j)=xn(j)/ddn
                  xt(j)=xt(j)/ddt
                enddo
              else
                do j=1,3
                  xn(j)=xn(j)/ddn
                  xt(j)=0.d0
                enddo
              endif
!     
!     calculate shear stress parallel to wall
!     
              tau=xt(1)*(stn(1,k)*xn(1)+stn(4,k)*xn(2)+stn(5,k)*xn(3))+
     &             xt(2)*(stn(4,k)*xn(1)+stn(2,k)*xn(2)+stn(6,k)*xn(3))+
     &             xt(3)*(stn(5,k)*xn(1)+stn(6,k)*xn(2)+stn(3,k)*xn(3))
              tau=dabs(tau)
!     
!     calculate density and dynamic viscosity
!     
              if(compressible.eq.1) then
                rho=vold(4,k)/(shcon(3,1,imat)*temp)
              else
                call materialdata_rho(rhcon,nrhcon,imat,rho,
     &               temp,ntmat_,ithermal)
              endif
              call materialdata_dvifem(imat,ntmat_,temp,shcon,nshcon,
     &             dvi)
!     
              utau=dsqrt(tau/rho)
              yplus=rho*utau*yy(i)/dvi
              uplus=ddt/utau
!     
c     write(*,*) 'frdfluidfem ',nkold(i),k,tau,yy(i),yplus
            else
!     
!     node at wall
!     
              yplus=0.d0
              uplus=0.d0
            endif
            write(13,100) m1,nkold(i),vold(5,i),vold(6,i),
     &           vold(5,i)/vold(6,i),yplus,uplus
          enddo
!     
          write(13,'(a3)') m3
        endif
      endif
!     
 100  format(a3,i10,1p,6e12.5)
!     
      return
      end
