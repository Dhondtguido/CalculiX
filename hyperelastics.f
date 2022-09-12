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
      subroutine hyperelastics(inpc,textpart,elcon,nelcon,
     &  nmat,ntmat_,ncmat_,irstrt,istep,istat,n,iperturb,iline,ipol,
     &  inl,ipoinp,inp,ipoinpc,ier)
!
!     reading the input deck: *HYPERELASTIC
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nelcon(2,*),nmat,ntmat,ntmat_,istep,istat,ipoinpc(0:*),
     &  n,key,i,j,k,ityp,iperturb(*),iend,jcoef(3,14),ncmat_,irstrt(*),
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),ier
!
      real*8 elcon(0:ncmat_,ntmat_,*),um
!
!     jcoef indicates for each hyperelastic model the position of
!     the compressibility coefficients in the field elcon (max. 3
!     positions per model)
!
      data jcoef /3,0,0,3,0,0,2,0,0,3,0,0,5,6,0,7,8,9,3,0,0,
     &            6,7,0,10,11,12,2,0,0,3,4,0,4,5,6,5,0,0,4,5,6/
!
      ntmat=0
      iperturb(1)=3
      iperturb(2)=1
      write(*,*) '*INFO reading *HYPERELASTIC: nonlinear geometric'
      write(*,*) '      effects are turned on'
      write(*,*)
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) 
     &       '*ERROR reading *HYPERELASTIC: *HYPERELASTIC should be'
         write(*,*) '  placed before all step definitions'
         ier=1
         return
      endif
!
      if(nmat.eq.0) then
         write(*,*) 
     &        '*ERROR reading *HYPERELASTIC: *HYPERELASTIC should be'
         write(*,*) '  preceded by a *MATERIAL card'
         ier=1
         return
      endif
!
      ityp=-7
!
      do i=2,n
         if(textpart(i)(1:12).eq.'ARRUDA-BOYCE') then
            ityp=-1
         elseif(textpart(i)(1:13).eq.'MOONEY-RIVLIN') then
            ityp=-2
         elseif(textpart(i)(1:8).eq.'NEOHOOKE') then
            ityp=-3
         elseif(textpart(i)(1:5).eq.'OGDEN') then
            ityp=-4
         elseif(textpart(i)(1:10).eq.'POLYNOMIAL') then
            ityp=-7
         elseif(textpart(i)(1:17).eq.'REDUCEDPOLYNOMIAL') then
            ityp=-10
         elseif(textpart(i)(1:11).eq.'VANDERWAALS') then
            ityp=-13
         elseif(textpart(i)(1:4).eq.'YEOH') then
            ityp=-14
         elseif(textpart(i)(1:2).eq.'N=') then
            if(textpart(i)(3:3).eq.'1') then
            elseif(textpart(i)(3:3).eq.'2') then
               if(ityp.eq.-4) then
                  ityp=-5
               elseif(ityp.eq.-7) then
                  ityp=-8
               elseif(ityp.eq.-10) then
                  ityp=-11
               else
                  write(*,*) '*WARNING reading *HYPERELASTIC: N=2 is not
     & applicable for this material type; '
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*HYPERELASTIC%",ier)
                  return
               endif
            elseif(textpart(i)(3:3).eq.'3') then
               if(ityp.eq.-4) then
                  ityp=-6
               elseif(ityp.eq.-7) then
                  ityp=-9
               elseif(ityp.eq.-10) then
                  ityp=-12
               else
                  write(*,*) '*WARNING reading *HYPERELASTIC: N=3 is not
     & applicable for this material type; '
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*HYPERELASTIC%",ier)
                  return
               endif
            else
               write(*,*) '*WARNING reading *HYPERELASTIC: only N=1, N=2  
     &, or N=3 are allowed; '
               call inputerror(inpc,ipoinpc,iline,
     &              "*HYPERELASTIC%",ier)
               return
            endif
         else
            write(*,*) 
     &       '*WARNING reading *HYPERELASTIC: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*HYPERELASTIC%")
         endif
      enddo
!
      nelcon(1,nmat)=ityp
!
      if((ityp.ne.-6).and.(ityp.ne.-9)) then
         if((ityp.eq.-3).or.(ityp.eq.-10)) then
            iend=2
         elseif((ityp.eq.-1).or.(ityp.eq.-2).or.(ityp.eq.-4).or.
     &      (ityp.eq.-7)) then
            iend=3
         elseif(ityp.eq.-11) then
            iend=4
         elseif(ityp.eq.-13) then
            iend=5
         elseif((ityp.eq.-5).or.(ityp.eq.-12).or.(ityp.eq.-14)) then
            iend=6
         elseif(ityp.eq.-8) then
            iend=7
         endif
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            ntmat=ntmat+1
            nelcon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*)'*ERROR reading *HYPERELASTIC: increase ntmat_'
               ier=1
               return
            endif
            do i=1,iend
               read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &                  elcon(i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*HYPERELASTIC%",ier)
                  return
               endif
            enddo
            read(textpart(iend+1)(1:20),'(f20.0)',iostat=istat) 
     &                  elcon(0,ntmat,nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*HYPERELASTIC%",ier)
               return
            endif
         enddo
      else
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            ntmat=ntmat+1
            nelcon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*)'*ERROR reading *HYPERELASTIC: increase ntmat_'
               ier=1
               return
            endif
            do i=1,8
               read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &                     elcon(i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*HYPERELASTIC%",ier)
                  return
               endif
            enddo
!
            if(ityp.eq.-6) then
               iend=1
            else
               iend=4
            endif
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
               write(*,*) 
     &           '*ERROR reading *HYPERELASTIC: hyperelastic definition'
               write(*,*) '  is not complete. '
               call inputerror(inpc,ipoinpc,iline,
     &              "*HYPERELASTIC%",ier)
               return
            endif
            do i=1,iend
               read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &                  elcon(8+i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*HYPERELASTIC%",ier)
                  return
               endif
            enddo
            read(textpart(iend+1)(1:20),'(f20.0)',iostat=istat) 
     &                  elcon(0,ntmat,nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*HYPERELASTIC%",ier)
               return
            endif
         enddo
      endif
!
!     if any of the compressibility coefficients is zero (incompressible
!     material), it is replaced. The lowest order coefficient is replaced
!     such that it corresponds to a Poisson coeffient of 0.475, the 
!     following ones are replaced by a power of the first one
!
      do j=1,ntmat
!
!        calculating the shear coefficient in the undeformed state
!
         if(ityp.eq.-1) then
            um=elcon(1,j,nmat)
         elseif(ityp.eq.-2) then
            um=2.d0*(elcon(1,j,nmat)+elcon(2,j,nmat))
         elseif(ityp.eq.-3) then
            um=2.d0*elcon(1,j,nmat)
         elseif(ityp.eq.-4) then
            um=elcon(1,j,nmat)
         elseif(ityp.eq.-5) then
            um=elcon(1,j,nmat)+elcon(3,j,nmat)
         elseif(ityp.eq.-6) then
            um=elcon(1,j,nmat)+elcon(3,j,nmat)+elcon(5,j,nmat)
         elseif((ityp.eq.-7).or.(ityp.eq.-8).or.(ityp.eq.-9)) then
            um=2.d0*(elcon(1,j,nmat)+elcon(2,j,nmat))
         elseif((ityp.eq.-10).or.(ityp.eq.-11).or.(ityp.eq.-12)) then
            um=2.d0*elcon(1,j,nmat)
         elseif(ityp.eq.-13) then
            um=elcon(1,j,nmat)
         elseif(ityp.eq.-14) then
            um=2.d0*elcon(1,j,nmat)
         endif
!
         do i=1,3
            k=jcoef(i,abs(ityp))
            if(k.eq.0) exit
            if(dabs(elcon(k,j,nmat)).lt.1.d-10) then
               elcon(k,j,nmat)=(0.1d0/um)**i
               write(*,*) 
     &             '*WARNING reading *HYPERELASTIC: default value was'
               write(*,*) '         used for compressibility coefficient
     &s'
               write(*,100) i,elcon(k,j,nmat)
            endif
         enddo
      enddo
!
 100  format('   D',i1,' = ',e11.4)
!
      return
      end

