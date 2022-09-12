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
      subroutine cyclichardenings(inpc,textpart,nelcon,nmat,ntmat_,
     &        npmat_,plicon,nplicon,ncmat_,elcon,matname,
     &        irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &        ier)
!
!     reading the input deck: *CYCLIC HARDENING
!
      implicit none
!
      character*1 inpc(*)
      character*80 matname(*)
      character*132 textpart(16)
!
      integer nelcon(2,*),nmat,ntmat_,ntmat,npmat_,npmat,istep,ier,
     &  n,key,i,nplicon(0:ntmat_,*),istat,ncmat_,itemp,id,ipoinpc(0:*),
     &  ndata,ndatamax,kin,irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*)
!
      real*8 plicon(0:2*npmat_,ntmat_,*),temperature,
     &  elcon(0:ncmat_,ntmat_,*),plconloc(802),t1l
!
      ntmat=0
      npmat=0
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) 
     &      '*ERROR reading *CYCLIC HARDENING: *CYCLIC HARDENING'
         write(*,*) '       should be placed before all step'
         write(*,*) '       definitions'
         ier=1
         return
      endif
!
      if(nmat.eq.0) then
         write(*,*) 
     &      '*ERROR reading *CYCLIC HARDENING: *CYCLIC HARDENING'
         write(*,*) '       should be preceded'
         write(*,*) '       by a *MATERIAL card'
         ier=1
         return
      endif
!
      if(((nelcon(1,nmat).ne.-51).and.(nelcon(1,nmat).ne.-114)).or.
     &   (nplicon(0,nmat).ne.0)) then
         write(*,*) 
     &       '*ERROR reading *CYCLIC HARDENING: *CYCLIC HARDENING'
         write(*,*) '       should be preceded'
         write(*,*) '  by an *PLASTIC,HARDENING=COMBINED card'
         ier=1
         return
      endif
!
      do i=2,n
         write(*,*) 
     &   '*WARNING reading *CYCLIC HARDENING: parameter not recognized:'
         write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
         call inputwarning(inpc,ipoinpc,iline,
     &"*CYCLIC HARDENING%")
      enddo
!
!        isotropic hardening coefficients
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) exit
         read(textpart(3)(1:20),'(f20.0)',iostat=istat) temperature
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*CYCLIC HARDENING%",ier)
            return
         endif
!
!           first temperature
!
         if(ntmat.eq.0) then
            npmat=0
            ntmat=ntmat+1
            if(ntmat.gt.ntmat_) then
               write(*,*) 
     &           '*ERROR reading *CYCLIC HARDENING: increase ntmat_'
               ier=1
               return
            endif
            nplicon(0,nmat)=ntmat
            plicon(0,ntmat,nmat)=temperature
!
!           new temperature
!
         elseif(plicon(0,ntmat,nmat).ne.temperature) then
            npmat=0
            ntmat=ntmat+1
            if(ntmat.gt.ntmat_) then
               write(*,*) 
     &           '*ERROR reading *CYCLIC HARDENING: increase ntmat_'
               ier=1
               return
            endif
            nplicon(0,nmat)=ntmat
            plicon(0,ntmat,nmat)=temperature
         endif
         do i=1,2
            read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &           plicon(2*npmat+i,ntmat,nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*CYCLIC HARDENING%",ier)
               return
            endif
         enddo
         npmat=npmat+1
         if(npmat.gt.npmat_) then
            write(*,*) 
     &         '*ERROR reading *CYCLIC HARDENING: increase npmat_'
            ier=1
            return
         endif
         nplicon(ntmat,nmat)=npmat
      enddo
!
      if(ntmat.eq.0) then
         write(*,*) 
     &       '*ERROR reading *CYCLIC HARDENING: *CYCLIC HARDENING card'
         write(*,*) '       without data encountered'
         ier=1
         return
      endif
!
!     elastically anisotropic materials: recasting the input data
!     in a format conform to the user routine umat_aniso_plas.f
!
      if(nelcon(1,nmat).eq.-114) then
!
!        isotropic hardening
!
!        interpolating the plastic data at the elastic temperature
!        data points
!
         ndatamax=0
         do i=1,nelcon(2,nmat)
            t1l=elcon(0,i,nmat)
!     
            if(nplicon(0,nmat).eq.1) then
               id=-1
            else
               call ident2(plicon(0,1,nmat),t1l,nplicon(0,nmat),
     &              2*npmat_+1,id)
            endif
!     
            if(nplicon(0,nmat).eq.0) then
               continue
            elseif((nplicon(0,nmat).eq.1).or.(id.eq.0).or.
     &              (id.eq.nplicon(0,nmat))) then
               if(id.le.0) then
                  itemp=1
               else
                  itemp=id
               endif
               kin=0
               call plcopy(plicon,nplicon,plconloc,npmat_,ntmat_,
     &              nmat,itemp,i,kin)
               if((id.eq.0).or.(id.eq.nplicon(0,nmat))) then
               endif
            else
               kin=0
               call plmix(plicon,nplicon,plconloc,npmat_,ntmat_,
     &              nmat,id+1,t1l,i,kin)
            endif
!     
            ndata=int(plconloc(801))
            if(ndata.eq.1) then
               elcon(10,i,nmat)=plconloc(2)
               elcon(11,i,nmat)=0.d0
            else
               elcon(10,i,nmat)=plconloc(2)
               elcon(11,i,nmat)=(plconloc(4)-plconloc(2))/
     &              (plconloc(3)-plconloc(1))
            endif
            ndatamax=max(ndata,ndatamax)
         enddo
         if(ndatamax.gt.2) then
            write(*,*) '*WARNING in plastics: isotropic hardening'
            write(*,*) '         curve is possibly nonlinear for'
            write(*,*) '         the elastically anisotropic'
            write(*,*) '         material ',matname(nmat)(71:80)
         endif
      endif
!
      return
      end

