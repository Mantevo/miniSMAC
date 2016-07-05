c copyright notice
c                      Copyright 2013 Sandia Corporation             
c                                                                    
c     Under the terms of Contract DE-AC04-94AL85000 with Sandia      
c     Corporation, the U.S. Government retains certain rights in this
c     software.                                                      
c                                                                    
c     Redistribution and use in source and binary forms, with or     
c     without modification, are permitted provided that the following
c     conditions are met:                                            
c                                                                    
c      1. Redistributions of source code must retain the above       
c         copyright notice, this list of conditions and the following
c         disclaimer.                                                
c                                                                    
c      2. Redistributions in binary form must reproduce the above    
c         copyright notice, this list of conditions and the following
c         disclaimer in the documentation and/or other materials     
c         provided with the distribution.                            
c                                                                    
c      3. Neither the name of the Corporation nor the names of the   
c         contributors may be used to endorse or promote products    
c         derived from this software without specific prior written  
c         permission.                                                
c                                                                    
c     THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
c     EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,  
c     THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A    
c     PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA    
c     CORPORATION OR THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT,      
c     INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL     
c     DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF         
c     SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;   
c     OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  
c     LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      
c     (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF  
c     THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
c     SUCH DAMAGE.                                                   
c     
c
c************************************************************
c Subroutines in this file:
c  1. visc4
c  2. turvis
c
c************************************************************
c
c
c**************************************************************
      subroutine viscterms(jmax,kmax,q,rtxy,dj,s,vnut,
     &                 btc,bjm,bjp,bkm,bkp)
c**************************************************************
c  Add the viscous terms to the right hand side array
c  and the implicit matrices
c
c  Assumes constant viscosity and an orthgonal grid for the
c  implicit terms.
c--------------------------------------------------------------
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"

      __REAL q(jmax,kmax,3),rtxy(jmax,kmax,2,3),
     &          dj(jmax,kmax),s(jmax,kmax,3),vnut(jmax,kmax),
     &          bjm(jmax,kmax,3,3), btc(jmax,kmax,3,3),
     &          bjp(jmax,kmax,3,3), bkm(jmax,kmax,3,3),
     &          bkp(jmax,kmax,3,3)
      real, allocatable:: udif(:),vdif(:),vf(:)
      logical DEBUG,DEBUG1

      DEBUG = .false.
c      DEBUG = .true.
 
      if (DEBUG) then
       if (nodeid.eq.0) then
        print*
        print*,' Entering visc.f/viscterms'
        print*, ' s = '
        do 10 j=1,60
        do 10 k=2,2
         write(*,400) j,k,s(j,k,1),s(j,k,2),s(j,k,3)
400      format(2x,i3,1x,i3,1x,1pe13.5,1x,1pe13.5,1x,1pe13.5)
10      continue
        print*
       endif
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c       stop 'stop: entering visc.f/viscterms'
      endif

c-----
c  J-sweep
c-----
c
      allocate(udif(jmax),vdif(jmax),vf(jmax))
      do 100 k=2,kmax-1
       do 60 j=1,jmax-1
        djv = 0.5*( dj(j,k)*vnu + dj(j+1,k)*vnu )
        udif(j) = q(j+1,k,2) - q(j,k,2)
        vdif(j) = q(j+1,k,3) - q(j,k,3)
        rx = 0.5*( rtxy(j,k,1,2) + rtxy(j+1,k,1,2) )
        ry = 0.5*( rtxy(j,k,1,3) + rtxy(j+1,k,1,3) )
        vf(j) = djv*( rx*rx + ry*ry )
60     continue
       do 80 j=2,jmax-1
        bjm(j,k,2,2) = bjm(j,k,2,2) - vf(j-1)
        bjm(j,k,3,3) = bjm(j,k,3,3) - vf(j-1)
        btc(j,k,2,2) = btc(j,k,2,2) + vf(j-1) + vf(j)
        btc(j,k,3,3) = btc(j,k,3,3) + vf(j-1) + vf(j)
        bjp(j,k,2,2) = bjp(j,k,2,2) - vf(j)
        bjp(j,k,3,3) = bjp(j,k,3,3) - vf(j)
        s(j,k,2) = s(j,k,2) + (vf(j)*udif(j)-vf(j-1)*udif(j-1))
        s(j,k,3) = s(j,k,3) + (vf(j)*vdif(j)-vf(j-1)*vdif(j-1))
80     continue
100   continue
      deallocate(udif,vdif,vf)

c-----
c  K sweep non-periodic
c-----
      allocate(udif(kmax),vdif(kmax),vf(kmax))
      do 300 j=2,jmax-1
       do 260 k=1,kmax-1
        djv = 0.5*( dj(j,k)*vnu + dj(j,k+1)*vnu )
        udif(k) = q(j,k+1,2) - q(j,k,2)
        vdif(k) = q(j,k+1,3) - q(j,k,3)
        rx = 0.5*( rtxy(j,k,2,2) + rtxy(j,k+1,2,2) )
        ry = 0.5*( rtxy(j,k,2,3) + rtxy(j,k+1,2,3) )
        vf(k) = djv*( rx*rx + ry*ry )
260    continue
       do 280 k=2,kmax-1
        bkm(j,k,2,2) = bkm(j,k,2,2) - vf(k-1)
        bkm(j,k,3,3) = bkm(j,k,3,3) - vf(k-1)
        btc(j,k,2,2) = btc(j,k,2,2) + vf(k-1) + vf(k)
        btc(j,k,3,3) = btc(j,k,3,3) + vf(k-1) + vf(k)
        bkp(j,k,2,2) = bkp(j,k,2,2) - vf(k)
        bkp(j,k,3,3) = bkp(j,k,3,3) - vf(k)
        s(j,k,2) = s(j,k,2)+(vf(k)*udif(k)-vf(k-1)*udif(k-1))
        s(j,k,3) = s(j,k,3)+(vf(k)*vdif(k)-vf(k-1)*vdif(k-1))
280    continue
300   continue

      deallocate(udif,vdif,vf)

      if (DEBUG) then
       if (nodeid.eq.0) then
        print*
        print*,' Exiting visc.f/viscterms'
        print*, ' s = '
        do 11 j=1,60
        do 11 k=2,2
         write(*,400) j,k,s(j,k,1),s(j,k,2),s(j,k,3)
11      continue
        print*
       endif
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       stop 'stop: end of visc.f/viscterms'
      endif

      return
      end
c
c**************************************************************
      subroutine viscterms0(jmax,kmax,q,rtxy,dj,s,vnut,
     &                 btc,bjm,bjp,bkm,bkp)
c**************************************************************
c  Add the viscous terms to the right hand side array
c  and the implicit matrices
c
c  Assumes non-orthgonal grid and non-constant viscosity
c  for the explicit terms.
c  Assumes constant viscosity and an orthgonal grid for the
c  implicit terms.
c--------------------------------------------------------------
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"

      __REAL q(jmax,kmax,3),rtxy(jmax,kmax,2,3),
     &          dj(jmax,kmax),s(jmax,kmax,3),vnut(jmax,kmax),
     &          bjm(jmax,kmax,3,3), btc(jmax,kmax,3,3),
     &          bjp(jmax,kmax,3,3), bkm(jmax,kmax,3,3),
     &          bkp(jmax,kmax,3,3)
      __REAL, allocatable:: uxi(:),vxi(:),ueta(:),veta(:),vf(:,:)
      logical DEBUG,DEBUG1

      DEBUG = .false.
c      DEBUG = .true.
      DEBUG1 = .false.
      DEBUG1 = .true.
 
      if (DEBUG1) then
       if (nodeid.eq.0) then
        print*
        print*,' Entering visc.f/viscterms'
        print*, ' s = '
        do 10 j=1,60
        do 10 k=2,2
         write(*,400) j,k,s(j,k,1),s(j,k,2),s(j,k,3)
400      format(2x,i3,1x,i3,1x,1pe13.5,1x,1pe13.5,1x,1pe13.5)
10      continue
        print*
       endif
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       stop 'stop: entering visc.f/viscterms'
      endif

c-----
c  J-sweep
c-----
      if (jmax.ge.kmax) then
       jk_max = jmax
      else
       jk_max = kmax
      endif

      allocate (uxi(jk_max),vxi(jk_max),ueta(jk_max),veta(jk_max))
      allocate (vf(jk_max,8))
c
      do 100 k=2,kmax-1
         if(k .eq. 1) then
            kp = k + 1
            kr = kmax
         elseif(k .eq. kmax) then
            kp = 1
            kr = k - 1
         else
            kp = k + 1
            kr = k - 1
         endif
         do 60 j=1,jmax-1
            djave = 0.5*( dj(j,k) + dj(j+1,k) )
            vnuave = vnu + 0.5*( vnut(j,k) + vnut(j+1,k) )
            djv = djave*vnuave
            uxi(j) = q(j+1,k,2) - q(j,k,2)
            vxi(j) = q(j+1,k,3) - q(j,k,3)
            ueta(j) = 0.25*( q(j+1,kp,2) + q(j,kp,2)
     &                     - q(j+1,kr,2) - q(j,kr,2) )
            veta(j) = 0.25*( q(j+1,kp,3) + q(j,kp,3)
     &                     - q(j+1,kr,3) - q(j,kr,3) )
            xix  = 0.5*( rtxy(j,k,1,2) + rtxy(j+1,k,1,2) )
            xiy  = 0.5*( rtxy(j,k,1,3) + rtxy(j+1,k,1,3) )
            etax = 0.5*( rtxy(j,k,2,2) + rtxy(j+1,k,2,2) )
            etay = 0.5*( rtxy(j,k,2,3) + rtxy(j+1,k,2,3) )
            vf(j,1) = djv*( xix*xix + xiy*xiy )
            vf(j,2) = djv*( xix*xix ) + vf(j,1)
            vf(j,3) = djv*( xiy*xiy ) + vf(j,1)
            vf(j,4) = djv*( 2.*xix*etax + xiy*etay )
            vf(j,5) = djv*( xix*etax + 2.*xiy*etay )
            vf(j,6) = djv*( xix*xiy )
            vf(j,7) = djv*( xiy*etax )
            vf(j,8) = djv*( xix*etay )
c            if (k.eq.2.and.nodeid.eq.4.and.j.eq.8) then
            if (DEBUG.and.nodeid.eq.0) then
              print*
              print*,' DELETE30: j,k,djave,vnuave = ',
     &  j,k,djave,vnuave
              print*,' uxi,vxi,ueta,veta(j) = ',
     &  uxi(j),vxi(j),ueta(j),veta(j)
              print*,' xix,xiy,etax,etay = ',
     &  xix,xiy,etax,etay
              print*,' vf(j,1-8) =',
     &  vf(j,1),vf(j,2),vf(j,3),vf(j,4),vf(j,5),vf(j,6),vf(j,7),vf(j,8)
              print*,'==========================================='
            endif
60       continue

        if (DEBUG1) then
         if (nodeid.eq.0) then
          print*
          print*, '===== before do 80 ======='
          print*,' s = ',s
         endif
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         call flush(6)
         call flush(istdout)
c         stop 'stop: before do 80 in visc.f/viscterms'
        endif 

         do 80 j=2,jmax-1
            bjm(j,k,2,2) = bjm(j,k,2,2) - vf(j-1,2)
            bjm(j,k,3,3) = bjm(j,k,3,3) - vf(j-1,3)
            btc(j,k,2,2) = btc(j,k,2,2) + vf(j-1,2) + vf(j,2)
            btc(j,k,3,3) = btc(j,k,3,3) + vf(j-1,3) + vf(j,3)
            bjp(j,k,2,2) = bjp(j,k,2,2) - vf(j,2)
            bjp(j,k,3,3) = bjp(j,k,3,3) - vf(j,3)
c            print*, ' DWB j,k,s2 = ',j,k
            s(j,k,2) = s(j,k,2) + vf(j,2)*uxi (j) - vf(j-1,2)*uxi (j-1)
     &                          + vf(j,4)*ueta(j) - vf(j-1,4)*ueta(j-1)
     &                          + vf(j,6)*vxi (j) - vf(j-1,6)*vxi (j-1)
     &                          + vf(j,7)*veta(j) - vf(j-1,7)*veta(j-1)

c            print*,' ... s(j,k,2) = ',s(j,k,2)
c            print*, 'DWB j,k,s3 = ',j,k,s(j,k,3)
            s(j,k,3) = s(j,k,3) + vf(j,3)*vxi (j) - vf(j-1,3)*vxi (j-1)
     &                          + vf(j,5)*veta(j) - vf(j-1,5)*veta(j-1)
     &                          + vf(j,6)*uxi (j) - vf(j-1,6)*uxi (j-1)
     &                          + vf(j,8)*ueta(j) - vf(j-1,8)*ueta(j-1)
c            print*, ' ... s(j,k,3) = ',s(j,k,3)

      if (DEBUG.and.nodeid.eq.0) then
       print*
c       print*, '-- in visc4, for j,k = ',j,k
       print*,'DELETE50: j,k,btc(j,k,2,2) = ',
     &  j,k,btc(j,k,2,2)
       endif
       
80       continue

100   continue

      deallocate (uxi,vxi,ueta,veta)
      deallocate (vf)

      if (DEBUG.and.nodeid.eq.0) then
c       stop 'stop: in visc.f/visc4 #1' 
      endif

c-----
c  K sweep non-periodic
c-----
      allocate (uxi(kmax),vxi(kmax),ueta(kmax),veta(kmax))
      allocate (vf(kmax,8))

         do 300 j=2,jmax-1
            do 260 k=1,kmax-1
               djave = 0.5*( dj(j,k) + dj(j,k+1) )
               vnuave = vnu + 0.5*( vnut(j,k) + vnut(j,k+1) )
               djv = djave*vnuave
               uxi(k) = 0.25*( q(j+1,k+1,2) + q(j+1,k,2)
     &                       - q(j-1,k+1,2) - q(j-1,k,2) )
               vxi(k) = 0.25*( q(j+1,k+1,3) + q(j+1,k,3)
     &                       - q(j-1,k+1,3) - q(j-1,k,3) )
               ueta(k) = q(j,k+1,2) - q(j,k,2)
               veta(k) = q(j,k+1,3) - q(j,k,3)
               xix  = 0.5*( rtxy(j,k,1,2) + rtxy(j,k+1,1,2) )
               xiy  = 0.5*( rtxy(j,k,1,3) + rtxy(j,k+1,1,3) )
               etax = 0.5*( rtxy(j,k,2,2) + rtxy(j,k+1,2,2) )
               etay = 0.5*( rtxy(j,k,2,3) + rtxy(j,k+1,2,3) )
               vf(k,1) = djv*( etax*etax + etay*etay )
               vf(k,2) = djv*( etax*etax ) + vf(k,1)
               vf(k,3) = djv*( etay*etay ) + vf(k,1)
               vf(k,4) = djv*( 2.*xix*etax + xiy*etay )
               vf(k,5) = djv*( xix*etax + 2.*xiy*etay )
               vf(k,6) = djv*( xix*etay )
               vf(k,7) = djv*( xiy*etax )
               vf(k,8) = djv*( etax*etay )
260         continue
            do 280 k=2,kmax-1
               bkm(j,k,2,2) = bkm(j,k,2,2) - vf(k-1,2)
               bkm(j,k,3,3) = bkm(j,k,3,3) - vf(k-1,3)
               btc(j,k,2,2) = btc(j,k,2,2) + vf(k-1,2) + vf(k,2)
               btc(j,k,3,3) = btc(j,k,3,3) + vf(k-1,3) + vf(k,3)
               bkp(j,k,2,2) = bkp(j,k,2,2) - vf(k,2)
               bkp(j,k,3,3) = bkp(j,k,3,3) - vf(k,3)
               s(j,k,2)=s(j,k,2) + vf(k,2)*ueta(k) - vf(k-1,2)*ueta(k-1)
     &                           + vf(k,4)*uxi (k) - vf(k-1,4)*uxi (k-1)
     &                           + vf(k,6)*vxi (k) - vf(k-1,6)*vxi (k-1)
     &                           + vf(k,8)*veta(k) - vf(k-1,8)*veta(k-1)
               s(j,k,3)=s(j,k,3) + vf(k,3)*veta(k) - vf(k-1,3)*veta(k-1)
     &                           + vf(k,5)*vxi (k) - vf(k-1,5)*vxi (k-1)
     &                           + vf(k,7)*uxi (k) - vf(k-1,7)*uxi (k-1)
     &                           + vf(k,8)*ueta(k) - vf(k-1,8)*ueta(k-1)

      if (DEBUG.and.nodeid.eq.0) then
       print*
c       print*, ' DELETE60: j,k,btc(j,k,2,2) = ',j,k,btc(j,k,2,2)
       print*, ' DELETE60: j,k,s(j,k,2) = ',j,k,s(j,k,2)
      endif

280         continue
300      continue

      deallocate (uxi,vxi,ueta,veta)
      deallocate (vf)

      if (DEBUG.and.nodeid.eq.0) then
c       stop 'stop: in visc.f/visc4 #2'
      endif

c DWB 082213: add some smoothing
c      if (dcoef2.gt.0.0.or.dcoef4.gt.0.0) then
c      do 400 k=2,kmax-1
c      do 400 j=2,jmax-1
c      do 400 nvar=1,3
c       if (k.eq.2.or.j.eq.2.or.k.eq.kmax-1.or.j.eq.jmax-1) then
cc 2nd order smoothing near the boundaries
c         dq_damp = dq(j+1,k,nvar) + dq(j-1,k,nvar) 
c     &      + dq(j,k+1,nvar) + dq(j,k-1,nvar) - 4.0d0*dq(j,k,nvar) 
c         s(j,k,nvar) = s(j,k,nvar) + dcoef2*dq_damp 
c       else
cc 4th order smoothing in the interior
c         dq_damp_j = dq(j+2,k,nvar) - 4.0d0*dq(j+1,k,nvar) 
c     &      + 6.0d0*dq(j,k,nvar) - 4.0d0*dq(j-1,k,nvar)
c     &      + dq(j-2,k,nvar)
c         dq_damp_k = dq(j,k+2,nvar) - 4.0d0*dq(j,k+1,nvar) 
c     &      + 6.0d0*dq(j,k,nvar) - 4.0d0*dq(j,k-1,nvar)
c     &      + dq(j,k-2,nvar)  
c         s(j,k,nvar) = s(j,k,nvar)
c     &    + dcoef4/8.0d0*(dq_damp_j + dq_damp_k)
c       endif
c400   continue
c      endif

c-----
c  End of visc4
c-----

      return
      end
c
c
c****************************************************************
      subroutine turvis(jmax,kmax,x,y,q,rtxy,
     &                  dj,vnut,turvar,
     &                  scr1,scr2,scr3)
c****************************************************************
c Called from: step.f/step

c  Call the Spalart-Allmaras turbulence model.
c     vnut:   eddy-viscosity
c     turvar: turbulent dependent variable(s)
c     scr*:   scratch memory
c
c----------------------------------------------------------------
#include "common.f"
      __REAL x(jmax,kmax),y(jmax,kmax),q(jmax,kmax,3),
     &   rtxy(jmax,kmax,2,3),dj(jmax,kmax),vnut(jmax,kmax),
     &   turvar(jmax,kmax,2)
      __REAL scr1(jmax*kmax,9), scr2(jmax*kmax,9), 
     &       scr3(jmax*kmax,9)  
c     &       scr4(jmax,kmax,9),scr5(jmax,kmax,9), scr6(jmax,kmax,9),

      __INTEGER jmax, kmax, jmxdum, kmxdum
c
      logical first
      data first /.true./
c
      logical DEBUG

#include "mpif.h"
#include "mpi_params.f"
      __INTEGER nodeid,numprocs,ierr,stat(MPI_STATUS_SIZE)

      DEBUG = .false.
c      DEBUG = .true.

      if (DEBUG.and.nodeid.eq.0) then
        print*
        print*, '>>> Node 0 - MODULE: visc.f, SUB: turvis'
        print*
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c-----
c  Call Spalart-Allmaras model.  Pass in different segments of the
c  scratch arrays for use in 1 equation solution.
c-----

         if (DEBUG.and.(nodeid.eq.8.or.nodeid.eq.4)) then
          print*, ' DWB0.3: nodeid,q(2,1,1) = ',nodeid,q(2,1,1)
          print*, ' DWB0.3: nodeid,q(10,1,2) = ',nodeid,q(10,1,2)
         endif

       call tursa(jmax,kmax,x,y,q,
     &  rtxy,dj,vnut,turvar,
     &  scr1(1,1),scr1(1,2),scr1(1,3),scr1(1,4),scr1(1,5),
     &  scr1(1,6),scr1(1,7),scr1(1,8),scr1(1,9),scr2(1,1),
     &  scr2(1,2),scr2(1,3),scr2(1,4),scr2(1,5),
     &  scr3)

      if (DEBUG) then
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       if (nodeid.eq.0) then
        print*
        print*,' Node 0: All nodes reached end of visc.f/turvis, 200'
       endif
      call flush(6)
c      stop ' stop: stopping in visc.f/turvis, 200'
      endif

c-----
c  End of turvis
c-----
      return
      end
