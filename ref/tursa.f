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
c  1. tursa
c  2. tsaupd
c  3. dfunsa
c  4. tsafill
c  5. vortsa
c  6. curvature
c  7. tsabc
c  8. tsaic
c
c************************************************************
c
c************************************************************************
c
c               Spalart-Allmaras Turbulence Model 
c
c   The last section of this file contains the geometry dependent
c   subroutines for the turbulence model, turbcexp, turic,
c   turinit.   These may need to be modified to implement special
c   situations.
c************************************************************************
c
c
c****************************************************************
      subroutine tursa(jmax,kmax,x,y,q,
     &                 rtxy,dj,vnut,anut,
     &                 vort,diff,conv,prod,dest,
     &                 fr1,tran,sjm,sjp,skm,
     &                 skp,stc,ds,sf,
     &                 scr1)
c
c 
#include "common.f"
      common/tsacoma/akarman,dtm,cb1,cb2,sigma,cw1,cw2,cw3,cv1,cv2,
     &               ct1,ct2,ct3,ct4,cr1,cr2,cr3,anutinf,savers
      __REAL x(jmax,kmax), y(jmax,kmax), q(jmax,kmax,3), 
     & rtxy(jmax,kmax,2,3),dj(jmax,kmax), vnut(jmax,kmax), 
c     & anut(jmax,kmax,2),  ! not sure why 2 is in original code 
     & anut(jmax,kmax), 
     & vort(jmax,kmax),diff(jmax,kmax), conv(jmax,kmax), 
     & prod(jmax,kmax),dest(jmax,kmax),fr1(jmax,kmax),tran(jmax,kmax), 
     & smin(jmax,kmax),
     & sjm(jmax,kmax),sjp(jmax,kmax),skm(jmax,kmax),skp(jmax,kmax),
     & stc(jmax,kmax),ds(jmax,kmax),sf(jmax,kmax),scr1(jmax,kmax,9)
      __INTEGER jmax, kmax
c
      logical first
      data first /.true./

      logical DEBUG

#include "mpif.h"
#include "mpi_params.f"
      __INTEGER stat(MPI_STATUS_SIZE)

      DEBUG = .false.
c      DEBUG = .true.

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, '>>> MODULE: tursa.f, SUB: tursa'
       print*
      endif

c-----
c  Version number
c-----
      savers = 1.00
      anutinf = 0.001

c-----
c  Initialize implicit arrays
c-----
      do 2 k=1,kmax
      do 2 j=1,jmax
         sjm(j,k) = 0.0
         sjp(j,k) = 0.0
         skm(j,k) = 0.0
         skp(j,k) = 0.0
         stc(j,k) = 0.0
         sf(j,k) = 0.0
         ds(j,k) = 0.0
         smin(j,k) = 0.0
2     continue

c-----
c  First time through: compute the generalized distance function;
c    Only needs to be done once
c-----
      if (first) then
         first = .false.

c compute generalized distance function for all points
         call dfunsa(jmax,kmax,x,y,smin)

c-----
c  Initial conditions
c-----
c   start code from initial conditions  
c ... computes anut(j,k); sets anut(j,k)=0 at walls
          call tsaic(jmax,kmax,x,y,q,rtxy,dj,anut,
     &                    vnut,vort,smin)
c-----
c  End of first-time-through initializations
c-----
      endif

c-----
c  For each zone, fill arrays with system of equations
c-----
      if (DEBUG) then
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (nodeid.eq.0) then
          print*
          print*, '>>> Node 0 - before tsafill in tursa.f/tursa' 
          print*
        endif
c        if (nodeid.eq.4) then
c          print*
c          print*, ' STC VALUES'
c          do 133 k=1,kmax
c          do 133 j=1,jmax
c           write(*,234) j,k,stc(j,k)
c134        format(i3,5x,i3,5x,1pe12.5)
c133       continue
c        endif
        call flush(6)
c        stop ' stop: 10 stopping in tursa.f/tursa before call tsafill'
      endif

       call tsafill(jmax,kmax,x,y,q,rtxy,
     &  dj,vnut,anut,vort,diff,conv,
     &  prod,dest,fr1,tran,
     &  sjm,sjp,skm,skp,stc,ds,
     &  sf,smin
     &  )

c add identity matrix to lhs
      call tsaibt(jmax,kmax,sjm,sjp,skm,skp,stc,sf)

      if (DEBUG) then
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (nodeid.eq.0) then
          print*
          print*, '>>> Node 0 - before tsabc in tursa.f/tursa' 
          print*
        endif
        call flush(6)
        call flush(istdout)
        if (nodeid.eq.4) then
          print*
          print*, ' STC VALUES after call tsaibt'
          do 233 k=1,kmax
          do 233 j=1,jmax
           write(*,234) j,k,stc(j,k)
234        format(i3,5x,i3,5x,1pe12.5)
233       continue
        endif
c        call MPI_FINALIZE(ierr)
c        stop ' stop: 20 stopping in tursa.f/tursa'
      endif

c-----
c  Boundary conditions
c-----
         call tsabc(jmax,kmax,anut(1,1),
     &              sjm(1,1),sjp(1,1),
     &              skm(1,1),skp(1,1),
     &              stc(1,1),sf(1,1),nodeid+1)

      if (DEBUG) then
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (nodeid.eq.0) then
          print*
          print*, '>>> Node 0 - before lreltur in tursa.f/tursa' 
          print*
        endif
        call flush(6)
        if (nodeid.eq.4) then
          print*
          print*, ' STC VALUES'
          do 333 k=1,kmax
          do 333 j=1,jmax
           write(*,334) j,k,stc(j,k)
334        format(i3,5x,i3,5x,1pe12.5)
333       continue
        endif
c        call MPI_FINALIZE(ierr)
c        stop ' stop: 22 stopping in tursa.f/tursa'
      endif


c-----
c  Solve system of equations
c    line relaxation
c-----
      if(impsch .eq. 1) then
         if (DEBUG.and.nodeid.eq.0) then
          print*
          print*, '> call lreltur'
c          stop 'stop: 23 stopping in tursa.f/tursa'
         endif
         call lreltur(jmax,kmax,anut,ds,sf,
     &                sjm,sjp,skm,skp,stc,nodeid+1)

c-----
c  lusgs -- symmetric gauss-seidel
c-----
      elseif(impsch .eq. 2) then
         neqs = 1
         call lusgs(jmax,kmax,neqs,anut,ds,sf,
     &              stc,sjm,sjp,skm,skp,scr1)

c-----
c  error if mean flow equations are using something else
c-----
      else
        write(istdout,*) ' ERROR in tursa: impsch out of range.'
        call exit(1)
      endif

      if (DEBUG.and.nodeid.eq.0) then
        print*
        print*, '>>> Node 0 - before tsaupd in tursa.f/tursa' 
        print*
        stop ' stop: 24 stopping in tursa.f/tursa'
      endif


c-----
c  Update anut array with new solution
c-----
         call tsaupd(jmax,kmax,anut(1,1),ds(1,1),1,nodeid)

      if (DEBUG.and.nodeid.eq.0) then
        print*
        print*, '>>> Node 0 - before eddy visc in tursa.f/tursa' 
        print*
        call flush(6)
c        stop ' stop: stopping in tursa.f/tursa'
      endif

c-----
c  Load eddy viscosity
c-----
      do 200 k=1,kmax
      do 200 j=1,jmax
         chi = anut(j,k)
         fv1 = chi**3/(chi**3 + cv1**3)
         vnut(j,k) = fv1*anut(j,k)*vnu
200   continue

      if (DEBUG.and.nodeid.eq.0) then
        print*
        print*,' All nodes have reached end of tursa.f/tursa, 200'
c       stop 'stopping in sub tursa.f/tursa, 200'
      endif 

       if (DEBUG.and.(nodeid.eq.8.or.nodeid.eq.4)) then
        print*, 'tursa DWB0.4: nodeid,q(2,1,1) = ',nodeid,q(2,1,1)
        print*, 'tursa DWB0.4: nodeid,q(10,1,2) = ',nodeid,q(10,1,2)
c        stop 'stopping in sub tursa.f/tursa, 210'
       endif

c-----
c  End of tursa
c-----
      return
      end
c
c
c************************************************************************
      subroutine tsaupd(jmax,kmax,anut,ds,nz,nodeid)
c************************************************************************
c  Update anut and compute rms change in anut
c------------------------------------------------------------------------
#include "common.f"
      common/tsacoma/akarman,dtm,cb1,cb2,sigma,cw1,cw2,cw3,cv1,cv2,
     &               ct1,ct2,ct3,ct4,cr1,cr2,cr3,anutinf,savers
      __REAL anut(jmax,kmax), ds(jmax,kmax)
      double precision, parameter :: zero=0.0
      logical DEBUG
     
      DEBUG = .false.
c      DEBUG = .true.
c
      sumn = 0.0
      dsmax = 0.0
      do 50 k=1,kmax
      do 50 j=1,jmax
         anut(j,k) = max(zero, anut(j,k) + ds(j,k))
         dsmax = max( ds(j,k), dsmax )
50    continue
      turres = dsmax

      if (DEBUG) then
       print*
       print*, ' sub tsaupd: nodeid, turres = ',nodeid, turres
       print*
      endif

c-----
c  End of tsaupd
c-----
      return
      end
c
c
c****************************************************************
      subroutine dfunsa(jmax,kmax,x,y,smin)
c****************************************************************
c  Compute generalized distance function for all point in all zones
c  and load into array smin.
c
c NOTES:
c  For static grids, this subroutine can be moved out from under 
c   subroutine step, which calls it for every iteration.
c  For moving grids, leave this subroutine in place to be called
c   by subroutine tursa.
c  Since for miniapps we may want to examine moving grids in the
c   future, this subroutine was left in place, called by tursa.
c----------------------------------------------------------------
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      INTEGER nodeid,numprocs,ierr,stat(MPI_STATUS_SIZE)

      __REAL smin(jmax,kmax), x(jmax,kmax), y(jmax,kmax)
      __INTEGER jmax, kmax
c      common/tsacoma/akarman,dtm,cb1,cb2,sigma,cw1,cw2,cw3,cv1,cv2,
c     &               ct1,ct2,ct3,ct4,cr1,cr2,cr3,anutinf,savers

      logical DEBUG
c-----

      DEBUG = .false.
c      DEBUG = .true.


      do 5 k=1,kmax
       do 5 j=1,jmax 
         smin(j,k) = 1.e6
5     continue

      if (DEBUG) then
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      endif

c
c  For every point on a wall, loop through all points in all grids
c
      do 30 nw=1,nwall
        if (jkwall(nw).eq.1) then 
c ... j=constant
         numpointswall = kwall2(nw) - kwall1(nw) + 1
        elseif (jkwall(nw).eq.2) then 
c ... k=constant
         numpointswall = jwall2(nw) - jwall1(nw) + 1
        endif
        do 20 npw=1,numpointswall
            xwall = xwallv(npw,nw)
            ywall = ywallv(npw,nw)
             do 10 k=1,kmax
              do 10 j=1,jmax
               dxx = x(j,k) - xwall
               dyy = y(j,k) - ywall
               smin(j,k) = min(smin(j,k), sqrt(dxx*dxx + dyy*dyy))
10          continue
20       continue
c        endif

30    continue
      
      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' j    k      smin(j,k)'
       do 50 j=1,jmax
       do 50 k=1,kmax
        write(*,60) j,k,smin(j,k)
60      format(2i5,1pe13.5)
50     continue
c       print*, 'stopping in tursa.f/dfunsa'
       call flush(istdout)
       call flush(6)
c       stop 'stop: stopping in tursa.f/dfunsa'
      endif

c-----
c  End of dfunsa
c-----
      return
      end
c
c
c****************************************************************
      subroutine tsafill(jmax,kmax,x,y,q,rtxy,dj,
     &                   vnut,anut,vort,diff,conv,prod,dest,
     &                   fr1,tran,sjm,sjp,skm,skp,stc,ds,sf,smin)
c****************************************************************
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      INTEGER nodeid,numprocs,ierr,stat(MPI_STATUS_SIZE)
c
      __REAL q(jmax,kmax,3), rtxy(jmax,kmax,2,3), dj(jmax,kmax),
     &  vnut(jmax,kmax), x(jmax,kmax),y(jmax,kmax),anut(jmax,kmax),
     &  vort(jmax,kmax), diff(jmax,kmax), conv(jmax,kmax),
     &  prod(jmax,kmax), dest(jmax,kmax),
     &  tran(jmax,kmax),  fr1(jmax,kmax),
     &   sjm(jmax,kmax),sjp(jmax,kmax),ds(jmax,kmax),
     &   skm(jmax,kmax),stc(jmax,kmax),skp(jmax,kmax),sf(jmax,kmax),
     &  smin(jmax,kmax)
      common/tsacoma/akarman,dtm,cb1,cb2,sigma,cw1,cw2,cw3,cv1,cv2,
     &               ct1,ct2,ct3,ct4,cr1,cr2,cr3,anutinf,savers
      logical first      
      data first /.true. /
      logical DEBUG
      double precision, parameter :: ten=10.0, biglim=1.e30
      double precision, parameter :: epss=1.e-12, epss2=0.00001, epss3=1.e-20

      DEBUG = .false.
c      DEBUG = .true.

c-----
c  Input data
c-----
      akarman = 0.41
      dtm     = dt
      cb1 = 0.1355
      cb2 = 0.622
      two = 2.0
      three = 3.0
      sigma = two/three
      cw1 = cb1/akarman**2 + (1. + cb2)/sigma
      cw2 = 0.3
      cw3 = 2.0
      cv1 = 7.1
      cv2 = 5.0
      ct1 = 1.0
      ct2 = 2.0
      ct3 = 1.2
      ct4 = 0.5
c
c   Curvature constants. They can be disabled by setting cr1 = -1.
c
      cr1 = 1
      cr2 = 12
      cr3 = 1
c
      cw3r6 = cw3**6
c-----
c  Output the input variables
c-----
      if(first)then
         first = .false.
         if (nodeid.eq.0) then
         write(istdout,*)' '
         write(istdout,*)'*********************************************'
         write(istdout,*)'  Node 0:  Spalart-Allmaras Turbulence Model'
         write(istdout,*)'               Version ',savers
         write(istdout,*)' '
         write(istdout,*)' karman constant = ',akarman
         write(istdout,*)' dtm = ',dtm
         write(istdout,*)' sigma = ',sigma
         write(istdout,*)' cb1 = ',cb1
         write(istdout,*)' cb2 = ',cb2
         write(istdout,*)' cw1 = ',cw1
         write(istdout,*)' cw2 = ',cw2
         write(istdout,*)' cw3 = ',cw3
         write(istdout,*)' cv1 = ',cv1
         write(istdout,*)' cv2 = ',cv2
         write(istdout,*)' ct1 = ',ct1
         write(istdout,*)' ct2 = ',ct2
         write(istdout,*)' ct3 = ',ct3
         write(istdout,*)' ct4 = ',ct4
         write(istdout,*)' cr1 = ',cr1
         write(istdout,*)' cr2 = ',cr2
         write(istdout,*)' cr3 = ',cr3
         write(istdout,*)' anutinf = ',anutinf
         write(istdout,*)' Transition locations: '
         write(istdout,*)'   nz,itrans = ',nodeid,itrans
         write(istdout,*)'     jtran1,jtran2 = ',
     &                            jtran1,jtran2
         write(istdout,*)'     xtran1,xtran2 = ',
     &                            xtran1,xtran2
         write(istdout,*)' Line-relaxation sweeps: '
         write(istdout,*)'    nz,ntjsp,ntksp = ',
     &                     nodeid,ntjsp,ntksp
         write(istdout,*)' Walls for original grid: nwall = ',nwall
         do 10 nw=1,nwall
           write(istdout,333) nw,nzwall(nw),jwall1(nw),jwall2(nw),
     &                kwall1(nw),kwall2(nw)
333   format(' nw = ',i3,'   nz = ',i5,'   j: ',i6,i6,'   k: ',i6,i6)
10       continue
         write(istdout,*)'*********************************************'
         write(istdout,*)' '
c-----
c  End of first-time through initialization
c-----
      endif
      endif
c-----

      if (DEBUG.and.nodeid.eq.0) then
       print*
       call flush(6)
       call flush(istdout)
c       stop 'stop: after S-A turb model output in tursa.f/tsafill'
      endif
c
c            Now the model    
c
c  Compute the vorticity magnitude.
c-----
      if (DEBUG.and.nodeid.eq.0) then
        print*
        print*, ' >>> Node 0 - before vortsa in tursa.f/tsafill'
        print*
c DWB: stop
c        stop ' stop: stopping before vortsa in tursa.f/tsafill'
      endif

      call vortsa(jmax,kmax,rtxy,dj,q,
     &   sjm,sjp,skm,skp,vort)

c-----
c  Compute curvature terms
c-----
      if (DEBUG.and.nodeid.eq.0) then
        print*
        print*, ' >>> Node 0 - before curvature in tursa.f/tsafill'
        print*
        call flush(6)
c        stop ' stop: stopping before curvature in tursa.f/tsafill'
      endif

      call curvature(jmax,kmax,x,y,rtxy,dj,q,
     &               sjm,sjp,skm,skp,fr1)

c-----
c  Now solve the real equations
c-----

c-----
c  f_eta_eta viscous terms
c-----
      if (DEBUG.and.nodeid.eq.0) then
        print*
        print*, ' >>> Node 0 - before f_eta_eta in tursa.f/tsafill'
        print*
        call flush(6)
c        stop ' stop: stopping before f_eta_eta terms in tursa.f/tsafill'
      endif

c            do 180 k=kend2,kendm
            do 180 k=2,kmax-1
               km1  = nint( mod( float(k-1+kmax), float(kmax)+.1 ) )
               kp1  = nint( mod( float(k+1), float(kmax)+.1 ) )
            do 180 j=2,jmax-1
                  xy3p     = .5*( rtxy(j,k  ,2,2)*dj(j,k  )
     &                          + rtxy(j,kp1,2,2)*dj(j,kp1))
                  xy4p     = .5*( rtxy(j,k  ,2,3)*dj(j,k  )
     &                          + rtxy(j,kp1,2,3)*dj(j,kp1))
                  ttp      = (xy3p*rtxy(j,k,2,2)*dj(j,k)
     &                     +  xy4p*rtxy(j,k,2,3)*dj(j,k))
c                  
                  xy3m     = .5*( rtxy(j,k  ,2,2)*dj(j,k  )
     &                          + rtxy(j,km1,2,2)*dj(j,km1))
                  xy4m     = .5*( rtxy(j,k  ,2,3)*dj(j,k  )
     &                          + rtxy(j,km1,2,3)*dj(j,km1))
                  ttm      = (xy3m*rtxy(j,k,2,2)*dj(j,k)
     &                     +  xy4m*rtxy(j,k,2,3)*dj(j,k))
c
                  coef = -cb2*vnu*(1. + anut(j,k))/sigma
                  diff(j,k) = 
     &                        coef*( ttp*(anut(j,kp1) - anut(j,k))
     &                             - ttm*(anut(j,k) - anut(j,km1)) )
                  skm(j,k) = -coef*ttm
                  skp(j,k) = -coef*ttp
                  stc(j,k) =  coef*(ttm + ttp)
c      if (DEBUG.and.nodeid.eq.0) then
c       print*, '1 tursa.f/tsafill: nodeid,j,k,stc(j,k)',
c     & nodeid,j,k,stc(j,k)
c       stop 'stop: after first print of stc'
c      endif
c                  
                  coef = vnu*(1.+cb2)/sigma
                  vnup = 1.0 + 0.5*(anut(j,kp1) + anut(j,k))
                  vnum = 1.0 + 0.5*(anut(j,km1) + anut(j,k))
                  ttp = ttp*vnup
                  ttm = ttm*vnum
                  diff(j,k) = diff(j,k)
     &                      + coef*( ttp*(anut(j,kp1) - anut(j,k))
     &                             - ttm*(anut(j,k) - anut(j,km1)) )
                  skm(j,k) = skm(j,k) - coef*ttm
                  skp(j,k) = skp(j,k) - coef*ttp
                  stc(j,k) = stc(j,k) + coef*(ttm + ttp)
c      if (DEBUG.and.nodeid.eq.0) then
c       print*, '2 tursa.f/tsafill: nodeid,j,k,stc(j,k)',
c     & nodeid,j,k,stc(j,k)
c       stop 'stop: after second print of stc'
c      endif
c                  
180       continue

      if (DEBUG.and.nodeid.eq.0) then
        print*
        print*, '2 tursa.f/tsafill'
        print*, '   j    k    stc(j,k)'
       do 999 j=2,jmax-1
       do 999 k=2,kmax-1
        write(*,995) j,k,stc(j,k)
995     format(2i5,1pe13.5)
999    continue
       call flush(istdout)
       call flush(6)
c       stop 'stop: after second print of stc'
      endif

      if (DEBUG) then
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      endif

c-----
c  Advective terms in eta
c-----
      if (DEBUG.and.nodeid.eq.0) then
        print*
        print*, ' >>> Node 0 - before advect in tursa.f/tsafill'
        print*
        call flush(6)
c        stop ' stop: stopping before advect terms in tursa.f/tsafill'
      endif

         one = 1.0
c         do 190 k=kend2,kendm
         do 190 k=2,kmax-1
            km1  = nint( mod( float(k-1+kmax), float(kmax)+.1 ) )
            kp1  = nint( mod( float(k+1), float(kmax)+.1 ) )
         do 190 j=2,jmax-1
            uu = (rtxy(j,k,2,2)*q(j,k,2)+rtxy(j,k,2,3)*q(j,k,3))
     &           *dj(j,k)
            sgnu = sign(one,uu)
            app  = .5*(1.+sgnu)
            apm  = .5*(1.-sgnu)
            conv(j,k) =      + uu*( app*(anut(j,k)-anut(j,km1))
     &                             +apm*(anut(j,kp1)-anut(j,k)) )
            skm(j,k)   = skm(j,k)   - uu*app
            stc(j,k)   = stc(j,k)   + uu*(app-apm)
      if (DEBUG.and.nodeid.eq.4) then
       print*, '3 tursa.f/tsafill: nodeid,j,k,stc(j,k)',
     & nodeid,j,k,stc(j,k)
c       stop 'stop: after second print of stc'
      endif
            skp(j,k)   = skp(j,k)   + uu*apm
190      continue

c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      call MPI_FINALIZE(ierr)
c      call flush(6)
c      call flush(istdout)
c      stop 'stop: after 190 continue'

c-----
c   e_xi_xi viscous terms
c-----
      if (DEBUG.and.nodeid.eq.0) then
        print*
        print*, ' >>> Node 0 - before e_xi_xi in tursa.f/tsafill'
        print*
        call flush(6)
c        stop ' stop: stopping before e_xi_xi terms in tursa.f/tsafill'
      endif

            do 220 k=2,kmax-1
            do 220 j=2,jmax-1
                  jp1 = j+1
                  jm1 = j-1
                  xy1p     = .5*(rtxy(j  ,k,1,2)*dj(j  ,k)
     &                         + rtxy(jp1,k,1,2)*dj(jp1,k))
                  xy2p     = .5*(rtxy(j  ,k,1,3)*dj(j  ,k)
     &                         + rtxy(jp1,k,1,3)*dj(jp1,k))
                  ttp      = (xy1p*rtxy(j,k,1,2)*dj(j,k)
     &                     +  xy2p*rtxy(j,k,1,3)*dj(j,k))
                  
                  xy1m     = .5*(rtxy(j  ,k,1,2)*dj(j  ,k)
     &                         + rtxy(jm1,k,1,2)*dj(jm1,k))
                  xy2m     = .5*(rtxy(j  ,k,1,3)*dj(j  ,k)
     &                         + rtxy(jm1,k,1,3)*dj(jm1,k))
                  ttm      = (xy1m*rtxy(j,k,1,2)*dj(j,k)
     &                     +  xy2m*rtxy(j,k,1,3)*dj(j,k))
c
                  coef = -cb2*vnu*(1. + anut(j,k))/sigma
                  diff(j,k) = diff(j,k)
     &                      + coef*( ttp*(anut(j+1,k) - anut(j,k))
     &                             - ttm*(anut(j,k) - anut(j-1,k)) )
                  sjm(j,k) = -coef*ttm
                  sjp(j,k) = -coef*ttp
                  stc(j,k) =  coef*(ttm + ttp) + stc(j,k)
c                  
                  coef = vnu*(1.+cb2)/sigma
                  vnup = 1.0 + 0.5*(anut(j+1,k) + anut(j,k))
                  vnum = 1.0 + 0.5*(anut(j-1,k) + anut(j,k))
                  ttp = ttp*vnup
                  ttm = ttm*vnum
                  diff(j,k) = diff(j,k)
     &                      + coef*( ttp*(anut(j+1,k) - anut(j,k))
     &                             - ttm*(anut(j,k) - anut(j-1,k)) )
                  sjm(j,k) = sjm(j,k) - coef*ttm
                  sjp(j,k) = sjp(j,k) - coef*ttp
                  stc(j,k) = stc(j,k) + coef*(ttm + ttp)
c      if (DEBUG.and.nodeid.eq.4) then
c       print*, '4 tursa.f/tsafill: nodeid,j,k,stc(j,k)',
c     & nodeid,j,k,stc(j,k)
c       stop 'stop: after second print of stc'
c      endif
220      continue

c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      call MPI_FINALIZE(ierr)
c      call flush(6)
c      call flush(istdout)
c      stop 'stop: after 220 continue'

c-----
c  Advective terms in xi
c-----
      if (DEBUG.and.nodeid.eq.0) then
        print*
        print*, ' >>> Node 0 - before xi adv terms in tursa.f/tsafill'
        print*
        call flush(6)
c        stop ' stop: stopping before xi adv terms in tursa.f/tsafill'
      endif

           one = 1.0
           do 230 k=2,kmax-1
           do 230 j=2,jmax-1
                 uu = (rtxy(j,k,1,2)*q(j,k,2)+rtxy(j,k,1,3)*q(j,k,3))
     &                *dj(j,k)
                 sgnu = sign(one,uu)
                 app  = .5*(1.+sgnu)
                 apm  = .5*(1.-sgnu)
                 conv(j,k)=conv(j,k)+uu*(app*(anut(j,k)-anut(j-1,k))
     >                                  +apm*(anut(j+1,k)-anut(j,k)) )
                 sjm(j,k)   = sjm(j,k)   - uu*app
                 stc(j,k)   = stc(j,k)   + uu*(app-apm)
      if (DEBUG.and.nodeid.eq.0) then
       print*, '5 tursa.f/tsafill: nodeid,j,k,stc(j,k)',
     & nodeid,j,k,stc(j,k)
c       stop 'stop: after second print of stc'
      endif
                 sjp(j,k)   = sjp(j,k)   + uu*apm
230        continue

c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      call MPI_FINALIZE(ierr)
c      call flush(6)
c      call flush(istdout)
c      stop 'stop: after 230 continue'

c           do 235 k=kend2,kendm
           do 235 k=2,kmax-1
           do 235 j=2,jmax-1
              sf(j,k) = -conv(j,k) + diff(j,k)
235        continue

c-----
c  Production and destruction terms
c-----
      if (DEBUG.and.nodeid.eq.0) then
        print*
        print*, ' >>> Node 0 - before prod terms in tursa.f/tsafill'
        print*, ' itrans = ',itrans
        print*
        call flush(6)
c        call MPI_FINALIZE(ierr)
c        stop ' stop: stopping before prod terms in tursa.f/tsafill'
      endif

         if (itrans .eq. 0) then
            ct3tmp= 0.0
         else
            ct3tmp = ct3
         endif
c
c         do 260 k=kend2,kendm
         do 260 k=2,kmax-1
         do 260 j=2,jmax-1
          if (DEBUG) then
           print*
           print*, ' in tursa.f/tsafill: j,k,epss,smin(j,k) =',
     &      j,k,epss,smin(j,k)
          endif
            fib = 1 
            sminlim = max( smin(j,k), epss )
            rkd2 = 1./( akarman*sminlim )**2
            ss = vort(j,k)
c  Spalart mod. Use limited version of chi.
c  All the functions are well-behaved as chi -> 0.
            chi = max ( anut(j,k), epss2 )
            fv1 = chi**3/(chi**3 + cv1**3)
c  Spalart mod. New damping functions so Stilde is always positive.
            fv2 = 1 / ( 1. + chi / cv2 ) ** 3
            fv3 = ( 1 + chi * fv1 ) / chi * ( 1 - fv2 )
c  Old fv2  fv2 = 1. - chi/(1. + chi*fv1)
            ft2 = ct3tmp*exp( -ct4*chi**2 )
            sst = fv3 * ss + vnu * chi * fv2 * rkd2
c  Old sst  sst = ss + vnu*chi*fv2*rkd2
c  Spalart mod. Protection for zero values.
c            rr = 100
c            if ( sst .ne. 0. ) rr = vnu*anut(j,k)*rkd2/sst
c            rr = min(rr, 10.)
            sst = max(sst, epss3)
            rr = min(ten, vnu*anut(j,k)*rkd2/sst)
c
            gg = rr + cw2*(rr**6 - rr)
c      if (DEBUG.and.nodeid.eq.4) then
c      if (DEBUG) then
c       print*
c       print*, 'In tsafill #1: nodeid,rr,cw2,gg,j,k,anut_jk = ',
c     &   nodeid,rr,cw2,gg,j,k,anut(j,k)
c       if (gg.gt.300007.) print*, ' DANGER: gg = ',gg 
c       call flush(6)
c       call flush(istdout)
c       print*, ' XXX  nodeid, gg, gg**6 = ',nodeid,gg,gg**6
c       stop 'stop: in tursa.f/tsafill #1'
c      endif
            gg6 = gg**6
            fg = (1. + cw3r6)/(gg6 + cw3r6)
            fwdg = fg**(1./6.)
            fw = gg*fwdg
            pp = cb1*sst*(1. - ft2)
            dd = vnu*(cw1*fw-cb1*ft2/(akarman**2))
     &              *anut(j,k)/(sminlim**2)
c
            ft2p = -2.*ct4*ct3tmp*chi*exp( -ct4*chi**2 )
            fv1p = 3.*chi**2*cv1**3/( (chi**3 + cv1**3)**2 )
            fv2p = - 3 / cv2 / ( 1. + chi / cv2 ) ** 4
            fv3p = ( ( fv1 + chi * fv1p ) * ( 1 - fv2 )
     &      - ( 1 + chi * fv1 ) * ( ( 1 - fv2 ) / chi + fv2p ) ) / chi
c  Old fv2  fv2p = ( chi**2*fv1p - 1. )/(1. + chi*fv1)**2
            sstp = fv3p * ss + vnu*(chi*fv2p+fv2)*rkd2
c  Old sst  sstp = vnu*rkd2*( fv2 + chi*fv2p )
c  Spalart mod. Protection for zero values.
c            rrp = 0
c            if ( sst .ne. 0. ) rrp = vnu*rkd2*(sst - chi*sstp)/(sst**2)
            rrp = vnu*rkd2*(sst - chi*sstp)/(sst**2)
            ggp = rrp*(1. + cw2*(6.*rr**5 - 1.))
c
c            fgpdenom = ((gg**6 + cw3r6)**2)
c            fgp = -6.*(1. + cw3r6)*gg**5*ggp/fgpdenom
c            fwp = ggp*fg**(1./6.) + gg*fg**(-5./6.)*fgp/6.
c
            fwp = ggp*fwdg*(1. - gg6/(gg6 + cw3r6))
c
            ppp = cb1*(sstp - ft2p*sst - ft2*sstp)
            ddp = vnu*( cw1*fw - cb1*ft2/akarman**2
     &                + (cw1*fwp-cb1*ft2p/akarman**2)*chi )/sminlim**2
c
            prod(j,k) = fib*pp*anut(j,k)*fr1(j,k)
            dest(j,k) = fib*dd*anut(j,k)
c
            sf(j,k) = sf(j,k) + prod(j,k) - dest(j,k)
c
            zero = 0.
            addpp = max( dd-pp, zero )
            addpppp = max( ddp-ppp, zero )
            stc(j,k) = stc(j,k)  + fib*addpp
            stc(j,k) = stc(j,k) + fib*addpppp*anut(j,k)
      if (DEBUG.and.nodeid.eq.0) then
       print*, '6 tursa.f/tsafill: nodeid,j,k,stc(j,k)',
     & nodeid,j,k,stc(j,k)
c       stop 'stop: after second print of stc'
      endif
260      continue
 
      if (DEBUG.and.nodeid.eq.0) then
       print*
c       print*, 'stop: in tursa.f/tsafill #2 ins'
       call flush(6)
       call flush(istdout)
c       stop 'stop: in tursa.f/tsafill #2 ins'
      endif

c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      call MPI_FINALIZE(ierr)
      if (DEBUG.and.nodeid.eq.4) then
c       call flush(6)
c       call flush(istdout)
c       stop 'stop: after 260 continue'
      endif

      if (DEBUG.and.nodeid.eq.0) then
c        print*
c        print*, ' >>> Node 0 - end tursa.f/tsafill'
c        print*
c        call flush(6)
c        call MPI_FINALIZE(ierr)
        stop ' stop: stopping end  tursa.f/tsafill'
      endif

      return
      end
c
c
c********************************************************************
      subroutine vortsa(jmax,kmax,rtxy,dj,q,
     &                   uj,vj,uk,vk,vort)
c********************************************************************
#include "common.f"
      __REAL rtxy(jmax,kmax,2,3),
     &    dj(jmax,kmax), q(jmax,kmax,3), vort(jmax,kmax),
     &    uj(jmax,kmax), vj(jmax,kmax), uk(jmax,kmax), vk(jmax,kmax)
c-----
c  j-direction differences
c-----
      do 10 j=2,jmax-1
         jmm = max( j-2, 1 )
         jpp = min( j+2, jmax )
      do 10 k=1,kmax
         ujc = 0.5*( q(j+1,k,2) - q(j-1,k,2) )
         vjc = 0.5*( q(j+1,k,3) - q(j-1,k,3) )
         ujp = 0.5*( -q(jpp,k,2) + 4.*q(j+1,k,2) - 3.*q(j,k,2) )
         vjp = 0.5*( -q(jpp,k,3) + 4.*q(j+1,k,3) - 3.*q(j,k,3) )
         ujm = 0.5*( 3.*q(j,k,2) - 4.*q(j-1,k,2) + q(jmm,k,2) )
         vjm = 0.5*( 3.*q(j,k,3) - 4.*q(j-1,k,3) + q(jmm,k,3) )
c
         uj(j,k) = ujc 
         vj(j,k) = vjc
10    continue

      j = 1
      do 20 k=1,kmax
         uj(j,k) = 0.5*( -q(j+2,k,2) + 4.*q(j+1,k,2) - 3.*q(j,k,2) )
         vj(j,k) = 0.5*( -q(j+2,k,3) + 4.*q(j+1,k,3) - 3.*q(j,k,3) )
20    continue

      j = jmax
      do 30 k=1,kmax
         uj(j,k) = 0.5*( 3.*q(j,k,2) - 4.*q(j-1,k,2) + q(j-2,k,2) )
         vj(j,k) = 0.5*( 3.*q(j,k,3) - 4.*q(j-1,k,3) + q(j-2,k,3) )
30    continue
c-----
c  k-direction differences
c-----
      do 40 k=2,kmax-1
         kmm = max( k-2, 1 )
         kpp = min( k+2, kmax )
      do 40 j=1,jmax
         ukc = 0.5*( q(j,k+1,2) - q(j,k-1,2) )
         vkc = 0.5*( q(j,k+1,3) - q(j,k-1,3) )
         ukp = 0.5*( -q(j,kpp,2) + 4.*q(j,k+1,2) - 3.*q(j,k,2) )
         vkp = 0.5*( -q(j,kpp,3) + 4.*q(j,k+1,3) - 3.*q(j,k,3) )
         ukm = 0.5*( 3.*q(j,k,2) - 4.*q(j,k-1,2) + q(j,kmm,2) )
         vkm = 0.5*( 3.*q(j,k,3) - 4.*q(j,k-1,3) + q(j,kmm,3) )
c
         uk(j,k) = ukc 
         vk(j,k) = vkc 
40    continue

       k = 1
       do 70 j=1,jmax
         uk(j,k) = 0.5*( -q(j,k+2,2) + 4.*q(j,k+1,2) - 3.*q(j,k,2) )
         vk(j,k) = 0.5*( -q(j,k+2,3) + 4.*q(j,k+1,3) - 3.*q(j,k,3) )
70     continue
       k = kmax
       do 80 j=1,jmax
          uk(j,k) = 0.5*( 3.*q(j,k,2) - 4.*q(j,k-1,2) + q(j,k-2,2) )
          vk(j,k) = 0.5*( 3.*q(j,k,3) - 4.*q(j,k-1,3) + q(j,k-2,3) )
80     continue

c-----
c  Compute vorticity
c-----
      do 100 k=1,kmax
      do 100 j=1,jmax
         vort(j,k) = abs( dj(j,k)*(rtxy(j,k,1,2)*vj(j,k)
     &                           + rtxy(j,k,2,2)*vk(j,k)
     &                           - rtxy(j,k,1,3)*uj(j,k)
     &                           - rtxy(j,k,2,3)*uk(j,k)) )
100   continue
c-----
c  End of vortsa
c-----
      return
      end
c
c
c********************************************************************
      subroutine curvature(jmax,kmax,x,y,rtxy,dj,q,
     &                     uj,vj,uk,vk,fr1)
c********************************************************************
#include "common.f"
      __REAL x(jmax,kmax), y(jmax,kmax), rtxy(jmax,kmax,2,3),
     &    dj(jmax,kmax), q(jmax,kmax,3), fr1(jmax,kmax), 
     &    uj(jmax,kmax), vj(jmax,kmax), uk(jmax,kmax), vk(jmax,kmax)
      common/tsacoma/akarman,dtm,cb1,cb2,sigma,cw1,cw2,cw3,cv1,cv2,
     &               ct1,ct2,ct3,ct4,cr1,cr2,cr3,anutinf,savers
c-----
c  j-direction differences
c-----
      do 10 j=2,jmax-1
         jmm = max( j-2, 1 )
         jpp = min( j+2, jmax )
      do 10 k=1,kmax
         ujc = q(j+1,k,2) - q(j-1,k,2)
         vjc = q(j+1,k,3) - q(j-1,k,3)
         ujp = -q(jpp,k,2) + 4*q(j+1,k,2) - 3*q(j,k,2)
         vjp = -q(jpp,k,3) + 4*q(j+1,k,3) - 3*q(j,k,3)
         ujm = 3*q(j,k,2) - 4*q(j-1,k,2) + q(jmm,k,2)
         vjm = 3*q(j,k,3) - 4*q(j-1,k,3) + q(jmm,k,3)
c
         uj(j,k) = .5 *  ujc 
         vj(j,k) = .5 * vjc 
10    continue

      j = 1
      do 20 k=1,kmax
         uj(j,k) = 0.5*( -q(j+2,k,2) + 4.*q(j+1,k,2) - 3.*q(j,k,2) )
         vj(j,k) = 0.5*( -q(j+2,k,3) + 4.*q(j+1,k,3) - 3.*q(j,k,3) )
20    continue

      j = jmax
      do 30 k=1,kmax
         uj(j,k) = 0.5*( 3.*q(j,k,2) - 4.*q(j-1,k,2) + q(j-2,k,2) )
         vj(j,k) = 0.5*( 3.*q(j,k,3) - 4.*q(j-1,k,3) + q(j-2,k,3) )
30    continue

c-----
c  k-direction differences
c-----
      do 40 k=2,kmax-1
         kmm = max( k-2, 1 )
         kpp = min( k+2, kmax )
      do 40 j=1,jmax
         ukc = q(j,k+1,2) - q(j,k-1,2)
         vkc = q(j,k+1,3) - q(j,k-1,3)
         ukp = -q(j,kpp,2) + 4*q(j,k+1,2) - 3*q(j,k,2)
         vkp = -q(j,kpp,3) + 4*q(j,k+1,3) - 3*q(j,k,3)
         ukm = 3*q(j,k,2) - 4*q(j,k-1,2) + q(j,kmm,2)
         vkm = 3*q(j,k,3) - 4*q(j,k-1,3) + q(j,kmm,3)
c
         uk(j,k) = .5*ukc 
         vk(j,k) = .5*vkc 
40    continue

       k = 1
       do 70 j=1,jmax
         uk(j,k) = 0.5*( -q(j,k+2,2) + 4.*q(j,k+1,2) - 3.*q(j,k,2) )
         vk(j,k) = 0.5*( -q(j,k+2,3) + 4.*q(j,k+1,3) - 3.*q(j,k,3) )
70     continue
       k = kmax
       do 80 j=1,jmax
         uk(j,k) = 0.5*( 3.*q(j,k,2) - 4.*q(j,k-1,2) + q(j,k-2,2) )
         vk(j,k) = 0.5*( 3.*q(j,k,3) - 4.*q(j,k-1,3) + q(j,k-2,3) )
80     continue

c-----
c  Compute gradient terms, and strain rate
c-----
      do 100 k=1,kmax
      do 100 j=1,jmax
         ux = dj(j,k)*(rtxy(j,k,1,2)*uj(j,k)+rtxy(j,k,2,2)*uk(j,k))
         vx = dj(j,k)*(rtxy(j,k,1,2)*vj(j,k)+rtxy(j,k,2,2)*vk(j,k))
         uy = dj(j,k)*(rtxy(j,k,1,3)*uj(j,k)+rtxy(j,k,2,3)*uk(j,k))
         vy = dj(j,k)*(rtxy(j,k,1,3)*vj(j,k)+rtxy(j,k,2,3)*vk(j,k))
c  Start computing curvature term. Eq. 4 in the paper.
c  To save storage, store strain terms in available arrays.
c  rstar will be in vj:
         vj(j,k) = sqrt(2.*(ux**2+vx*uy+vy**2)+uy**2+vx**2)
c  if the vorticity is zero, production is zero, so don't bother to divide.
         if ( vx - uy .ne. 0. ) vj(j,k) = vj(j,k) / abs ( vx - uy)
c  s11 = - s22 will be in uj:
         uj(j,k) = ( ux - vy ) / 2
c  s12 will be in uk:
         uk(j,k) = ( uy + vx ) / 2
c  fr1 will hold the beginning of the rotation-rate formula:
         fr1(j,k) = 2 * ( vx - uy )
         denom = (2.*(ux**2+vx*uy+vy**2)+uy**2+vx**2) *
     &       ( ux ** 2 + vx ** 2 + uy ** 2 + vy ** 2 )
c  if the denominator is 0, so will the numerator be, so don't bother to divide.
        if( denom .ne. 0. ) fr1(j,k) = fr1(j,k) / denom
100   continue
c
c  Use centered differences for interior points
      do 101 k=1,kmax
        km = max (    1, k - 1 )
        kp = min ( kmax, k + 1 )
        do 102 j=2,jmax-1
          uu = (rtxy(j,k,2,2)*q(j,k,2)+rtxy(j,k,2,3)*q(j,k,3))*dj(j,k)
          vv = (rtxy(j,k,1,2)*q(j,k,2)+rtxy(j,k,1,3)*q(j,k,3))*dj(j,k)
          ds11dt = uu * ( uj(j, kp) - uj(j, km) ) / ( kp - km )
     &           + vv * ( uj(j+1,k) - uj(j-1,k) ) / 2
          ds12dt = uu * ( uk(j, kp) - uk(j, km) ) / ( kp - km )
     &           + vv * ( uk(j+1,k) - uk(j-1,k) ) / 2
          rtilde = ( uj(j,k) * ds12dt - uk(j,k) * ds11dt ) * fr1(j,k)
          fr1(j,k) = ( 1 + cr1 ) * 2 * vj(j,k) / ( 1 + vj(j,k) ) *
     &    ( 1 - cr3 * atan ( cr2 * rtilde ) ) - cr1
102     continue
c  End points
        uu = (rtxy(1,k,2,2)*q(1,k,2)+rtxy(1,k,2,3)*q(1,k,3))*dj(1,k)
        vv = (rtxy(1,k,1,2)*q(1,k,2)+rtxy(1,k,1,3)*q(1,k,3))*dj(1,k)
        ds11dt = uu * ( uj(1, kp) - uj(1, km) ) / ( kp - km )
     &         + vv * ( uj(2,k) - uj(1,k) )
        ds12dt = uu * ( uk(1, kp) - uk(1, km) ) / ( kp - km )
     &         + vv * ( uk(2,k) - uk(1,k) )
        rtilde = ( uj(1,k) * ds12dt - uk(1,k) * ds11dt ) * fr1(1,k)
        fr1(1,k) = ( 1 + cr1 ) * 2 * vj(1,k) / ( 1 + vj(1,k) ) *
     &  ( 1 - cr3 * atan ( cr2 * rtilde ) ) - cr1
        uu = (rtxy(jmax,k,2,2)*q(jmax,k,2)+
     &        rtxy(jmax,k,2,3)*q(jmax,k,3))*dj(jmax,k)
        vv = (rtxy(jmax,k,1,2)*q(jmax,k,2)+
     &        rtxy(jmax,k,1,3)*q(jmax,k,3))*dj(jmax,k)
        ds11dt = uu * ( uj(jmax, kp) - uj(jmax, km) ) / ( kp - km )
     &         + vv * ( uj(jmax,k) - uj(jmax-1,k) )
        ds12dt = uu * ( uk(jmax, kp) - uk(jmax, km) ) / ( kp - km )
     &         + vv * ( uk(jmax,k) - uk(jmax-1,k) )
        rtilde = ( uj(jmax,k) * ds12dt - uk(jmax,k) * ds11dt ) *
     &                                                    fr1(jmax,k)
        fr1(jmax,k) = ( 1 + cr1 ) * 2 * vj(jmax,k) / ( 1 + vj(jmax,k) )*
     &  ( 1 - cr3 * atan ( cr2 * rtilde ) ) - cr1
101   continue
c-----
c  End of curvature
c-----
      return
      end
c
c
c************************************************************************
      subroutine tsabc(jmax,kmax,anut,sjm,sjp,skm,skp,stc,sf,nz)
c************************************************************************
c  This routine handles boundary conditions for the turbulence model
c  using info from the flow code's boundary condition surface arrays.
c  It is called once for each zone.
c
c           solid walls -    anut(j,k) = 0.0 (no change)
c           inflow      -    no change
c           outflow     -    extrapolate
c------------------------------------------------------------------------
#include "common.f"
      __REAL anut(jmax,kmax), sjm(jmax,kmax), sjp(jmax,kmax),
     &       skm(jmax,kmax), skp(jmax,kmax), stc(jmax,kmax),
     &       sf(jmax,kmax)
      INTEGER nz 
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

c-----
c  Process the boundary condition surface if it belongs
c  to the current zone.
c-----
      do 100 nr=1,nbcreg
c         if(nzbc(nr) .eq. nz) then
         if(nzbc(nr) .eq. nz) then
c-----
c  Outflow or symmetry
c-----
            if(ibcval(nr) .eq. 30 .or. ibcval(nr) .eq. 31
     &        .or. ibcval(nr) .eq. 32 .or. ibcval(nr) .eq. 33
     &        .or. ibcval(nr) .eq. 40) then
               if (DEBUG) then
                print*
                write(*,50) ibcval(nr),nz
50    format(/'Applying turb bc # ',i3,' for zone ',i4)
               endif
               if(jkbc(nr) .eq. 1 .and. jbcb(nr) .eq. 1) then
                  if (DEBUG) then
                   print*, ' jkbc(nr), jbcb(nr) = ',jkbc(nr),jbcb(nr)
                  endif
                  j = jbcb(nr)
                  do 20 k=kbcb(nr),kbce(nr)
                     stc(j,k) = 1.0
                     sjp(j,k) =-1.0
                     sf(j,k) = anut(j+1,k) - anut(j,k)
20                continue
               elseif(jkbc(nr) .eq. 1 .and. jbcb(nr) .eq. jmax) then
                  if (DEBUG) then
                   print*, ' jkbc(nr), jbcb(nr) = ',jkbc(nr),jbcb(nr)
                  endif
                  j = jbcb(nr)
                  do 22 k=kbcb(nr),kbce(nr)
                     stc(j,k) = 1.0
                     sjm(j,k) =-1.0
                     sf(j,k) = anut(j-1,k) - anut(j,k)
22                continue
               elseif(jkbc(nr) .eq. 2 .and. kbcb(nr) .eq. 1) then
                  if (DEBUG) then
                   print*, ' jkbc(nr), kbcb(nr) = ',jkbc(nr),kbcb(nr)
                  endif
                  k = kbcb(nr)
                  do 24 j=jbcb(nr),jbce(nr)
                     stc(j,k) = 1.0
                     skp(j,k) =-1.0
                     sf(j,k) = anut(j,k+1) - anut(j,k)
24                continue
               elseif(jkbc(nr) .eq. 2 .and. kbcb(nr) .eq. kmax) then
                  if (DEBUG) then
                   print*, ' jkbc(nr), kbcb(nr) = ',jkbc(nr),kbcb(nr)
                  endif
                  k = kbcb(nr)
                  do 26 j=jbcb(nr),jbce(nr)
                     stc(j,k) = 1.0
                     skm(j,k) =-1.0
                     sf(j,k) = anut(j,k-1) - anut(j,k)
26                continue
               endif
            endif
         endif
100   continue
c-----
c  End of tsabc
c-----
      return
      end
c
c
c*****************************************************************
      subroutine tsaic(jmax,kmax,x,y,q,rtxy,dj,anut,vnut,vort,
     &                 smin)
c*****************************************************************
c  Initializes turbulence quantites for Spalart-Allmaras turbulence
c  model
c-----------------------------------------------------------------
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      INTEGER nodeid,numprocs
      __REAL x(jmax,kmax), y(jmax,kmax), q(jmax,kmax,3),
     &     rtxy(jmax,kmax,2,3), dj(jmax,kmax), vort(jmax,kmax),
     &     anut(jmax,kmax), vnut(jmax,kmax), smin(jmax,kmax)
      common/tsacoma/akarman,dtm,cb1,cb2,sigma,cw1,cw2,cw3,cv1,cv2,
     &               ct1,ct2,ct3,ct4,cr1,cr2,cr3,anutinf,savers

      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

c----
c  Initialize edge points to free-stream value, and inner points to a higher
c  value, in an attempt to help the initial transient with more glue.
c  Change by Spalart, 4/15/93
c
c  Modify this to only have an effect close to the wall, using smin array.
c-----
      do 90 k=1,kmax
      do 90 j=1,jmax
         anut(j,k) = anutinf + 10.0*exp( -10.*smin(j,k) )
90    continue
c-----                                              
c  Initial conditions: zero anut on walls
c-----
      do 110 nw=1,nwall
         if (nzwall(nw).eq.nodeid+1) then
            do 100 k=kwall1(nw),kwall2(nw)
            do 100 j=jwall1(nw),jwall2(nw)
               anut(j,k) = 0.0
100         continue
         endif
110   continue

c      if (DEBUG.and.nodeid.eq.0) then
      if (DEBUG.and.nodeid.eq.4) then
       print*
       print*, ' nodeid  j   k    smin(j,k)     anut(j,k)'
       do 200 k=1,kmax
       do 200 j=1,jmax
         write(*,205) nodeid,j,k,smin(j,k),anut(j,k)
205      format(i5,2x,2i4,1p2e13.5)
200    continue
       stop ' stop: in tsaic at end'
      endif

c-----
c  End of tsaic
c-----
      return
      end
c
c
c********************************************************************
      subroutine tsaibt(jmax,kmax,sjm,sjp,skm,skp,stc,sf)
c********************************************************************
c  Adds identity matrix terms to lhs; blanks out boundary points
c--------------------------------------------------------------------
#include "common.f"
#include "mpi_params.f"
      __REAL sjm(jmax,kmax),sjp(jmax,kmax),
     &       skm(jmax,kmax),skp(jmax,kmax),
     &       stc(jmax,kmax), sf(jmax,kmax)
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

c-----
c  add identity matrix
c-----
      
      if (DEBUG.and.nodeid.eq.4) then
       print*, ' nodeid, dt = ',nodeid,dt
       call flush(6)
       call flush(istdout)
       call MPI_FINALIZE(ierr)
       stop 'stop: in tursa.f/tsaibt' 
      endif

      rterm = 1./dt
      do 20 j=1,jmax
      do 20 k=1,kmax
       if (j.eq.1.or.j.eq.jmax.or.k.eq.1.or.k.eq.kmax) then
        fib = 0
       else
        fib = 1
       endif
         sjm(j,k) = fib*sjm(j,k)
         sjp(j,k) = fib*sjp(j,k)
         skm(j,k) = fib*skm(j,k)
         skp(j,k) = fib*skp(j,k)
         stc(j,k) = fib*(stc(j,k) + rterm) + (1.-fib)
         sf(j,k)  = fib*sf(j,k)
20    continue
c-----
c  End of tsaibt
c-----
      return
      end

