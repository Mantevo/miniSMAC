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
c  1. step
c  2. resid
c  3. update
c  4. blankit
c
c************************************************************
c#include "blasdef.h"
c
c
c********************************************************************
      subroutine step(jmax,kmax,x,y,rtxy,dj,q,dq,s,
     &                vnut,turvar)
c********************************************************************
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      __REAL q(jmax,kmax,3),dq(jmax,kmax,3),
     &   x(jmax,kmax), y(jmax,kmax), rtxy(jmax,kmax,2,3),
     &   dj(jmax,kmax), s(jmax,kmax,3),
     &   bjm(jmax,kmax,3,3), btc(jmax,kmax,3,3), bjp(jmax,kmax,3,3),
     &   bkm(jmax,kmax,3,3), bkp(jmax,kmax,3,3),
     &   scr1(jmax,kmax,3,3), vnut(jmax,kmax), turvar(jmax,kmax,2)

      __INTEGER jmax, kmax
      __INTEGER stat(MPI_STATUS_SIZE)

      logical DEBUG,DEBUG1,DEBUG2

      DEBUG = .false.
c      DEBUG = .true.
      DEBUG1 = .false.
c      DEBUG1 = .true.
      DEBUG2 = .false.
c      DEBUG2 = .true.

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' --------- in step.f/step --------------'
       print*, ' beta = ',beta
c       stop 'stop: AT START OF step.f/step'
      endif

c output boundary values at j=jmaxb for Node 0
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntref1 = 1 
      ntref2 = ntmax 
      if (DEBUG1.and.numprocs.eq.2.and.nodeid.eq.0.and.
     &  (nt.eq.ntref1.or.nt.eq.ntref2)) then 
       print*
       print*, 'START of step.f/step -- Overlap BCs for nodeid:',nodeid
       print*, '  jmax,kmax = ',jmax,kmax
       print*, '  time step nt =',nt
       j=jmax
       print*, ' nodeid  j  k     q_jk1       q_jk2        q_jk3'
       do 901 k=1,kmax
        write(*,810) nodeid,j,k,q(j,k,1),q(j,k,2),q(j,k,3)
810     format(3i5,1p3e15.7)
901    continue
      endif

c zero arrays btc, bjp, bjm, bkm, bkp, scr1
      do 5000 m=1,3
      do 5000 n=1,3
      do 5000 k=1,kmax
      do 5000 j=1,jmax
       btc(j,k,n,m) = 0.0
       bjp(j,k,n,m) = 0.0
       bjm(j,k,n,m) = 0.0
       bkm(j,k,n,m) = 0.0
       bkp(j,k,n,m) = 0.0
       scr1(j,k,n,m) = 0.0
5000  continue

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*,' After do 5000 in step.f/step:'
c       print*, ' s = ',s
c       stop ' stop: after do 5000 in step.f/step'
      endif

      if (DEBUG.and.(nodeid.eq.0)) then
      icount = 0
      print*
      print*, '  node  i   j   k   dq1        dq2        dq3'
      print*, ' This is test 0 in step.f/step'
      itest = 0 
c      do 4000 j=1,jmax
c      do 4000 k=1,kmax
c       icount = icount + 1
c       write(*,4005) itest,nodeid,icount,j,k,dq(j,k,1),dq(j,k,2),
c    &  dq(j,k,3)
c4005   format(i4,3x,i4,3x,i3,2x,i3,2x,i3,1x,1pe10.3,1x,1pe10.3,
c     &  1x,1pe10.3)
c4000  continue
      print*
      print*
      endif

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' DELETE1: nodeid, q(32,2,1) = ',nodeid,q(32,2,1)
      endif
c      stop 'stop: entering step.f/step'

      if (DEBUG) then
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (nodeid.eq.0) then
         print*
         print*, ' All nodes have reached step.f/step, 100'
        endif
      endif

      if (DEBUG.and.(nodeid.eq.0)) then
        print*, 'DWB0: nodeid,q(2,1,1) = ',nodeid,q(2,1,1)
        print*, 'DWB0: nodeid,q(10,1,2) = ',nodeid,q(10,1,2)
      endif

c-----
c  Compute eddy viscosity if this is a turbulent calculation
c  for all zones.
c-----

         if (DEBUG.and.(nodeid.eq.0)) then
          print*, ' DWB0.2: nodeid,q(2,1,1) = ',nodeid,q(2,1,1)
          print*, ' DWB0.2: nodeid,q(10,1,2) = ',nodeid,q(10,1,2)
          print*, ' DWB0.2: nodeid, btc(1,1,1,1) = ',nodeid,
     &      btc(1,1,1,1)
          print*, ' DWB0.2: nodeid, btc(2,2,1,1) = ',nodeid,
     &      btc(2,2,1,1)
         endif

      if (DEBUG.and.nodeid.eq.0) then
      print*
      print*, ' DELETE2: nodeid, q(2,2,2) = ',nodeid,q(2,2,2)
      endif
c      stop 'stop: before call turvis in step.f/step'

c if viscous turbulent flow...
      if (ivis.eq.1.and.iturb.eq.1) then
         call turvis(jmax,kmax,x,y,q,rtxy,
     &               dj,vnut,turvar,
     &               scr1,btc,bjm)

      endif

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' DELETE3: nodeid, q(32,2,1) = ',nodeid,q(32,2,1)
      endif
c      stop 'stop: after call turvis in step.f/step'

         if (DEBUG.and.(nodeid.eq.0)) then
          print*, ' DWB0.51: nodeid,q(2,1,1) = ',nodeid,q(2,1,1)
          print*, ' DWB0.51: nodeid,q(10,1,2) = ',nodeid,q(10,1,2)
          print*, ' DWB0.51: nodeid, btc(1,1,1,1) = ',nodeid,
     &      btc(1,1,1,1)
          print*, ' DWB0.51: nodeid, btc(2,2,1,1) = ',nodeid,
     &      btc(2,2,1,1)
         endif
c-----
c  Compute entire system of equations for all zones
c-----
       if (DEBUG.and.nodeid.eq.0) then
        print*
        print*, 'XXXXXXX sub step XXXXXXXXXXXXXXXX'
        print*, ' ... impsch = ',impsch
        print*, ' '
       endif

c-----
c  Upwind differencing of inviscid terms
c-----
      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' node 0 dwb: call conup'
       print*
      endif

c compute convective terms using upwind differencing;
c  this is where btc, bjm, bjp, bkm, and bkp are defined.
c      do 100 m=1,3
c      do 100 n=1,3
c      do 100 k=1,kmax
c      do 100 j=1,jmax
c       btc(j,k,m,n) = 0.
c       bjm(j,k,m,n) = 0.
c       bjp(j,k,m,n) = 0.
c       bkm(j,k,m,n) = 0.
c       bkp(j,k,m,n) = 0.
c100    continue

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call flush(6)
      call flush(istdout)
      icount = 0
      print*
      if (DEBUG.and.nodeid.eq.0) then
       print*, '  node  i   j   k   dq1        dq2        dq3'
       print*,' This is test 1 in step.f/step'
       itest = 1
       do 2000 j=1,jmax
       do 2000 k=1,kmax
        icount = icount + 1
        write(*,2005) itest,nodeid,icount,j,k,dq(j,k,1),dq(j,k,2),
     &   dq(j,k,3)
2005    format(i4,3x,i4,3x,i4,2x,i3,2x,i3,1x,1pe10.3,1x,1pe10.3,
     &  1x,1pe10.3)
2000   continue
       print*
       print*, 'BEFORE call conup in step.f/step 7:'
       print*, ' btc(2,2,1,1), btc(2,2,2,2) =',
     &   btc(2,2,1,1),btc(2,2,2,2)
       print*
c       stop 'stop: before call conup in step.f/step 7'
      endif

c upwind difference using flux-difference splitting
c       call conup(jmax,kmax,q(1,1,1),rtxy(1,1,1,1),dj(1,1),s(1,1,1),
c     & btc(1,1,1,1),bjm(1,1,1,1),bjp(1,1,1,1),bkm(1,1,1,1),bkp(1,1,1,1)
c     & )

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call flush(6)
      call flush(istdout)

      call conup(jmax,kmax,q,rtxy,dj,s,btc,bjm,bjp,bkm,bkp)

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' AFTER call conup in step.f/step 8:'
       print*, '  btc(2,2,1,1), btc(2,2,2,2) = ',
     &  btc(2,2,1,1),btc(2,2,2,2)
       print*
c       stop 'stop: AFTER call conup in step.f/step 8:'
      endif

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' step.f/step DELETE4: nodeid, q(2,2,2),dq(2,2,2) = ',
     &  nodeid,q(2,2,2),dq(2,2,2)
      endif
c      stop 'stop: before call viscterms in step.f/step'

      if (DEBUG.and.(nodeid.eq.0)) then
       print*, ' DWB0.52: nodeid,q(2,1,1) = ',nodeid,q(2,1,1)
       print*, ' DWB0.52: nodeid,q(10,1,2) = ',nodeid,q(10,1,2)
      endif

c      stop 'stop: before call viscterms in step.f/step'

c-----
c  Viscous terms
c-----
      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' node 0 dwb: call viscterms'
       print*
      endif

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, 'btc in step.f/step, before call viscterms:'
       do 811 k=1,3
       do 801 j=1,jmax
       print*, ' DELETE20.1: j,k,btc(j,k,1,1),dj(j,k),dq(j,k,2) = ',
     &   j,k,btc(j,k,1,1),dj(j,k),dq(j,k,2)
801    continue
       print*, '=================================================='
811    continue
       print*
      endif


      if (DEBUG) then
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call flush(6)
      call flush(istdout)
      icount = 0
      print*
      if (DEBUG.and.nodeid.eq.0) then
      print*, 'itst  node  i   j   k   dq1        dq2        dq3'
      print*,' This is test 2'
      itest = 2
      do 1000 j=1,jmax
      do 1000 k=1,kmax
       icount = icount + 1
       write(*,1005) itest,nodeid,icount,j,k,dq(j,k,1),dq(j,k,2),
     &  dq(j,k,3)
1005   format(i4,3x,i4,3x,i3,2x,i3,2x,i3,1x,1pe10.3,1x,1pe10.3,1x,
     &  1pe10.3)
1000  continue
      endif
c      stop 'stop: before call viscterms in step.f/step 250'
      endif

      if (DEBUG2) then
       if (nodeid.eq.0) then
        print*
        print*,' Before call viscterms in step.f/step:'
        print*, ' s = ',s
       endif
c       stop ' stop: Before call viscterms in step.f/step'
      endif

      if (ivis.eq.1) then
            call viscterms(jmax,kmax,q,rtxy,dj,s,vnut,btc,bjm,
     &         bjp,bkm,bkp
     &         )
      endif

      if (DEBUG) then
      call flush(6)
      call flush(istdout)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (nodeid.eq.0) then
      print*
      print*, ' DELETE5: nodeid, q(32,2,1) = ',nodeid,q(32,2,1)
      endif
c      stop 'stop: before call bcmain in step.f/step'
      endif

       if (DEBUG.and.(nodeid.eq.0)) then
        print*, ' DWB0.53: nodeid,q(2,1,1) = ',nodeid,q(2,1,1)
        print*, ' DWB0.53: nodeid,q(10,1,2) = ',nodeid,q(10,1,2)
       endif

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, 'btc before blankit, node 0:'
c       k=2
c       do 802 j=1,jmax
c       print*, ' DELETE20.2: j,k,btc(j,k,1,1),dj(j,k),dq(j,k,2) = ',
c     &   j,k,btc(j,k,1,1),dj(j,k),dq(j,k,2)
c802    continue
       print*, ' btc(2,2,1,1), btc(2,2,2,2) =',
     &   btc(2,2,1,1),btc(2,2,2,2)
       print*
c       stop 'stop: before call blankit in step.f/step 9'
      endif

c-----
c  at boundary points where ib<=0, values are set to zero.
c-----
      do 380 m=1,3
      do 380 k=1,kmax
      do 380 j=1,jmax
       if (j.eq.1.or.j.eq.jmax.or.k.eq.1.or.k.eq.kmax) then
c set boundary values to zero
         fib = 0
       else
c otherwise, no change
         fib = 1
       endif
c         fib = float( max(0, min(1, ib(i)) ) )
         btc(j,k,m,1) = fib*btc(j,k,m,1)
         btc(j,k,m,2) = fib*btc(j,k,m,2)
         btc(j,k,m,3) = fib*btc(j,k,m,3)
         bjp(j,k,m,1) = fib*bjp(j,k,m,1)
         bjp(j,k,m,2) = fib*bjp(j,k,m,2)
         bjp(j,k,m,3) = fib*bjp(j,k,m,3)
         bjm(j,k,m,1) = fib*bjm(j,k,m,1)
         bjm(j,k,m,2) = fib*bjm(j,k,m,2)
         bjm(j,k,m,3) = fib*bjm(j,k,m,3)
         bkp(j,k,m,1) = fib*bkp(j,k,m,1)
         bkp(j,k,m,2) = fib*bkp(j,k,m,2)
         bkp(j,k,m,3) = fib*bkp(j,k,m,3)
         bkm(j,k,m,1) = fib*bkm(j,k,m,1)
         bkm(j,k,m,2) = fib*bkm(j,k,m,2)
         bkm(j,k,m,3) = fib*bkm(j,k,m,3)
         s(j,k,m) = fib*s(j,k,m)
380   continue

c-----
c add identity matrix to lhs
c-----
c form identity matrix for steady-state calculations
      do 385 k=1,kmax
      do 385 j=1,jmax
       if (j.eq.1.or.j.eq.jmax.or.k.eq.1.or.k.eq.kmax) then
c set boundary values to zero
         fib = 0
       else
c otherwise, no change
         fib = 1
       endif
       rdjdtau = fib/(dj(j,k)*dtau) + 1. - fib
       btc(j,k,1,1) = btc(j,k,1,1) + rdjdtau
       btc(j,k,2,2) = btc(j,k,2,2) + rdjdtau
       btc(j,k,3,3) = btc(j,k,3,3) + rdjdtau
      if (DEBUG.and.nodeid.eq.0.and.j.eq.2.and.k.eq.2) then
        print*
        print*, ' In blankit 2: btc(j,k,2,2) = ',
     &   j,k,btc(j,k,2,2),btc(j,k,1,1),rdjdtau,dtau,dj(j,k),fib
      endif
385   continue

      if (DEBUG.and.nodeid.eq.0) then
c        print*, ' stop after blankit loops'
c        stop 'stop: after blankit loops'
      endif

c-----
c  Boundary conditions
c-----
      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' node 0 dwb: before call bcmain'
       print*, ' ... node 0: btc(1,1,1,1) = ',btc(1,1,1,1)
       print*
      endif

      if (DEBUG) then 
       call flush(6)
       call flush(istdout)
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      stop 'stop: before call bcmain in step.f/step'
      endif

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, 'btc before call bcmain:'
       print*, ' node 0: btc(2,2,1,1), btc(2,2,2,2):',
     &   btc(2,2,1,1),btc(2,2,2,2)
c       stop 'stop: before call bcmain 15'
      endif

         call bcmain(jmax,kmax,x,y,q,rtxy,
     &        dj,nodeid+1,s,btc,bjm,
     &        bjp,bkm,bkp
     &        )

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' DELETE6: nodeid, q(32,2,1) = ',nodeid,q(32,2,1)
      endif
c      stop 'stop: before call resid in step.f/step'

      if (DEBUG.and.nodeid.eq.0) then
       print*, ' DWB5: nodeid,q(2,1,1) = ',nodeid,q(2,1,1)
       print*, ' DWB5: nodeid,q(10,1,2) = ',nodeid,q(10,1,2)
       print*, ' DWB5: nodeid,btc(1,1,1,1) = ',nodeid,btc(1,1,1,1)
      endif
c-----
c  Find max residual, max divergence, and rmsrhs
c-----
      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' node 0 dwb: call resid'
       print*
      endif

         call resid(jmax,kmax,s,scr1,dj,nodeid+1)

      if (DEBUG.and.(nodeid.eq.0)) then
      print*
      print*, ' DELETE7: nodeid, q(32,2,1) = ',nodeid,q(32,2,1)
      endif
     
      if(DEBUG) then
      call flush(6)
      call flush(istdout)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      stop 'stop: before call lusgs in step.f/step'
      endif

      if (DEBUG.and.(nodeid.eq.0)) then
       print*, ' DWB5.1: nodeid, q(2,1,1) = ',nodeid,q(2,1,1)
       print*, ' DWB5.1: nodeid, q(10,1,2) = ',nodeid,q(10,1,2)
      endif

      rhsnorm = 0.0
      do 112 i=1,3
      do 112 k=1,kmax
      do 112 j=1,jmax
         rhsnorm = rhsnorm + s(j,k,i)**2
112   continue
      rhsnorm = sqrt( rhsnorm )

c-----
c  Solve system of equations
c    Line relaxation
c-----

c-----
c  Line relaxation
c-----
      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' node 0 dwb: BEFORE calling lrelax'
       print*, ' node 0 dwb: btc(2,2,1,1), btc(2,2,2,2) =',btc(2,2,1,1),
     &    btc(2,2,2,2)
c       stop 'stop: before call lrelax in step.f/step 10' 
      endif

      if(impsch .eq. 1) then
         call lrelax(jmax,kmax,q,dq,s,
     &               btc,bjm,bjp,bkm,bkp)

      if (DEBUG.and.nodeid.eq.0) then
       print*, ' DWB7: nodeid, q(2,1,1) = ',nodeid,q(2,1,1)
       print*, ' DWB7: nodeid, q(10,1,2) = ',nodeid,q(10,1,2)
       print*, ' DWB7: nodeid, q(32,2,1) = ',nodeid,q(32,2,1)
      endif

c-----
c  lusgs
c-----
      elseif(impsch .eq. 2) then
         neqs = 3
         call lusgs(jmax,kmax,neqs,
     &              q,dq,s,btc,bjm,bjp,bkm,bkp,scr1)
c-----
c  Error in impsch
c-----
      else
         write(istdout,*) ' ERROR in step: impsch out of range.'
         call exit(1)
      endif

      if (DEBUG.and.nodeid.eq.0) then
      print*
      print*, ' DELETE8: nodeid, q,dq(32,2,1) = ',
     &  nodeid,q(32,2,1),dq(32,2,1)
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call flush(istdout)
      call flush(6)
c      stop 'stop: after lusgs but before update in step.f/step'

c-----
c  Update flow variables to next time step 
c-----
      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' node 0 dwb: call update'
       print*
      endif

         call update(jmax,kmax,q,dq,x,y,nt,dcoef2,dcoef4)

      if (DEBUG.and.nodeid.eq.0) then
      print*
      print*, ' DELETE9: nodeid, q,dq(32,2,1) = ',
     & nodeid,q(32,2,1),dq(32,2,1)
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call flush(istdout)
c      stop 'stop: at end of step.f/step'

      if (DEBUG.and.(nodeid.eq.8)) then
       print*, ' DWB8: nodeid, q(2,1,1), nt = ',nodeid,q(2,1,1),nt
       print*, ' DWB8: nodeid, q(2,2,1), nt = ',nodeid,q(2,2,1),nt
       print*, ' DWB8: nodeid, q(10,1,2), nt = ',nodeid,q(10,1,2),nt
       print*, ' DWB8: nodeid, q(10,2,2), nt = ',nodeid,q(10,2,2),nt
       print*, ' DWB8: nodeid, q(32,2,1), nt = ',nodeid,q(32,2,1),nt
       call flush(6)
       call flush(istdout)
c       print*, ' stop: at end of step'
c       stop 'stop: at end of step'
      endif

c-----
c  End of step
c-----

      return
      end
c
c
c***************************************************************
      subroutine resid(jmax,kmax,s,dq,dj,nz)
c***************************************************************
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"

c-----
c  Find max residual for all points and max 
c  divergence of velocity for all interior points
c-----
      __REAL s(jmax,kmax,3), dq(jmax,kmax,3), dj(jmax,kmax)
      __INTEGER nodeid,numprocs,ierr,stat(MPI_STATUS_SIZE)
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, '----------- sub resid -----------'
       print*, 'nt,beta =',nt,beta
       print*, ' numprocs =',numprocs
       print*
       print*, '   dq in step.f/resid'
       print*, '   j   k        dq(1)        dq(2)        dq(3)'
       print*, '  --- --- -------------- -------------- -------------'
      endif

      do 40 k=1,kmax
      do 40 j=1,jmax
         dq(j,k,1) = s(j,k,1)*dj(j,k)
         dq(j,k,2) = s(j,k,2)*dj(j,k)
         dq(j,k,3) = s(j,k,3)*dj(j,k)
c-----
         if (DEBUG.and.nodeid.eq.0) then
          if (j.eq.1.and.k.eq.1) then
           print*, ' -- for nodeid =',nodeid
          endif
          write(*,400) j,k,dq(j,k,1),dq(j,k,2),dq(j,k,3)
c          write(*,400) j,k,s(j,k,1),s(j,k,2),s(j,k,3)
c          write(*,400) j,k,dj(j,k)
400       format(2x,i3,1x,i3,1x,1pe13.5,1x,1pe13.5,1x,1pe13.5)
         endif
c-----
40    continue

      if (DEBUG) then
       stop 'stop in step.f/resid'
      endif

c      irmax = BLAS_IAMAX( 3*jmax*kmax, dq, 1 )
      irmax = isamax( 3*jmax*kmax, dq, 1 )
      jk = (irmax-1)/(jmax*kmax) + 1
      irmax = irmax - (jk-1)*jmax*kmax
      kres = (irmax-1)/jmax + 1
      jres = irmax - (kres-1)*jmax
      resmax = abs( dq(jres,kres,jk) )
      if(nt.eq.1) then
        resmax0 = resmax
        if( resmax0 .lt. 1. ) resmax0 = 1.0
        write(*,*) ' Initial resmax0 for zone ',nz,' : ',resmax0
c send all resmax0 (initial resmax) to node 0; node 0 will determine the
c  max of the resmax0 values to use for normalizing the residuals at 
c  each time step and then MPI_BCAST the result

c------------------------------------
      if (numprocs.gt.1) then
c------------------------------------
c------------------------------------
       if (nodeid.eq.0) then
c------------------------------------
        do i=2,numprocs
         node_recv = i-1

#ifdef D_PRECISION
c         call MPI_RECV(resmax0_temp,1,MPI_REAL8,node_recv,node_recv,
         call MPI_RECV(resmax0_temp,1,MPI_REAL8,node_recv,MPI_ANY_TAG,
     &   MPI_COMM_WORLD,stat,ierr)
#else
c         call MPI_RECV(resmax0_temp,1,MPI_REAL,node_recv,node_recv,
         call MPI_RECV(resmax0_temp,1,MPI_REAL,node_recv,MPI_ANY_TAG,
     &   MPI_COMM_WORLD,stat,ierr)
#endif
         if (resmax0_temp.gt.resmax0) then
          resmax0 = resmax0_temp
         endif
        enddo
c------------------------------------
       else
c------------------------------------

#ifdef D_PRECISION
        call MPI_SEND(resmax0,1,MPI_REAL8,0,nodeid,
     &   MPI_COMM_WORLD,ierr)
#else
        call MPI_SEND(resmax0,1,MPI_REAL,0,nodeid,
     &   MPI_COMM_WORLD,ierr)
#endif

c------------------------------------
       endif
c------------------------------------
c------------------------------------
      endif
c------------------------------------

      endif
      
c send resmax0 result to all nodes
      if (numprocs.gt.1) then

#ifdef D_PRECISION
       call MPI_BCAST(resmax0,1,MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr)
#else
       call MPI_BCAST(resmax0,1,MPI_REAL,0,
     &  MPI_COMM_WORLD,ierr)
#endif

      endif
c
c      idmax = BLAS_IAMAX( jmax*kmax, dq, 1 )
      idmax = isamax( jmax*kmax, dq, 1 )
      k = (idmax-1)/jmax + 1
      j = idmax - (k-1)*jmax
      divmax = abs( dq(j,k,1)/beta )

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' In step.f: j,k,dq(j,k,1),beta = ',j,k,dq(j,k,1),beta
       print*, 'nodeid,irmax,jk,kres,jres,resmax0,resmax,idmax,divmax=',
     &  nodeid,irmax,jk,kres,jres,resmax0,resmax,idmax,divmax
       print*
       call flush(6)
       call flush(istdout)
c       print*, 'stop at end of step.f/resid'
c       stop 'stop at end of step.f/resid'
      endif

c-----
c  End of resid
c-----
      return
      end
c
c
c*************************************************************
      subroutine update(jmax,kmax,q,dq,x,y,nt,dcoef2,dcoef4)
c*************************************************************
#include "precis.h"
#include "mpif.h"
#include "mpi_params.f"
      __REAL q(jmax,kmax,3), dq(jmax,kmax,3)
      __REAL x(jmax,kmax),y(jmax,kmax)
      logical DEBUG,DEBUG1
c-----
c  Update variables to new time step
c-----

      DEBUG = .false.
c      DEBUG = .true.
      DEBUG1 = .false.
c      DEBUG1 = .true.

      if (DEBUG.and.(nodeid.eq.0)) then
       print*
       print*, ' B4 UPDATE: nodeid, q(2,1,1), dq(2,1,1) = ',nodeid,
     &  q(2,1,1),dq(2,1,1)
       print*, ' B4 UPDATE:: nodeid, q(10,1,2), dq(10,1,2) = ',nodeid,
     &  q(10,1,2),dq(10,1,2)
       print*, ' x(2,1),y(2,1) = ',x(2,1),y(2,1)
       print*, ' x(10,1),y(10,1) = ',x(10,1),y(10,1)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntref1 = 1
      ntref2 = 3
      if (DEBUG1.and.numprocs.eq.2.and.nodeid.eq.0.and.
     &  (nt.eq.ntref1.or.nt.eq.ntref2)) then
       print*
       print*, 'START of step.f/update -- updated q for nodeid:',nodeid
       print*, '   jmax,kmax =',jmax,kmax
       print*, '   time step nt =',nt
       j=jmax
       print*, ' nodeid   j   k      q_jk1        q_jk2        q_jk3'
       do 901 k=1,kmax
        write(*,810) nodeid,j,k,q(j,k,1),q(j,k,2),q(j,k,3)
810     format(3i5,1p3e15.7)
901    continue
      endif

      do 100 k=1,kmax
      do 100 j=1,jmax
         q(j,k,1) = q(j,k,1) + dq(j,k,1)
         q(j,k,2) = q(j,k,2) + dq(j,k,2)
         q(j,k,3) = q(j,k,3) + dq(j,k,3)
100   continue

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' AFTR UPDATE: nodeid, q(2,1,1),dq(2,1,1) = ',nodeid,
     &  q(2,1,1),dq(2,1,1)
       print*, ' AFTR UPDATE:: nodeid, q(10,1,2), dq(10,1,2) = ',nodeid,
     &  q(10,1,2),dq(10,1,2)
      endif

      if (DEBUG.and.nt.eq.2) then
       call flush(6)
       call flush(istdout)
       print*, ' time step nt = ',nt
       stop 'stop: in subroutine udpate 10'
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG1.and.numprocs.eq.2.and.nodeid.eq.0.and.
     &  (nt.eq.ntref1.or.nt.eq.ntref2)) then
       print*
       print*, 'END of step.f/update -- updated q for nodeid:',nodeid
       print*, '   jmax,kmax =',jmax,kmax
       print*, '   time step nt =',nt
       j=jmax
       print*, ' nodeid   j   k      q_jk1        q_jk2        q_jk3'
       do 902 k=1,kmax
        write(*,810) nodeid,j,k,q(j,k,1),q(j,k,2),q(j,k,3)
902    continue
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG1.and.numprocs.eq.2.and.
     &  (nt.eq.ntref1.or.nt.eq.ntref2)) then
       call flush(6)
       call flush(istdout)
       if (nt.eq.ntref2) then
        call MPI_FINALIZE(ierr)
        stop 'stop: stopping at end step.f/update'
       endif
      endif 

c-----
c  End of update
c-----

      return
      end
c
c
c********************************************************************
      subroutine print_subgrid(node,layers,x,y,q,jmax,kmax)
c********************************************************************
c Purpose:
c  print grid points and flow variables in outer grid point layers to 
c  check if boundary condition updates are working on the subgrids
c
#include "mpif.h"
#include "mpi_params.f"
c
      __REAL x(jmax,kmax),y(jmax,kmax),q(jmax,kmax,3)
c
      if (node.ne.nodeid) return
c
      print*
      write(*,5) layers,node,node+1
5     format(' Printing ',i3,' outer layers (outermost first) for node '
     &,i4,', subgrid #',i4)

      do 100 layer=1,layers
c start at botton left corner of subrid and work in counter-clockwise 
c  direction

c ... bottom of grid
      print*, '> Bottom of grid, layer #',layer
      write(*,8)
8     format('  j    k    x(j,k)     y(j,k)    q(j,k,1)    q(j,k,2)    q
     &(j,k,3)',/,
     &       ' ---  ---   ------     ------    --------    --------    -
     &-------'
     &)
      k = layer
      do 10 j=layer,jmax - layer + 1
       write(*,15) j,k,x(j,k),y(j,k),q(j,k,1),q(j,k,2),q(j,k,3) 
15     format(1x,i3,2x,i3,f8.3,3x,f8.3,3x,1pe10.3,3x,1pe10.3,3x,1pe10.3)
10    continue

c ... right side of grid
      print*, '> Right side of grid, layer #',layer
      write(*,8)
      j = jmax - layer + 1
      do 20 k=layer,kmax - layer + 1
       write(*,15) j,k,x(j,k),y(j,k),q(j,k,1),q(j,k,2),q(j,k,3)
20    continue

c ... top of grid
      print*, '> Top of grid, layer #',layer
      write(*,8)
      k = kmax - layer + 1
      do 30 j=layer,jmax - layer + 1
       write(*,15) j,k,x(j,k),y(j,k),q(j,k,1),q(j,k,2),q(j,k,3)
30    continue

c ... left side of grid
      print*, '> Left side of grid, layer #',layer
      write(*,8)
      j = layer
      do 40 k=layer,kmax - layer + 1
       write(*,15) j,k,x(j,k),y(j,k),q(j,k,1),q(j,k,2),q(j,k,3)
40    continue

100   continue

      return
      end 
