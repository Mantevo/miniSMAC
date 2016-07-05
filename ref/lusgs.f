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
c  1. lusgs
c  2. lusini3
c  3. lusini1
c
c************************************************************
c
c
c************************************************************************
      subroutine lusgs(jmax,kmax,neqs,q,dq,rhs,
     &                 btc,bjm,bjp,bkm,bkp,scr1)
c************************************************************************
c
c Purpose: solve system of equations using symmetric gauss-seidel
c
c
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      INTEGER nodeid,numprocs,ierr,stat(MPI_STATUS_SIZE)

      __REAL q(jmax,kmax,neqs), dq(jmax,kmax,neqs), rhs(jmax,kmax,neqs),
     &  bjm(jmax,kmax,neqs,neqs), btc(jmax,kmax,neqs,neqs),
     &  bjp(jmax,kmax,neqs,neqs), bkm(jmax,kmax,neqs,neqs), 
     &  bkp(jmax,kmax,neqs,neqs), scr1(jmax,kmax,neqs,neqs)
      __INTEGER jmax, kmax
      logical DEBUG
 
      DEBUG = .false.
c      DEBUG = .true.

c-----
c Zero out dq array
c-----
c      do 9 i=1,neqs
c      do 9 k=1,kmax
c      do 9 j=1,jmax
c        dq(j,k,i) = 0.0
c        rhs(j,k,i) = 0.0
c9     continue

      if (DEBUG.and.nodeid.eq.0.and.neqs.eq.3) then
      print*
      print*, ' btc before lusini3: nodeid,neqs = ',nodeid,neqs
      k=2
      do 800 j=1,jmax
      print*,' DELETE14.1: j,k,btc(j,2,2,2),dq(j,2,2) = ',
     & j,k,btc(j,k,2,2),dq(j,2,2)
800   continue
      print*
      endif

      if (DEBUG.and.neqs.eq.3) then
c      print*, 'nodeid: neqs = ',nodeid,neqs
c      stop 'stop: before calling lusini3'
      endif

c-----
c  Perform lu decomposition of main diagonal block or scalar
c   (generates the modified block array, stored in lu factored form,
c    with diagonals stored as reciprocals)
c-----
      if(neqs .eq. 3) then
            call lusini3(jmax,kmax,btc,scr1)
      elseif(neqs .eq. 1) then
            call lusini1(jmax,kmax,btc,scr1)
      else
        print*
        print*, ' ERROR: neqs value should be 1 or 3'
        print*, '    neqs = ',neqs
        print*, ' stopping'
        print*
        stop ' stopping: neqs value should be 1 or 3'
      endif

      if (DEBUG.and.nodeid.eq.4) then
       print*
       print*, ' DELETE8.1: nodeid,neqs,q,dq(2,2,1),btc(2,2,1,1) = ',
     &  nodeid,neqs,q(2,2,1),dq(2,2,1),btc(2,2,1,1)
      endif

      if (DEBUG.and.neqs.eq.3) then
      print*, 'In lusgs.f -- nodeid: neqs = ',nodeid,neqs
c      stop 'stop: before calling ilusol3'
      endif

c-----
c  Perform lu solution
c-----
      if (neqs.eq.3) then
c NOTE: this is where dq(j,k,n) are calculated
            call ilusol3(jmax,kmax,dq,rhs,scr1,
     &       bjm,bjp,bkm,bkp
     &       )
      else if (neqs.eq.1) then
            call ilusol1(jmax,kmax,dq,rhs,scr1,
     &       bjm,bjp,bkm,bkp
     &       )
      endif

      if (DEBUG.and.nodeid.eq.4) then
       print*
       print*, ' DELETE8.2: nodeid,neqs,q,dq(2,2,1),btc(2,2,1,1) = ',
     &  nodeid,neqs,q(2,2,1),dq(2,2,1),btc(2,2,1,1)
      endif

      if (DEBUG.and.neqs.eq.3) then
      call flush(6)
      call flush(istdout)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      print*, ' Nodeid: neqs = ',nodeid,neqs
c      stop 'stop: stopping before bcimpds'
      endif

c-----
c  Zonal Boundary conditions
c-----
      call bcimpds(jmax,kmax,neqs,q,dq,rhs,
     &        nodeid+1)
c     &                   jmax,kmax,q,dq,rhs)

      if (DEBUG.and.neqs.eq.3) then
      print*, '10 Nodeid: neqs = ',nodeid,neqs
      call flush(istdout)
c      stop 'stop: stopping after bcimpds'
      endif

      if (DEBUG.and.nodeid.eq.4.and.neqs.eq.3) then
       print*
       print*, ' DELETE8.3: nodeid, q,dq(2,2,2),rhs(2,2,2) = ',
     &  nodeid,q(2,2,2),dq(2,2,2),rhs(2,2,2)
       call flush(istdout)
c       stop 'stop: at end of lusgs.f/lusgs'
      endif

c-----
c  End of lusgs
c-----
      return
      end
c
c
c************************************************************************
      subroutine lusini3(jmax,kmax,btc,bprm)
c************************************************************************
c  Vectorized version of this routine.
c
c  Performs ilu(0) decomposition: compute the modified main block
c  diagonal.  Store these blocks in an lu decomposition form, with
c  the diagonal elements of the block stored as reciprocals.
c------------------------------------------------------------------------
#include "precis.h"
      __REAL btc(jmax,kmax,3,3), bprm(jmax,kmax,3,3)
c-----
      do 10 k=1,kmax
      do 10 j=1,jmax
         bprm(j,k,1,1) = 1./ btc(j,k,1,1)
         bprm(j,k,1,2) =     btc(j,k,1,2)*bprm(j,k,1,1)
         bprm(j,k,1,3) =     btc(j,k,1,3)*bprm(j,k,1,1)
         bprm(j,k,2,1) =     btc(j,k,2,1)
         bprm(j,k,2,2) = 1./(btc(j,k,2,2) - bprm(j,k,2,1)*bprm(j,k,1,2))
         bprm(j,k,2,3) =    (btc(j,k,2,3) - bprm(j,k,2,1)*bprm(j,k,1,3))
     &                      *bprm(j,k,2,2)
         bprm(j,k,3,1) =     btc(j,k,3,1)
         bprm(j,k,3,2) =     btc(j,k,3,2) - bprm(j,k,3,1)*bprm(j,k,1,2)
         bprm(j,k,3,3) = 1./(btc(j,k,3,3) - bprm(j,k,3,1)*bprm(j,k,1,3)
     &                                - bprm(j,k,3,2)*bprm(j,k,2,3) )
10    continue
c-----
c  End of lusini3
c-----
      return
      end
c
c
c************************************************************************
      subroutine lusini1(jmax,kmax,btc,bprm)
c************************************************************************
c  Vectorized version of this routine.
c
c  Performs ilu(0) decomposition: compute the modified main
c  diagonal.  Store recipricol of diagonal elements.
c------------------------------------------------------------------------
#include "precis.h"
      __REAL btc(jmax,kmax), bprm(jmax,kmax)
c-----
      do 10 k=1,kmax
      do 10 j=1,jmax
         bprm(j,k) = 1./btc(j,k)
10    continue
c-----
c  End of lusini1
c-----
      return
      end
