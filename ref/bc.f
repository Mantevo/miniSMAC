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
c  1. bcmain 
c  2. bcimpds 
c  3. bcintds 
c  4. bcnoslj
c  5. bcnoslk
c  6. bccgridj
c  7. bccgridk
c  8. bcoutj
c  9. bcoutk
c
c************************************************************
c#include "blasdef.h"

c--------------------------------------------------------------
c     Boundary condition subroutines for SMAC2D 
c
c***********************************************************************
      subroutine bcmain(jmax,kmax,x,y,q,rtxy,dj,nz,s,
     &                  btc,bjm,bjp,bkm,bkp)
c***********************************************************************
c  This routine has master control over all boundary condition
c  routines.  It is called once for each zone at the beginning
c  of each iteration.
c-------------------------------------------------------------
#include "common.f"

      __REAL q(jmax,kmax,3), x(jmax,kmax), y(jmax,kmax),
     & rtxy(jmax,kmax,2,3), dj(jmax,kmax), s(jmax,kmax,3),
     & btc(jmax,kmax,3,3),bjm(jmax,kmax,3,3),bjp(jmax,kmax,3,3),
     & bkm(jmax,kmax,3,3),bkp(jmax,kmax,3,3)
      character*9 bcfile
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

      node = nz - 1

      if (DEBUG.and.node.eq.0) then
       print*
       print*, ' >>> In sub. bc.f/bcmain, nodeid, nbcreg =',
     &   node,nbcreg
       print*

       print*, '  nr nzbc ibcval jkbc jbcb jbce kbcb kbce'
       do 250 nr = 1,nbcreg
        print 240,nr,nzbc(nr),ibcval(nr),jkbc(nr),jbcb(nr),
     &    jbce(nr),kbcb(nr),kbce(nr)
240     format(8i5)
250    continue

      endif

      if (DEBUG) then
       print 290,node,jmax,kmax
290    format('+++ Node, jmax,kmax =',3i4)
c       stop 'stop: stopping near beginning of bc.f/bcmain'
      endif

c-----
c  Process the boundary condition surface if it belongs to
c  the current zone
c-----
      do 100 nr=1,nbcreg
         if (DEBUG.and.node.eq.5) then
          print*,' nr nzbc(nr) node ibcval(nr) ... for node',node
          print 200, nr,nzbc(nr),node,ibcval(nr)
200       format(' nr, nzbc(nr), nz, ibcval(nr) = ',4i4)
         endif

         if(nzbc(nr) .eq. nz) then
c-----
c  No-slip walls
c-----
            if (ibcval(nr) .eq. 0)  then
               if (DEBUG.and.nz.eq.5) then
                print*,'5 no-slip walls: ibcval(nr),jkbc(nr),nzbc,nodeid
     &  =',
     &            ibcval(nr),jkbc(nr),nzbc(nr),nodeid
               endif
               if (jkbc(nr) .eq. 1) then
                  if (DEBUG.and.nz.eq.5) then
                   print*, ' ...... call bcnoslj (j=constant)'
                  endif
                  call bcnoslj(jmax,kmax,q,rtxy,dj,nz,s,btc,bjm,bjp,
     &                         ibcval(nr),jbcb(nr),kbcb(nr),kbce(nr))
               elseif (jkbc(nr) .eq. 2) then
                  if (DEBUG.and.(nz.eq.5.or.nz.eq.9)) then
                   print*, ' ...... 5: call bcnoslk (k=constant)'
                   print*, ' ...... 5: kbeg,jbeg,jend = ',kbcb(nr),
     &                jbcb(nr),jbce(nr)
                  endif
                  call bcnoslk(jmax,kmax,q,rtxy,dj,nz,s,btc,bkm,bkp,
     &                         ibcval(nr),kbcb(nr),jbcb(nr),jbce(nr))
                  if (DEBUG.and.(nz.eq.5.or.nz.eq.9)) then
                   print*,'DWB101: after bcnoslk, node,btc(1,1,1,1) = ',
     &               node,btc(1,1,1,1)
                  endif
               endif

c-----
c  c-grid outer boundary
c-----
            elseif (ibcval(nr) .eq. 25) then
               if (DEBUG.and.nz.eq.5) then
                print*, ' ... 5 c-grid outer boundary: ibcval(nr), jkbc(
     &nr) =',ibcval(nr),jkbc(nr)
               endif
               if (jkbc(nr) .eq. 1) then
c ... j = constant
                  if (DEBUG.and.nz.eq.5) then
                   print*, ' ......call bccgridj'
                  endif
                  call bccgridj(jmax,kmax,x,y,q,rtxy,nz,s,btc,bjm,bjp,
     &                        ibcval(nr),pin,jbcb(nr),kbcb(nr),kbce(nr))
               elseif (jkbc(nr) .eq. 2) then
c ... k = constant
                  if (DEBUG.and.nz.eq.5) then
                   print*, ' ......call bccgridk'
                  endif
                  call bccgridk(jmax,kmax,x,y,q,rtxy,nz,s,btc,bkm,bkp,
     &                        ibcval(nr),pin,kbcb(nr),jbcb(nr),jbce(nr))
               else
                print*
                print*, 'ERROR: jkbc(nr) out of bounds in bc.f'
                print*, '  jkbc(nr) should be 1 or 2'
                print*, '  currently, nr, jkbc(nr) = ',
     &            nr,jkbc(nr)
                print*
                print*, ' ... program terminating'
                stop 'stop: jkbc(nr) out of bounds in bc.f'  
               endif

c-----
c  Outflow
c-----
            elseif (ibcval(nr) .eq. 31) then
               if (DEBUG.and.nz.eq.5) then
                print*, ' ... 5 outflow: ibcval(nr), jkbc(nr) = ',
     &           ibcval(nr),jkbc(nr)
               endif
               if (jkbc(nr) .eq. 1) then
                  if (DEBUG.and.nz.eq.5) then
                   print*, ' ......call bcoutj'
                  endif

                  call bcoutj(jmax,kmax,q,rtxy,dj,nz,s,btc,bjm,bjp,
     &                       ibcval(nr),pout,jbcb(nr),kbcb(nr),kbce(nr))
               elseif (jkbc(nr) .eq. 2) then
                  if (DEBUG.and.nz.eq.5) then
                   print*, ' ......call bcoutk'
                  endif
                  call bcoutk(jmax,kmax,q,rtxy,dj,nz,s,btc,bkm,bkp,
     &                       ibcval(nr),pout,kbcb(nr),jbcb(nr),jbce(nr))
               endif

c-----
c  Error in bc's
c-----
            elseif (ibcval(nr).ne.60) then
              print*
              print*, ' ERROR: invalid ibcval(nr) in sub. bc.f/bcmain'
              print*, '    node, ibcval(',nr,') = ',nz,ibcval(nr)
              print*, ' Program stopping'
              print*
              stop ' stop: error in ibcval in sub. bc.f/bcmain'

            endif
         endif
100   continue

       call flush(6)
       call flush(istdout)
      if (DEBUG) then
c       print*, ' Nodeid = ',nodeid,' stopping in bc.f/bcmain at END'
c       stop 'stop: stopping in bc.f/bcmain at END'
      endif
c-----
c  End of bcmain
c-----
      return
      end
c
c
c**************************************************************
      subroutine bcimpds(jmaxb,kmaxb,neqs,s,ds,rhs,nz)
c Equivalence for lrelax.f, lusgs.f:  sb=q, dsb=dq, rhsb=s, nz=nodeid+1
c Equivalence for lreltur.f: sb=turvar, dsb=ds, rhsb=rhs
c**************************************************************
c  Imposes boundary conditions during implicit schemes.
c  All patched interfaces contained in input files are handled
c  automatically here.
c--------------------------------------------------------------
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      INTEGER nodeid,numprocs,ierr,stat(MPI_STATUS_SIZE)
c
      __REAL s(jmaxb,kmaxb,neqs),ds(jmaxb,kmaxb,neqs),
     &  rhs(jmaxb,kmaxb,neqs)
      __INTEGER jmaxb, kmaxb
      logical DEBUG,DEBUG1,DEBUG2
      __REAL, allocatable:: sb_line(:,:), dsb_line(:,:)
      __REAL, allocatable:: st_line(:,:), dst_line(:,:)
      __REAL, allocatable:: s_b(:,:), ds_b(:,:)
      __REAL, allocatable:: s_t(:,:), ds_t(:,:)
      integer isize_of_message
      logical flag_probe

      DEBUG = .false.
c      DEBUG = .true.
      DEBUG1 = .false.
c      DEBUG1 = .true.
      DEBUG2 = .false.
c      DEBUG2 = .true.

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' >> Entering bc.f/bcimpds'
       print*, '    neqs, nt, numprocs, nodeid = ',neqs,nt,numprocs,
     1  nodeid
       print*
      endif

c-----
c  Zonal boundary interface for patched grids:
c  Directly inject s base values located in the zone which was just
c  updated into the target locations.
c-----

c-----
c  1. Apply interpolation boundary conditions for wake-cut boundary 
c     condition. 
c-----

      if (DEBUG) then
       if (nodeid.eq.0.and.neqs.eq.3) then
        print*
        print*
        print*, ' Entering sub bcimpds, nodeid=0, neqs=3'
        print*, '  nbcreg_wake =',nbcreg_wake
        print*, '  s(10,2,2) =',s(10,2,2)
        print*, '  ds(10,2,2) =',ds(10,2,2)
        call flush(6)
        call flush(istdout)
c        stop 'stop: after entering bc.f/bcimpds'
       endif
      endif
c
c reference time steps used later in this subroutine
      ntref1 = 1 
      ntref2 = 3
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if (DEBUG) then 
       if (nodeid.eq.0) then 
        print*
        print*, ' Node 0: 444 all nodes reached this in bc.f/bcimpds'
        print*
       endif
c use probe to determine status and size of any incoming data;
c  used to check if there are any outstanding messages upto now
c       call MPI_IProbe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,
c     1  flag_probe,stat,ierr)
c       if (flag_probe) then
c         call MPI_Get_Count(stat,MPI_REAL8,isize_of_message,ierr)
c         print*
c         print*,' After top MPI_BARRIER in bc.f/bcimpds:'
c         print*, '  number of elements to RECV = ',isize_of_message
c         call flush(istdout)
c       else
c         print*
c         print*, ' 123 -- No data to receive per MPI_IProbe #1'
c         print*
c       endif
      endif

      if (DEBUG2.and.nodeid.eq.0.and.neqs.eq.3.and.
     &  (nt.eq.ntmax)) then
c check some variables at end of run when nt=ntmax
       print*
       print*, 'START of bcimpds -- Overlap BCs for nodeid:',nodeid
       print*, '  jmaxb,kmaxb =',jmaxb,kmaxb
       print*, '  time step nt =',nt
       print*, '  neqs = ',neqs
       j=jmaxb
       print*, ' nodeid  j  k     q_jk1       q_jk2        q_jk3'
       do 901 k=1,kmaxb
        write(*,810) nodeid,j,k,s(j,k,1),s(j,k,2),s(j,k,3)
901    continue
      endif

c==================================
      if (numprocs.gt.1) then
c==================================

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' For nodeid = ',nodeid
       print*, '   #      nzbc_t_wake     nzbc_b_wake'
       do 38 nireg = 1,nbcreg_wake
        write(6,39) nireg,nzbc_t_wake(nireg),nzbc_b_wake(nireg)
39      format(i5,5x,i8,10x,i8)
38     continue
       print*
c       print*, ' stop: before do 40 in bc.f/bcimpds'
c       stop 'stop: before do 40 in bc.f/bcimpds'
      endif

      do 40 nireg = 1,nbcreg_wake

c--------------------------------------------
       if ((nzbc_t_wake(nireg).eq.nzbc_b_wake(nireg)).and.
     1   nzbc_t_wake(nireg).eq.nz) then
c--------------------------------------------
c no need for mpi calls, since all data is located in same zone
c
c calculate number of interpolation points to use from this grid
          num_pts_base = abs(jbce_b_wake(nireg) - jbcb_b_wake(nireg))
     &      + 1

         jnew = 0

         do 45 jt = jbcb_t_wake(nireg),jbce_t_wake(nireg)
          jnew = jnew + 1
          j_index_b = jbcb_b_wake(nireg) - jnew + 1
         do 45 neq = 1,neqs
          snew = 0.5*(s(jt,2,neq) + s(j_index_b,2,neq))
          dsnew = 0.5*(ds(jt,2,neq) + ds(j_index_b,2,neq))
          ds(jt,1,neq) = snew + dsnew - s(jt,1,neq)
          rhs(jt,1,neq) = ds(jt,1,neq)
45       continue

c--------------------------------------------
       else if (nzbc_t_wake(nireg).eq.nz) then
c--------------------------------------------
c receive base values from corresponding subgrid
         node_from = nzbc_b_wake(nireg) - 1
c calculate number of interpolation points to receive in this subgrid
         num_pts_target = abs(jbce_t_wake(nireg) - jbcb_t_wake(nireg)) 
     &    + 1
         if (DEBUG) then
          print*
          print*, '** B4 RECV: Nodeid,node_from,num_pts_target,neqs = ',
     1       nodeid,node_from,num_pts_target,neqs
          print*, ' -- B4 RECV: nireg, nbcreg_wake = ',nireg,nbcreg_wake
          print*, ' -- B4 RECV: nzbc_t_wake,nz = ',
     1     nzbc_t_wake(nireg),nz
          call flush(istdout)
          print*
c          print*, ' stop in bc.f/bcimpds before r1 mpi_recv'
c          stop 'stop in bc.f/bcimpds before r1 mpi_recv'
         endif

c allocate variables now that number of points are known
         allocate (s_t(num_pts_target,neqs),ds_t(num_pts_target,neqs))

c use probe to determine status and size of incoming data
c        if (DEBUG) then
c         call MPI_IProbe(node_from,MPI_ANY_TAG,MPI_COMM_WORLD,
c     & flag_probe,stat,ierr)
c        if (flag_probe) then
c         call MPI_Get_Count(stat,MPI_REAL8,isize_of_message,ierr)
c         print*
c         print*,' After MPI_PROBE and MPI_GET_COUNT b4 MPI_RECV:'
c         print*, '  number of elements to RECV = ',isize_of_message
c         call flush(istdout)
c        else
c         print*
c         print*, ' 456 -- No data to receive per MPI_IProbe #2'
c         print*
c        endif
c        endif

c receive s_b (q) values at k=2 from 'node_from' for s_t target values
#ifdef D_PRECISION
          call MPI_RECV(s_t,neqs*num_pts_target,MPI_REAL8,node_from,
     &     node_from,MPI_COMM_WORLD,stat,ierr)
#else
          call MPI_RECV(s_t,neqs*num_pts_target,MPI_REAL,node_from,
     &     node_from,MPI_COMM_WORLD,stat,ierr)
#endif

      if (DEBUG) then
       print*
       print*, ' **R1**: nodeid =',nodeid
       print*, ' s_t for R1 = ',s_t
       call flush(istdout)
       print*
      endif

c receive ds_b (dq) values at k=2 from 'node_from' for ds_t target values
#ifdef D_PRECISION
         call MPI_RECV(ds_t,neqs*num_pts_target,MPI_REAL8,node_from,
     &    node_from + 1,MPI_COMM_WORLD,stat,ierr)
#else
         call MPI_RECV(ds_t,neqs*num_pts_target,MPI_REAL,node_from,
     &    node_from + 1,MPI_COMM_WORLD,stat,ierr)
#endif

      if (DEBUG) then
       print*
       print*, ' **R2**: nodeid =',nodeid
       print*, ' ds_t for R2 = ',ds_t
       call flush(istdout)
       print*
      endif

c compute bcs on target surface using target and base values
c         do 41 ni = 1,num_pts_target
        if (DEBUG1.and.neqs.eq.3) then
         if (nodeid.eq.0) then
          print*, ' 0:ds(10,2,1,2,3) =',ds(10,2,1),ds(10,2,2),ds(10,2,3)
         endif
         if (nodeid.eq.1) then
          print*, ' 1:ds(51,2,1,2,3) =',ds(51,2,1),ds(51,2,2),ds(51,2,3)
         endif
        endif
c
        if (DEBUG.and.nodeid.eq.0) then
         print*
         print*, ' nodeid, ntime = ',nodeid,ntime
         print*, 'j  k neq nireg    snew   dsnew   s    rhs'
        endif
c
         jnew = 0

         do 41 jt = jbcb_t_wake(nireg),jbce_t_wake(nireg)
          jnew = jnew + 1
         do 41 neq = 1,neqs
          snew = 0.5*(s(jt,2,neq) + s_t(jnew,neq))
          dsnew = 0.5*(ds(jt,2,neq) + ds_t(jnew,neq))
          ds(jt,1,neq) = snew + dsnew - s(jt,1,neq)
          rhs(jt,1,neq) = ds(jt,1,neq)
41       continue

       if (DEBUG1.and.neqs.eq.3) then
        if (nodeid.eq.0) then
         print*, ' 0: dsinterp(10,1,123) =',ds(10,1,1),ds(10,1,2),
     &    ds(10,1,3)
        endif
        if (nodeid.eq.1) then
         print*, ' 1: dsinterp(51,1,123) =',ds(51,1,1),ds(51,1,2),
     &    ds(51,1,3)
        endif
       endif

c deallocate s_t, ds_t
         deallocate (s_t,ds_t)

c--------------------------------------------
        else if (nzbc_b_wake(nireg).eq.nz) then
c--------------------------------------------
c send base values to target grid
          node_to = nzbc_t_wake(nireg) - 1
c calculate number of interpolation points to send from this subgrid
          num_pts_base = abs(jbce_b_wake(nireg) - jbcb_b_wake(nireg))
     &      + 1

c allocate variables now that number of points are known
          allocate (s_b(num_pts_base,neqs),ds_b(num_pts_base,neqs))

c extract s and ds values at k=2
          jnew = 0

          do 44 jb=jbcb_b_wake(nireg),jbce_b_wake(nireg),-1
          jnew = jnew + 1
          do 44 neq = 1,neqs
           s_b(jnew,neq) = s(jb,2,neq)
           ds_b(jnew,neq) = ds(jb,2,neq)
44        continue
           
         if (DEBUG) then
          print*
          print*, '** B4 SEND: Nodeid, node_to, num_pts_base, neqs = ',
     1     nodeid,node_to,num_pts_base,neqs
          print*, ' -- B4 SEND: nireg, nbcreg_wake = ',nireg,nbcreg_wake
          print*
         endif

c send s_b values at k=2 in this base subgrid
          if (DEBUG) then
           print*, ' **S1** neqs,num_pts_base,node_to,nodeid = ',
     &      neqs,num_pts_base,node_to,nodeid
           print*, ' s_b for S1 = ',s_b
           call flush(istdout)
          endif

#ifdef D_PRECISION
          call MPI_SEND(s_b,neqs*num_pts_base,MPI_REAL8,node_to,
     &     nodeid,MPI_COMM_WORLD,ierr)
#else
          call MPI_SEND(s_b,neqs*num_pts_base,MPI_REAL,node_to,
     &     nodeid,MPI_COMM_WORLD,ierr)
#endif

c send ds_b values at k=2 in this base subgrid 
          if (DEBUG) then 
           print*, ' **S2** neqs,num_pts_base,node_to,nodeid+1 =',
     &      neqs,num_pts_base,node_to,nodeid+1
           print*, ' ds_b for S2 = ',ds_b
           call flush(istdout)
          endif

#ifdef D_PRECISION
          call MPI_SEND(ds_b,neqs*num_pts_base,MPI_REAL8,node_to,
     &     nodeid+1,MPI_COMM_WORLD,ierr)
#else
          call MPI_SEND(ds_b,neqs*num_pts_base,MPI_REAL,node_to,
     &     nodeid+1,MPI_COMM_WORLD,ierr)
#endif

c deallocate s_b, ds_b
          deallocate (s_b,ds_b)

c--------------------------------------------
       endif
c--------------------------------------------

40    continue


c==================================
      else
c==================================
c for 1 processor; no MPI calls needed


      do 42 nireg = 1,nbcreg_wake
c calculate number of interpolation points to use from this grid
          num_pts_base = abs(jbce_b_wake(nireg) - jbcb_b_wake(nireg))
     &      + 1

         jnew = 0

         do 43 jt = jbcb_t_wake(nireg),jbce_t_wake(nireg)
          jnew = jnew + 1
          j_index_b = jbcb_b_wake(nireg) - jnew + 1
         do 43 neq = 1,neqs
c          snew = 0.5*(s(jt,2,neq) + s_t(jnew,neq))
          snew = 0.5*(s(jt,2,neq) + s(j_index_b,2,neq))

c          dsnew = 0.5*(ds(jt,2,neq) + ds_t(jnew,neq))
          dsnew = 0.5*(ds(jt,2,neq) + ds(j_index_b,2,neq))
          ds(jt,1,neq) = snew + dsnew - s(jt,1,neq)
          rhs(jt,1,neq) = ds(jt,1,neq)
43       continue

42      continue

c==================================
      endif
c==================================

      if (DEBUG) then 
       print*, ' > MPI_BARRIER for nodeid = ',nodeid
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      endif

      if (DEBUG) then 
       print*, ' All nodes reached this 500 , nodeid = ',nodeid
       stop 'stop: after MPI_BARRIER in bc.f/bcimpds 500'
      endif

      if (DEBUG.and.neqs.eq.3.and.nodeid.eq.0) then
       print*
       print*, ' After equivalent to call bcimpds:'
       print*, '  s(10,2,2) = ',s(10,2,2)
       print*, '  ds(10,2,2) = ',ds(10,2,2)
       print*
       print*, ' In bc.f/bcimpds, rhs from bcintds = '
       print*, '   j    k      rhs_jk1         rhs_jk2        rhs_jk3
     &       dsb_jk2'
       k = 2
       do 805 j=1,17
        write(*,905) j,k,rhs(j,k,1),rhs(j,k,2),rhs(j,k,3),
     &   ds(j,k,2)
905     format(2i5,1p4e13.5)
805    continue
c       stop 'stop: in bc.f/bcintds - neqs = 3'
      endif


      if (DEBUG.and.nodeid.eq.0) then
c        print*, ' Nodeid, nintreg = ',nodeid,nintreg
        print*, ' After imposing wake bcs for Nodeid, nbcreg_wake = ',
     &    nodeid,nbcreg_wake
        print*, ' ... rhs(10,2,1) = ',rhs(10,2,1)
      endif


      if (DEBUG.and.nodeid.eq.0.and.neqs.eq.3) then
       print*
       print*,' In bc.f/bcimpds, rhs from bcintds = '
       print*, ' nid  j   k    rhs_jk1      rhs_jk2      rhs_jk3'
       k = 2
c       do 800 j=1,jmaxb
       do 800 j=1,17
        write(*,900) nz,j,k,rhs(j,k,1),rhs(j,k,2),rhs(j,k,3),
     &   ds(j,k,1)
900     format(3i5,1p4e13.5)
800    continue
       call flush(6)
       call flush(istdout)
c       stop 'stop: In bc.f/bcimpds, after call bcintds'
      endif


      if (DEBUG.and.neqs.eq.3) then
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call flush(6)
       print*, '20 Nodeid: neqs = ',nodeid,neqs
       print*, ' nreg,nz = ',nreg,nz
       if (nodeid.eq.0) then
        print*
        print*, '  nr   nzt(nr)    nzb(nr)'
        do 1000 nr=1,nreg
        write(*,1010) nr,nzt(nr),nzb(nr)
1010    format(1x,i4,3x,i5,5x,i5)
1000    continue
       endif
       call flush(istdout)
       call flush(6)
c       stop 'stop: stopping after bcintds in bc.f/bcimpds'
      endif

c=====================================
      if (numprocs.gt.1) then
c=====================================
c following not needed for 1 processor
c
c for overlapped interior subgrids
c pre-post receives
      do 50 nr=1,nreg
c       if (nzt(nr).eq.nz) then
c++++++++++++++++++++++++++++++++++++++
       if (nodeid.eq.nzt(nr)-1) then
c++++++++++++++++++++++++++++++++++++++
        node_in = nzb(nr)-1  ! node_in corresponds to MPI rank
        if (DEBUG) then
         print*
         print*, ' In bc.f/bcimpds: nr, node_in = ',nr,node_in
        endif
c  note that either jmax or kmax, or possibly both, will be same in 
c   target and base grids
c  note that either jeb=jbb or keb=kbb
        if (jeb(nr).eq.jbb(nr)) then
c ... constant j
          ipoints=kmaxb
          iflag=1
        elseif (keb(nr).eq.kbb(nr)) then
c ... constant k
          ipoints=jmaxb
          iflag=2
        else
c ... could not find constant j or constant k ==> error
       write(*,500) nodeid,node_out,nr,jeb(nr),jbb(nr),keb(nr),kbb(nr)
        endif

      if (DEBUG) then
       print*, ' --- For nodeid, ipoints,iflag,neqs = ',nodeid,ipoints,
     &  iflag,neqs
      endif

c now that we know how many points involved, allocate memory
        allocate ( st_line(ipoints,neqs), dst_line(ipoints,neqs))

c get values from line of interior points in base grid
        if (DEBUG) then
         print*
         print*,' In bc.f/bcimpds: calling first MPI_RECV, nodeid= ',
     1    nodeid 
         call flush(6)
        endif

#ifdef D_PRECISION
        call MPI_RECV(st_line,ipoints*neqs,MPI_REAL8,node_in,node_in+1,
     &        MPI_COMM_WORLD,stat,ierr)
#else
        call MPI_RECV(st_line,ipoints*neqs,MPI_REAL,node_in,node_in+1,
     &        MPI_COMM_WORLD,stat,ierr)
#endif

c get delta values from line of interior points in base grid
        if (DEBUG) then
         print*
         print*,' In bc.f/bcimpds: calling second MPI_RECV ',ipoints
         call flush(6)
        endif

#ifdef D_PRECISION
        call MPI_RECV(dst_line,ipoints*neqs,MPI_REAL8,node_in,node_in+2,
     &        MPI_COMM_WORLD,stat,ierr)
#else
        call MPI_RECV(dst_line,ipoints*neqs,MPI_REAL,node_in,node_in+2,
     &        MPI_COMM_WORLD,stat,ierr)
#endif

c
c update s and ds values

!$omp parallel do private(ip)
        do 600 n=1,neqs
         ip = 0
        do 600 kt=kbt(nr),ket(nr)
        do 600 jt=jbt(nr),jet(nr)
          ip = ip + 1
          s(jt,kt,n) = st_line(ip,n)
          ds(jt,kt,n) = dst_line(ip,n)
          rhs(jt,kt,n) = ds(jt,kt,n)
600      continue

!$omp end parallel do

        deallocate (st_line,dst_line)

c       if (nzb(nr).eq.nz) then
c++++++++++++++++++++++++++++++++++++++
       elseif (nodeid.eq.nzb(nr)-1) then
c++++++++++++++++++++++++++++++++++++++
c sends
        node_out = nzt(nr)-1  ! node_out corresponds to MPI rank
        if (DEBUG) then
         print*
         print*, ' In bc.f/bcimpds: nr, node_out = ',nr,node_out
        endif
c note that either jmax or kmax, or possibly both, will be same in 
c  target and base grids
c note that either jeb=jbb or keb=kbb
        if (jeb(nr).eq.jbb(nr)) then
c ... constant j
         ipoints=kmaxb
        elseif (keb(nr).eq.kbb(nr)) then
c ... constant k
         ipoints=jmaxb
        else
c ... could not find constant j or constant k ==> error
       write(*,500) nodeid,node_out,nr,jeb(nr),jbb(nr),keb(nr),kbb(nr)
500    format('//Nodeid: ',i4,' to node ',i4,' for nr = ',i4,
     & '/  indexing error: jeb <> jbb and keb <> kbb',
     & '/    jeb = ',i4,'  jbb = ',i4,'  keb = ',i4,'  kbb = ',i4,
     & '/stopping execution.//'
     & )
        endif

        if (DEBUG) then
         print*, ' For node_out = ',node_out,' jmaxb,kmaxb = ',jmaxb,
     &     kmaxb,' ipoints = ',ipoints
         call flush(6)
        endif

c now that we know how many points are involved, allocate memory
        allocate ( sb_line(ipoints,neqs), dsb_line(ipoints,neqs))

        if (DEBUG) then
         print*,'999 For node_out =', node_out,', after allocate for sb'
         print*
         print*, ' nr,kbb(nr),keb(nr) =',
     &             nr,kbb(nr),keb(nr)
         print*, ' nr,jbb(nr),jeb(nr) =',
     &             nr,jbb(nr),jeb(nr)
        endif


!$omp parallel do private(ip)

        do 100 n=1,neqs
         ip = 0
        do 100 kb=kbb(nr),keb(nr)
        do 100 jb=jbb(nr),jeb(nr)
         ip = ip + 1
         sb_line(ip,n) = s(jb,kb,n)
         dsb_line(ip,n) = ds(jb,kb,n)
100     continue

!$omp end parallel do

        if (DEBUG) then
         print*
         do i=1,ipoints
          print*, ' sb_line(',i,',1) = ',sb_line(i,1)
         enddo
        endif

c send values from line of interior points in base grid
        if (DEBUG) then
         print*
         print*, ' In node_out = ',node_out,' ipoints = ',ipoints
         print*
         print*,' In bc.f/bcimpds: calling first MPI_SEND ',ipoints
         call flush(6)
        endif
#ifdef D_PRECISION
       call MPI_SEND(sb_line,ipoints*neqs,MPI_REAL8,node_out,nodeid+1,
     &       MPI_COMM_WORLD,ierr)
#else
       call MPI_SEND(sb_line,ipoints*neqs,MPI_REAL,node_out,nodeid+1,
     &       MPI_COMM_WORLD,ierr)
#endif

c send delta values from line of interior points in base grid
        if (DEBUG) then
         print*
         print*,' In bc.f/bcimpds: calling second MPI_SEND ',ipoints
         call flush(6)
        endif

#ifdef D_PRECISION
       call MPI_SEND(dsb_line,ipoints*neqs,MPI_REAL8,node_out,nodeid+2,
     &       MPI_COMM_WORLD,ierr)
#else
       call MPI_SEND(dsb_line,ipoints*neqs,MPI_REAL,node_out,nodeid+2,
     &       MPI_COMM_WORLD,ierr)
#endif

        deallocate (sb_line, dsb_line)

c++++++++++++++++++++++++++++++++++++++
       endif
c++++++++++++++++++++++++++++++++++++++

50    continue

c check overlapped bcs for the airfoil two node case; keep MPI_BARRIERs
c  as they make sure printout is ordered as shown below

c ... print j=jmax values for node 0 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG) then
        print*
        print*, ' 88 In bc.f/bcimpds, This is nodeid,neqs = ',nodeid,neqs
c        stop 'stop: after last MPI_BARRIER in bc.f/bcimpds'
      endif

c ... print j=jmax values for node 0
      if (DEBUG2.and.numprocs.eq.2.and.nodeid.eq.0.and.neqs.eq.3.and.
     &  (nt.eq.ntref1.or.nt.eq.ntref2)) then
       print*
       print*, 'End of bc.f/bcimpds -- Overlap BCs for nodeid:',nodeid
       print*, '  jmaxb,kmaxb = ',jmaxb,kmaxb
       print*, ' time step nt = ',nt
       j=jmaxb
       print*, ' nodeid  j  k     q_jk1       q_jk2        q_jk3'
       do 910 k=1,kmaxb
        write(*,810) nodeid,j,k,s(j,k,1),s(j,k,2),s(j,k,3)
810     format(3i5,1p3e15.7)
910    continue
      endif
 
c ... print j=2 values for node 1, which should match j=jmax for node 0
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG2.and.numprocs.eq.2.and.nodeid.eq.1.and.neqs.eq.3.and.
     &  (nt.eq.ntref1.or.nt.eq.ntref2)) then
       print*
       print*, 'End of bc.f/bcimpds -- Overlap BCs for nodeid:',nodeid
       print*, '  jmaxb,kmaxb = ',jmaxb,kmaxb
       print*, ' time step nt = ',nt
       j=29
       print*, ' nodeid  j  k     q_jk1       q_jk2        q_jk3'
       do 915 k=1,kmaxb
        write(*,810) nodeid,j,k,s(j,k,1),s(j,k,2),s(j,k,3)
915    continue
      endif

c ... print j=jmax-1 values for node 0 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG2.and.numprocs.eq.2.and.nodeid.eq.0.and.neqs.eq.3.and.
     &  (nt.eq.ntref1.or.nt.eq.ntref2)) then
       print*
       print*, 'End of bc.f/bcimpds -- Overlap BCs for nodeid:',nodeid
       print*, '  jmaxb,kmaxb = ',jmaxb,kmaxb
       print*, ' time step nt = ',nt
       j=jmaxb-1
       print*, ' nodeid  j  k     q_jk1       q_jk2        q_jk3'
       do 920 k=1,kmaxb
        write(*,810) nodeid,j,k,s(j,k,1),s(j,k,2),s(j,k,3)
920    continue
      endif

c ... print j=1 values for node 1, which should be same for j=jmax-1 for node 0 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG2.and.numprocs.eq.2.and.nodeid.eq.1.and.neqs.eq.3.and.
     &  (nt.eq.ntref1.or.nt.eq.ntref2)) then
       print*
       print*, 'End of bc.f/bcimpds -- Overlap BCs for nodeid:',nodeid
       print*, '  jmaxb,kmaxb = ',jmaxb,kmaxb
       print*, ' time step nt = ',nt
       j=1
       print*, ' nodeid  j  k     q_jk1       q_jk2        q_jk3'
       do 925 k=1,kmaxb
        write(*,810) nodeid,j,k,s(j,k,1),s(j,k,2),s(j,k,3)
925    continue
      endif

c ... print j=1 values for node 0 to check boundary conditions 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG2.and.numprocs.eq.2.and.nodeid.eq.0.and.neqs.eq.3.and.
     &  (nt.eq.ntref1.or.nt.eq.ntref2)) then
       print*
       print*, 'End of bc.f/bcimpds -- Overlap BCs for nodeid:',nodeid
       print*, '  jmaxb,kmaxb = ',jmaxb,kmaxb
       print*, ' time step nt = ',nt
       j=1
       print*, ' nodeid  j  k     q_jk1       q_jk2        q_jk3'
       do 930 k=1,kmaxb
        write(*,810) nodeid,j,k,s(j,k,1),s(j,k,2),s(j,k,3)
930    continue
      endif

c ... print j=jmaxb/2 values for node 0 for sanity check
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG2.and.numprocs.eq.2.and.nodeid.eq.0.and.neqs.eq.3.and.
     &  (nt.eq.ntref1.or.nt.eq.ntref2)) then
       print*
       print*, 'End of bc.f/bcimpds -- Overlap BCs for nodeid:',nodeid
       print*, '  jmaxb,kmaxb = ',jmaxb,kmaxb
       print*, '  time step nt = ',nt
       j=jmaxb/2
       print*, ' nodeid  j  k     q_jk1       q_jk2        q_jk3'
       do 935 k=1,kmaxb
        write(*,810) nodeid,j,k,s(j,k,1),s(j,k,2),s(j,k,3)
935    continue
      endif

c ... print j=jmaxb/2 values for node 1 for sanity check 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG2.and.numprocs.eq.2.and.nodeid.eq.1.and.neqs.eq.3.and.
     &  (nt.eq.ntref1.or.nt.eq.ntref2)) then
       print*
       print*, 'End of bc.f/bcimpds -- Overlap BCs for nodeid:',nodeid
       print*, '  jmaxb,kmaxb = ',jmaxb,kmaxb
       print*, ' time step nt = ',nt
       j=jmaxb/2
       print*, ' nodeid  j  k     q_jk1       q_jk2        q_jk3'
       do 940 k=1,kmaxb
        write(*,810) nodeid,j,k,s(j,k,1),s(j,k,2),s(j,k,3)
940    continue
      endif

c ... print j=jmaxb values for node 1 for check for boundary conditions
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG2.and.numprocs.eq.2.and.nodeid.eq.1.and.neqs.eq.3.and.
     &  (nt.eq.ntref1.or.nt.eq.ntref2)) then
       print*
       print*, 'END of bc.f/bcimpds -- Overlap BCs for nodeid:',nodeid
       print*, '  jmaxb,kmaxb =',jmaxb,kmaxb
       print*, '  time step nt =',nt
       j=jmaxb
       print*, ' nodeid  j  k     q_jk1       q_jk2        q_jk3'
       do 945 k=1,kmaxb
        write(*,810) nodeid,j,k,s(j,k,1),s(j,k,2),s(j,k,3)
945    continue
      endif

c=====================================
      endif
c=====================================

      if (DEBUG2.and.neqs.eq.3.and.
     & (nt.eq.ntref1.or.nt.eq.ntref2)) then
       call flush(6)
       call flush(istdout)
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c       if (nt.eq.ntref2) then
c       call MPI_FINALIZE(ierr)
c       print*, '30 Nodeid: neqs = ',nodeid,neqs
c       stop 'stop: stopping at end bc.f/bcimpds'
c       endif
      elseif (DEBUG2) then
c       print*, ' stop: at end of bc.f/bcimpds'
c       stop ' stop: at end of bc.f/bcimpds'
      endif

c-----
c  End of bcimpds
c-----

      return
      end
c
c
c***********************************************************************
      subroutine bcintds(jmax,kmax,neqs,s,ds,rhs,nireg,nz,
     &                   jintt,kintt,jintb1,kintb1,jintb2,kintb2)
c***********************************************************************
c  Handles boundary conditions at computational boundaries where the
c  values may be interpolated from any 3 different points in the zone.
c  Values of s & ds interpolated from the 3 base points are injected
c  into the target point for the solution vector.
c-----------------------------------------------------------------------
#include "common.f"
      __REAL s(jmax,kmax,neqs), ds(jmax,kmax,neqs), rhs(jmax,kmax,neqs)
      __INTEGER jintt(*), kintt(*),
     &        jintb1(*), kintb1(*),
     &        jintb2(*), kintb2(*)
      logical DEBUG
      
      DEBUG = .false.
c      DEBUG = .true.

      nodeid = nz - 1
      if (DEBUG.and.nodeid.eq.0) then
       print*
c       print*, 'j  k  n  ni nireg'
       print*, 'j  k  n  ni nireg'
      endif
c

!$omp parallel do private(j,k,dsnew,snew)
      do 10 n=1,neqs
      do 10 ni=iintbeg(nireg),iintend(nireg)
         j = jintt(ni)
         k = kintt(ni)
c for interpolation, just average points above and below target point
c   on C-grid wake cut
         dsnew = 0.5*(ds(jintb1(ni),kintb1(ni),n) +
     &                 ds(jintb2(ni),kintb2(ni),n))
         snew = 0.5*(s(jintb1(ni),kintb1(ni),n) +
     &                s(jintb2(ni),kintb2(ni),n))
         ds(j,k,n) = snew + dsnew - s(j,k,n)
         rhs(j,k,n) = ds(j,k,n)
10    continue

!$omp end parallel do

      if (DEBUG.and.neqs.eq.3) then
       stop 'stop: in bc.f/bcintds - neqs = 3'
      endif

c-----
c  End of bcintds
c-----
      return
      end  
c
c
c**********************************************************************
      subroutine bcnoslj(jmax,kmax,q,rtxy,dj,nz,s,btc,bjm,bjp,
     &                   ibval,j,kbeg,kend)
c**********************************************************************
c imposes no-slip conditions along j=constant, in j direction
c-------------------------------------------------------------
#include "common.f"
      __REAL q(jmax,kmax,3), rtxy(jmax,kmax,2,3), dj(jmax,kmax),
     &       s(jmax,kmax,3), btc(jmax,kmax,3,3),
     &       bjm(jmax,kmax,3,3), bjp(jmax,kmax,3,3)
c-----
c  No-slip with flow at j+1
c-----

      if(ibval .eq. 0) then

!$omp parallel do 
         do 10 k=kbeg,kend
            btc(j,k,1,1) = 1.0
            btc(j,k,2,2) = 1.0
            btc(j,k,3,3) = 1.0
            bjp(j,k,1,1) =-1.0
            s(j,k,1) = q(j+1,k,1) - q(j,k,1)
            s(j,k,2) = -q(j,k,2)
            s(j,k,3) = -q(j,k,3)
10       continue
!$omp end parallel do

c-----
c  No-slip with flow at j-1
c-----
      else

!$omp parallel do 
         do 20 k=kbeg,kend
            btc(j,k,1,1) = 1.0
            btc(j,k,2,2) = 1.0
            btc(j,k,3,3) = 1.0
            bjm(j,k,1,1) =-1.0
            s(j,k,1) = q(j-1,k,1) - q(j,k,1)
            s(j,k,2) = -q(j,k,2)
            s(j,k,3) = -q(j,k,3)
20       continue
!$omp end parallel do

      endif

c-----
c  End of bcnoslj
c-----
      return
      end
c
c
c**********************************************************************
      subroutine bcnoslk(jmax,kmax,q,rtxy,dj,nz,s,btc,bkm,bkp,
     &                   ibval,k,jbeg,jend)
c**********************************************************************
c imposes no-slip conditions along k=constant, in k direction
c-------------------------------------------------------------
#include "common.f"
      __REAL q(jmax,kmax,3), rtxy(jmax,kmax,2,3), dj(jmax,kmax),
     &       s(jmax,kmax,3), btc(jmax,kmax,3,3),
     &       bkm(jmax,kmax,3,3), bkp(jmax,kmax,3,3)
c-----
c  No-slip with flow at k+1
c-----

      if(ibval .eq. 0) then

!$omp parallel do
         do 10 j=jbeg,jend
            btc(j,k,1,1) = 1.0
            btc(j,k,2,2) = 1.0
            btc(j,k,3,3) = 1.0
            bkp(j,k,1,1) =-1.0
            s(j,k,1) = q(j,k+1,1) - q(j,k,1)
            s(j,k,2) = -q(j,k,2)
            s(j,k,3) = -q(j,k,3)
10       continue
!$omp end parallel do

c-----
c  No-slip with flow at k-1
c-----
      else

!$omp parallel do 
         do 20 j=jbeg,jend
            btc(j,k,1,1) = 1.0
            btc(j,k,2,2) = 1.0
            btc(j,k,3,3) = 1.0
            bkm(j,k,1,1) =-1.0
            s(j,k,1) = q(j,k-1,1) - q(j,k,1)
            s(j,k,2) = -q(j,k,2)
            s(j,k,3) = -q(j,k,3)
20       continue
!$omp end parallel do

      endif

c-----
c  End of bcnoslk
c-----
      return
      end
c
c
c**********************************************************************
      subroutine bccgridj(jmax,kmax,x,y,q,rtxy,nz,s,btc,bjm,bjp,
     &                    ibval,pres,j,kbeg,kend)
c**********************************************************************
#include "common.f"
      __REAL q(jmax,kmax,3), x(jmax,kmax), y(jmax,kmax),
     &     s(jmax,kmax,3), btc(jmax,kmax,3,3), rtxy(jmax,kmax,2,3),
     &   bjm(jmax,kmax,3,3), bjp(jmax,kmax,3,3)
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.
c-----
c  Freestream velocities at angle of attack
c-----
      jk = 1
      pi = 4.d0*datan( 1.d0 )
      uinf = cos( alpha*pi/180. )
      vinf = sin( alpha*pi/180. )
c-----
c  Point vortex velocity field
c-----

!$omp parallel

!$omp do private(dx,dy,ucirc,vcirc)
      do 10 k=kbeg,kend
         dx = x(j,k) - 0.25
         dy = y(j,k)
         ucirc = clift*dy/( 4.*pi*(dx*dx + dy*dy) )
         vcirc =-clift*dx/( 4.*pi*(dx*dx + dy*dy) )
c
         btc(j,k,2,2) = 1.0
         btc(j,k,3,3) = 1.0
         s(j,k,2) = uinf+ucirc - q(j,k,2)
         s(j,k,3) = vinf+vcirc - q(j,k,3)
10    continue
!$omp end do

c-----
c  j = 1: characteristic relation for pressure
c-----
      if(ibval .eq. 25) then
         if(j .eq. 1) then

!$omp do private(qq,cc,dq1,dq2,dq3)
            do 20 k=kbeg,kend
               qq = rtxy(j,k,jk,2)*q(j,k,2) + rtxy(j,k,jk,3)*q(j,k,3)
               cc = sqrt( qq*qq + beta*(rtxy(j,k,jk,2)**2
     &                                + rtxy(j,k,jk,3)**2) )
               bjp(j,k,1,1) = rtxy(j,k,jk,2)**2 + rtxy(j,k,jk,3)**2
               bjp(j,k,1,2) = rtxy(j,k,jk,2)*(qq - cc)
               bjp(j,k,1,3) = rtxy(j,k,jk,3)*(qq - cc)
c
               btc(j,k,1,1) =-bjp(j,k,1,1)
               btc(j,k,1,2) =-bjp(j,k,1,2)
               btc(j,k,1,3) =-bjp(j,k,1,3)
c
               dq1 = q(j+1,k,1) - q(j,k,1)
               dq2 = q(j+1,k,2) - q(j,k,2)
               dq3 = q(j+1,k,3) - q(j,k,3)
               s(j,k,1) = - bjp(j,k,1,1)*dq1 - bjp(j,k,1,2)*dq2
     &                    - bjp(j,k,1,3)*dq3
20          continue
!$omp end do

c-----
c  j = jmax: characteristic relation for pressure
c-----
         elseif(j .eq. jmax) then

!$omp do private(qq,cc,dq1,dq2,dq3)
            do 30 k=kbeg,kend
               qq = rtxy(j,k,jk,2)*q(j,k,2) + rtxy(j,k,jk,3)*q(j,k,3)
               cc = sqrt( qq*qq + beta*(rtxy(j,k,jk,2)**2
     &                                + rtxy(j,k,jk,3)**2) )
               bjm(j,k,1,1) =-(rtxy(j,k,jk,2)**2 + rtxy(j,k,jk,3)**2)
               bjm(j,k,1,2) =-rtxy(j,k,jk,2)*(qq + cc)
               bjm(j,k,1,3) =-rtxy(j,k,jk,3)*(qq + cc)
c
               btc(j,k,1,1) =-bjm(j,k,1,1)
               btc(j,k,1,2) =-bjm(j,k,1,2)
               btc(j,k,1,3) =-bjm(j,k,1,3)
c
               dq1 = q(j,k,1) - q(j-1,k,1)
               dq2 = q(j,k,2) - q(j-1,k,2)
               dq3 = q(j,k,3) - q(j-1,k,3)
               s(j,k,1) = bjm(j,k,1,1)*dq1 + bjm(j,k,1,2)*dq2
     &                  + bjm(j,k,1,3)*dq3
30          continue
!$omp end do

c-----
c  Error
c-----
         else
            write(istdout,40)
            stop
40          format(' ERROR: in bccgridj, j must be 1 or jmax ')
         endif
c-----
c  j = 1: extrapolate pressure
c-----
      elseif(ibval .eq. 26) then
         if(j .eq. 1) then

!$omp do 
            do 50 k=kbeg,kend
               bjp(j,k,1,1) =-1.0
               btc(j,k,1,1) = 1.0
               s(j,k,1) = q(j+1,k,1) - q(j,k,1)
50          continue
!$omp end do

c-----
c  j = jmax: characteristic relation for pressure
c-----
         elseif(j .eq. jmax) then

!$omp do 
            do 60 k=kbeg,kend
               btc(j,k,1,1) = 1.0
               bjm(j,k,1,1) =-1.0
               s(j,k,1) = q(j-1,k,1) - q(j,k,1)
60          continue
!$omp end do

c-----
c  Error
c-----
         else
            write(istdout,115)
            stop
115         format(' ERROR: in bccgridj, j must be 1 or jmax ')
         endif
      endif

!$omp end parallel

c-----
c  End of bccgridj
c-----
      return
      end
c
c
c**********************************************************************
      subroutine bccgridk(jmax,kmax,x,y,q,rtxy,nz,s,btc,bkm,bkp,
     &                    ibval,pres,k,jbeg,jend)
c**********************************************************************
#include "common.f"
      __REAL q(jmax,kmax,3), x(jmax,kmax), y(jmax,kmax),
     &     s(jmax,kmax,3), btc(jmax,kmax,3,3), rtxy(jmax,kmax,2,3),
     &   bkm(jmax,kmax,3,3), bkp(jmax,kmax,3,3)
      logical DEBUG,DEBUG1

      DEBUG = .false.
c      DEBUG = .true.
      DEBUG1 = .false.
c      DEBUG1 = .true.

c-----
c  Freestream velocities at angle of attack
c-----
      jk = 2
      pi = 4.d0*datan( 1.d0 )
      uinf = cos( alpha*pi/180. )
      vinf = sin( alpha*pi/180. )

!$omp parallel

c-----
c  Point vortex velocity field
c-----

!$omp do private(dx,dy,ucirc,vcirc)

      do 10 j=jbeg,jend
         dx = x(j,k) - 0.25
         dy = y(j,k)
         ucirc = clift*dy/( 4.*pi*(dx*dx + dy*dy) )
         vcirc =-clift*dx/( 4.*pi*(dx*dx + dy*dy) )
c
         btc(j,k,2,2) = 1.0
         btc(j,k,3,3) = 1.0
         s(j,k,2) = uinf+ucirc - q(j,k,2)
         s(j,k,3) = vinf+vcirc - q(j,k,3)
10    continue

!$omp end do

c-----
c  k = 1: Characteristic relation for pressure
c-----
      if(ibval .eq. 25) then
c C-grid outer boundary
         if(k .eq. 1) then

!$omp do private(qq,cc,dq1,dq2,dq3)

            do 20 j=jbeg,jend
               qq = rtxy(j,k,jk,2)*q(j,k,2) + rtxy(j,k,jk,3)*q(j,k,3)
               cc = sqrt( qq*qq + beta*(rtxy(j,k,jk,2)**2
     &                                + rtxy(j,k,jk,3)**2) )
               bkp(j,k,1,1) = rtxy(j,k,jk,2)**2 + rtxy(j,k,jk,3)**2
               bkp(j,k,1,2) = rtxy(j,k,jk,2)*(qq - cc)
               bkp(j,k,1,3) = rtxy(j,k,jk,3)*(qq - cc)
c
               btc(j,k,1,1) =-bkp(j,k,1,1)
               btc(j,k,1,2) =-bkp(j,k,1,2)
               btc(j,k,1,3) =-bkp(j,k,1,3)
c
               dq1 = q(j,k+1,1) - q(j,k,1)
               dq2 = q(j,k+1,2) - q(j,k,2)
               dq3 = q(j,k+1,3) - q(j,k,3)
               s(j,k,1) = - bkp(j,k,1,1)*dq1 - bkp(j,k,1,2)*dq2
     &                    - bkp(j,k,1,3)*dq3
20          continue

!$omp end do

c-----
c  k = kmax: characteristic relation for pressure
c-----
         elseif(k .eq. kmax) then

!$omp do private(qq,cc,dq1,dq2,dq3)
            do 30 j=jbeg,jend
               qq = rtxy(j,k,jk,2)*q(j,k,2) + rtxy(j,k,jk,3)*q(j,k,3)
               cc = sqrt( qq*qq + beta*(rtxy(j,k,jk,2)**2
     &                                + rtxy(j,k,jk,3)**2) )
               bkm(j,k,1,1) =-(rtxy(j,k,jk,2)**2 + rtxy(j,k,jk,3)**2)
               bkm(j,k,1,2) =-rtxy(j,k,jk,2)*(qq + cc)
               bkm(j,k,1,3) =-rtxy(j,k,jk,3)*(qq + cc)
c
               btc(j,k,1,1) =-bkm(j,k,1,1)
               btc(j,k,1,2) =-bkm(j,k,1,2)
               btc(j,k,1,3) =-bkm(j,k,1,3)
c
               dq1 = q(j,k,1) - q(j,k-1,1)
               dq2 = q(j,k,2) - q(j,k-1,2)
               dq3 = q(j,k,3) - q(j,k-1,3)
               s(j,k,1) = bkm(j,k,1,1)*dq1 + bkm(j,k,1,2)*dq2
     &                  + bkm(j,k,1,3)*dq3
30          continue
!$omp end do

c-----
c  Error
c-----
         else
            write(istdout,40)
            print 40
            print*, '>bc.f  nz  ibval kmax   k  jbeg jend'
            print 45,nz,ibval,kmax,k,jbeg,jend
45          format(5x,6i5)
            stop 'stop: stopping in bc.f/bccgridk'
40          format(' ERROR: in bccgridk, k must be 1 or kmax ')
         endif

      elseif(ibval .eq. 26) then
c C-grid outer boundary with pressure extrapolated
c-----
c  k = 1: Characteristic relation for pressure
c-----
         if(k .eq. 1) then
!$omp do
            do 50 j=jbeg,jend
               btc(j,k,1,1) = 1.0
               bkp(j,k,1,1) =-1.0
               s(j,k,1) = q(j,k+1,1) - q(j,k,1)
50          continue
!$omp end do

c-----
c  k = kmax: characteristic relation for pressure
c-----
         elseif(k .eq. kmax) then

!$omp do
            do 60 j=jbeg,jend
               btc(j,k,1,1) = 1.0
               bkm(j,k,1,1) =-1.0
               s(j,k,1) = q(j,k-1,1) - q(j,k,1)
60          continue
!$omp end do

c-----
c  Error
c-----
         else
            write(istdout,70)
            stop
70          format(' ERROR: in bccgridk, k must be 1 or kmax ')
         endif
      endif

!$omp end parallel

c-----
c  End of bccgridk
c-----
      return
      end
c
c
c**********************************************************************
      subroutine bcoutj(jmax,kmax,q,rtxy,dj,nz,s,btc,bjm,bjp,
     &                  ibval,pres,j,kbeg,kend)
c**********************************************************************
c-------------------------------------------------------------
#include "common.f"
      __REAL q(jmax,kmax,3), rtxy(jmax,kmax,2,3), dj(jmax,kmax),
     &       s(jmax,kmax,3), btc(jmax,kmax,3,3),
     &       bjm(jmax,kmax,3,3), bjp(jmax,kmax,3,3)
      __REAL el(2)
c
      jk = 1
c-----
c  Characteristic relation for 2 variables at j=1
c  Assumes that 1/dtau = 0.
c-----

!$omp parallel

      if(j .eq. 1 .and. (ibval.eq.30.or.ibval.eq.32)) then

!$omp do private(qq,cc,a11,a12,a13,a32,a33,maxel1,maxel2,
!$omp+ dq1,dq2,dq3)

         do 100 k=kbeg,kend
            qq = rtxy(j,k,jk,2)*q(j,k,2) + rtxy(j,k,jk,3)*q(j,k,3)
            cc = sqrt( qq*qq + beta*(rtxy(j,k,jk,2)**2
     &                             + rtxy(j,k,jk,3)**2) )
            a11 = rtxy(j,k,jk,3)*q(j,k,2) - rtxy(j,k,jk,2)*q(j,k,3)
            a12 = -qq*q(j,k,3) - beta*rtxy(j,k,jk,3)
            a13 =  qq*q(j,k,2) + beta*rtxy(j,k,jk,2)
            a31 = -cc - qq
            a32 = beta*rtxy(j,k,jk,2)
            a33 = beta*rtxy(j,k,jk,3)
c-----
c  Pivot 2nd and 3rd rows to put max elements on diagonal
c-----
            el(1) = a12
            el(2) = a13
c            maxel1 = BLAS_IAMAX(2, el, 1)
            maxel1 = isamax(2, el, 1)
            bjp(j,k,maxel1+1,1) = a11
            bjp(j,k,maxel1+1,2) = a12
            bjp(j,k,maxel1+1,3) = a13
c
            el(1) = a32
            el(2) = a33
            el(maxel1) = 0.0
c            maxel2 = BLAS_IAMAX(2, el, 1)
            maxel2 = isamax(2, el, 1)
            bjp(j,k,maxel2+1,1) = a31
            bjp(j,k,maxel2+1,2) = a32
            bjp(j,k,maxel2+1,3) = a33
c
            btc(j,k,2,1) = -bjp(j,k,2,1)
            btc(j,k,2,2) = -bjp(j,k,2,2)
            btc(j,k,2,3) = -bjp(j,k,2,3)
            btc(j,k,3,1) = -bjp(j,k,3,1)
            btc(j,k,3,2) = -bjp(j,k,3,2)
            btc(j,k,3,3) = -bjp(j,k,3,3)
c
            dq1 = q(j+1,k,1) - q(j,k,1)
            dq2 = q(j+1,k,2) - q(j,k,2)
            dq3 = q(j+1,k,3) - q(j,k,3)
            s(j,k,2) = - bjp(j,k,2,1)*dq1 - bjp(j,k,2,2)*dq2
     &                 - bjp(j,k,2,3)*dq3
            s(j,k,3) = - bjp(j,k,3,1)*dq1 - bjp(j,k,3,2)*dq2
     &                 - bjp(j,k,3,3)*dq3
100      continue

!$omp end do

c-----
c  Characteristic relation for 2 variables at j=jmax
c-----
      elseif(j .eq. jmax .and. (ibval.eq.30.or.ibval.eq.32)) then

!$omp do private(qq,cc,a11,a12,a13,a31,a32,a33,maxel1,maxel2,
!$omp+ dq1,dq2,dq3)
         do 110 k=kbeg,kend
            qq = rtxy(j,k,jk,2)*q(j,k,2) + rtxy(j,k,jk,3)*q(j,k,3)
            cc = sqrt( qq*qq + beta*(rtxy(j,k,jk,2)**2
     &                             + rtxy(j,k,jk,3)**2))
            a11 = rtxy(j,k,jk,3)*q(j,k,2) - rtxy(j,k,jk,2)*q(j,k,3)
            a12 = -qq*q(j,k,3) - beta*rtxy(j,k,jk,3)
            a13 =  qq*q(j,k,2) + beta*rtxy(j,k,jk,2)
            a31 = cc - qq
            a32 = beta*rtxy(j,k,jk,2)
            a33 = beta*rtxy(j,k,jk,3)
c-----
c  Pivot 2nd and 3rd rows to put max elements on diagonal
c-----
            el(1) = a12
            el(2) = a13
c            maxel1 = BLAS_IAMAX(2, el, 1)
            maxel1 = isamax(2, el, 1)
            bjm(j,k,maxel1+1,1) = -a11
            bjm(j,k,maxel1+1,2) = -a12
            bjm(j,k,maxel1+1,3) = -a13
c
            el(1) = a32
            el(2) = a33
            el(maxel1) = 0.0
c            maxel2 = BLAS_IAMAX(2, el, 1)
            maxel2 = isamax(2, el, 1)
            bjm(j,k,maxel2+1,1) = -a31
            bjm(j,k,maxel2+1,2) = -a32
            bjm(j,k,maxel2+1,3) = -a33
c
            btc(j,k,2,1) = -bjm(j,k,2,1)
            btc(j,k,2,2) = -bjm(j,k,2,2)
            btc(j,k,2,3) = -bjm(j,k,2,3)
            btc(j,k,3,1) = -bjm(j,k,3,1)
            btc(j,k,3,2) = -bjm(j,k,3,2)
            btc(j,k,3,3) = -bjm(j,k,3,3)
c
            dq1 = q(j,k,1) - q(j-1,k,1)
            dq2 = q(j,k,2) - q(j-1,k,2)
            dq3 = q(j,k,3) - q(j-1,k,3)
            s(j,k,2) = bjm(j,k,2,1)*dq1 + bjm(j,k,2,2)*dq2
     &               + bjm(j,k,2,3)*dq3
            s(j,k,3) = bjm(j,k,3,1)*dq1 + bjm(j,k,3,2)*dq2
     &               + bjm(j,k,3,3)*dq3
110      continue
!$omp end do

c-----
c  Extrapolate velocity, j=1
c-----
      elseif((ibval.eq.31.or.ibval.eq.33) .and. j .eq. 1) then

!$omp do 
         do 120 k=kbeg,kend
            btc(j,k,2,2) = 1.0
            btc(j,k,3,3) = 1.0
            bjp(j,k,2,2) =-1.0
            bjp(j,k,3,3) =-1.0
            s(j,k,2) = q(j+1,k,2) - q(j,k,2)
            s(j,k,3) = q(j+1,k,3) - q(j,k,3)
120      continue
!$omp end do

c-----
c  Extrapolate velocity, j=jmax
c-----
      elseif((ibval.eq.31.or.ibval.eq.33) .and. j .eq. jmax) then

!$omp do 
         do 130 k=kbeg,kend
            btc(j,k,2,2) = 1.0
            btc(j,k,3,3) = 1.0
            bjm(j,k,2,2) =-1.0
            bjm(j,k,3,3) =-1.0
            s(j,k,2) = q(j-1,k,2) - q(j,k,2)
            s(j,k,3) = q(j-1,k,3) - q(j,k,3)
130      continue

!$omp end do

c-----
c  Error 
c-----
      else
         write(istdout,*) ' ERROR: in bcoutj - ibval,j,jmax,nz: = ',
     &                                        ibval,j,jmax,nz
         call exit(1)
      endif
c-----
c  Constant static pressure: pres
c-----
      if(ibval.eq.30 .or. ibval.eq.31) then

!$omp do 
         do 200 k=kbeg,kend
            btc(j,k,1,1) = 1.0
            s(j,k,1) = pres - q(j,k,1)
200      continue
!$omp end do

      endif

!$omp end parallel

c-----
c  End of bcoutj
c-----
      return
      end
c
c
c**********************************************************************
      subroutine bcoutk(jmax,kmax,q,rtxy,dj,nz,s,btc,bkm,bkp,
     &                  ibval,pres,k,jbeg,jend)
c**********************************************************************
c-------------------------------------------------------------
#include "common.f"
      __REAL q(jmax,kmax,3), rtxy(jmax,kmax,2,3), dj(jmax,kmax),
     &       s(jmax,kmax,3), btc(jmax,kmax,3,3),
     &       bkm(jmax,kmax,3,3), bkp(jmax,kmax,3,3)
      __REAL el(2)
c
      jk = 2
c-----
c  Characteristic relation for 2 variables at k=1
c  Assumes that 1/dtau = 0.
c-----

!$omp parallel

      if((ibval.eq.30.or.ibval.eq.32) .and. k .eq. 1) then

!$omp do private(qq,cc,a11,a12,a13,a31,a32,a33,maxel1,maxel2,
!$omp+ dq1,dq2,dq3)
         do 100 j=jbeg,jend
            qq = rtxy(j,k,jk,2)*q(j,k,2) + rtxy(j,k,jk,3)*q(j,k,3)
            cc = sqrt( qq*qq + beta*(rtxy(j,k,jk,2)**2
     &                             + rtxy(j,k,jk,3)**2))
            a11 = rtxy(j,k,jk,3)*q(j,k,2) - rtxy(j,k,jk,2)*q(j,k,3)
            a12 = -qq*q(j,k,3) - beta*rtxy(j,k,jk,3)
            a13 =  qq*q(j,k,2) + beta*rtxy(j,k,jk,2)
            a31 = -cc - qq
            a32 = beta*rtxy(j,k,jk,2)
            a33 = beta*rtxy(j,k,jk,3)
c-----
c  Pivot 2nd and 3rd rows to put max elements on diagonal
c-----
            el(1) = a12
            el(2) = a13
c            maxel1 = BLAS_IAMAX(2, el, 1)
            maxel1 = isamax(2, el, 1)
            bkp(j,k,maxel1+1,1) = a11
            bkp(j,k,maxel1+1,2) = a12
            bkp(j,k,maxel1+1,3) = a13
c
            el(1) = a32
            el(2) = a33
            el(maxel1) = 0.0
c            maxel2 = BLAS_IAMAX(2, el, 1)
            maxel2 = isamax(2, el, 1)
            bkp(j,k,maxel2+1,1) = a31
            bkp(j,k,maxel2+1,2) = a32
            bkp(j,k,maxel2+1,3) = a33
c
            btc(j,k,2,1) = -bkp(j,k,2,1)
            btc(j,k,2,2) = -bkp(j,k,2,2)
            btc(j,k,2,3) = -bkp(j,k,2,3)
            btc(j,k,3,1) = -bkp(j,k,3,1)
            btc(j,k,3,2) = -bkp(j,k,3,2)
            btc(j,k,3,3) = -bkp(j,k,3,3)
c
            dq1 = q(j,k+1,1) - q(j,k,1)
            dq2 = q(j,k+1,2) - q(j,k,2)
            dq3 = q(j,k+1,3) - q(j,k,3)
            s(j,k,2) = - bkp(j,k,2,1)*dq1 - bkp(j,k,2,2)*dq2
     &                 - bkp(j,k,2,3)*dq3
            s(j,k,3) = - bkp(j,k,3,1)*dq1 - bkp(j,k,3,2)*dq2
     &                 - bkp(j,k,3,3)*dq3
100      continue

!$omp end do

c-----
c  Characteristic relation for 2 variables at k=kmax
c-----
      elseif((ibval.eq.30.or.ibval.eq.32) .and. k .eq. kmax) then

!$omp do private(qq,cc,a11,a12,a13,a31,a32,a33,maxel1,maxel2,
!$omp+ dq1,dq2,dq3)

         do 110 j=jbeg,jend
            qq = rtxy(j,k,jk,2)*q(j,k,2) + rtxy(j,k,jk,3)*q(j,k,3)
            cc = sqrt( qq*qq + beta*(rtxy(j,k,jk,2)**2
     &                             + rtxy(j,k,jk,3)**2) )
            a11 = rtxy(j,k,jk,3)*q(j,k,2) - rtxy(j,k,jk,2)*q(j,k,3)
            a12 = -qq*q(j,k,3) - beta*rtxy(j,k,jk,3)
            a13 =  qq*q(j,k,2) + beta*rtxy(j,k,jk,2)
            a31 = cc - qq
            a32 = beta*rtxy(j,k,jk,2)
            a33 = beta*rtxy(j,k,jk,3)
c
c-----
c  Pivot 2nd and 3rd rows to put max elements on diagonal
c-----
            el(1) = a12
            el(2) = a13
c            maxel1 = BLAS_IAMAX(2, el, 1)
            maxel1 = isamax(2, el, 1)
            bkm(j,k,maxel1+1,1) = -a11
            bkm(j,k,maxel1+1,2) = -a12
            bkm(j,k,maxel1+1,3) = -a13
c
            el(1) = a32
            el(2) = a33
            el(maxel1) = 0.0
c            maxel2 = BLAS_IAMAX(2, el, 1)
            maxel2 = isamax(2, el, 1)
            bkm(j,k,maxel2+1,1) = -a31
            bkm(j,k,maxel2+1,2) = -a32
            bkm(j,k,maxel2+1,3) = -a33
c
            btc(j,k,2,1) = -bkm(j,k,2,1)
            btc(j,k,2,2) = -bkm(j,k,2,2)
            btc(j,k,2,3) = -bkm(j,k,2,3)
            btc(j,k,3,1) = -bkm(j,k,3,1)
            btc(j,k,3,2) = -bkm(j,k,3,2)
            btc(j,k,3,3) = -bkm(j,k,3,3)
c
            dq1 = q(j,k,1) - q(j,k-1,1)
            dq2 = q(j,k,2) - q(j,k-1,2)
            dq3 = q(j,k,3) - q(j,k-1,3)
            s(j,k,2) = bkm(j,k,2,1)*dq1 + bkm(j,k,2,2)*dq2
     &               + bkm(j,k,2,3)*dq3
            s(j,k,3) = bkm(j,k,3,1)*dq1 + bkm(j,k,3,2)*dq2
     &               + bkm(j,k,3,3)*dq3
110      continue
!$omp end do

c-----
c  Extrapolate velocity, k=1
c-----
      elseif((ibval.eq.31.or.ibval.eq.33) .and. k .eq. 1) then

!$omp do
         do 120 j=jbeg,jend
            btc(j,k,2,2) = 1.0
            btc(j,k,3,3) = 1.0
            bkp(j,k,2,2) =-1.0
            bkp(j,k,3,3) =-1.0
            s(j,k,2) = q(j,k+1,2) - q(j,k,2)
            s(j,k,3) = q(j,k+1,3) - q(j,k,3)
120      continue
!$omp end do

c-----
c  Extrapolate velocity, k=kmax
c-----
      elseif((ibval.eq.31.or.ibval.eq.33) .and. k .eq. kmax) then

!$omp do 
         do 130 j=jbeg,jend
            btc(j,k,2,2) = 1.0
            btc(j,k,3,3) = 1.0
            bkm(j,k,2,2) =-1.0
            bkm(j,k,3,3) =-1.0
            s(j,k,2) = q(j,k-1,2) - q(j,k,2)
            s(j,k,3) = q(j,k-1,3) - q(j,k,3)
130      continue
!$omp end do

c-----
c  Error 
c-----
      else
         write(istdout,*)' ERROR: in bcoutk, ibval,k = ',ibval,k
         call exit(1)
      endif
c-----
c  Constant static pressure: pres
c-----
      if(ibval .eq. 30 .or. ibval .eq. 31) then

!$omp do 
         do 140 j=jbeg,jend
            btc(j,k,1,1) = 1.0
            s(j,k,1) = pres - q(j,k,1)
140      continue
!$omp end do

      endif

!$omp end parallel

c-----
c  End of bcoutk
c-----
      return
      end
c
c============================================
