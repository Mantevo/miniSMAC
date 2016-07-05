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
c  1. rbcfile
c  2. rbcmain
c  3. rbcpzn
c  4. rbccgrid
c  5. loadint
c  6. testbcmain
c  7. testbcflag
c
c************************************************************
c
c
c************************************************************************
      subroutine rbcfile(jmax,kmax,x,y,
     & idimj,idimk,
     & jindex_global,kindex_global,ioverlap
     & ) 
c************************************************************************
c  This routine calls all the routines which read in data files containing
c  info used to update boundary points.
c
c-------------------------------------------------------------------------
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"

      __REAL x(jmax,kmax),y(jmax,kmax)
      __REAL bcflag(jmax,kmax)

      __INTEGER stat(MPI_STATUS_SIZE)

      dimension jindex_global(numprocs),kindex_global(numprocs)

      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

      if (nodeid.eq.0) then
      write(istdout,10)
10    format(//'----------------- Node 0 -----------------------------',
     &        /'| bc bc bc bc bc bc bc bc bc bc bc bc bc bc bc bc bc |',
     &        /'------------------------------------------------------',
     &        /'               Boundary Conditions: '/)
      endif

c-----
c  read file bcmain.dat:
c-----

      call rbcmain(jmax,kmax,x,y,idimj,idimk,
     & jindex_global,kindex_global,ioverlap
     & )

      if (DEBUG) then
       print*, '** Node ',nodeid,': After call rbcmain, jmax,kmax = ',
     1 jmax,kmax
       call flush(6)
       print*
c       print*, ' stop in rbc.f/rbcfile #1 -- examine bcmain.dat'
c       stop 'stop in rbc.f/rbcfile #1 -- examine bcmain.dat'
      endif

c-----
c  read file bcpzn.dat; not needed for a single processor run
c-----
      if (numprocs.gt.1) then
       call rbcpzn(jmax,kmax)
      endif
      if (DEBUG) then
       print*, '** Node ',nodeid,': After call rbcpzn, jmax,kmax = ',
     1 jmax,kmax
       call flush(6)
      endif

c-----
c  take care of c-grid wake-cut boundaries
c-----
c      call rbccgrid(jmax,kmax)

      goto 500

c check boundary conditions for multiple values

c +++++++++++++++++++++++
      if (numprocs.gt.1) then
c +++++++++++++++++++++++
      if (nodeid.eq.0) then
c +++++++++++++++++++++++
c-----
c  ... Initialize bcflag array
c-----
      do 34 k=1,kmax
      do 34 j=1,jmax
         bcflag(j,k) = 0.0
34    continue

c-----
c  ... Test bcmain.dat boundary conditions
c-----
      if (DEBUG) then
       print*
       print*, ' Node ',nodeid,': testing ',nbcreg,' bcmain.dat bcs'
       print*
      endif
      do 36 nr=1,nbcreg
       if(ibcval(nr) .ge. 0 .and. ibcval(nr) .ne. 60 .and.
     &   ibcval(nr) .ne. 61) then
         call testbcmain(bcflag,jmax,kmax)
       endif
36    continue

c-----
c  ... Test for multiple boundary conditions applied to the same point
c-----
      if (DEBUG) then
       print*
       print*, ' Node ', nodeid,': testing for multiple bcs'
      endif
      call testbcflag(bcflag,jmax,kmax)
c-----

c +++++++++++++++++++++++
      endif
c +++++++++++++++++++++++
      endif
c +++++++++++++++++++++++

500   continue

      if (nodeid.eq.0) then
      write(istdout,100)
100   format(/'------------------------------------------------------',
     &      /'| bc bc bc bc bc bc bc bc bc bc bc bc bc bc bc bc bc |',
     &      /'----------------- End Node 0  -------------------------')
      endif

c-----
c  End of rbcfile
c-----
      return
      end
c
c
c************************************************************************
      subroutine rbcmain(jmax,kmax,x,y,idimj,idimk,
     & jindex_global,kindex_global,ioverlap
     & )
c************************************************************************
c  This routine reads in the following mandatory boundary condition file:
c
c Input:
c     bcmain.dat_original   bc file for original grid; one zone only
c     idimj,idimk           index of blocks (not zones)
c Output
c     bcmain.dat:           Contains boundary conditions such as walls, 
c                             inflow/outflow, symmetry, etc. for all subgrids
c     bcinterp.dat:         Contains indices for interpolating data across
c                             C-grid wake cuts
c     bcwake.dat            File containing indices for generating wake
c                             boundary conditions
c     jindex_global:        global gridpoint indices of blocks
c     kindex_global
c 
c------------------------------------------------------------------------
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      __INTEGER jmax, kmax
      character*110 line, filename, filename_new
c      character*110 file_home_dir, file_name
      character*40 bc_string
      character*1 tabchar, space
      logical DEBUG,DEBUG1,DEBUG_150,DEBUG2
c#include "mpif.h"
c#include "mpi_params.f"
      __INTEGER stat(MPI_STATUS_SIZE)
      __INTEGER, allocatable:: ivarbc(:,:),ivarwall(:,:),
     &                         ivarbc_wake(:,:)
      __INTEGER, allocatable:: ibcval_o(:),nzbc_o(:),jbcb_o(:),
     &                         jbce_o(:),kbcb_o(:),kbce_o(:)
      __INTEGER, allocatable:: ibcval_wake_1(:),nzbc_wake_1(:),
     &     jbcb_wake_1(:),jbce_wake_1(:),kbcb_wake_1(:),kbce_wake_1(:),
     &     jmax_wake_1(:),kmax_wake_1(:)
      __INTEGER, allocatable:: ibcval_wake_2(:),nzbc_wake_2(:),
     &     jbcb_wake_2(:),jbce_wake_2(:),kbcb_wake_2(:),kbce_wake_2(:),
     &     jmax_wake_2(:),kmax_wake_2(:)

c      __INTEGER jkbc(:)
c      INTEGER, allocatable:: jmaxx(:),kmaxx(:)
      __INTEGER jmaxx(numprocs),kmaxx(numprocs)
      integer ivarjk(2)
      dimension jindex_global(numprocs),kindex_global(numprocs)
      __REAL x(jmax,kmax),y(jmax,kmax)
      parameter (isurface_bc_max=20)

      DEBUG = .false.
c      DEBUG = .true.
      DEBUG1 = .false.
c      DEBUG1 = .true.
      DEBUG2 = .false.
c      DEBUG2 = .true.
      DEBUG_150 = .false.
c      DEBUG_150 = .true.

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c
c++++++++++++++++++++++++
      if (nodeid.eq.0) then
c++++++++++++++++++++++++
      open(99,status='scratch')
c      open(99,status='unknown')
2     format(a80)
      tabchar = '	'
      space = ' '

c-----
c  Read in main patched boundary condition info for original grid
c-----
c      file_home_dir =
c     1 '/gscratch1/dwbarne/ms2d_2d_only_partitioned_grids/src/'
c      file_name = 'bcmain.dat_original'
c      filename = trim(file_home_dir)//trim(file_name) 
      filename = trim(bcmain_directory)//'/'//trim(bcmain_filename) 
      print*
      print*,' filename = ',filename
      print*, '  string length = ',len_trim(filename),'/110 characters'
      print*
c
      open(49,file=filename,form='formatted',status='old',err=10)
      goto 20

10    continue
      write(istdout,*)' ERROR: Could not find boundary condition '
      write(istdout,*)'        input file = ',filename
      write(istdout,*)' Program terminated.'
      call MPI_ABORT(MPI_COMM_WORLD,1)
      stop ' -- STOP: could not find boundary condition input file'

20    continue

c find out how many original bc's there are, to dimension variables
      nbcreg_o = 0 
300   continue
      read(49,2,end=400) line
      line = trim(line)
      if( line(1:1) .eq. 'c' .or. line(1:1) .eq. 'C'.or.line.eq.'') 
c skip comment lines in file
     &  goto 300
      nbcreg_o = nbcreg_o + 1
      goto 300
400   continue

      if (nbcreg_o.eq.0) then
       print*
c       print*, 'ERROR: no bcs found in file bcmain.dat_original'
       print*, 'ERROR: no bcs found in file ',trim(bcmain_filename)
       print*, '  This file must contain boundary conditions pertaining'
       print*, '  to the original grid.'
       print*, '  Program must terminate.'
       print*
       stop 'stop: no bcs in main bc file'
      else
       if (DEBUG) then
        print*
        write(6,405) trim(filename), nbcreg_o
405     format(' Number of bcs in file ',a,' = ',i5)
       endif
      endif

c allocate memory for the following temporary values
      allocate (ibcval_o(nbcreg_o),nzbc_o(nbcreg_o),jbcb_o(nbcreg_o),
     &          jbce_o(nbcreg_o),kbcb_o(nbcreg_o),kbce_o(nbcreg_o))
c
c rewind and read in original data
      rewind (49)
      do 40 nr = 1,nbcreg_o
30    continue
      read(49,2) line
c skip processing if line is a comment
      if( line(1:1) .eq. 'c' .or. line(1:1) .eq. 'C') goto 30
        do 35 icol=1,80
c convert any tabs into spaces
          if( line(icol:icol) .eq. tabchar) line(icol:icol) = space
35      continue
      rewind(99)
      print*, ' Before write(99,2), line = ',line
      write(99,2) line
      rewind(99)
      read(99,*,end=30,err=31) ibcval_o(nr),nzbc_o(nr),jbcb_o(nr),
     &   jbce_o(nr),kbcb_o(nr),kbce_o(nr)
      if (DEBUG) then
       print*
       print*, 'After read(99:'
       print*, ' nr = ',nr
       print*, ' nbcreg_o = ',nbcreg_o
       print*, ' line = ',line
       print*, ' ibcval_o(nr) = ',ibcval_o(nr)
       print*, ' nzbc_o(nr) = ',nzbc_o(nr)
       print*, ' jbcb_o(nr) = ',jbcb_o(nr)
       print*, ' jbce_o(nr) = ',jbce_o(nr)
       print*, ' kbcb_o(nr) = ',kbcb_o(nr)
       print*, ' kbce_o(nr) = ',kbce_o(nr)
       print*
      endif
      goto 40
31    continue
c read 99 err
      print*
      print*, ' Error reading scratch file unit 99 in rbc.f/rbcmain'
      print*, ' nr = ',nr
      print*, ' nbcreg_o = ',nbcreg_o
      print*, ' line = ',line
      print*, ' ibcval_o(nr) = ',ibcval_o(nr)
      print*, ' nzbc_o(nr) = ',nzbc_o(nr)
      print*, ' jbcb_o(nr) = ',jbcb_o(nr)
      print*, ' jbce_o(nr) = ',jbce_o(nr)
      print*, ' kbcb_o(nr) = ',kbcb_o(nr)
      print*, ' kbce_o(nr) = ',kbce_o(nr)
      print*
      print*, '  terminating program'
      print*
      stop 'stop: in rbc.f/rbcmain'
c
40    continue


      print*
      write(istdout,46) filename,nbcreg_o
46    format(' Number of BC regions read in from ',a20,' = ',i3)
      print*,'#   ibcval_o  nzbc_o  jbcb_o  jbce_o  kbcb_o  kbce_o'
      print*,'-   --------  ------  ------  ------  ------  ------'
      do 74 nr=1,nbcreg_o
       write(*,77) nr,ibcval_o(nr),nzbc_o(nr),jbcb_o(nr),jbce_o(nr),
     &             kbcb_o(nr),kbce_o(nr)
77     format(i2,'. ',i7,i8,1x,4i8)
74    continue
      print*

c close main bc file
      close(49)

c take care of negative indices on original grid
c ... negative index '-n' indicates, for example, jmax + 1 - n
c ... made it easier for users not to have to keep track of jmax
      do 81 nr=1,nbcreg_o
c if -1
       if (jbcb_o(nr).eq.-1) jbcb_o(nr) = jmax_o 
       if (jbce_o(nr).eq.-1) jbce_o(nr) = jmax_o 
       if (kbcb_o(nr).eq.-1) kbcb_o(nr) = kmax_o 
       if (kbce_o(nr).eq.-1) kbce_o(nr) = kmax_o
c if -2
       if (jbcb_o(nr).eq.-2) jbcb_o(nr) = jmax_o - 1 
       if (jbce_o(nr).eq.-2) jbce_o(nr) = jmax_o - 1
       if (kbcb_o(nr).eq.-2) kbcb_o(nr) = kmax_o - 1
       if (kbce_o(nr).eq.-2) kbce_o(nr) = kmax_o - 1
81    continue

c write original boundary conditions with actual indices
      print*
      write(istdout,47) filename,nbcreg_o
47    format(' Re-indexed (no neg indices): Number of BC regions from ',
     &a20,' = ',i3)
      print*,'#   ibcval_o  nzbc_o  jbcb_o  jbce_o  kbcb_o  kbce_o'
      print*,'-   --------  ------  ------  ------  ------  ------'
      do 75 nr=1,nbcreg_o
       write(*,77) nr,ibcval_o(nr),nzbc_o(nr),jbcb_o(nr),jbce_o(nr),
     &             kbcb_o(nr),kbce_o(nr)
75    continue
      print*

c+++++++++++++++++++++++++
      endif
c++++++++++++++++++++++++


c get jmax,kmax from all other nodes

c++++++++++++++++++++++++
      if (numprocs.gt.1) then
c++++++++++++++++++++++++
      if (nodeid.eq.0) then
c+++++++++++++++++++++++g

c get jmax,kmax from each node
c ... first for node 0
      jmaxx(1) = jmax
      kmaxx(1) = kmax
c ... now the other nodes
      do i=2,numprocs
        node_in = i-1
        call MPI_RECV(ivarjk,2,MPI_INTEGER,node_in,node_in,
     &   MPI_COMM_WORLD,stat,ierr)
c unpack ivarjk
        jmaxx(i) = ivarjk(1)
        kmaxx(i) = ivarjk(2)
      enddo

c++++++++++++++++++++++++
      else
c++++++++++++++++++++++++

c send jmax,kmax from all other nodes to node 0
      if (DEBUG) then
        write(*,73) nodeid,jmax,kmax
73      format(' rbcmain Node 0: MPI_SEND to Node: ',i5,
     1   ' jmax = ',i6,'  kmax = ',i6)
      endif
      ivarjk(1) = jmax
      ivarjk(2) = kmax
      call MPI_SEND(ivarjk,2,MPI_INTEGER,0,nodeid,
     & MPI_COMM_WORLD,ierr)

c++++++++++++++++++++++++
      endif 
c++++++++++++++++++++++++
c for numprocs = 1
      else
c++++++++++++++++++++++++
       jmaxx(1) = jmax
       kmaxx(1) = kmax
c++++++++++++++++++++++++
      endif
c++++++++++++++++++++++++

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c construct bcmain.dat for use with subgrids, constructed from
c  the original main bc file 

c++++++++++++++++++++++++
      if (DEBUG.and.nodeid.eq.0) then
c++++++++++++++++++++++++
      print*
      print*, '  i  jmaxx(i)  kmaxx(i)'
      do 250 i=1,numprocs
       write(*,252) i,jmaxx(i),kmaxx(i)
252    format(i3,i6,i10)
250   continue
c++++++++++++++++++++++++
      endif
c++++++++++++++++++++++++

c++++++++++++++++++++++++
      if (DEBUG.and.nodeid.eq.0) then
c++++++++++++++++++++++++
      print*
      print*,' In rbcmain: idimj,idimk = ',idimj,idimk
      print*
      print*,' Global indices without overlap'
      print*,'  j   k   index   jindex_g  kindex_g'
      print*,'--- ---   -----   --------  --------'
      do 304 j=1,idimj
      do 304 k=1,idimk
      index = (j-1)*idimk + k
      write(*,303) j,k,index,jindex_global(index),kindex_global(index)
303   format(i3,i4,i6,2i10)
304   continue
c++++++++++++++++++++++++
      endif
c++++++++++++++++++++++++

c determine indices of subgrids over which original indices span
c    jmax_o, kmax_o: dimensions of original grid
c original variables are
c   from nr = 1 to nbcreg_o
c  nr,ibcval_o(nr),nzbc_o(nr),jbcb_o(nr),jbce_o(nr),kbcb_o(nr),kbce_o(nr)
c nbcreg = number of subgrid boundary condition regions
c j1:  min j index value for subgrid j,k 
c j2:  max j index
c k1:  min k index value for subgrid j,k
c k2:  max k index

c j,k = global index array
c++++++++++++++++++++++++
      if (nodeid.eq.0) then
c++++++++++++++++++++++++
c open new bc file
      filename_new = 'bcmain.dat'
      open(46,file=filename_new,form='formatted',status='unknown',
     & err=22)
      goto 24
22    continue
      write(istdout,*)' ERROR: Could not open boundary condition'
      write(istdout,*)'        input file = ',filename_new
      write(istdout,*)' Program must terminate.'
      call MPI_FINALIZE(ierr)
      call exit(0)
24    continue

c write header
      write(46,5) filename_new

      if (DEBUG) then
       print*
       print*, ' jmax_o, kmax_o = ',jmax_o,kmax_o
       print*, '  j   k    index   j1   j2   k1   k2'
       print*, ' --- ---  ------  ---- ---- ---- ----'
      endif

c set counters and flags 
      nbcreg = 0    ! total boundary condition count for all subgrids
      iflag_wake_cut_surface = 0
      isurface_1_bc_count = 0    ! top surface, or surface 1
      isurface_2_bc_count = 0    ! bottom surface, or surface 2

      do 306 nr=1,nbcreg_o
      ibc_count = 0   ! refreshes for each new original bc
      ibcval_temp = ibcval_o(nr)
      jfound = 0
      kfound = 0

c loop over global blocks
      do 305 j=1,idimj   
      do 305 k=1,idimk   

c If more than 1 processor, calculate where boundary conditions begin
c  and end on the subgrids, which gets pretty messy
c==========
      if (numprocs.gt.1) then
c==========

c use index for blocks to compute subgrid indices
       index = (j-1)*idimk + k
       index_jm1 = (j-2)*idimk + k
       index_km1 = (j-1)*idimk + (k-1)
c for j index
       if (jmax_o.eq.1) then
        j2 = jindex_global(index)
        j1 = 1
       else if (j.eq.1) then
        j2 = jindex_global(index) + ioverlap
        j1 = 1
       else if (j.eq.idimj) then
        j2 = jindex_global(index)
        j1 = jindex_global(index_jm1)
       else
        j2 = jindex_global(index) + ioverlap
        j1 = jindex_global(index_jm1)
       endif
c for k index
       if (kmax_o.eq.1) then
        k2 = kindex_global(index)
        k1 = 1
       else if (k.eq.1) then
        k2 = kindex_global(index) + ioverlap
        k1 = 1
       else if (k.eq.idimk) then
        k2 = kindex_global(index)
        k1 = kindex_global(index_km1)
       else
        k2 = kindex_global(index) + ioverlap
        k1 = kindex_global(index_km1)
       endif
c print some output, but only for the first bc (nr=1)
       if (DEBUG.and.nr.eq.1) then
        write(*,257) j,k,index,j1,j2,k1,k2
257     format(2i4,i6,i7,3i5)
       endif

c find if jbcb_o, jbce_o, kbcb_o, or kbce_o is bounded by j1,j2 and k1,k2
c ... if jbcb_o(nr) and jbce_o(nr) lie outside current subgrid dimensions
       jvalid = 0
       jsuccess = 0
       if (jbcb_o(nr).lt.j1.and.jbce_o(nr).gt.j2) then
        jb = 1 
        je = j2 - j1 + 1
        jvalid = 1
        jfound = 1 
        jsuccess = 1

c ... if jbcb_o(nr) lies within and jbce_o(nr) is outside
       elseif (jbcb_o(nr).ge.j1.and.jbcb_o(nr).lt.j2
     & .and.jbce_o(nr).gt.j2) then
        jb = jbcb_o(nr) - j1 + 1
        je = j2 - j1 + 1
        jvalid = 1
        jfound = 1
        jsuccess = 2

c ... if jbcb_o(nr) lies outside and jbce_o(nr) is inside
       elseif (jbcb_o(nr).lt.j1.and.jbce_o(nr).le.j2.
     & and.jbce_o(nr).gt.j1) then
        jb = 1
        je = jbce_o(nr) - j1 + 1 
        jvalid = 1
        jfound = 1
        jsuccess = 3

c ... if jbcb_o(nr) and jbce_o(nr) lies within current subgrid
       elseif (jbcb_o(nr).ge.j1.and.jbce_o(nr).le.j2) then
        jb = jbcb_o(nr) - j1 + 1
        je = jbce_o(nr) - j1 + 1
        jvalid = 1
        jfound = 1
        jsuccess = 4
c
       else
        goto 305
c
       endif

c
c if j is found above, next find k

      if (jvalid.eq.1) then
       kvalid = 0
       ksuccess = 0
c ... if kbcb_o(nr) and kbce_o(nr) lie outside current subgrid dimensions
       if (kbcb_o(nr).lt.k1.and.kbce_o(nr).gt.k2) then
        kb = 1
        ke = k2 - k1 + 1
        kvalid = 1
        kfound = 1 
        ksuccess = 1

c ... if kbcb_o(nr) lies within and kbce_o(nr) is outside
       elseif (kbcb_o(nr).ge.k1.and.kbcb_o(nr).lt.k2
     & .and.kbce_o(nr).gt.k2) then
        kb = kbcb_o(nr) - k1 + 1
        ke = k2 - k1 + 1 
        kvalid = 1
        kfound = 1 
        ksuccess = 2

c ... if (kbcb_o(nr) lies outside and kbce_o(nr) is inside
       elseif (kbcb_o(nr).lt.k1.and.kbce_o(nr).le.k2
     & .and.kbce_o(nr).gt.k1) then
        kb = 1
        ke = kbce_o(nr) - k1 + 1
        kvalid = 1
        kfound = 1 
        ksuccess = 3

c ... if kbcb_o(nr) and kbce_o(nr) lies within current subgrid
       elseif (kbcb_o(nr).ge.k1.and.kbce_o(nr).le.k2) then
        kb = kbcb_o(nr) - k1 + 1
        ke = kbce_o(nr) - k1 + 1
        kvalid = 1
        kfound = 1 
        ksuccess = 4
c ... 
      else
        goto 305
      endif
      endif

      if (DEBUG) then
       print*
       write(istdout,470) nr,jbcb_o(nr),jbce_o(nr),kbcb_o(nr),kbce_o(nr)
470    format(' === nr,jbcb_o,jbce_o,kbcb_o,kbce_o: ',5i5)
       write(istdout,471) j1,j2,k1,k2
471    format('   j1,j2,k1,k2 = ',4i5)
       write(istdout,472) jsuccess,ksuccess
472    format(' jsuccess, ksuccess = ',2i5)
       write(istdout,*) ' '
      endif
c
c==========
      else
c=========
c for just 1 processor...
       jb = jbcb_o(nr)
       je = jbce_o(nr)
       kb = kbcb_o(nr)
       ke = kbce_o(nr)
       j1 = jb
       j2 = je
       k1 = kb
       k2 = ke
       index = 1
c
       if (DEBUG) then
        print*
        write(istdout,473) 1,1,index,j1,j2,k1,k2
473     format(1x,i3,1x,i3,2x,i6,2x,4i5)
       endif
c==========
      endif
c==========

      ibc_count = ibc_count + 1
      nbcreg = nbcreg + 1
c
      if (ibc_count.eq.1) then
       write(46,490) 
490    format('c')
       if(ibcval_temp.eq.0) then
         bc_string = 'No-slip wall'
       elseif (ibcval_temp.eq.25) then
         bc_string = 'C-grid outer boundary w/ point vortex'
       elseif (ibcval_temp.eq.31) then
         bc_string = 'Outflow boundary'
       elseif (ibcval_temp.eq.60) then
         bc_string = 'C-grid wake cut surface'
c set flag for first (top) surface (1) or second (bottom) surface (2)
         iflag_wake_cut_surface = iflag_wake_cut_surface + 1 
       else
         bc_string = 'ERROR'
         print*
         print*, ' ERROR: invalid value for bc_string in preproc.f'
         print*, '   Examine main bc file ',trim(bcmain_filename)
         print*, '   Currently, valid values for ibcval are'
         print*, '     0 - no-slip wall'
         print*, '    25 - c-grid outer boundary w/ point vortex'
         print*, '    31 - outflow boundary'
         print*, '    60 - c-grid wake cut'
         print*
         print*, ' Program is terminating'
         print*
         call exit(1)
       endif
       write(46,495) bc_string
495    format('c ',a)
 
      endif

      if (bc_string .eq. 'C-grid wake cut surface') then
       print*, ' Wake index,jb,je,jmaxx(index) =',index,jb,je,
     1  jmaxx(index)
c stay away from j=1 or j=jmax; these are updated elsewhere via point-to-point
c  match up with neighboring grids
       if (jb.eq.1) then
        jb = 2
       endif
       if (je.eq.jmaxx(index)) then
        je = je - 1
       endif
      endif

c write new boundary conditions to file 'bcmain.dat'
      write(46,7) ibcval_temp,index,jb,je,kb,ke
7     format(6x,i2,4x,i5,1x,i5,1x,i5,1x,i5,1x,i5)

c if wake cut boundary, save zones to which they apply 
      if (bc_string.eq.'C-grid wake cut surface') then
        if (iflag_wake_cut_surface.eq.1) then
         isurface_1_bc_count = isurface_1_bc_count + 1
         if (isurface_1_bc_count.gt.isurface_bc_max) then
          print*
          print*, ' ERROR: isurface_1_bc_count > isurface_bc_max in rbc.
     &f/rbcmain'
          print*, ' isurface_1_bc_count =',isurface_1_bc_count
          print*, ' isurface_bc_max =',isurface_bc_max
          print*, '  Increase the value for parameter isurface_bc_max in
     & rbc.f/rbcmain'
          print*
          print*, ' Program terminating.'
          print*
          stop 'stop: isurface_1_bc_count > isurface_bc_max in rbc.f/rbc
     &main'
         endif
        else if (iflag_wake_cut_surface.eq.2) then
         isurface_2_bc_count = isurface_2_bc_count + 1
         if (isurface_2_bc_count.gt.isurface_bc_max) then
          print*
          print*, ' ERROR: isurface_2_bc_count > isurface_bc_max in rbc.
     &f/rbcmain'
          print*, ' isurface_2_bc_count =',isurface_2_bc_count
          print*, ' isurface_bc_max =',isurface_bc_max
          print*, '  Increase the value for parameter isurface_bc_max in
     & rbc.f/rbcmain'
          print*
          print*, ' Program terminating.'
          print*
          stop 'stop: isurface_2_bc_count > isurface_bc_max in rbc.f/rbc
     &main'
         endif
        endif
      endif

305   continue

c let user know if nothing found for this boundary condition

c+++++
      if (numprocs.gt.1) then
c+++++
      if (jfound.eq.0.or.kfound.eq.0) then
       write(istdout,*) ' '
       write(istdout,*) 'ERROR: boundary condition not found in grid dom
     &ain!'
       write(istdout,450) nr,ibcval(nr),jbcb_o(nr),jbce_o(nr),
     & kbcb_o(nr),kbce_o(nr)
450    format(i4,'. ',5i5)
       if (jfound.eq.0) then
        write(istdout,*) jbcb_o(nr),jbce_o(nr)
455     format(' Nothing found for j between ',i5,' and ',i5)
       endif
       if (kfound.eq.0) then
        write(istdout,460) kbcb_o(nr),kbce_o(nr)
460     format(' Nothing found for k between ',i5,' and ',i5)
       endif
      endif
c+++++
      endif
c+++++

306   continue

c---------------------------
      endif
c---------------------------

c+++++++++++++++++++++++++++
      if (nodeid.eq.0) then
c+++++++++++++++++++++++++++
c deallocate unneeded variables
      deallocate (ibcval_o,nzbc_o,jbcb_o,jbce_o,kbcb_o,kbce_o)

      write(istdout,465) nbcreg
465   format(/'Total number of subgrid boundary conditions: ',i5)
      write(istdout,*) ' '
c+++++++++++++++++++++++++++
      endif
c+++++++++++++++++++++++++++

c      if (DEBUG) then
c       stop 'stop: after 306 continue in rbc.f/rbcmain'
c      endif

c now that file bcmain.dat contains boundary conditions for the
c  subgrids, read these boundary conditions into appropriate variables

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if (DEBUG2) then
       print*
       call flush(6)
       call flush(istdout)
       print*, '======= node responding:',nodeid
c       stop 'stop: after node responding'
      endif

c++++++++++++++++++++++++
      if (nodeid.eq.0) then
c++++++++++++++++++++++++

      rewind(46)

      linecount = 0
      lines_skipped = 0
      do 4000 nr=1,nbcreg
3000  continue

c read line from new bc file bcmain.dat
      read(46,2,end=4000) line 
      linecount = linecount + 1

c skip over comment lines 
      if (line(1:1).eq.'c'.or.line(1:1).eq.'C') then 
       if (DEBUG) then
        write(*,3200) filename_new,linecount
3200  format('>> In ',a20,', skipping line ',i4)
c        print*, ' >> In bcmain.dat, skipping line ',linecount
       endif
       lines_skipped = lines_skipped + 1
       goto 3000
      endif

c check for tabs; if tabs exist, replace with spaces
      do 3500 icol=1,80
       if (line(icol:icol).eq.tabchar) line(icol:icol) = space
3500  continue

c write line to scratch file 99
      rewind(99)
      write(99,2) line 
      rewind(99)
c assign values to variables
      read(99,*,end=3000) ibcval(nr),nzbc(nr),
     &                    jbcb(nr),jbce(nr),kbcb(nr),kbce(nr) 

4000  continue

c close scratch file
c      close(99)

c close bcmain.dat file
      close(46)

c-----
c  Process boundaries 
c-----
      do 50 nr=1,nbcreg
c ... check if beginning and end j values are identical
         if (jbcb(nr) .eq. jbce(nr)) then
            jkbc(nr) = 1
c ... check if beginning and end k values are identical
         elseif (kbcb(nr) .eq. kbce(nr)) then
            jkbc(nr) = 2
c ... if none of above is true, then we have a problem with grid or bc's
         else
            write(istdout,55) nr,jbcb(nr),jbce(nr),kbcb(nr),kbce(nr)
            stop 'stop: error in surface number in bcmain.dat'
         endif
50    continue
55    format(/' ERROR: surface number ',i5,' in bcmain.dat',
     &        ' does not describe a surface.',/,
     &        '     jbeg,jend: ',2i5,/,
     &        '     kbeg,kend: ',2i5)

c ... apply some human-readable description to the bc's 
c      and print to output file
      if (DEBUG1) then
       bc_string = ''
       print*
       print*, '        bc_string       ibcval nzbc jkbc jbcb jbce kbcb 
     & kbce'
       do 60 nr=1,nbcreg
         if(ibcval(nr) .eq. 0) then
            bc_string = 'No-slip wall'
         elseif(ibcval(nr) .eq. 25) then
            bc_string = 'C-grid outer boundary w/ point vortex'
         elseif(ibcval(nr) .eq. 31) then
            bc_string = 'Outflow boundary'
         elseif(ibcval(nr) .eq. 60) then
             bc_string = 'C-grid wake cut surface'
         else
            bc_string = 'ERROR'
         endif
         write(istdout,70) bc_string,ibcval(nr),nzbc(nr),jkbc(nr),
     &               jbcb(nr),jbce(nr),kbcb(nr),kbce(nr)
60     continue
70     format(5x,a23,i3,i4,3i5,2i5)
      endif

c-----
c  Process walls
c-----
      nwall = 0
      do 80 nr=1,nbcreg
         if(ibcval(nr) .eq. 0 .or. ibcval(nr) .eq. 1
     & .or. ibcval(nr) .eq. 4 .or. ibcval(nr) .eq. 5) then
            nwall = nwall + 1
            jkwall(nwall) = jkbc(nr)
            nzwall(nwall) = nzbc(nr)
            jwall1(nwall) = jbcb(nr)
            jwall2(nwall) = jbce(nr)
            kwall1(nwall) = kbcb(nr)
            kwall2(nwall) = kbce(nr)
            if (ibcval(nr) .eq. 0 .or. ibcval(nr) .eq. 4) then
              jkinc(nwall) = 1
            else
              jkinc(nwall) = -1
            endif
         endif
80    continue

c check nwall with value for nwallmax in common.f
      if (nwall.gt.nwallmax) then
       print*
       print*, ' ERROR: nwall > nwallmax'
       print*, '  nwall =',nwall,', nwallmax = ',nwallmax
       print*, ' nwall must be less than or equal to nwallmax'
       print*
       print*, ' Increase the value for nwallmax in common.f and rerun.'
       print*
       print*, ' Program is terminating.'
       call exit(0)
      endif

      write(istdout,90) nwall
90    format(/' Total of ',i3,' no-slip walls processed.')

      print*
      print*, ' No-slip walls for subgrids:'
      print*, '    #  nzbc  jkbc  jbcb  jbce  kbcb  kbce'
      do nw = 1,nwall
      write(istdout,71) nw,nzwall(nw),jkwall(nw),jwall1(nw),jwall2(nw),
     1  kwall1(nw),kwall2(nw) 
71    format(3x,i3,6i6)
      enddo
      print*

c separate the c-grid wake cut bc's from the main bcs since wake cut bcs
c  are handled separately via interpolation than the others 
c ... find out how many 
      num_wake_bcs_total = isurface_1_bc_count + isurface_2_bc_count 

c skip if no C-grid wake bcs
      if (num_wake_bcs_total.eq.0) then 
       print*
       print*, 'No wake boundary conditions exist.'
       goto 600
      else
       print*
       print*, 'Number of wake boundary conditions:',num_wake_bcs_total
      endif

c allocate memory
c ... surface 1
      is1 = isurface_1_bc_count
      allocate(
     &  ibcval_wake_1(is1),nzbc_wake_1(is1),
     &  jbcb_wake_1(is1),jbce_wake_1(is1),
     &  kbcb_wake_1(is1),kbce_wake_1(is1),
     &  jmax_wake_1(is1),kmax_wake_1(is1)
     & )
c ... surface 2
      is2 = isurface_2_bc_count
      allocate(
     &  ibcval_wake_2(is2),nzbc_wake_2(is2),
     &  jbcb_wake_2(is2),jbce_wake_2(is2),
     &  kbcb_wake_2(is2),kbce_wake_2(is2),
     &  jmax_wake_2(is2),kmax_wake_2(is2)
     & )

c assign wake values to new variables for convenience, and
c  separate surface 1 bcs from surface 2 bcs
      icount = 0
      do 1460 nw = 1,nbcreg
       if (ibcval(nw).eq.60) then
        icount = icount + 1
        if (icount.le.isurface_1_bc_count) then
c ... 1st surface
         if (DEBUG.and.icount.eq.1) then
           print*
           print*, ' Wake bcs separated into surface 1 and surface 2:'
           print*, ' surf icnt  # ibcval nzbc jbcb jbce kbcb kbce jmx  k
     &max'
          endif
         icount1 = icount
         ibcval_wake_1(icount1) = ibcval(nw)
         nzbc_wake_1(icount1) = nzbc(nw)
         jbcb_wake_1(icount1) = jbcb(nw)
         jbce_wake_1(icount1) = jbce(nw)
         kbcb_wake_1(icount1) = kbcb(nw)
         kbce_wake_1(icount1) = kbce(nw)
         jmax_wake_1(icount1) = jmaxx(nzbc(nw)) - 1
         kmax_wake_1(icount1) = kmaxx(nzbc(nw)) - 1
         if (DEBUG) then
          write(*,2020) 1,icount,icount1,ibcval_wake_1(icount1),
     &     nzbc_wake_1(icount1),
     &     jbcb_wake_1(icount1),jbce_wake_1(icount1),
     &     kbcb_wake_1(icount1),kbce_wake_1(icount1),
     &     jmax_wake_1(icount1),kmax_wake_1(icount1)
2020      format(11i5)
         endif
        else if (icount.gt.isurface_1_bc_count.and.
     &           icount.le.num_wake_bcs_total) then
c ... 2nd surface
         icount2 = icount - isurface_1_bc_count
         ibcval_wake_2(icount2) = ibcval(nw)
         nzbc_wake_2(icount2) = nzbc(nw)
         jbcb_wake_2(icount2) = jbcb(nw)
         jbce_wake_2(icount2) = jbce(nw)
         kbcb_wake_2(icount2) = kbcb(nw)
         kbce_wake_2(icount2) = kbce(nw)
         jmax_wake_2(icount2) = jmaxx(nzbc(nw)) - 1
         kmax_wake_2(icount2) = kmaxx(nzbc(nw)) - 1
         if (DEBUG) then
          write(*,2030) 2,icount,icount2,ibcval_wake_2(icount2),
     &     nzbc_wake_2(icount2),
     &     jbcb_wake_2(icount2),jbce_wake_2(icount2),
     &     kbcb_wake_2(icount2),kbce_wake_2(icount2),
     &     jmax_wake_2(icount2),kmax_wake_2(icount2)
2030      format(11i5)
         endif
        else
         print*
         print*, ' ERROR: icount out of bounds in rbc.f/rbcmain'
         print*, '   icount should be less than or equal to',
     &    num_wake_bcs_total  
         print*, '   current value of icount:',icount
         print*
         print*, ' isurface_1_bc_count =',isurface_1_bc_count
         print*, ' isurface_2_bc_count =',isurface_2_bc_count
         print*
         print*, ' Program terminating.'
         stop ' stop: icount out of bounds in rbc.f/rbcmain'
        endif
       endif
1460  continue

c      stop 'stop: after 1460 in rbc.f/rbcmain'

      call flush(6)
      call flush(istdout)
c      stop ' stop: after do 1470 in rbc.f'

c open file for wake bc zone links
      filename = 'bcwake.dat'
      open(150,file=filename,form='formatted',status='unknown',err=1910)
      goto 1920
1910  continue
      print*
      print*, ' ERROR: In rbc.f/rbcmain, cannot open file =',filename
      print*
      print*, '   Program must terminate.'
      print*
      stop 'stop: in rbc.f/rbcmain, cannot open bc file for wake'
1920  continue
      write(150,*) 'c These wake boundary conditions come in row pairs:'
      write(150,*) 'c - first row lists surface 1 target points at k=1'
      write(150,*) 'c - second row lists surface 2 base points at k=2'
      write(150,*) 'c - base points are extracted from first row at k=2'
      write(150,*) 'c    as well' 
      write(150,*) 'c q values at k=2 on surface 1 and surface 2 are'
      write(150,*) 'c averaged to calculate the wake cut bcs'
      write(150,*) ' '
      write(150,*) 'c zone jbeg jend jinc kbeg kend kinc'
      write(150,*) 'c ---- ---- ---- ---- ---- ---- ----'

c determine wake 'surface 1' zone links with 'surface 2'
c assumes j conforms to walls, k is perpendicular to walls, 
c  and k=1 surface is essentially along the wake centerline.
c Target points on surface 1 at k=1 use base points from surface 1
c  at k=2 as well as base points from surface 2 at k=2.
c Target and base points on surface 2 are then taken directly from
c  the determination of points used for surface 1 -- there is no need
c  to re-calculate target and base points for surface 2!
c Also assumes that no more than 2 base grids come within one interval
c  of the target grid; if grid partition is load balanced properly,
c  this constraint should not occur.
      i1 = 1
      i2 = isurface_2_bc_count
      istep = 0
      nbcreg_wake = 0
c search in backward j-direction along surface 2 for points matching surface 1
      jbcb_2 = jbcb_wake_2(i2)
      jbce_2 = jbce_wake_2(i2)
      jbcb_2_actual = jbcb_2
      jbce_2_actual = jbce_2
      idiff_2 = jbce_2 - jbcb_2 
c search forward along surface 1
c have to assume that there is point-to-point match between surfaces 1 and 2
      jbcb_1 = jbcb_wake_1(i1)
      jbce_1 = jbce_wake_1(i1)
      jbcb_1_actual = jbcb_1
      jbce_1_actual = jbce_1
      idiff_1 = jbce_1 - jbcb_1

      if (idiff_2.lt.idiff_1) then
c actual end of surface 1 is end of surface 2 
       jbce_1_actual = jbce_2 
      else if (idiff_2.eq.idiff_1) then
c nothing to do here
      else if (idiff_2.gt.idiff_1) then
c actual end of surface 2 is end of surface 1
       jbce_2_actual = jbce_1
      endif

c begin looping thru wake boundary conditions

1930  continue
      istep = istep + 1

      if (DEBUG) then
       print*
       print*, ' -------------------------------------'
       print*, ' istep =',istep
       print*, ' i1, i2 =',i1,i2
       print*, ' isurface_1_bc_count, isurface_2_bc_count = ',
     &   isurface_1_bc_count,isurface_2_bc_count
       print*, ' jbcb_1, jbce_1 =',jbcb_1,jbce_1
       print*, ' idiff_1 =',idiff_1
       print*, ' jbcb_2, jbce_2 =',jbcb_2,jbce_2
       print*, ' idiff_1, idiff_2 =',idiff_1,idiff_2
       if (idiff_2.lt.idiff_1) then
        print*, '    idiff_2 < idiff_1' 
       else if (idiff_2.eq.idiff_1) then
        print*, '    idiff_2 = idiff_1'
       else if (idiff_2.gt.idiff_1) then
        print*, '    idiff_2 > idiff_1'
       else
        print*, ' ERROR: wrong values for idiff_2 and/or idiff_1'
        print*, '   Program terminating'
        print*
        stop ' stop: wrong values for idiff_2 and/or idiff_1'
       endif
      endif


c  idiff_2 < idiff_1
c-------------------
      if (idiff_2.lt.idiff_1) then
c-------------------
c write to wake bc file
       if (DEBUG_150) then
        write(150,*) ' '
        write(150,*) 'istep = ',istep
        write(150,*) 'i1,i2 = ',i1,i2
        write(150,*) 'idiff_1,idiff_2 = ',idiff_1,idiff_2
       endif
       write(150,1990) nzbc_wake_1(i1),jbcb_1_actual,jbce_1_actual,
     &  1,1,1,1
       write(150,1990) nzbc_wake_2(i2),jbce_2_actual,jbcb_2_actual,
     &  -1,2,2,1
c increment nbcreg_wake
       nbcreg_wake = nbcreg_wake + 1 
c check limits on number of wake bc's
       if (nbcreg_wake.gt.ibcmax) then
        write(*,1560) nbcreg_wake,ibcmax
       endif
c ... target and base zones affected
       nzbc_t_wake(nbcreg_wake) = nzbc_wake_1(i1)
       nzbc_b_wake(nbcreg_wake) = nzbc_wake_2(i2)
c ... target indices and increment
       jbcb_t_wake(nbcreg_wake) = jbcb_1_actual
       jbce_t_wake(nbcreg_wake) = jbce_1_actual
       jinc_t_wake(nbcreg_wake) = 1
c ... base indices and increment
       jbcb_b_wake(nbcreg_wake) = jbce_2_actual
       jbce_b_wake(nbcreg_wake) = jbcb_2_actual
       jinc_b_wake(nbcreg_wake) = -1
       write(150,1990) nzbc_wake_2(i2),jbcb_2_actual,jbce_2_actual,
     &  1,1,1,1
       write(150,1990) nzbc_wake_1(i1),jbce_1_actual,jbcb_1_actual,
     &  -1,2,2,1
c increment nbcreg_wake
       nbcreg_wake = nbcreg_wake + 1 
c check limits on number of wake bc's
       if (nbcreg_wake.gt.ibcmax) then
        write(*,1560) nbcreg_wake,ibcmax
       endif
c ... target and base zones affected
       nzbc_t_wake(nbcreg_wake) = nzbc_wake_2(i2)
       nzbc_b_wake(nbcreg_wake) = nzbc_wake_1(i2)
c ... target indices and incremet
       jbcb_t_wake(nbcreg_wake) = jbcb_2_actual
       jbce_t_wake(nbcreg_wake) = jbce_2_actual
       jinc_t_wake(nbcreg_wake) = 1
c ... base indices and increment
       jbcb_b_wake(nbcreg_wake) = jbce_1_actual
       jbce_b_wake(nbcreg_wake) = jbcb_1_actual
       jinc_b_wake(nbcreg_wake) = -1
c decrement surface 2
       i2 = i2 - 1
c check if done
       if (i2.lt.1) then
        goto 1485
       endif
c calculate beginning indices for next pair of surfaces 1 and 2
c ... surface 1
       jbce_1_actual = jbce_wake_1(i1)
       jbcb_1_actual = jbce_wake_1(i1) - 
     &                 (jbce_wake_1(i1) - jbce_1_actual)
c ... surface 2
       jbce_2_actual = jbce_wake_2(i2)
       jbcb_2_actual = jbce_wake_2(i2) - 
     &                 (jbce_wake_1(i1) - jbce_1_actual)
c calculate differences
       idiff_1 = jbce_1_actual - jbcb_1_actual
       idiff_2 = jbce_2_actual - jbcb_2_actual

       if (DEBUG) then
        print*, 'For next output:'
        print*, ' i1,i2 =',i1,i2
        print*, ' nz1,nz2 =',nzbc_wake_1(i1),nzbc_wake_2(i2)
        print*, ' idiff_1,idiff_2 =',idiff_1,idiff_2
        print*, ' actual jbcb_1,jbce_1,jbcb_2,jbce_2 =',
     &   jbcb_1_actual,jbce_1_actual,jbcb_2_actual,jbce_2_actual
c        stop 'stop: 100 in rbc.f/rbcmain'
       endif

        
c  idiff_2 = idiff_1
c-------------------
      else if (idiff_2.eq.idiff_1) then
c-------------------
       if (i1.eq.isurface_1_bc_count) then
c last subgrid for surface 1, so recalculate jbcb_2_actual and jbce_2_actual 
c   for surface 2 since surface 1 length is probably shorter than 
c         jbcb_2_actual = jbce_2_actual - idiff_1
       endif

       if (DEBUG_150) then
        write(150,*) ' '
        write(150,*) 'istep = ',istep
        write(150,*) 'i1,i2 = ',i1,i2
        write(150,*) 'idiff_1,idiff_2 = ',idiff_1,idiff_2
       endif
       write(150,1990) nzbc_wake_1(i1),jbcb_1_actual,jbce_1_actual,
     &  1,1,1,1
       write(150,1990) nzbc_wake_2(i2),jbce_2_actual,jbcb_2_actual,
     &  -1,2,2,1
c increment nbcreg_wake
       nbcreg_wake = nbcreg_wake + 1 
c check limits on number of wake bc's
       if (nbcreg_wake.gt.ibcmax) then
        write(*,1560) nbcreg_wake,ibcmax
       endif
c ... target and base zones affected
       nzbc_t_wake(nbcreg_wake) = nzbc_wake_1(i1)
       nzbc_b_wake(nbcreg_wake) = nzbc_wake_2(i2)
c ... target indices and increment
       jbcb_t_wake(nbcreg_wake) = jbcb_1_actual
       jbce_t_wake(nbcreg_wake) = jbce_1_actual
       jinc_t_wake(nbcreg_wake) = 1
c ... base indices and increment
       jbcb_b_wake(nbcreg_wake) = jbce_2_actual
       jbce_b_wake(nbcreg_wake) = jbcb_2_actual
       jinc_b_wake(nbcreg_wake) = -1
       write(150,1990) nzbc_wake_2(i2),jbcb_2_actual,jbce_2_actual,
     &  1,1,1,1
       write(150,1990) nzbc_wake_1(i1),jbce_1_actual,jbcb_1_actual,
     &  -1,2,2,1
c increment nbcreg_wake
       nbcreg_wake = nbcreg_wake + 1 
c check limits on number of wake bc's
       if (nbcreg_wake.gt.ibcmax) then
        write(*,1560) nbcreg_wake,ibcmax
       endif
c ... target and base zones affected
       nzbc_t_wake(nbcreg_wake) = nzbc_wake_2(i2)
       nzbc_b_wake(nbcreg_wake) = nzbc_wake_1(i2)
c ... target indices and incremet
       jbcb_t_wake(nbcreg_wake) = jbcb_2_actual
       jbce_t_wake(nbcreg_wake) = jbce_2_actual
       jinc_t_wake(nbcreg_wake) = 1
c ... base indices and increment
       jbcb_b_wake(nbcreg_wake) = jbce_1_actual
       jbce_b_wake(nbcreg_wake) = jbcb_1_actual
       jinc_b_wake(nbcreg_wake) = -1
c move to next subgrids of either or both surfaces 1 and 2
        if (DEBUG) then
         print*
         print*, ' idiff_2 = idiff_1: ',idiff_2,idiff_1
         print*, ' jbce_1_actual, jbce_2_actual = ',
     &      jbce_1_actual, jbce_2_actual
         print*, ' nz1,nz2,jmaxx1,jmaxx2 = ',
     &      nzbc_wake_1(i1),nzbc_wake_2(i2),
     &      jmaxx(nzbc_wake_1(i1)),jmaxx(nzbc_wake_2(i2))
         print*, ' i1,i2, = ',i1,i2
        endif
        if (jbce_1_actual.eq.jmaxx(nzbc_wake_1(i1)) - 1) then
         i1 = i1 + 1
         if (i1.gt.isurface_1_bc_count) goto 1485
        endif
        if (jbce_2_actual.eq.1) then
         i2 = i2 - 1 
         if (i2.lt.1) goto 1485
        endif
c calculate beginning and end indices for next pair of surfaces 1 and 2
c increment surface 1
        i1 = i1 + 1
c decrement surface 2
        i2 = i2 - 1
c check if done
        if (i1.gt.isurface_1_bc_count.or.i2.lt.1) then
         goto 1485
        endif
c ... surface 1
        jbcb_1_actual = jbcb_wake_1(i1)
        jbce_1_actual = jbce_wake_1(i1)
c ... surface 2
        jbcb_2_actual = jbcb_wake_2(i2)
        jbce_2_actual = jbce_wake_2(i2)
c calculate differences
        idiff_1 = jbce_1_actual - jbcb_1_actual
        idiff_2 = jbce_2_actual - jbcb_2_actual
c make sure no endpoints are involved in surface 1
        if ((jbce_1_actual - idiff_2).eq.1) then
         jbce_2_actual = jbce_2_actual - 1
         idiff_2 = idiff_2 - 1
        endif

        if (DEBUG) then
         print*, 'For next output:'
         print*, ' i1,i2: ',i1,i2
         print*, ' nz1,nz2 =',nzbc_wake_1(i1),nzbc_wake_2(i2)
         print*, ' idiff_1,idiff_2 =',idiff_1,idiff_2
         print*, ' actual jbcb_1,jbce_1,jbcb_2,jbce_2 =',
     &    jbcb_1_actual,jbce_1_actual,jbcb_2_actual,jbce_2_actual
         if (istep.gt.2) then
          stop ' stop: 200 in rbc.f/rbcmain'
         endif 
        endif
         

c  idiff_2 > idiff_1
c-------------------
       else if (idiff_2.gt.idiff_1) then
c-------------------
c write to wake bc file
        if (DEBUG_150) then
         write(150,*) ' '
         write(150,*) 'istep = ',istep
         write(150,*) 'i1,i2 = ',i1,i2
         write(150,*) 'idiff_1,idiff_2 = ',idiff_1,idiff_2
        endif  
        write(150,1990) nzbc_wake_1(i1),jbcb_1_actual,jbce_1_actual,
     &   1,1,1,1
        write(150,1990) nzbc_wake_2(i2),jbce_2_actual,jbcb_2_actual,
     &   -1,2,2,1
c increment nbcreg_wake
       nbcreg_wake = nbcreg_wake + 1 
c check limits on number of wake bc's
       if (nbcreg_wake.gt.ibcmax) then
        write(*,1560) nbcreg_wake,ibcmax
       endif
c ... target and base zones affected
       nzbc_t_wake(nbcreg_wake) = nzbc_wake_1(i1)
       nzbc_b_wake(nbcreg_wake) = nzbc_wake_2(i2)
c ... target indices and increment
       jbcb_t_wake(nbcreg_wake) = jbcb_1_actual
       jbce_t_wake(nbcreg_wake) = jbce_1_actual
       jinc_t_wake(nbcreg_wake) = 1
c ... base indices and increment
       jbcb_b_wake(nbcreg_wake) = jbce_2_actual
       jbce_b_wake(nbcreg_wake) = jbcb_2_actual
       jinc_b_wake(nbcreg_wake) = -1
        write(150,1990) nzbc_wake_2(i2),jbcb_2_actual,jbce_2,actual,
     &   1,1,1,1
        write(150,1990) nzbc_wake_1(i1),jbce_1_actual,jbcb_1_actual,
     &   -1,2,2,1
c increment nbcreg_wake
       nbcreg_wake = nbcreg_wake + 1 
c check limits on number of wake bc's
       if (nbcreg_wake.gt.ibcmax) then
        write(*,1560) nbcreg_wake,ibcmax
       endif
c ... target and base zones affected
       nzbc_t_wake(nbcreg_wake) = nzbc_wake_2(i2)
       nzbc_b_wake(nbcreg_wake) = nzbc_wake_1(i2)
c ... target indices and incremet
       jbcb_t_wake(nbcreg_wake) = jbcb_2_actual
       jbce_t_wake(nbcreg_wake) = jbce_2_actual
       jinc_t_wake(nbcreg_wake) = 1
c ... base indices and increment
       jbcb_b_wake(nbcreg_wake) = jbce_1_actual
       jbce_b_wake(nbcreg_wake) = jbcb_1_actual
       jinc_b_wake(nbcreg_wake) = -1
c increment surface 1
        i1 = i1 + 1
c finish if i1 has run out of surfaces 
        if (i1.gt.isurface_1_bc_count) then
         goto 1485
        endif
c calculate beginning indices for next pair of surfaces 1 and 2
c ... surface 1
        jbce_1_actual = jbce_wake_2(i2) - jbce_wake_1(i1-1) + 1
        jbcb_1_actual = jbcb_wake_1(i1)
c ... surface 2
        jbce_2_actual = jbce_wake_2(i2)
        jbcb_2_actual = jbce_wake_2(i2) - 
     &                  (jbce_1_actual - jbcb_1_actual)
c calculate differences
        idiff_1 = jbce_1_actual - jbcb_1_actual
        idiff_2 = jbce_2_actual - jbcb_2_actual

        if (DEBUG) then
         print*, 'For next output:'
         print*, ' i1,i2 =',i1,i2
         print*, ' nz1,nz2 =',nzbc_wake_1(i1),nzbc_wake_2(i2)
         print*, ' idiff_1,idiff_2 =',idiff_1,idiff_2
         print*, ' actual jbcb_1,jbce_1,jbcb_2,jbce_2 =',
     &    jbcb_1_actual,jbce_1_actual,jbcb_2_actual,jbce_2_actual
c         stop 'stop: 300 in rbc.f/rbcmain'
        endif
         
c move to next subgrids of either or both surfaces 1 and 2

c        stop 'stop: 400 in rbc.f/rbcmain'

c end of 'idiff' comparisons
c------------------
       endif
c------------------

       call flush(6)
       call flush(istdout)

c      stop 'stop: after first pair of wake bc indices'
      goto 1930


1485   continue

      close(150)

c      stop 'stop: after 1485 in rbc.f/rbcmain'

c deallocate ibcval_wake, etc.
c ... surface 1
      deallocate (ibcval_wake_1,nzbc_wake_1,
     &  jbcb_wake_1,jbce_wake_1,
     &  kbcb_wake_1,kbce_wake_1,
     &  jmax_wake_1, kmax_wake_1
     & )
c ... surface 2
      deallocate (ibcval_wake_2,nzbc_wake_2,
     &  jbcb_wake_2,jbce_wake_2,
     &  kbcb_wake_2,kbce_wake_2,
     &  jmax_wake_2,kmax_wake_2
     & )

600   continue

c      stop 'stop: after 600 continue in rbc.f/rbcmain'


c-----
c  Error checking on original main bc file input
c-----
c ... check if any bc zone is greater than number of processors`
      print*, ' --- in rbc.f ---- '
      print*, '  nr   nzbc   jbcb   kbcb   jbce   kbce   jmaxx   kmaxx'
      print*, ' ----  ----  ------ ------ ------ ------ ------- -------'
      do 92 nr=1,nbcreg
         write(*,920) nr,nzbc(nr),jbcb(nr),kbcb(nr),jbce(nr),kbce(nr),
     1    jmaxx(nzbc(nr)),kmaxx(nzbc(nr))
920      format(1x,i4,3x,i4,2x,i6,2x,i6,2x,i6,2x,i6,2x,i7,2x,i7)
         if( nzbc(nr) .gt. numprocs ) then
            write(istdout,93)
            write(istdout,94) nr
            print*, ' ERROR: nzbc(nr) is greater than numprocs'
            print*, '   Check file bcmain.dat'
            print*, '   Stopping: execution halted on Node 0'
            stop ' stopping: nzbc(nr) is greater than numprocs'
         endif

c ... check if any k or j values for boundary conditions exceed
c      the subgrid dimensions
         if( jbce(nr) .gt. jmaxx(nzbc(nr)) .or.
     &       kbce(nr) .gt. kmaxx(nzbc(nr)) ) then
            write(istdout,93)
            write(istdout,95) jbce(nr), kbce(nr), nzbc(nr)
            print*, ' ERROR: k or j bc value exceeds max subgrid dim'
            print*, '   Check main bc file ', trim(bcmain_filename)
            print*, '   Stopping: execution halted on Node 0'
            stop ' stopping: k or j value for bc exceeds subgrid dim'
         endif

c ... check if any beginning j or k values for bc's exceed end values
         if( jbcb(nr) .gt. jbce(nr) .or.
     &       kbcb(nr) .gt. kbce(nr) ) then
            write(istdout,93)
            write(istdout,96) nr
            print*, ' ERROR: beginning j or k value exceed end value'
            print*, '   Check file bcmain.dat'
            print*, '   Stopping: execution halted on Node 0'
            stop ' stopping: beginning j or k value exceeds end value'
         endif
92    continue

      close(99)

       print*
       print*,' NODE 0: All boundary conditions have been checked.'
       print*, '  ... number of regular boundary conditions: ',nbcreg
       print*, '  ... number of wall boundary conditions: ',nwall

c================================
      if (numprocs.gt.1) then
c================================
c send boundary conditions to all nodes
c general bcs
c ... pack up variables

      allocate(ivarbc(6,nbcreg))

c ... pack indices to all subgrid boundary conditions, to be sent
c  to all nodes
      do i=1,nbcreg
        ivarbc(1,i) = ibcval(i)
        ivarbc(2,i) = nzbc(i)
        ivarbc(3,i) = jbcb(i)
        ivarbc(4,i) = jbce(i)
        ivarbc(5,i) = kbcb(i)
        ivarbc(6,i) = kbce(i)
      enddo

c wall bcs
      allocate (ivarwall(7,nwall))

c ... pack indices to subgrid wall boundary conditions, to be sent
c  to all nodes
      do i=1,nwall
        ivarwall(1,i) = jkwall(i)
        ivarwall(2,i) = nzwall(i)
        ivarwall(3,i) = jwall1(i)
        ivarwall(4,i) = jwall2(i)
        ivarwall(5,i) = kwall1(i)
        ivarwall(6,i) = kwall2(i)
        ivarwall(7,i) = jkinc(i)
      enddo

c wake bcs
      allocate (ivarbc_wake(8,nbcreg_wake))

c ... pack indices to sugrid wake cut boundary conditions, to be sent
c  to all nodes
      do i=1,nbcreg_wake
       ivarbc_wake(1,i) = nzbc_t_wake(i)
       ivarbc_wake(2,i) = nzbc_b_wake(i)
       ivarbc_wake(3,i) = jbcb_t_wake(i)
       ivarbc_wake(4,i) = jbce_t_wake(i)
       ivarbc_wake(5,i) = jinc_t_wake(i)
       ivarbc_wake(6,i) = jbcb_b_wake(i)
       ivarbc_wake(7,i) = jbce_b_wake(i)
       ivarbc_wake(8,i) = jinc_b_wake(i)
      enddo

c end 'if (numprocs.gt.1) then'
c================================
      endif
c================================

c end node 0
c------------------
      endif
c------------------


c send variables
c ... nbcreg

c================================
      if (numprocs.gt.1) then
c================================
      if (nodeid.eq.0) then
       if (DEBUG2) then
        print*
        print*, 'nodeid,numprocs,nbcreg,MPI_SEND(nbcreg,)',
     &   nodeid,numprocs,nbcreg
       endif 
       do i=2,numprocs
        node_out = i-1
        call MPI_SEND(nbcreg,1,MPI_INTEGER,node_out,node_out,
     &    MPI_COMM_WORLD,ierr)
       enddo
      else
       call MPI_RECV(nbcreg,1,MPI_INTEGER,0,nodeid,
     &  MPI_COMM_WORLD,stat,ierr)
       if (DEBUG2) then
        print*
        print*, 'nodeid,numprocs,nbcreg,MPI_SEND(nbcreg,)',
     &   nodeid,numprocs,nbcreg
       endif 
      endif
c================================
      endif
c================================

c===================================
      if (numprocs.gt.1) then
c===================================

c ... boundary conditions from bcmain.dat
      if (nodeid.eq.0) then
       if (DEBUG2) then
        print*
        print*, 'nodeid,numprocs,MPI_SEND(ivarbc,)',nodeid,numprocs
       endif 
       do i=2,numprocs
        node_out = i-1 
        call MPI_SEND(ivarbc,6*nbcreg,MPI_INTEGER,node_out,node_out+1,
     &    MPI_COMM_WORLD,ierr)
       enddo
      else
c allocate ivarbc
       allocate(ivarbc(6,nbcreg))
c receive boundary conditions from node 0
       call MPI_RECV(ivarbc,6*nbcreg,MPI_INTEGER,0,nodeid+1,
     &  MPI_COMM_WORLD,stat,ierr)
       if (DEBUG2) then
        print*
        print*, 'nodeid,numprocs,MPI_RECV(ivarbc,)',nodeid,numprocs
       endif 
c unpack variables
       do i=1,nbcreg
        ibcval(i) = ivarbc(1,i)
        nzbc(i) = ivarbc(2,i)
        jbcb(i) = ivarbc(3,i)
        jbce(i) = ivarbc(4,i)
        kbcb(i) = ivarbc(5,i)
        kbce(i) = ivarbc(6,i)
       enddo
      endif

c ... nwall
      if (nodeid.eq.0) then
       if (DEBUG2) then
        print*
        print*, 'nodeid,numprocs,nwall,MPI_SEND(nwall,)',
     &   nodeid,numprocs,nwall
       endif 
       do i=2,numprocs
        node_out = i-1
        call MPI_SEND(nwall,1,MPI_INTEGER,node_out,node_out+2,
     &    MPI_COMM_WORLD,ierr)
       enddo
      else
c receive number of wall boundary conditions
       call MPI_RECV(nwall,1,MPI_INTEGER,0,nodeid+2,
     &  MPI_COMM_WORLD,stat,ierr)
       if (DEBUG2) then
        print*
        print*, 'nodeid,numprocs,nwall,MPI_RECV(nwall,)',
     &   nodeid,numprocs,nwall
       endif 
      endif

c ... nwall boundary conditions
      if (nodeid.eq.0) then
       if (DEBUG2) then
        print*
        print*, 'nodeid,numprocs,MPI_SEND(ivarwall,)',nodeid,numprocs
       endif 
       do i=2,numprocs
        node_out = i-1
        call MPI_SEND(ivarwall,7*nwall,MPI_INTEGER,node_out,node_out+3,
     &   MPI_COMM_WORLD,ierr)
       enddo
      else
c allocate ivarwall
       allocate (ivarwall(7,nwall))
c receive wall boundary conditions
       call MPI_RECV(ivarwall,7*nwall,MPI_INTEGER,0,nodeid+3,
     &  MPI_COMM_WORLD,stat,ierr)
       if (DEBUG2) then
        print*
        print*, 'nodeid,numprocs,MPI_RECV(ivarwall,)',nodeid,numprocs
       endif 
c unpack variables
       do i=1,nwall
        jkwall(i) = ivarwall(1,i)
        nzwall(i) = ivarwall(2,i)
        jwall1(i) = ivarwall(3,i)
        jwall2(i) = ivarwall(4,i)
        kwall1(i) = ivarwall(5,i)
        kwall2(i) = ivarwall(6,i)
         jkinc(i) = ivarwall(7,i)
       enddo
      endif

c ... jkbc values
      if (nodeid.eq.0) then
       if (DEBUG2) then
        print*
        print*, 'nodeid,numprocs,MPI_SEND(jkbc,)',nodeid,numprocs
       endif 
       do i=2,numprocs
        node_out = i-1
        call MPI_SEND(jkbc,nbcreg,MPI_INTEGER,node_out,node_out+4,
     &   MPI_COMM_WORLD,ierr)
       enddo
      else
c receive jkbc values
       call MPI_RECV(jkbc,nbcreg,MPI_INTEGER,0,nodeid+4,
     &  MPI_COMM_WORLD,stat,ierr)
       if (DEBUG2) then
        print*
        print*, 'nodeid,numprocs,MPI_RECV(jkbc,)',nodeid,numprocs
       endif 
      endif

c ... send number of wake boundary conditions
      if (nodeid.eq.0) then
       if (DEBUG2) then
        print*
        print*, 'nodeid,numprocs,nbcreg_wake,MPI_SEND(nbcreg_wake,)',
     &   nodeid,numprocs,nbcreg_wake
       endif 
       do i=2,numprocs
        node_out = i - 1
        call MPI_SEND(nbcreg_wake,1,MPI_INTEGER,node_out,node_out+5,
     &   MPI_COMM_WORLD,ierr)
       enddo
      else
c ... receive number of wake boundary conditions
       call MPI_RECV(nbcreg_wake,1,MPI_INTEGER,0,nodeid+5,
     &  MPI_COMM_WORLD,stat,ierr)
       if (DEBUG2) then
        print*
        print*, 'nodeid,numprocs,nbcreg_wake,MPI_RECV(nbcreg_wake,)',
     &   nodeid,numprocs,nbcreg_wake
       endif 
      endif

      if (nbcreg_wake.gt.0) then
c ... send wake boundary conditions
       if (nodeid.eq.0) then
        do i=2,numprocs
         node_out = i - 1
         call MPI_SEND(ivarbc_wake,8*nbcreg_wake,MPI_INTEGER,node_out,
     &                node_out+6,MPI_COMM_WORLD,ierr)
        enddo
       else
c allocate ivarbc_wake
        allocate (ivarbc_wake(8,nbcreg_wake))
c receive wake bcs
        call MPI_RECV(ivarbc_wake,8*nbcreg_wake,MPI_INTEGER,0,
     &   nodeid+6,MPI_COMM_WORLD,stat,ierr)
        if (DEBUG2) then
         print*
         print*, 'nodeid,numprocs,MPI_RECV(ivarbc_wake,)',nodeid,numprocs
        endif 
c unpack variables
        do i=1,nbcreg_wake
         nzbc_t_wake(i) = ivarbc_wake(1,i)
         nzbc_b_wake(i) = ivarbc_wake(2,i)
         jbcb_t_wake(i) = ivarbc_wake(3,i)
         jbce_t_wake(i) = ivarbc_wake(4,i)
         jinc_t_wake(i) = ivarbc_wake(5,i)
         jbcb_b_wake(i) = ivarbc_wake(6,i)
         jbce_b_wake(i) = ivarbc_wake(7,i)
         jinc_b_wake(i) = ivarbc_wake(8,i)
        enddo
       endif
      endif

      deallocate (ivarbc,ivarwall)
      if (nbcreg_wake.gt.0) then
       deallocate (ivarbc_wake)
      endif

c===================================
c end of 'if (numprocs.gt.1) then'
      endif
c===================================

      if (DEBUG2) then
       print*
       print*, 'After deallocating ivarbc,ivarwall:'
       print*, ' nodeid,nbcreg_wake,numprocs =',
     &   nodeid, nbcreg_wake, numprocs
       call flush(6)
       call flush(istdout)
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call flush(6)
       call flush(istdout)
       stop 'stop: after deallocating ivarbc,ivarwall in rbc.f/rbcmain'
      endif

c calculate local x,y values for wall and send to each node for use in
c  subroutine dfunsa in tursa.f
      nodecountwall = 0
      do 5000 nw=1,nwall
       nodewall = nzwall(nw) - 1
       if (DEBUG) then
        write(istdout,6080) nw,nwall,nodewall,nodeid,jwall1(nw),
     1   kwall1(nw),jwall2(nw),kwall2(nw),jkwall(nw)
6080    format('Loop 5000: nw,nwall,nodewall,nodeid,jwall1(nw),kwall1(nw
     1),jwall2(nw),kwall2(nw),jkwall(nw)  = ',9i4)
       endif
c1-------------------------------------
       if (nodeid.eq.nodewall) then
c1-------------------------------------
        nodecountwall = nodecountwall + 1
c calculate xwallval, ywallval which are local grid point values at the wall
        if (jkwall(nw).eq.1) then
c ... j=constant
         jmaxwall = jwall2(nw)  ! = jwall1(nw)
         numpointswall = kwall2(nw) - kwall1(nw) + 1
         kcount = 0
         do 5010 k=kwall1(nw),kwall2(nw),jkinc(nw)
          kcount = kcount + 1
          xwallval(kcount) = x(jmaxwall,k)
          ywallval(kcount) = y(jmaxwall,k)
5010     continue
        else
c ... k=constant
         numpointswall = jwall2(nw) - jwall1(nw) + 1
         kmaxwall = kwall2(nw)  ! = kwall1(nw)
         jcount = 0
         do 5020 j=jwall1(nw),jwall2(nw),jkinc(nw)
          jcount = jcount + 1
          xwallval(jcount) = x(j,kmaxwall)
          ywallval(jcount) = y(j,kmaxwall)
          if (DEBUG) then
           print*, 'do 5020: j,jcount,xwallval,ywallval = ',
     1      j,jcount,xwallval(jcount),ywallval(jcount)
          endif
5020     continue
        endif


c collect all results on node 0
c2----------------------------
       if (nodewall.eq.0) then
c2----------------------------
c ... we are already on node 0 if we come here
        if (jkwall(nw).eq.1) then
         do 5030 k=1,kcount
          xwallv(k,nw) = xwallval(k)
          ywallv(k,nw) = ywallval(k)
5030     continue
        elseif (jkwall(nw).eq.2) then
         do 5040 j=1,jcount
          xwallv(j,nw) = xwallval(j)
          ywallv(j,nw) = ywallval(j)
5040     continue
        endif
c skip mpi_send and mpi_recv
        goto 5000

c2----------------------------
       elseif (numprocs.gt.1) then
c2----------------------------
c send to node 0
         if (jkwall(nw).eq.1) then
          call MPI_SEND(kcount,1,MPI_INTEGER,0,nodewall,
     1     MPI_COMM_WORLD,ierr)

#ifdef D_PRECISION
          call MPI_SEND(xwallval,kcount,MPI_REAL8,0,nodewall+1,
     1     MPI_COMM_WORLD,ierr)
          call MPI_SEND(ywallval,kcount,MPI_REAL8,0,nodewall+2,
     1     MPI_COMM_WORLD,ierr)
#else
          call MPI_SEND(xwallval,kcount,MPI_REAL,0,nodewall+1,
     1     MPI_COMM_WORLD,ierr)
          call MPI_SEND(ywallval,kcount,MPI_REAL,0,nodewall+2,
     1     MPI_COMM_WORLD,ierr)
#endif

         elseif (jkwall(nw).eq.2) then

          call MPI_SEND(jcount,1,MPI_INTEGER,0,nodewall,
     1     MPI_COMM_WORLD,ierr)
#ifdef D_PRECISION
          call MPI_SEND(xwallval,jcount,MPI_REAL8,0,nodewall+1,
     1     MPI_COMM_WORLD,ierr)
          call MPI_SEND(ywallval,jcount,MPI_REAL8,0,nodewall+2,
     1     MPI_COMM_WORLD,ierr)
#else
          call MPI_SEND(xwallval,jcount,MPI_REAL,0,nodewall+1,
     1     MPI_COMM_WORLD,ierr)
          call MPI_SEND(ywallval,jcount,MPI_REAL,0,nodewall+2,
     1     MPI_COMM_WORLD,ierr)
#endif

         endif 
c2---------------------------
      endif
c2---------------------------

c1---------------------------
      elseif (nodeid.eq.0.and.nodewall.ne.0) then
c1---------------------------

         call MPI_RECV(icount,1,MPI_INTEGER,nodewall,nodewall,
     1    MPI_COMM_WORLD,stat,ierr)
#ifdef D_PRECISION
         call MPI_RECV(xwallval,icount,MPI_REAL8,nodewall,nodewall+1,
     1    MPI_COMM_WORLD,stat,ierr)
         call MPI_RECV(ywallval,icount,MPI_REAL8,nodewall,nodewall+2,
     1    MPI_COMM_WORLD,stat,ierr)
#else
         call MPI_RECV(xwallval,icount,MPI_REAL,nodewall,nodewall+1,
     1    MPI_COMM_WORLD,stat,ierr)
         call MPI_RECV(ywallval,icount,MPI_REAL,nodewall,nodewall+2,
     1    MPI_COMM_WORLD,stat,ierr)
#endif

         do 5050 i=1,icount
          xwallv(i,nw) = xwallval(i)
          ywallv(i,nw) = ywallval(i)
5050     continue

c1---------------------------
      endif
c1---------------------------

5000  continue

c print x,y wall values
      if (DEBUG.and.nodeid.eq.0) then
      do 5500 nw=1,nwall
       if (jkwall(nw).eq.1) then
c ... j=constant
        numpointswall = kwall2(nw) - kwall1(nw) + 1
       elseif (jkwall(nw).eq.2) then
c ... k=constant
        numpointswall = jwall2(nw) - jwall1(nw) + 1
       endif
       print*
       print*, ' nw  npw        xwallv         ywallv'
       do 5500 npw=1,numpointswall
        write(*,5510) nw,npw,xwallv(npw,nw),ywallv(npw,nw)
5510    format(2i5,1p2e13.5)
5500   continue
c      stop 'stop: in rbc.f/rbcmain 100' 
      endif

c send x,y wall values to all mpi ranks, to be used in subroutine dfunsa
      if (numprocs.gt.1) then

#ifdef D_PRECISION
c ... xwallv values
       call MPI_BCAST(xwallv,jkmax*nwall,MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr)
c ... ywallv values
       call MPI_BCAST(ywallv,jkmax*nwall,MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr)
#else
c ... xwallv values
       call MPI_BCAST(xwallv,jkmax*nwall,MPI_REAL,0,
     &  MPI_COMM_WORLD,ierr)
c ... ywallv values
       call MPI_BCAST(ywallv,jkmax*nwall,MPI_REAL,0,
     &  MPI_COMM_WORLD,ierr)
#endif

      endif

c------------------------------------------------------------
c   format statements
c------------------------------------------------------------
5     format('c filename: ',a50,
     1/,'c Flow over an airfoil, C-grid',
     2/,'c - subgrids are overlapped',
     3/,'c - NOTE: All bcs use j=1 to jmax, or k=1 to kmax'
     4/,'c',
     5/,'c------------------------------------------------------------',
     6/,'c   ibcval   zone   jbeg   jend   kbeg   kend',
     7/,'c------------------------------------------------------------'
     8 )
93    format(//,' ERROR in main bc file input:')
94    format('      Zone number for boundary condition region number',
     &              i4,/,
     &       '      exceeds number of zones in grid.  Stopping.')
95    format('      One of the two ending indices = ',2i5,/
     &       '      exceeds maximum dimensions in zone number ',i4,/
     &       '      Stopping.')
96    format('      One of the beginning indices exceeds the ending ',/
     &       '      index for boundary condition region number ',i4,/
     &       '      Stopping.')
1560  format(//,' BOUNDS ERROR: nbcreg_wake > ibcmax',//,
     &          '   nbcreg_wake = ',i5,/,
     &          '   ibcmax = ',i5,/,
     &          ' nbcreg_wake must be less than ibcmax!',/,
     &          ' adjust ibcmax in the file common.f and re-run.',//,
     &          ' Program terminated.',/
     &       )       
1990   format(7i5)

c-----
c  End of rbcmain
c-----
c      stop 'stop in rbc.f/rbcmain at end'

      return
      end
c
c
c************************************************************************
      subroutine rbcpzn(jmax,kmax)
c************************************************************************
c  This routine reads in the following boundary condition file:
c     bcpzn.dat:   Patched zones with one-to-one grid point connectivity
c------------------------------------------------------------------------
#include "common.f"

#include "mpif.h"
#include "mpi_params.f"
      __INTEGER stat(MPI_STATUS_SIZE)
      __INTEGER, allocatable:: ivari(:,:)
      __INTEGER jmaxx(numprocs),kmaxx(numprocs)
      __INTEGER ivarjk(2)

      __INTEGER jmax, kmax
      character*80 filename
      character*80 line
      character*20 string
      character*1 tabchar, space
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.
c
c+++++++++++++++++++++++++++
      if (nodeid.eq.0) then
c+++++++++++++++++++++++++++

      open(99,status='scratch')
2     format(a80)
      tabchar = '	'
      space = ' '
c-----
c  File name
c-----
      filename = 'bcpzn.dat'
      open(51,file=filename,form='formatted',status='old',err=10)
      goto 20
10    continue
         write(istdout,11) filename
11       format(' ERROR: patch zonal interface file not found: ',a,/,
     &          '        Stopping.')
      call MPI_FINALIZE(ierr)
      stop ' stopping: file bcpzn.dat not found'

20    continue
c
      nr = 1
30    continue

c read target (recipient) zone and points

c-----
c  Read in patched zonal interface info
c   b = base; t = target in index names
c   for example, jeb is the ending jth index for the base grid
c-----

         read(51,2,end=80) line
c
c  Memory check
c
         if(nr .gt. ibcmax) then
            write(istdout,35) ibcmax
            call exit(1)
         endif
35       format(/,' ERROR: dimensioning parameter ibcmax is less than',
     &            '  number of input lines in file bcpzn.dat.',
     &            '  ibcmax: ',i5)
c
         if( line(1:1) .eq. 'c' .or. line(1:1) .eq. 'C') goto 30
         do 40 icol=1,80
            if( line(icol:icol) .eq. tabchar) line(icol:icol) = space
40       continue
         rewind(99)
         write(99,2) line
         rewind(99)
         read(99,*,end=30) nzt(nr),jbt(nr),jet(nr),
     &                      kbt(nr),ket(nr)
50       continue

c read base (donor) zone and points

         read(51,2) line
         if( line(1:1) .eq. 'c' .or. line(1:1) .eq. 'C') goto 50
         do 60 icol=1,80
            if( line(icol:icol) .eq. tabchar) line(icol:icol) = space
60       continue
         rewind(99)
         write(99,2) line
         rewind(99)
         read(99,*,end=50) nzb(nr),jbb(nr),jeb(nr),
     &                      kbb(nr),keb(nr)

         nr = nr + 1
         goto 30
80    continue
c-----
c  End of file reached
c-----
      nreg = nr - 1
      close(51)

      if (DEBUG) then
       print*
       print*, ' TARGET POINTS:'
       print*, ' nreg = ',nreg
       print*, '        #          nzt         jbt         jet          
     1kbt          ket'
       do nrt=1,nreg
        print*,nrt,nzt(nrt),jbt(nrt),jet(nrt),kbt(nrt),ket(nrt)
       enddo
       print*
       print*, ' BASE POINTS:'
       print*, '        #          nzb         jbb         jeb          
     1kbb          keb'
       do nrb=1,nreg
        print*,nrb,nzb(nrb),jbb(nrb),jeb(nrb),kbb(nrb),keb(nrb)
       enddo
       print*
      endif


c check if any target bc's are listed more than once
      print*
      print*, ' File bcpzn.dat -- Checking if target bcs are listed more 
     1 than once' 
      ierrflagi1 = 0

      do nrt=1,nreg-1
        nztarget = nzt(nrt)
        do nrt2 = nrt+1,nreg
          if (nzt(nrt2).eq.nztarget) then
            if (jbt(nrt2).eq.jbt(nrt).and.kbt(nrt2).eq.kbt(nrt).and.
     1          jet(nrt2).eq.jet(nrt).and.ket(nrt2).eq.ket(nrt)) then
             ierrflag1 = ierrflag1 + 1
             print*, ' > Identical boundary condition, file bcpzn.dat:'
             print*, '      target subgrid: ',nztarget
             print*, '      identical target subgrid: ',nzt(nrt2)
             print*, '      indices nrt, nrt2: ',nrt,nrt2
             print*, '      jbt, kbt: ',jbt(nrt),kbt(nrt)
             print*, '      jet, ket: ',jet(nrt),ket(nrt)
             print*, '      jbt2, kbt2: ',jbt(nrt2),kbt(nrt2)
             print*, '      jet2, ket2: ',jet(nrt2),ket(nrt2)
            endif
          endif
        enddo
       enddo

c check if any donor bc's are listed more than once
      print*
      print*, ' File bcpzn.dat -- Checking if donor bcs are listed more 
     1than once'
      ierrflag2 = 0

      do nrb=1,nreg-1
        nzbase = nzb(nrb)
        do nrb2 = nrb+1,nreg
          if (nzb(nrb2).eq.nzbase) then
            if (jbb(nrb2).eq.jbb(nrb).and.kbb(nrb2).eq.kbb(nrb).and.
     1          jeb(nrb2).eq.jeb(nrb).and.keb(nrb2).eq.ket(nrb)) then
             ierrflag2 = ierrflag2 + 1
             print*, ' > Identical boundary condition, file bcpzn.dat:'
             print*, '      donor subgrid: ',nzbase
             print*, '      identical donor subgrid: ',nzb(nrb2)
             print*, '      indices nrb, nrb2: ',nrb,nrb2
             print*, '      jbb, kbb: ',jbb(nrb),kbb(nrb)
             print*, '      jeb, keb: ',jeb(nrb),keb(nrb)
             print*, '      jbb2, kbb2: ',jbb(nrb2),kbb(nrb2)
             print*, '      jeb2, keb2: ',jeb(nrb2),keb(nrb2)
            endif
          endif
        enddo
       enddo

c ... check all endpoints for multiple bcs 

c ... for kbt = ket (constant k) 
        print*
        print*, ' File bcpzn.dat -- Checking j endpoints, k=constant'
        ierrflag3 = 0

c ... simply flags points that can be updated from multiple base points,
c ...   for j endpoints, k=constant
c DWB: SKIP THIS FOR NOW

        goto 950

        do nrt=1,nreg
          if (kbt(nrt).eq.ket(nrt)) then
           nztsave = nzt(nrt)
           jbtsave = jbt(nrt)
           jetsave = jet(nrt)
c ...   now look thru target data for identical points that get updated
          do nrt2=1,nreg
            if (nzt(nrt2).eq.nztsave.and.nrt2.ne.nrt) then
             if (jbt(nrt2).eq.jbtsave) then
               ierrflag3 = ierrflag3 + 1
               print*, ' > Multiple update on j=1 k=constant grid points
     1:'
               print*, '     line number of original: ',nrt
               print*, '                 zone number: ',nztsave
               print*, '    line number of duplicate: ',nrt2
               print*, '                 zone number: ',nzt(nrt2)
               print*, '     kbt(nrt) = ket(nrt): ',ket(nrt)
               print*, '     duplicate point: jbt: ',jbtsave
             endif
             if (jet(nrt2).eq.jetsave) then
               ierrflag3 = ierrflag3 + 1
               print*, ' > Multiple update on j=end k=constant grid poin
     1ts:'
               print*, '     line number of original: ',nrt
               print*, '                 zone number: ',nztsave
               print*, '    line number of duplicate: ',nrt2
               print*, '                 zone number: ',nzt(nrt2)
               print*, '     kbt(nrt) = ket(nrt): ',ket(nrt)
               print*, '     duplicate point: jet: ',jetsave
             endif
           endif
         enddo

        endif
       enddo

950   continue

c ... for jbt = jet (constant j)
      print*
      print*, ' File bcpzn.dat -- Checking k endpoints, j=constant'
      ierrflag4 = 0

c ... simply flags points that can be updated from multiple base points,
c ...   for k endpoints, j=constant
c DWB: SKIP THIS FOR NOW

      goto 951

      do nrt=1,nreg
        if (jbt(nrt).eq.jet(nrt)) then
          nztsave = nzt(nrt)
          kbtsave = kbt(nrt)
          ketsave = ket(nrt)
c ...   now look thru target data for identical points that get updated
          do nrt2=1,nreg
            if (nzt(nrt2).eq.nztsave.and.nrt2.ne.nrt) then
              if (kbt(nrt2).eq.kbtsave) then
                ierrflag4 = ierrflag4 + 1
                print*, ' > Multiple udpate on k=1 j=constant grid point
     1s:'
                print*, '     line number of original: ',nrt
                print*, '                 zone number: ',nztsave
                print*, '    line number of duplicate: ',nrt2
                print*, '                 zone number: ',nzt(nrt2)
                print*, '     jbt(nrt) = jet(nrt): ',jet(nrt)
                print*, '     duplicate point: kbt: ',kbtsave
              endif
              if (ket(nrt2).eq.ketsave) then
                ierrflag4 = ierrflag4 + 1
                print*, ' > Multiple udpate on k=end j=constant grid poi
     1nts:'
                print*, '     line number of original: ',nrt
                print*, '                 zone number: ',nztsave
                print*, '    line number of duplicate: ',nrt2
                print*, '                 zone number: ',nzt(nrt2)
                print*, '     jbt(nrt) = jet(nrt): ',jet(nrt)
                print*, '     duplicate point: ket: ',ketsave
              endif
            endif
          enddo

        endif
      enddo

951   continue

      print*
      print*, ' Summary for boundary condition checks, file bcpzn.dat:'
      print*, '   (for info only; none of these should be harmful)'
      if (ierrflag1.eq.0) then
       print*, '   -- no multiply-defined target bcs'
      else
       print*, '   >> multiply-defined target bcs found:',ierrflag1
      endif
      if (ierrflag2.eq.0) then
       print*, '   -- no multiply-defined donor bcs'
      else
       print*, '   >> multiply-defined donor bcs found:',ierrflag2
      endif
      if (ierrflag3.eq.0) then
       print*, '   -- no multiply-defined bcs for constant k'
      else
       print*, '   >> multiply-defined bcs for constant k:',ierrflag3
      endif
      if (ierrflag4.eq.0) then
       print*, '   -- no muliply-defined bcs for constant j'
      else
       print*, '   >> multiply-defined bcs for constant j:',ierrflag4
      endif
      print*

      if (ierrflag1.ne.0.or.ierrflag2.ne.0.or.
     1 ierrflag3.ne.0.or.ierrflag4.ne.0) then
      if (DEBUG) then
       print*
       print*, ' ERRORS detected in file bcpzn.dat'
      endif
c       print*, '  stop: stopping execution'
c       stop ' stop: stopping execution due to multiply-defined bcs'
      endif
                   
c++++++++++++++++++++++++++
      endif
c++++++++++++++++++++++++++

c-----
c  Process negative indices
c-----

c++++++++++++++++++++++++++
      if (nodeid.eq.0) then
c++++++++++++++++++++++++++

c get jmax,kmax from each node
c ... first for node 0
      jmaxx(1) = jmax
      kmaxx(1) = kmax
      do i=2,numprocs
        node_in = i-1
        call MPI_RECV(ivarjk,2,MPI_INTEGER,node_in,node_in,
     &   MPI_COMM_WORLD,stat,ierr)
c unpack ivarjk
        jmaxx(i) = ivarjk(1)
        kmaxx(i) = ivarjk(2)
      enddo

c ... check for negative indices
      ierrflag = 0
      do 91 nr=1,nreg
       if (jbt(nr).lt.0) then
        ierrflag = 1
        print*, ' ERROR in patched bc: jbt for bc #',nr,' is negative!'
       endif
       if (jet(nr).lt.0) then
        ierrflag = 1
        print*, ' ERROR in patched bc: jet for bc #',nr,' is negative!'
       endif
       if (kbt(nr).lt.0) then
        ierrflag = 1
        print*, ' ERROR in patched bc: kbt for bc #',nr,' is negative!'
       endif
       if (ket(nr).lt.0) then
        ierrflag = 1
        print*, ' ERROR in patched bc: ket for bc #',nr,' is negative!'
       endif
       if (jbb(nr).lt.0) then
        ierrflag = 1
        print*, ' ERROR in patched bc: jbb for bc #',nr,' is negative!'
       endif
       if (jeb(nr).lt.0) then
        ierrflag = 1
        print*, ' ERROR in patched bc: jeb for bc #',nr,' is negative!'
       endif
       if (kbb(nr).lt.0) then
        ierrflag = 1
        print*, ' ERROR in patched bc: kbb for bc #',nr,' is negative!'
       endif
       if (keb(nr).lt.0) then
        ierrflag = 1
        print*, ' ERROR in patched bc: keb for bc #',nr,' is negative!'
       endif
91     continue

       if (ierrflag.gt.0) then
        print*
        print*, ' stop: negative indices -- stopping in rbc.f/rbcpzn'
        stop ' stop: negative indices -- stopping in rbc.f/rbcpzn'
       endif 

c-----
c  Error checking: indices inbounds
c-----
      do 100 nr=1,nreg
         if(nzt(nr) .gt. numprocs .or.
     &      jbt(nr) .gt. jmaxx(nzt(nr)) .or.
     &      kbt(nr) .gt. kmaxx(nzt(nr)) .or.
     &      jet(nr) .gt. jmaxx(nzt(nr)) .or.
     &      ket(nr) .gt. kmaxx(nzt(nr)) ) then
            write(istdout,110) 2*nr-1
            stop
         elseif (nzb(nr) .gt. numprocs .or.
     &      jbb(nr) .gt. jmaxx(nzb(nr)) .or.
     &      kbb(nr) .gt. kmaxx(nzb(nr)) .or.
     &      jeb(nr) .gt. jmaxx(nzb(nr)) .or.
     &      keb(nr) .gt. kmaxx(nzb(nr)) ) then
            write(istdout,110) 2*nr
            stop
         endif
100   continue

c+++++++++++++++++++++++++++++++++
      else
c+++++++++++++++++++++++++++++++++

c send jmax,kmax from all other nodes to node 0
      ivarjk(1) = jmax
      ivarjk(2) = kmax
      call MPI_SEND(ivarjk,2,MPI_INTEGER,0,nodeid,
     &  MPI_COMM_WORLD,ierr)

c+++++++++++++++++++++++++++++++++
      endif
c+++++++++++++++++++++++++++++++++
      
110   format(/,' ERROR in indices in file bcpzn.dat: zone number ',
     &       /,'       or index is greater than the dimensions of',
     &       /,'       the input grid on input line ',i5,'  Stopping.')
c-----
c  Error checking: indices describe lines with equal number of points
c     and ensure zero increment on constant index
c-----

c++++++++++++++++++++++++++
      if (nodeid.eq.0) then
c++++++++++++++++++++++++++

      do 120 nr=1,nreg
c
c target line
c
         if( jbt(nr) .eq. jet(nr) ) then
            if(kbt(nr) .eq. ket(nr)) then
               ntarpts = 1
            else
               ntarpts = (ket(nr) + 1 - kbt(nr))
            endif
         elseif( kbt(nr) .eq. ket(nr) ) then
            ntarpts = (jet(nr) + 1 - jbt(nr))
         else
            write(istdout,130) 2*nr-1,2*nr
            call exit(1)
         endif
c
c base line
c
         if( jbb(nr) .eq. jeb(nr) ) then
            if(kbb(nr) .eq. keb(nr)) then
               nbsepts = 1
            else
               nbsepts = (keb(nr) + 1 - kbb(nr))
            endif
         elseif( kbb(nr) .eq. keb(nr) ) then
               nbsepts = (jeb(nr) + 1 - jbb(nr))
         else
            write(istdout,130) 2*nr-1,2*nr
            call exit(1)
         endif
c
c check number of points
c
         if(ntarpts .le. 0 .or. nbsepts .le. 0) then
            print*,' ntarpts, nbsepts = ',ntarpts,nbsepts
            write(istdout,140) 2*nr-1,2*nr
            call exit(1)
         elseif( ntarpts .ne. nbsepts ) then
            write(istdout,150) 2*nr-1,2*nr,ntarpts,nbsepts
            call exit(1)
         endif
120   continue
c-----
130   format(/,' ERROR: set of indices input from file bcpzn.dat',
     &       /,'        do not form a single computational line.',
     &       /,'        Error is in input lines ',i5,' and ',i5,
     &       /,'        Stopping.')
135   format(/,' ERROR: in file bcpzn.dat: illegal zero increment',
     &       /,'        found in input indices in input line numbers',
     &       /,'        ',i5,' and ',i5,'  Stopping.')
140   format(/,' ERROR: in file bcpzn.dat: wrong sign on increment',
     &       /,'        found in input indices in input line numbers',
     &       /,'        ',i5,' and ',i5,'  Stopping.')
150   format(/,' ERROR: in file bcpzn.dat: different number of points',
     &       /,'        in the base and target lines.  On the input',
     &       /,'        lines number ',i5,' and ',i5,' there are ',i5,
     &       /,'        target points, and ',i5,' base points.',
     &       /,'        Stopping.')
c-----
c  Write out info
c-----
      if(ntime .eq. 0) then
         write(istdout,300) nreg
      endif
300   format(/' Patched grid zonal interface info read in for',
     &       i5,' regions.')
      do 350 nr=1,nreg
        if(ntime .eq. 0) then
         write(istdout,360) nr,nzt(nr),nzb(nr)
         write(istdout,365) jbt(nr),jet(nr),
     &                      kbt(nr),ket(nr)
         write(istdout,370) jbb(nr),jeb(nr),
     &                      kbb(nr),keb(nr)
        endif
350   continue
360   format('   Patch #',i4,' updates zone ',i4,' using zone ',i4)
365   format('        j-target: ',2i5,' k-target: ',2i5)
370   format('        j-base  : ',2i5,' k-base  : ',2i5)
      close(99)

c++++++++++++++++++++++++++++
      endif
c++++++++++++++++++++++++++++

c++++++++++++++++++++++++++++
      if (nodeid.eq.0) then
c++++++++++++++++++++++++++++
      allocate (ivari(10,nreg))
c pack up variables
      do i = 1,nreg
        ivari(1,i) = nzt(i)
        ivari(2,i) = jbt(i)
        ivari(3,i) = jet(i)
        ivari(4,i) = kbt(i)
        ivari(5,i) = ket(i)
        ivari(6,i) = nzb(i)
        ivari(7,i) = jbb(i)
        ivari(8,i) = jeb(i)
        ivari(9,i) = kbb(i)
        ivari(10,i) = keb(i)
      enddo
c send to other processors
      do inode = 2,numprocs
        node_out = inode - 1
        call MPI_SEND(nreg,1,MPI_INTEGER,node_out,node_out,
     &    MPI_COMM_WORLD,ierr)
        call MPI_SEND(ivari,10*nreg,MPI_INTEGER,node_out,node_out,
     &    MPI_COMM_WORLD,ierr)
      enddo

c++++++++++++++++++++++++++++
      else
c++++++++++++++++++++++++++++

      call MPI_RECV(nreg,1,MPI_INTEGER,0,nodeid,
     &  MPI_COMM_WORLD,stat,ierr)

      allocate (ivari(10,nreg))
      call MPI_RECV(ivari,10*nreg,MPI_INTEGER,0,nodeid,
     &  MPI_COMM_WORLD,stat,ierr)

c unpack variables
      do i = 1,nreg
        nzt(i) = ivari(1,i)
        jbt(i) = ivari(2,i)
        jet(i) = ivari(3,i)
        kbt(i) = ivari(4,i)
        ket(i) = ivari(5,i)
        nzb(i) = ivari(6,i)
        jbb(i) = ivari(7,i)
        jeb(i) = ivari(8,i)
        kbb(i) = ivari(9,i)
        keb(i) = ivari(10,i)
      enddo

c++++++++++++++++++++++++++++
      endif
c++++++++++++++++++++++++++++

      deallocate(ivari)


c sync all nodes
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG) then
c       print*
c       print *, ' stop in rbc.f/rbcpzn'
c       stop 'stop: in rbc.f/rbcpzn'
      endif

c-----
c  End of rbcpzn
c-----
      return
      end
c
c
c************************************************************************
      subroutine testbcmain(bcflag,jmax,kmax)
c************************************************************************
#include "precis.h"
      __REAL bcflag(jmax,kmax)
c
c      do 10 j=jbeg,jend
c      do 10 k=kbeg,kend
      do 10 j=1,jmax
      do 10 k=1,kmax
         bcflag(j,k) = bcflag(j,k) + 1.
10    continue
c-----
c  End of testbcmain
c-----
      return
      end
c
c
c************************************************************************
      subroutine testbcflag(bcflag,jmax,kmax)
c************************************************************************
c  Test for multiple boundary conditions applied to the same point
c     bcflag = 1           for each bcmain.dat  BC point
c            + 100         for each bcpzn.dat   BC point
c------------------------------------------------------------------------
#include "common.f"
      __REAL bcflag(jmax,kmax)
      save ibcerr
c
c      if(nz .eq. 1) ibcerr = 0
      ibcerr = 0
c
      do 70 j=1,jmax
      do 70 k=1,kmax
         if(bcflag(j,k) .ne. 0.0 .and.
     &      bcflag(j,k) .ne. 1.0 .and.
     &      bcflag(j,k) .ne. 100.0 ) then
c
            if(ibcerr .lt. 20) write(istdout,75) ibcerr,nodeid,j,k,bcflag(j,k)
            if(ibcerr .eq. 20) write(istdout,80)
            ibcerr = ibcerr + 1
         endif
70    continue
75    format(i3,'. WARNING: duplicate BC: node,j,k,bcflag = ',3i4,f12.1)
80    format('   Further duplicate BC warnings will not be reported...')
c-----
      if(ibcerr .gt. 0) then
         write(istdout,100) ibcerr
100      format(/,' Duplicate boundary conditions found at ',i6,
     &            ' points.  Applying more than one boundary',/,
     &            ' condition at a single point can cause problems.',/,
     &            ' The value of bcflag is computed from',/,
     &            ' bcflag = 1           for each bcmain.dat  BC point'/
     &            '        + 100         for each bcpzn.dat   BC point'/
     &   )
      endif
c-----
c  End of testbcflag
c-----
      return
      end
