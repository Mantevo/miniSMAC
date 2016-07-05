c
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
c  1. corners
c  2. func_index
c  3. solve
c  4. dateline
c  5. iterout
c
c************************************************************

c*****************************************************************
c
c                       smac2d
c                Sandia miniAero Code 2D
c                     Version 2.0 
c
c         Incompressible Navier-Stokes Flow Solver
c Supplied test case: NACA 4412 airfoil at user-specified angle of attack
c
c Required input files:
c   1. xy.dat_original    default grid file, 1 block only; formatted, plot3d format.
c                         copy one of the following files to xy.dat_original:
c                          > xy.dat_original_119x31 (3,689 grid points)
c                           - included in distribution tarball
c                          > xy.dat_original_5kx5k (25M grid points)
c                           - download, then untar, from:
c                           http://mantevo.org/downloads/Test_case_5kx5k.tar.gz
c                          > xy.dat_original_7kx7k (49M grid points)
c                           - download, then untar, from:
c                           http://mantevo.org/downloads/Test_case_7kx7k.tar.gz
c                         Default name can be changed in smac2d.in.
c   2. bcmain.dat_original default main boundary condition file corresponding to  
c          		  grid file xy.dat_original; text file.
c                         copy one of the following files to bcmain.dat_original
c                          > bcmain.dat_original_119x31
c                           - included in distribution tarball
c                          > bcmain.dat_original_5kx5k
c                           - included in Test_case_5kx5k.tar.gz
c                          > bcmain.dat_original_7kx7k
c                           - included in Test_case_7kx7k.tar.gz  
c                          Default name can be changed in smac2d.in.
c                 
c   All other necessary files are internally generated for the resulting
c   partitioned grid::
c   1. GRIDS2D.p3d
c   2. bcpzn.dat
c   3. bcmain.dat
c
c To compile:
c   1. edit 'Makefile' to make sure all commands are appropriate for
c      the machine on which you wish to run; check directory paths
c      and preprocessor commands for appropriateness; various  options 
c      are listed in the Makefile.
c   2. type 'make' to create the executable 'smac2d'
c
c To run 32-node case interactively on Sandia's redsky:
c   1. allocate 4 nodes, where FYxxxxxxx comes from Sandia's Workload
c      Characterization Tool, or WC Tool:
c       salloc -N4 --time=4:00:00 --account=FYxxxxxxx bash
c   2. compile code:
c       cd <code_directory>
c       make mpi  ! creates executable smac2d
c   3. to run 8 cores per node (for 32 cores total):
c       mpirun -npernode 8 smac2d > smac2d.out

c To run 16-node case with 1 core per node and 8 threads per rank:
c  1. copy 'my_batch_script_mpi_and_threads.bash' to 'my_batch_script'
c  2. edit my_batch_script and set the appropriate variables to run with
c     the above parameters; also check source code directory
c  3. sbatch my_batch_script (to run the code in batch mode)

c To monitor progress of your code when running batch:
c   If you want to monitor progress of your code, run the Python
c   script 'sqcheck.py' as follows:
c       python sqcheck.py
c  
c         
c
c
c   Description: 
c    1. This code solves the incompressible Navier-Stokes
c       equations in two-dimensional generalized coordinates for 
c       steady-state flow.  
c    2. The equations are formulated into a hyperbolic set of PDE using 
c       the method of artificial compressibility.  
c    3. Convective terms are differenced using an upwind biased 
c       flux-difference splitting.  
c    4. The equations are solved using a line relaxation or symmetric
c       Gauss-Seidel solver.
c    5. For parallel processing, a one-block grid is partitioned into
c       a user-specified number of subgrids. The code takes care of 
c       generating one grid cell overlap with point-to-point matching,
c       as well as ensuring the partitioning process minimizes the
c       surface-area-to-volume ratio of the subgrids based on the number
c       of subgrids chosen by the user..
c    6. The number of subgrids for the solution space is defined by the
c       number of MPI ranks assigned to the run. For example, if as shown
c       above, salloc allocates 4 nodes, and there are 8 cores per node
c       as indicated in the mpirun command, then the original 1-block grid
c       will be partitioned into 32 subgrids.
c
c---------------------------------------------------------------------------
c   This is an open-source code.
c
c  If you have comments or questions regarding this code, contact:
c       Daniel W. Barnette
c       Sandia National Laboratories
c       MS 1319
c       Albuquerque, NM 87185
c       dwbarne@sandia.gov 
c
c*****************************************************************
      program main
c*****************************************************************

#include "common.f"

c********************************************************
#include "mpif.h"
#include "mpi_params.f"

#ifdef _OPENMP
c for openmp runtime calls
#include "omp_lib.h"
#endif

      INTEGER nodeid,numprocs,ierr,stat(MPI_STATUS_SIZE)
c********************************************************

      INTEGER ivar(2),jmax,kmax
      character*100 chline
      character*110 filename
      character*60 title1,geometry1,geometry_yaml
      logical iblank,form,mzone
      logical DEBUG,DEBUG1,PRINT_CORNERS
      character*20 date_time
      character*20 my_username
      character*20 hostname
      character*10 date_calendar
      character*8  time_clock
      character*5 YAML_EXTENSION 
      character*110 filename_yaml
      character*5 numprocs_txt
      character*3 numthreads_txt
      character*6 jmax_o_txt,kmax_o_txt
      character*6 ntmax_txt
      character*30 impsch_txt
      character*6 fp_precision
      integer yaml_unit
      data yaml_unit/500/
      data YAML_EXTENSION/'.yaml'/
c Timings
      double precision startTimeInitialize,startTimeSendSubgrids,
     &  startTimeSolver
      double precision endTimeInitialize,endTimeSendSubgrids,
     &  endTimeSolver
      double precision timeInitialize,timeSendSubgrids,timeSolver
c      real timeInitialize,timeSendSubgrids,timeSolver
c allocatable arrays
      __REAL, allocatable:: xorig(:),yorig(:)
      __REAL, allocatable:: x(:,:),y(:,:)
      __REAL, allocatable:: j_subgrid(:),k_subgrid(:)
      __REAL, allocatable:: xgrid(:,:),ygrid(:,:)
      integer, allocatable:: jgrid(:),kgrid(:)
      integer, allocatable:: jindex_global(:),kindex_global(:)

#ifdef _OPENMP
      integer function omp_get_num_threads()
      integer function omp_get_thread_limit() 
#endif

      DEBUG = .false.
c      DEBUG = .true.
      DEBUG1 = .false.
c      DEBUG1 = .true.
      PRINT_CORNERS = .false.
c      PRINT_CORNERS = .true.

c-----
c  MPI initialization 
c-----

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, nodeid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

c Timings:
c ... start initialization timer
c+++++++++++++++++++++++++++++++++++
      if (nodeid.eq.0) then
c+++++++++++++++++++++++++++++++++++
        startTimeInitialize = MPI_WTIME()
c++++++++++++++++++++++++++++
      endif
c++++++++++++++++++++++++++++

c allocate memory for global values
      allocate (jindex_global(numprocs),kindex_global(numprocs))

c+++++++++++++++++++++++++++++++++++
      if (nodeid.eq.0) then
c+++++++++++++++++++++++++++++++++++

      if (DEBUG) then
        print*
        print*, '>>> MODULE: smac2d.f, SUBROUTINE: main'
        print*
      endif

      if (DEBUG) then
        print*
        print*, 'from "common.f" parameters:'
        print*, '  nzne = ',nzne
        print*, '  ibcmax = ',ibcmax
        print*, '  jkmax = ',jkmax
        print*, '  istdout = ',istdout
        print*
      endif

c-----
c  Standard input and output
c-----
      istdout = 6

c-----
c  Output version number, time, date
c-----

      call dateline(chline,date_calendar,time_clock,date_time,hostname)

c do NOT use ':' to separate hrs-mins-secs in 'time_clock'
c   since YAML interprets colons as something else
      time_clock = time_clock(1:2)//''//time_clock(4:5)//''//
     1 time_clock(7:8)

      if (DEBUG) then
        print*, ' '
        print*, '** In program main: from nodeid = ',nodeid
        print*
        print*, 'date_time = ',date_time
        print*, 'hostname = ',hostname
        print *, ' '
      endif

c title
      write(istdout,1) chline
1     format('********************************************************',
     &      /'    Sandia MiniAero Code SMAC2D-UP Version 1.0',/,10x,a,/,
     &       '********************************************************')
c copyright notice
      write(istdout,2)
2     format(/
     &'                 Copyright 2013 Sandia Corporation          ',/,
     &'                                                            ',/,
     &'Under the terms of Contract DE-AC04-94AL85000 with Sandia   ',/,
     &'Corporation, the U.S. Government retains certain rights in  ',/,
     &'this software.                                              ',/,
     &'                                                            ',/,
     &'Redistribution and use in source and binary forms, with or  ',/,
     &'without modification, are permitted provided that the       ',/,
     &'following conditions are met:                               ',/,
     &'                                                            ',/,
     &' 1. Redistributions of source code must retain the above    ',/,
     &'    copyright notice, this list of conditions and the       ',/,
     &'    following disclaimer.                                   ',/,
     &'                                                            ',/,
     &' 2. Redistributions in binary form must reproduce the above ',/,
     &'    copyright notice, this list of conditions and the       ',/,
     &'    following disclaimer in the documentation and/or other  ',/,
     &'    materials provided with the distribution.               ',/,
     &'                                                            ',/,
     &' 3. Neither the name of the Corporation nor the names of th ',/,
     &'    contributors may be used to endorse or promote products ',/,
     &'    derived from this software without specific prior       ',/,
     &'    written permission.                                     ',/,
     &'                                                            ',/,
     &'THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ',/,
     &'ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT       ',/,
     &'LIMITED TO THE IMPLIED WARRANTIES OF MERCHANTABILITY AND    ',/,
     &'FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT',/,
     &'SHALL SANDIA CORPORATION OR THE CONTRIBUTORS BE LIABLE FOR  ',/,
     &'ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR    ',/,
     &'CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       ',/,
     &'PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,   ',/,
     &'DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED  ',/,
     &'AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT ',/,
     &'LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)      ',/,
     &'ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ',/,
     &'ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                  '
     &)

c print number of MPI ranks or processes in this run
      print *
      write(*,85) numprocs
85    format('Total no. of MPI ranks for this run = ',i6)
      numthreads = 0

#ifdef _OPENMP
!$omp parallel shared(numthreads)
!$omp single
c print number of OPENMP threads per MPI rank for this run
      numthreads = omp_get_num_threads()
c      print*, 'numthreads per MPI rank within $omp parallel = ',
c     1 numthreads
c      print*, ' ... thread limit within $omp parallel = ',
c     1   omp_get_thread_limit()
      print*, '--------------------------'
!$omp end single
!$omp end parallel
#else
      numthreads = 0
#endif

      write(*,86) numthreads
86    format('No. of OPENMP threads per MPI rank for this run = ',i6)

#ifdef _OPENMP
      print*, '   ... thread limit = ',omp_get_thread_limit()
#endif

c      print*, 'stopped after MPI ranks and OPENMP threads/rank print'
c      call MPI_ABORT(MPI_COMM_WORLD,1)
c      stop 'stopped after MPI ranks and OPENMP threads/rank print'


c********************************************************
c end nodeid.eq.0
      endif
c********************************************************

c********************************************************
c sync all processes
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
c check node availability
      if (DEBUG) then
      if(nodeid.eq.0.and.numprocs.gt.1) then
        print*
        print*, 'Checking MPI node availability:'
        print*
        print*, ' ...numprocs = ',numprocs,', nodeid = ',nodeid
        do node = 1,numprocs-1
          call MPI_RECV(idum,1,MPI_INTEGER,node,MPI_ANY_TAG,
     &         MPI_COMM_WORLD,stat,ierr) 
        enddo
        print*          
        print*, ' -- ALL NODES HAVE RESPONDED --'
        print*

      else
        call MPI_SEND(1,1,MPI_INTEGER,0,0,
     &      MPI_COMM_WORLD,ierr )
      endif
      endif

c********************************************************
c-----
c  Initialize fixed variables in all zones
c-----

      call initia(geometry1,title1)

      if (DEBUG1.and.nodeid.eq.0) then
       print*
       print*, ' AFTER call initia in smac2d.f/smac2d'
       print*, '    beta = ',beta
       print*, '    dtau = ',dtau
       print*, '    dt   = ',dt
       print*, '    time = ',time
       print*, '    kadd = ',kadd
c      call MPI_ABORT(MPI_COMM_WORLD,1)
c       stop 'stop: AFTER call initia in smac2d.f/smac2d'
      endif


c-----
c  Read in original grid file header, determine grid file format, 

      iunit = 10
c++++++++++++++++++++++++++++++++
      if (nodeid.eq.0) then
c++++++++++++++++++++++++++++++++

       if (DEBUG) then
        print*
        print*, 'node 0 calling: rgridh in smac2d.f/main'
        print*
       endif

c read grid header; iselect=0 means don't send grid header to 
c  other nodes; usually don't want to send since grid needs
c  partitioning before sending to other nodes. Partitioning is 
c  a serial process and handled by node 0 alone.
      iselect=0
      call rgridh(iunit,jmax,kmax,form,mzone,iblank,
     1  nzne,nzone,jkmax,istdout,largest_dim,
     2  iselect)

c save jmax,kmax for original grid
      jmax_o = jmax
      kmax_o = kmax

c allocate arrays needed for breakup; 
c  include kadd for constructing larger grid size than the original
      jdim = jmax
      kdim = kmax + kadd 
      jk = jdim*kdim
      itotalOrigPoints = jdim*kdim

c allocate memory for dynamic grid points
      print*, ' Node 0: allocating memory for dynamic grid points x,y'
      print*, '   using dimensions j,k+kadd = ',jdim,kdim
      print*
c      allocate(x(jk,numprocs),y(jk,numprocs))

c allocate memory for static grid points (original grid)
      print*, ' allocating memory for static grid points xorig, etc.'
      allocate(xorig(jk),yorig(jk))

c allocate memory for local j,k,l values for each subgrid
      print*, ' allocating memory for local jk values j_subgrid, etc.'
      allocate(j_subgrid(numprocs),k_subgrid(numprocs))

c read original grid, using original dimensions 
      kindex = 1
      if (DEBUG) then
       print*, ' Node 0: calling readgo in smac2d.f/main'
       print*
      endif
      call readgo(iunit,xorig,yorig,jdim,kdim,kadd,mzone)
c
c close original grid file
      close(iunit)

c if kadd is not zero, add grid points to outer boundary in the radial
c  (k) direction. Number of grid points in j direction is not changed.
      if (kadd.ge.1) then
       call extend_grid_k(xorig,yorig,jdim,kdim,kadd)
       jmax_o = jdim
       kmax_o = kdim
      endif 
c++++++++++++++++++++++++++++++++
      endif
c++++++++++++++++++++++++++++++++

c++++++++++++++++++++++++++++++++
      if (nodeid.eq.0.and.numprocs.gt.1) then
c++++++++++++++++++++++++++++++++
        nghost = 1

        if (PRINT_CORNERS) then
         call corners(nodeid,jdim,kdim,jk,
     1    xorig,yorig)
        endif

        print*
        print*, ' Grid is to be broken up into ',numprocs,
     &    ' subgrids'
        print*, '   with ',nghost,' layer of ghost cells.'
        print*
c++++++++++++++++++++++++++++++++
      endif
c++++++++++++++++++++++++++++++++
      print*
c      print*, ' Node ',nodeid

c+++++++++++++++++++++++++++
      if (nodeid.eq.0.and.numprocs.gt.1) then
c+++++++++++++++++++++++++++
       if (DEBUG) then
        print*
        print*, ' ***** Node 0: entering breakup routines in smac2d.f/ma
     &in *****'
        print*
       endif

      if (PRINT_CORNERS) then
       call corners(nodeid,jdim,kdim,jk,
     1  xorig(1),yorig(1))
      endif

c write corner points of grid 
      print*
      print*, '+++++++++++++++++++++++++++++++++++++++++++++++'
      print*, ' calling: cube in smac2d.f/main'
      print*, '   jdim,kdim = ',jdim,kdim
      print*
      call cube(xorig(1),yorig(1),
     & jdim,kdim)

      if (PRINT_CORNERS) then
       call corners(nodeid,jdim,kdim,jk,
     1  xorig(1),yorig(1))
      endif
c
c partition grid into 'numprocs' partitions; output is a static grid
c   load balanced, ghost celled, and ready for parallel processing
       print*
       print*, '+++++++++++++++++++++++++++++++++++++++++++++++'
       print*, ' Node 0: calling partishn in smac2d.f/main'
       print*
      call partishn(form,mzone,xorig(1),yorig(1),
     & jdim,kdim,idimj,idimk,jk,nodeid,numprocs,
     & j_subgrid(1),k_subgrid(1),
     & jindex_global(1),kindex_global(1),ioverlap
     & )

c print grid corners to see if all is still ok
      if (PRINT_CORNERS) then
       call corners(nodeid,jdim,kdim,jk,
     1   xorig(1),yorig(1))
      endif

      print*
      print*, ' Node 0: finished calling partishn in smac2d.f/main'
      print*, '+++++++++++++++++++++++++++++++++++++++++++++++++++'
 
c+++++++++++++++++++++++++++
      endif
c+++++++++++++++++++++++++++

c send globals to all nodes
c+++++++++++++++++++++++++++
      if (numprocs.gt.1) then
c+++++++++++++++++++++++++++
       call MPI_BCAST(jindex_global,numprocs,MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr)
       call MPI_BCAST(kindex_global,numprocs,MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr)
       call MPI_BCAST(idimj,1,MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr)
       call MPI_BCAST(idimk,1,MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ioverlap,1,MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr)
c make sure all nodes have global values before continuing
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c+++++++++++++++++++++++++++
      else
c+++++++++++++++++++++++++++
       idimj = 1
       idimk = 1
c write out grid to GRIDS2D.p3d
       call write_grid_for_1_processor(jmax_o,kmax_o,xorig,yorig)
       if (DEBUG) then
        print*, ' stop: after call write_grid_for_1_processor'
        call MPI_ABORT(MPI_COMM_WORLD,1)
        stop ' stop: after call write_grid_for_1_processor'
       endif
c+++++++++++++++++++++++++++
      endif
c+++++++++++++++++++++++++++
      
c-----
c  Initialize fixed variables in all zones
c-----

c call inizone with all nodes to define zone-specific values
      call inizone

c-----
c  User-definable geometry dependent input
c  .... All nodes define alpha (angle of attack), etc.
c-----
c read geometry-dependent variables
      call geominit


c+++++++++++++++++++++++++++
      if (nodeid.eq.0) then
c+++++++++++++++++++++++++++

c 'getlog' is an intrinsic function to obtain username. not portable,
c   so it doesn't always work
      call getlog(my_username)
      if(my_username.eq.'') then
        my_username = 'unavailable'
      endif

      if (DEBUG) then
        print*, ' '
        print*, ' date_time = ',date_time
        print*, ' title = ',title1
        print*, ' geometry = ',geometry1
        print*, ' hostname = ',hostname
        print*, ' my_username = ',trim(my_username)
        print*, ' '
      endif

c+++++++++++++++++++++++++++
      endif
c+++++++++++++++++++++++++++



c++++++++++++++++++++++++++++++
      if (nodeid.eq.0) then
c++++++++++++++++++++++++++++++

c open yaml output file here, to capture output in other subroutines
c change grid dimensions to text
      write(jmax_o_txt,'(i6)') jmax_o
      jmax_o_txt = trim(adjustl(jmax_o_txt))
      write(kmax_o_txt,'(i6)') kmax_o
      kmax_o_txt = trim(adjustl(kmax_o_txt))

c change number of zones to text
      write(numprocs_txt,'(i5)') numprocs 
      numprocs_txt = trim(adjustl(numprocs_txt))
      write(numthreads_txt,'(i3)') numthreads
      numthreads_txt = trim(adjustl(numthreads_txt))
      write(ntmax_txt,'(i6)') ntmax
      ntmax_txt = trim(adjustl(ntmax_txt))

      if (DEBUG) then
        print*, ' jmax_o_txt = ',jmax_o_txt
        print*, ' kmax_o_txt = ',kmax_o_txt
        print*, ' numprocs = ',numprocs
        print*, ' numprocs_txt = ',numprocs_txt
        print*, ' numthreads_txt = ',numthreads_txt
        print*, ' ntmax = ',ntmax
        print*, ' ntmax_txt = ',ntmax_txt
      endif

c substitute spaces with underlines in the 'geometry1' variable
      geometry_yaml = trim(geometry1)
c      print*
c      print*, geomtry_yaml,' is the old geometry string'
c      print*, ' number of characters = ',len_trim(geometry_yaml)
      do icol=1,len_trim(geometry_yaml)
       if (geometry_yaml(icol:icol).eq.' ') geometry_yaml(icol:icol)='_'
      enddo
c      pcrint*, geometry_yaml,' is the new geometry string'
c      call MPI_ABORT(MPI_COMM_WORLD,1)
c      stop 'stop: printed geometry_yaml string'
      print*, 'geometry1 = ',trim(geometry1)
      print*, 'hostname = ',trim(hostname)
      print*, 'numprocs_txt = ',trim(numprocs_txt)
      print*, 'numthreads_txt = ',trim(numthreads_txt)
      print*, 'ntmax_txt = ',trim(ntmax_txt)
      print*, 'date_time = ',date_time
#ifdef D_PRECISION
      fp_precision = 'double'
#else
      fp_precision = 'single'
#endif
      print*, 'fp_precision = ',fp_precision
      print*, 'YAML_EXTENSION = ',YAML_EXTENSION

c trim trailing blanks from inputs to filename
c      filename_yaml = 'smac2d'//'_'//trim(geometry1)//'_'//
      filename_yaml = 'smac2d'//'_'//trim(geometry_yaml)//'_'//
     1 trim(hostname)//'_'//trim(jmax_o_txt)//'x'//trim(kmax_o_txt)//
     2 '_'//trim(numprocs_txt)//'ranks_'//trim(numthreads_txt)//
     3 'threads_'//trim(ntmax_txt)//'iters_'//fp_precision//'_'//
     4 date_time//YAML_EXTENSION
      print*
      print*, 'Write yaml-formatted output to file'
      print*, 
     1'   string length:',len_trim(filename_yaml),'/110 characters'
      print*, filename_yaml 
      print*, ' '

c      call MPI_ABORT(MPI_COMM_WORLD,1)
c      stop 'stop: after writing filename_yaml to output'

      open(yaml_unit,file=filename_yaml,form='formatted',status='new')

      write(yaml_unit,*) 'user: ',trim(my_username)
      write(yaml_unit,*) 'code: smac2d'
      write(yaml_unit,*) 'host: ',hostname
      write(yaml_unit,*) 'title: ',trim(title1)
      write(yaml_unit,*) 'geometry: ',trim(geometry1)
      write(yaml_unit,*) 'date_yyyy_mm_dd: ',date_calendar
      write(yaml_unit,*) 'time_of_day: ',time_clock
      write(yaml_unit,*) 'num_grids: ',numprocs
      write(yaml_unit,*) 'mpi_ranks: ',numprocs
      write(yaml_unit,*) 'num_threads_per_node: ',numthreads
      if (numthreads.ne.0) then
        write(yaml_unit,*) 'total_processes: ',numprocs*numthreads
      else
        write(yaml_unit,*) 'total_processes: ',numprocs
      endif
      write(yaml_unit,*) 'grid_file: ',trim(xygrid_filename)
      write(yaml_unit,*) 'largest_dim: ',largest_dim
      write(yaml_unit,*) 'j-dim: ',jdim
      write(yaml_unit,*) 'k-dim: ',kdim
      write(yaml_unit,*) 'total_points: ',itotalOrigPoints
      write(yaml_unit,*) 'beta: ',beta
#ifdef D_PRECISION
      write(yaml_unit,*) 'fp_precision: double'
#else
      write(yaml_unit,*) 'fp_precision: single'
#endif

      write(yaml_unit,900) dtau
900   format(' dtau: ',1pe12.4)
      write(yaml_unit,*) 'reynolds_number: ',reynum
      write(yaml_unit,*) 'angle_of_attack_degs: ',alpha
      if (impsch.eq.1) then
        impsch_txt = 'line relaxation'
      elseif (impsch.eq.2) then
        impsch_txt = 'lusgs factorization'
      else
        impsch_txt = 'unknown'
      endif
      write(yaml_unit,*) 'implicit_scheme: ',trim(impsch_txt)

c++++++++++++++++++++++++++++++
      endif
c++++++++++++++++++++++++++++++

c+++++++++++++++++++++++++++++++
      if (numprocs.gt.1) then
c+++++++++++++++++++++++++++++++
       if (nodeid.eq.0) then
c+++++++++++++++++++++++++++++++

c deallocate grid variables no longer needed
c       deallocate (x,y)
       deallocate (xorig,yorig)
       deallocate (j_subgrid,k_subgrid)

c read in new grid that's been divided into subgrids
      filename='GRIDS2D.p3d'
      open(22,file=trim(filename),iostat=ierr_unit22,form='formatted',
     & status='unknown')
      if(ierr_unit22.ne.0) then
        write(6,*) ' '
        write(6,*) ' ERROR: cannot open GRIDS2D.p3d file'
        write(6,*) '  Program stopped: problem with opening grid file GR
     &IDS2D.p3d'
        write(6,*) ' '
        call MPI_ABORT(MPI_COMM_WORLD,1)
        stop 'stopped: problem with opening grid file GRIDS2D.p3d'
      endif

c read number of sub-grids
      read(22,*) num_grids

      allocate (jgrid(num_grids),kgrid(num_grids))

      do 800 i=1,num_grids
       read(22,*) jgrid(i), kgrid(i)
800   continue

c find max grid values
      jgridmax = maxval(jgrid)
      kgridmax = maxval(kgrid)

c check jkmax value in common.f to see if it is big enough
      if (jgridmax.gt.jkmax.or.kgridmax.gt.jkmax) then
       if (jgridmax.gt.jkmax) then
        print*
        print*, ' ERROR: jgridmax > jkmax'
        print*, '   jgridmax =',jgridmax,', jkmax =',jkmax
        print*, ' jgridmax must be less than or equal to jkmax'
       elseif (kgridmax.gt.jkmax) then
        print*
        print*, ' ERROR: kgridmax > jkmax'
        print*, '   kgridmax =',kgridmax,', jkmax =',jkmax
        print*, ' kgridmax must be less than or equal to jkmax'
       endif
       print*
       print*, ' Increase the value for jkmax in common.f and rerun.' 
       print*
       print*, ' Program is terminating.'
       call exit(0)
      endif

c for node 0:
      jmax = jgrid(1)
      kmax = kgrid(1)

c allocate grid variables for Node 0 that depend on LOCAL jmax,kmax values
      allocate (x(jmax,kmax),y(jmax,kmax))

c send jgridmax,kgridmax to all other nodes and allocate memory for
c  xgrid, ygrid
c       ivar(1) = jgridmax
c       ivar(2) = kgridmax
       do 755 i = 2,num_grids
         node_out = i-1
         ivar(1) = jgrid(i)
         ivar(2) = kgrid(i)
         call MPI_SEND(ivar,2,MPI_INTEGER,node_out,node_out,
     1    MPI_COMM_WORLD,ierr)
755    continue

c+++++++++++++++++++++++++++++++
      else
c+++++++++++++++++++++++++++++++

      call MPI_RECV(ivar,2,MPI_INTEGER,0,nodeid,
     1 MPI_COMM_WORLD,stat,ierr)
      jmax = ivar(1)
      kmax = ivar(2) 
      allocate (x(jmax,kmax),y(jmax,kmax))

c+++++++++++++++++++++++++++++++
      endif
c+++++++++++++++++++++++++++++++


c+++++++++++++++++++++++++++++++
      else     ! numprocs = 1 
c+++++++++++++++++++++++++++++++
c allocate grid variables for Node 0 that depend on LOCAL jmax,kmax values
       allocate(jgrid(1),kgrid(1))
       allocate (x(jmax*kmax,1),y(jmax*kmax,1))

c!$omp parallel do shared(jk)
       do 565 j=1,jk
        x(j,1) = xorig(j)
        y(j,1) = yorig(j)
565    continue
c!$omp end parallel do

       jgrid(1) = jmax_o
       kgrid(1) = kmax_o
       jmax = jmax_o
       kmax = kmax_o
       deallocate(xorig,yorig)

       ioverlap = 1
       jindex_global = 1
       kindex_global = 1

c+++++++++++++++++++++++++++++++
      endif
c+++++++++++++++++++++++++++++++

c Timings
c ... end initialization timer
c ... start subgrid timer
c+++++++++++++++++++++++++++++++
      if (nodeid.eq.0) then
c+++++++++++++++++++++++++++++++
        endTimeInitialize = MPI_WTIME()
        startTimeSendSubgrids = MPI_WTIME()
c+++++++++++++++++++++++++++++++
      endif
c+++++++++++++++++++++++++++++++

c determine x and y for subgrids and send to proper mpi ranks
c+++++++++++++++++++++++++++++++
      if (numprocs.gt.1) then
c+++++++++++++++++++++++++++++++
      if (nodeid.eq.0) then
c+++++++++++++++++++++++++++++++
       print*
       print*, ' num_grids = ',num_grids
       print*
       print*, '          #       jgrid      kgrid'
       do 805 i=1,num_grids
       print*, i,jgrid(i),kgrid(i)
805    continue
       print*
       print*, ' max jgrid = ',jgridmax
       print*, ' max kgrid = ',kgridmax
       print*

c xgrid, ygrid are temporary arrays 
       do 810 i=1,num_grids
c allocate xgrid, ygrid
        allocate (xgrid(jgrid(i),kgrid(i)),ygrid(jgrid(i),kgrid(i)))
        print*
        print*, ' Reading subgrid number: ',i
        print*, '   jgrid(i) = ',jgrid(i)
        print*, '   kgrid(i) = ',kgrid(i)
c read from grid file
        read(22,*) ((xgrid(j,k),j=1,jgrid(i)),k=1,kgrid(i)),
     1             ((ygrid(j,k),j=1,jgrid(i)),k=1,kgrid(i))
        print*,'    Subgrid read: successful'

c DEBUG1 check
        if (DEBUG1.and.i.eq.2) then
         print*
         print*, ' In smac2d.f/main, just after read(22,*):'
         print*, '  i,jgrid(i),kgrid(i) = ',i,jgrid(i),kgrid(i)
         print*, ' nodeid_dest  j   k    xgrid(jk)     ygrid(jk)'
         do k=15,17
         do j=30,34
          write(*,555) i-1,j,k,xgrid(j,k),ygrid(j,k)
         enddo
         enddo
        endif

        if (i.eq.1) then
c for node 0
c          jmax = jgrid(i)
c          kmax = kgrid(i)
          print*, ' For nodeid ',nodeid,': jmax,kmax = ',jmax,kmax
c assign x and y for node 0 subgrid
c          allocate(x(jmax,kmax),y(jmax,kmax))
c          do 812 k=1,kgrid(i)
c          do 812 j=1,jgrid(i)

c!$omp parallel do shared(kmax,jmax)
          do 812 k=1,kmax
          do 813 j=1,jmax
            x(j,k) = xgrid(j,k)
            y(j,k) = ygrid(j,k)
813       continue
812       continue
c!$omp end parallel do

c deallocate xgrid, ygrid
          deallocate(xgrid,ygrid)
          print*, ' Node 0: done - xgrid, ygrid assigned'
          goto 810
        endif

c pack integers; no need to pack reals
c        ivar(1) = jgrid(i)
c        ivar(2) = kgrid(i)
c define node to which data will be sent
        node_out = i-1
c DEBUG1 CHECK
        if (DEBUG1.and.i.eq.2) then
         print*
         print*, ' In smac2d.f/main -- from node 0 -- SENDs:'
         print*, ' i,jgrid(i),kgrid(i) =',i,jgrid(i),kgrid(i)
         print*, ' nodeid_dest  j   k      xgrid(jk)   ygrid(jk)'
         do k=15,17
         do j=30,34
          write(*,555) node_out,j,k,xgrid(j,k),ygrid(j,k)
555       format(5x,i3,4x,2i3,1p2e15.7)
         enddo
         enddo
c         print*, ' stop: in smac2d.f/main 111'
c         call MPI_ABORT(MPI_COMM_WORLD,1)
c         stop 'stop: in smac2d.f/main 111'
        endif

#ifdef D_PRECISION
c send x grid
        call MPI_SEND(xgrid,jgrid(i)*kgrid(i),MPI_REAL8,node_out,
     1    node_out+1,MPI_COMM_WORLD,ierr)
c send y grid
        call MPI_SEND(ygrid,jgrid(i)*kgrid(i),MPI_REAL8,node_out,
     1    node_out+2,MPI_COMM_WORLD,ierr)
#else
c send x grid
        call MPI_SEND(xgrid,jgrid(i)*kgrid(i),MPI_REAL,node_out,
     1    node_out+1,MPI_COMM_WORLD,ierr)
c send y grid
        call MPI_SEND(ygrid,jgrid(i)*kgrid(i),MPI_REAL,node_out,
     1    node_out+2,MPI_COMM_WORLD,ierr)
#endif

      deallocate(xgrid,ygrid)

810     continue

c******************************       
      else
c******************************       

#ifdef D_PRECISION
c      call MPI_RECV(x,jmax*kmax,MPI_REAL8,0,nodeid+1,
      call MPI_RECV(x,jmax*kmax,MPI_REAL8,0,MPI_ANY_TAG,
     1  MPI_COMM_WORLD,stat,ierr)
c      call MPI_RECV(y,jmax*kmax,MPI_REAL8,0,nodeid+2,
      call MPI_RECV(y,jmax*kmax,MPI_REAL8,0,MPI_ANY_TAG,
     1  MPI_COMM_WORLD,stat,ierr)
#else
c      call MPI_RECV(x,jmax*kmax,MPI_REAL,0,nodeid+1,
      call MPI_RECV(x,jmax*kmax,MPI_REAL,0,MPI_ANY_TAG,
     1  MPI_COMM_WORLD,stat,ierr)
c      call MPI_RECV(y,jmax*kmax,MPI_REAL,0,nodeid+2,
      call MPI_RECV(y,jmax*kmax,MPI_REAL,0,MPI_ANY_TAG,
     1  MPI_COMM_WORLD,stat,ierr)
#endif

      if (DEBUG1.and.nodeid.eq.1) then
       print*
       print*, ' In smac2d.f/main -- RECVs:'
       print*, '   jmax,kmax = ',jmax,kmax
       print*, '   jgridmax,kgridmax = ',jgridmax,kgridmax
       print*, '  nodeid   j   k      x(j,k)        y(jk)'
c       print*, '  nodeid   j   k    xgrid(j,k)    ygrid(j,k)'
c       print*, '                    xgrid(j,k)    ygrid(j,k)'
       do k=15,17
       do j=30,34
        write(*,555) nodeid,j,k,x(j,k),y(j,k)
c        write(*,556) xgrid(j,k),ygrid(j,k)
556     format(16x,1p2e15.7)
       enddo
       enddo
      endif

c******************************       
      endif
c******************************
      endif
c******************************

c Timings:
c ... end subgrid timer
c ... start solver timer
c******************************
      if (nodeid.eq.0) then
c******************************
        endTimeSendSubgrids = MPI_Wtime()
        startTimeSolver = MPI_WTIME()
c******************************
      endif
c******************************

c steady state solution only

      iunsdim = 0

c for each subgrid pass following values: jmax,kmax,x,y

      call solve(jmax,kmax,x,y,yaml_unit,startTimeSolver,
     & idimj,idimk,
     & jindex_global,kindex_global,ioverlap
     & )

c******************************
      if (nodeid.eq.0) then
c******************************
c Timings
c ... end solver
        endTimeSolver = MPI_WTIME()
c ... time deltas
        timeInitialize = 
     &    endTimeInitialize - startTimeInitialize
        timeSendSubgrids = 
     &    endTimeSendSubgrids - startTimeSendSubgrids
        timeSolver = 
     &    endTimeSolver - startTimeSolver
        print*
        write(*,1000) timeInitialize
1000    format(' TIMER: initialization (sec): ',1pe12.5)
        write(*,1005) timeSendSubgrids
1005    format(' TIMER: send subgrids (sec): ',1pe12.5)
        write(*,1010) timeSolver
1010    format(' TIMER: final solver time (sec) = ',1pe12.5)
        write(*,1012) timeSolver/60.0
1012    format(' TIMER: final solver time (min) = ',1pe12.5)
        write(*,1014) timeSolver/3600.0
1014    format(' TIMER: final solver time (hr) = ',1pe12.5)
        print*
c write to yaml
        cputiter = timeSolver/float(ntmax)
        cpuitpt = cputiter/float(itotalOrigPoints)
        write(yaml_unit,2000) cputiter
2000    format(' sec_per_iteration: ',1pe12.5)
        write(yaml_unit,2010) cpuitpt*1.0e06
2010    format(' microsec_per_iter_per_point: ',1pe12.5)
        write(yaml_unit,*) 'number_of_iterations: ',ntmax
        write(yaml_unit,*) 'alpha_deg: ',alpha
        write(yaml_unit,2030) timeInitialize
2030    format(' initialize_sec: ',1pe12.5)
        write(yaml_unit,2040) timeSendSubgrids
2040    format(' send_subgrids_sec: ',1pe12.5)
        write(yaml_unit,2050) timeSolver
2050    format(' solve_sec: ',1pe12.5)
        write(yaml_unit,2060) (timeInitialize + timeSendSubgrids 
     &    + timeSolver)
2060    format(' total_time_sec: ',1pe12.5)
        write(yaml_unit,2070) filename_yaml
2070    format(' yaml_filename: ',a)
c inform user re yaml file
        print*, 'Yaml-formatted output written to file:',
     &   filename_yaml

c******************************
      endif
c******************************

      close( istdout )

c-----
c  MPI finish
c-----
c sync all processes and end
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
    
c-----
      
c  End of main
c-----
      end

c
c***********************************************************************`
      subroutine corners(nodeid,jdim,kdim,jk,
     1 xorig,yorig)
c
c     purpose:
c       print corners of the grid
c
      __REAL xorig(jk),yorig(jk)

      print*
      print*, '++++++++ In sub. CORNERS ++++++++++++++'
      print*, ' Node ',nodeid
      print*, ' jdim,kdim = ',jdim,kdim
      print*, ' jk (jdim*kdim) = ',jk

      call func_index(1,1,jdim,kdim,index1)
      print*, ' 1. xy(1,1): ',xorig(index1),yorig(index1)

      call func_index(jdim,1,jdim,kdim,index1)
      print*, ' 2. xy(jmx,1): ',xorig(index1),yorig(index1)

      call func_index(1,kdim,jdim,kdim,index1)
      print*, ' 3. xy(1,kmx): ',xorig(index1),yorig(index1)

      call func_index(jdim,kdim,jdim,kdim,index1)
      print*, ' 4. xy(jmx,kmx): ',xorig(index1),yorig(index1)

      print*
      print*, '+++++++++++++++++++++++++++++++++++++++++++++++'

      return
      end
c
c
c***********************************************************************`
      subroutine func_index(j,k,jdim,kdim,index_return)
c
c purpose: calculates the index for a one-dimensional vector representing
c          a 3-D array
c          Equation assumes column-major storage (Fortran)
      index_return = (k-1)*jdim + j

      return
      end
c
c
c***********************************************************************
      subroutine solve(jmax,kmax,x,y,yaml_unit,startTimeSolver,
     & idimj,idimk,
     & jindex_global,kindex_global,
     & ioverlap
     & )
c***********************************************************************
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      INTEGER nodeid,numprocs,ierr,stat(MPI_STATUS_SIZE)

      __REAL q(jmax,kmax,3),x(jmax,kmax),y(jmax,kmax),
     1 rtxy(jmax,kmax,2,3),dj(jmax,kmax),turvar(jmax,kmax,2),
     2 dq(jmax,kmax,3),s(jmax,kmax,3),vnut(jmax,kmax)
      dimension jindex_global(numprocs),kindex_global(numprocs)

      double precision startTimeSolver

      character*11 fname
      character*110 filename
      character*19 date_time
      character*20 my_username,hostname
      character*10 date_calendar
      character*8  time_clock
      integer yaml_unit
      logical exist, iblank
      logical DEBUG, DEBUG1

      DEBUG = .false.
c      DEBUG = .true.
      DEBUG1 = .false.
c      DEBUG1 = .true.

c
      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' In smac2d.f/solve:'
       print*, '  jmax,kmax,idimj,idimk = ',jmax,kmax,idimj,idimk
       print*, '  jindex_global,kindex_global,ioverlap = ',
     &   jindex_global,kindex_global,ioverlap
       print*
c       print*,' Stopping in smac2d.f/solve #1'
c       call MPI_ABORT(MPI_COMM_WORLD,1)
c       stop 'stop in smac2d.f/solve #1'
      endif

      if (DEBUG.and.nodeid.eq.0.and.numprocs.eq.1) then
       print*
       print*, ' nodeid  j   k      x         y'
       do 1103 k=1,1
       do 1103 j=21,99
        write(*,1107) nodeid,j,k,x(j,k),y(j,k)
1107    format(i5,1x,i5,1x,i3,1p2e13.5)
1103   continue
       print*
c       print*, ' stop: in smac2d.f/solve 10'
c       call MPI_ABORT(MPI_COMM_WORLD,1)
c       stop 'stop: in smac2d.f/solve 10'
      endif

c zero main arrays q, dq and s 
      do 1000 k=1,kmax
      do 1000 j=1,jmax
      do 1000 n=1,3
       q(j,k,n) = 0.d0
       dq(j,k,n) = 0.d0
       s(j,k,n) = 0.d0
1000  continue

      if (DEBUG) then
        print*, ' '
        print*, '** At start of sub solve: nodeid = ',nodeid
        print*, ' '
      endif

c-----
c  Read, process, and test boundary condition information
c-----

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if (DEBUG.and.nodeid.eq.0) then
      icount = 0 
      print*
      print*, '  node  i   j   k   dq1        dq2        dq3'
      print*,' This is test 130 in smac2d.f/solve'
      itest = 130 
      do 4030 j=1,jmax
      do 4030 k=1,kmax
       icount = icount + 1 
       write(*,4005) itest,nodeid,icount,j,k,dq(j,k,1),dq(j,k,2),
     &  dq(j,k,3)
4030  continue
      print*
c      print*, ' stop in smac2d.f/solve #2'
c      call MPI_ABORT(MPI_COMM_WORLD,1)
c      stop 'stop in smac2d.f/solve #2'
      endif

c-----


c read bcmain.dat for patched boundary conditions and conditions for
c  walls, inflow/outflow, wake-cut, etc.
      call rbcfile(jmax,kmax,x,y,
     & idimj,idimk,
     & jindex_global,kindex_global,ioverlap
     & )

      if (DEBUG.and.(nodeid.eq.0)) then
      icount = 0 
      print*
      print*, '  node  i   j   k   dq1        dq2        dq3'
      print*,' This is test 120 in smac2d.f/solve'
      itest = 120 
      do 4020 j=1,jmax
      do 4020 k=1,kmax
       icount = icount + 1 
       write(*,4005) itest,nodeid,icount,j,k,dq(j,k,1),dq(j,k,2),
     &  dq(j,k,3)
4020  continue
      print*
      print*
      endif

c-----

c-----
c  Compute metrics
c-----
c all nodes compute metrics (rtxy) separately

      call metric(jmax,kmax,x,y,rtxy,dj)

      if (DEBUG.and.(nodeid.eq.0)) then
      icount = 0 
      print*
      print*, '  node  i   j   k   dq1        dq2        dq3'
      print*,' This is test 110 in smac2d.f/solve'
      itest = 110 
      do 4010 j=1,jmax
      do 4010 k=1,kmax
       icount = icount + 1 
       write(*,4005) itest,nodeid,icount,j,k,dq(j,k,1),dq(j,k,2),
     &  dq(j,k,3)
4010  continue
      print*
      print*
      endif

c-----
c  set initial conditions for q variables
c  ... initializes  q(j,k,n) and sets velocities at walls to zero.
c-----
      call ic(jmax,kmax,q)

      if (DEBUG.and.(nodeid.eq.0)) then
        print*, 'DWB0: nodeid,q(2,1,1) = ',nodeid,q(2,1,1)
        print*, 'DWB0: nodeid,q(10,1,2) = ',nodeid,q(10,1,2)
      endif

      if (DEBUG.and.(nodeid.eq.0)) then
      icount = 0 
      print*
      print*, '  node  i   j   k   dq1        dq2        dq3'
      print*,' This is test 100 in smac2d.f/solve'
      itest = 100 
      do 4000 j=1,jmax
      do 4000 k=1,kmax
       icount = icount + 1 
       write(*,4005) itest,nodeid,icount,j,k,dq(j,k,1),dq(j,k,2),
     &  dq(j,k,3)
4005   format(i4,3x,i4,3x,i4,2x,i3,2x,i3,1x,1pe10.3,1x,1pe10.3,
     &  1x,1pe10.3)
4000  continue
      print*
      print*
      endif

c-----
c  CPU time spent in startup
c-----
      cputinit = MPI_WTIME()

c-----
c  Start of time iterations
c-----

c sync all nodes
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c start marching equations in time
      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, 'Number of time steps ntmax = ',ntmax
       print*
      endif

      do 100 ntime=1,ntmax

       if (DEBUG.and.nodeid.eq.0) then
        print*
        print*, ' **** Time step ',ntime,' of ',ntmax,' time steps'
        print*, ' q(32,1,X)=',q(32,1,1),q(32,1,2),q(32,1,3)
        print*, ' q(32,2,X)=',q(32,2,1),q(32,2,2),q(32,2,3)
        print*
       endif

       nt = nt + 1

c increment time using steady-state time step
       time = time + dtau

c-----
c  Advance solution
c-----

       call step(jmax,kmax,x,y,rtxy,dj,q,dq,s,
     &                vnut,turvar)
       if (DEBUG) then
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (nodeid.eq.0) then
         print*
         print*,' All nodes have finished with sub. step'
         print*
        endif
       endif

c calculate resmax and other values
         call iterout(jmax,kmax,q,x,y,rtxy,dj,vnut,
     &    yaml_unit,startTimeSolver,ioverlap)
      if (DEBUG) then
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       if (nodeid.eq.0) then
        print*
        print*,' All nodes have finished with sub. iterout'
        print*
       endif
      endif

c-----
c  Converged and Diverging solution check
c-----
         divmaxtot = 0.0
         if( resmax0 .eq. 0. ) resmax0 = 1.0

#ifdef D_PRECISION 
         divmaxtot = dmax1(divmaxtot, divmax )
#else
         divmaxtot = amax1(divmaxtot, divmax )
#endif

c end of iteration loop
100   continue

c print some boundaries to ensure bc updates are correct
c-----start-----
c node 0, j=1
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG.and.nodeid.eq.0) then
        print*
        do 7990 k=1,kmax
         print*, ' nodeid0,j=1,q123 = ',
     &    k,q(1,k,1),q(1,k,2),q(1,k,3)
7990    continue
      endif
c-----
c node 0, j=jmax
      call flush(6)
      call flush(istdout)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG.and.nodeid.eq.0) then
        print*
        do 8000 k=1,kmax
         print*, 'nodeid0,jmax,q123 = ',
     &    k,q(jmax,k,1),q(jmax,k,2),q(jmax,k,3)
8000    continue
      endif
c-----
c node 1, j=2
      call flush(6)
      call flush(istdout)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG.and.nodeid.eq.1) then
       print*
       do 8010 k=1,kmax
        print*, 'nodeid1,j_2,k_1tokmax = ',
     &   k,q(2,k,1),q(2,k,2),q(2,k,3)
8010   continue
      endif
c-----
c node 1, j=jmax
      call flush(6)
      call flush(istdout)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG.and.nodeid.eq.1) then
       print*
       do 8015 k=1,kmax
        print*, 'nodeid1, jmax, k_1tokmax = ',
     &   k,q(jmax,k,1),q(jmax,k,2),q(jmax,k,3)
8015   continue
      endif
c-----
      call flush(6)
      call flush(istdout)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG.and.nodeid.eq.0) then
       print*
       do 8020 j=1,jmax
        print*, 'nodeid0, j_1tojmax, k_1 = ',
     &   j,q(j,1,1),q(j,1,2),q(j,1,3)
8020   continue
      endif
c-----
      call flush(6)
      call flush(istdout)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG.and.nodeid.eq.0) then
       print*
       do 8025 j=1,jmax
        print*, 'nodeid0, j_1tojmax, k_1 = ',
     &   j,q(j,1,1),q(j,1,2),q(j,1,3)
8025   continue
      endif
c-----
      call flush(6)
      call flush(istdout)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG.and.nodeid.eq.1) then
       print*
       do 8030 j=1,jmax
        print*, 'nodeid1, j_1tojmax, k_1 = ',
     &   j,q(j,1,1),q(j,1,2),q(j,1,3)
8030   continue
      endif
c-----
      call flush(6)
      call flush(istdout)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (DEBUG.and.nodeid.eq.1) then
       print*
       do 8035 j=1,jmax
        print*, 'nodeid1, j_1tojmax, k_kmax = ',
     &   j,q(j,kmax,1),q(j,kmax,2),q(j,kmax,3)
8035   continue
      endif

c-----end-----


      ntime = ntime - 1

      return
      end
c
c
c************************************************************************
      subroutine dateline(chline,date_calendar,time_clock,date_time,
     & hostname)
c************************************************************************
#include "precis.h"
      character*(*) chline
      character*32 string
      character*80 host
      integer*4 ishost, hostnm
      character*20 date_time
      character*10 date_calendar
      character*8 time_clock
      character*20 hostname
c      character*20 username 
      character*2 day_num,month_num
      character*3 month_text
      character*4 year_num
      character*8 date, time
      character*2 time_hour,time_min,time_sec
      character*32 date_string 
      character*2 Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec
      logical DEBUG
c define months
      data Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec 
     &  /'01','02','03','04','05','06','07','08','09','10',
     &   '11','12'/

c
      DEBUG = .false.
c      DEBUG = .true.
c
      chline(:) = ' '

c---
c linux
c----	
      call fdate(date_string)
      chline(1:24) = date_string 
c 'hostnm' is an intrinsic function; = 0 on success
      icall = hostnm(host)
      if (icall.ne.0) then
        host = 'UNKNOWN'
        ihost = 7
      else 
         do i=2,80
            if(host(i:i) .eq. ' ' .and. host(i-1:i-1) .ne. ' ') then
               ihost = i-1
               goto 12
            else
               ihost = i
            endif
         enddo
12       continue
         chline(25:28) = ' on '
         chline(29:28+ihost) = host(1:ihost)
      endif
c---

      year_num = chline(21:24)
      month_text = chline(5:7)
      day_num = chline(9:10)
      time = chline(12:19)
      time_hour = chline(12:13)
      time_min = chline(15:16)
      time_sec = chline(18:19)
      if (DEBUG) then
        print*, ' ihost = ',ihost
        print*, ' host = ',host
      endif
      hostname = host(1:ihost) 
      if (DEBUG) then
        print*, 'chline = ',chline
      endif

c pad single-number days with zero
      if (day_num(1:1).eq.' ') then
        day_num = '0'//day_num(2:2)
      endif

c change month to numbers
      if (month_text .eq. 'Jan') then
        month_num = Jan
      else if (month_text .eq. 'Feb') then
        month_num = Feb
      else if (month_text .eq. 'Mar') then
        month_num = Mar
      else if (month_text .eq. 'Apr') then
        month_num = Apr
      else if (month_text .eq. 'May') then
        month_num = May 
      else if (month_text .eq. 'Jun') then
        month_num = Jun 
      else if (month_text .eq. 'Jul') then
        month_num = Jul 
      else if (month_text .eq. 'Aug') then
        month_num = Aug 
      else if (month_text .eq. 'Sep') then
        month_num = Sep 
      else if (month_text .eq. 'Oct') then
        month_num = Oct 
      else if (month_text .eq. 'Nov') then
        month_num = Nov 
      else if (month_text .eq. 'Dec') then
        month_num = Dec 
      endif

      if (DEBUG) then
        print*
        print*, 'Date parameters:'
        print*, 'year_num = ',year_num
        print*, 'month_text = ',month_text
        print*, 'month_num = ',month_num
        print*, 'day_num = ',day_num
        print*, 'time = ',time
        print*, 'hostname = ',hostname
        print*
      endif

      date_calendar = year_num//'_'//month_num//'_'//day_num
      time_clock = time

c new way
      date_time = year_num//'_'//month_num//'_'//day_num//
     &  '__'//time_hour//'_'//time_min//'_'//time_sec

c-----
c  End of dateline
c-----
      return
      end
c
c
c******************************************************************
      subroutine iterout(jmax,kmax,q,x,y,rtxy,dj,vnut,
     & yaml_unit,startTimeSolver,ioverlap)
c******************************************************************
c  Purpose:
c   prints each iteration step
c------------------------------------------------------------------
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      __REAL q(jmax,kmax,3), x(jmax,kmax), y(jmax,kmax), 
     &  rtxy(jmax,kmax,2,3), dj(jmax,kmax)
      __INTEGER jmax, kmax
c
      __REAL cdc, clc, cmc
      __REAL rvar(5)
      __INTEGER ivar(2)
      common/mulelr/chord,pinf,xrot,yrot,xmom,ymom
      common/write_p3d_q/write_plot3d_q_file
      character*1 nf_text
      character*1 nzslat_text
      character*110 filename
      integer yaml_unit
      logical DEBUG
      logical write_plot3d_q_file
      double precision startTimeSolver,del_cputime
      save cl,cd,cm,clold,cdold,cmold,ntimeold
      __INTEGER stat(MPI_STATUS_SIZE)
      __INTEGER jmax_zone(numprocs),kmax_zone(numprocs)
      __REAL,allocatable:: q_temp(:,:,:)

      DEBUG = .false.
c      DEBUG = .true.

      if (DEBUG.and.nodeid.eq.0) then
       print*
       print*, ' >>> In sub geom.f/iterout. Nodeid = ',nodeid 
       print*
      endif

c-----
c  Compute lift, drag and moment due to pressure and skin friction
c-----
      pi = 4.d0*datan( 1.d0 )
      ca = cos( alpha*pi/180.d0 )
      sa = sin( alpha*pi/180.d0 )
      xmoml = xmom
      ymoml = ymom
      if(ntimeold .ne. ntime) then
         cdold = cd
         clold = cl
         cmold = cm
         ntimeold = ntime
      endif
      cd = 0.0
      cl = 0.0
      cm = 0.0
      cdp = 0.0
      clp = 0.0
c
      fxp = 0.0
      fyp = 0.0
      fxs = 0.0
      fys = 0.0
      fmom = 0.0

      call force(jmax,kmax,x,y,q,rtxy,dj,
     &          xmoml,ymoml,fxp,fyp,fxs,fys,fmom,ioverlap)


c values for this zone
      cd = (2.*(fxp + fxs)*ca + 2.*(fyp + fys)*sa)/chord
      cl = (2.*(fyp + fys)*ca - 2.*(fxp + fxs)*sa)/chord
      cm = -2.*fmom/(chord*chord)

c total values
c ... cd due to pressure only
      cdp = (2.*fxp*ca + 2.*fyp*sa)/chord
c ... cl due to pressure only
      clp = (2.*fyp*ca - 2.*fxp*sa)/chord

      if (DEBUG) then
       print*
       print*, ' nodeid,fxp,fyp,fxs,fys = ',nodeid,fxp,fyp,fxs,fys
       print*, ' nodeid,ca,sa,chord = ',nodeid,ca,sa,chord
       print*, ' Node',nodeid,': cd,cl,cm = ',cd,cl,cm
       print*
       call flush(6)
       call flush(istdout)
c       print*, 'stop: in smac2d.f/iterout for forces'
c       call MPI_ABORT(MPI_COMM_WORLD,1)
c       stop 'stop: in smac2d.f/iterout for forces'
      endif

c gather coefficients from all zones and determine total values

c*********************************
      if (numprocs.gt.1) then
c*********************************
      if (nodeid.eq.0) then
c*********************************
        if (DEBUG.and.nt.eq.ntmax) then
         print*, ' node_in, cl_temp = ',nodeid,cl
        endif
        do 200 i=2,numprocs
          node_in = i-1
c receive
#ifdef D_PRECISION
          call MPI_RECV(rvar,5,MPI_REAL8,node_in,node_in,
     &      MPI_COMM_WORLD,stat,ierr)
#else
          call MPI_RECV(rvar,5,MPI_REAL,node_in,node_in,
     &      MPI_COMM_WORLD,stat,ierr)
#endif

c unpack reals
          cd_temp = rvar(1)
          cl_temp = rvar(2)
          cm_temp = rvar(3)
          cdp_temp = rvar(4)
          clp_temp = rvar(5)
c add to node 0 totals thus far
          cd = cd + cd_temp 
          cl = cl + cl_temp
          cm = cm + cm_temp
          cdp = cdp + cdp_temp
          clp = clp + clp_temp
          if (DEBUG.and.nt.eq.ntmax) then
           print*,' node_in, cl_temp,cl = ',node_in,cl_temp,cl
          endif
200     continue

c*********************************
      else
c*********************************

c pack variables
      rvar(1) = cd
      rvar(2) = cl
      rvar(3) = cm
      rvar(4) = cdp
      rvar(5) = clp

c send

#ifdef D_PRECISION
      call MPI_SEND(rvar,5,MPI_REAL8,0,nodeid,
     &  MPI_COMM_WORLD,ierr)
#else
      call MPI_SEND(rvar,5,MPI_REAL,0,nodeid,
     &  MPI_COMM_WORLD,ierr)
#endif

c*******************************t*
      endif
c*********************************
c*********************************
      endif
c*********************************
      clift = cl

c-----
c  Find max residual on each node
c-----
         rmax = 0.0
         tmax = 0.0
         dmax = 0.0
         nmxr = 0
         tmax = max( tmax, turres) 
         dmax = max( dmax, divmax)
         if (resmax0.eq.0.d0) resmax0 = 1.0
         rmax = resmax/resmax0
         jmxr = jres
         kmxr = kres

c*********************************
      if (numprocs.gt.1) then
c*********************************
c*********************************
      if (nodeid.eq.0) then
c*********************************
        if (DEBUG) then
         print*, ' nodeid,dmax,tmax = ',0,dmax,tmax
        endif
        if (resmax.eq.0.d0) then
         resmax = 1.d0
        endif
        rmax = resmax/resmax0      ! for node 0
        if (DEBUG) then
         print*, ' i, nodeid, turres = ',1,nodeid,turres
        endif
        do 275 i=2,numprocs
          node_in = i-1

c receive

#ifdef D_PRECISION
          call MPI_RECV(rvar,3,MPI_REAL8,node_in,node_in,
     &      MPI_COMM_WORLD,stat,ierr)
#else
          call MPI_RECV(rvar,3,MPI_REAL,node_in,node_in,
     &      MPI_COMM_WORLD,stat,ierr)
#endif

          call MPI_RECV(ivar,2,MPI_INTEGER,node_in,node_in+1,
     &      MPI_COMM_WORLD,stat,ierr)

c unpack variables; denote as temp unless greater than current
c  max values
          if (DEBUG) then
           print*, ' i, nodeid, turres = ',i,nodeid,rvar(1)
           print*, ' nodeid,dmax,tmax = ',node_in,dmax,tmax
           print*, ' nodeid,resmax,resmax0 = ',resmax,resmax0
           print*, ' -------------------------------------------'
          endif
          turres = max(tmax,rvar(1))
          divmax = max(dmax,rvar(2))
          resmax_temp = rvar(3)
          jres_temp = ivar(1)
          kres_temp = ivar(2)
          tmax = max(tmax,turres)
          dmax = max(dmax,divmax)
         if( resmax_temp/resmax0 .gt. rmax) then
            resmax = resmax_temp
            rmax = resmax/resmax0
            jres = jres_temp
            kres = kres_temp
            jmxr = jres
            kmxr = kres
            nmxr = node_in + 1   ! zone, not the nodeid
         endif
         if (DEBUG) then
          print*
          write(*,500) node_in,resmax,resmax0,rmax,tmax
500       format(/' In iterout: node_in,resmax,resmax0,rmax,tmax = ',
     &      i3,1p4e15.7)
         endif
275    continue
c

c*********************************
      else
c*********************************

c pack variables
c ... real
        rvar(1) = turres
        rvar(2) = divmax
        rvar(3) = resmax
c ... integer
        ivar(1) = jres
        ivar(2) = kres

c send

#ifdef D_PRECISION
        call MPI_SEND(rvar,3,MPI_REAL8,0,nodeid,
     &    MPI_COMM_WORLD,ierr)
#else
        call MPI_SEND(rvar,3,MPI_REAL,0,nodeid,
     &    MPI_COMM_WORLD,ierr)
#endif

        call MPI_SEND(ivar,2,MPI_INTEGER,0,nodeid+1,
     &    MPI_COMM_WORLD,ierr)


c*********************************
      endif
c*********************************
c*********************************
      else      ! numprocs = 1
c*********************************

      tmax = max(tmax,turres)
      dmax = max(dmax,divmax)
      if (resmax/resmax0 .gt. rmax) then
        rmax = resmax/resmax0
        jmxr = jres
        kmaxr = kres
        nmxr = 1
      endif
c*********************************
      endif
c*********************************

      if (DEBUG) then
       if (nodeid.eq.0) then
        print*
        print*, ' max turres = ',tmax
c       print*, ' stop in smac2d.f'
c       call MPI_ABORT(MPI_COMM_WORLD,1)
c       stop 'stop in smac2d.f'
       endif
      endif

c*********************************
       if (nodeid.eq.0) then
c*********************************
c-----
c  CPU time
c-----
         cputime = MPI_WTIME()
         del_cputime = cputime - startTimeSolver

c-----
c  Output at every niter iterations -- Steady-state
c-----
c 
c gather info from all nodes and print out max and total values

        niter = 1

        if (DEBUG) then
         print*
         print*, ' Values to print from  Node 0:'
         print*, ' niter,ifreq,nt = ',niter,ifreq,nt
        endif

         ifreq = 25*niter
         if(ifreq .eq. 0) ifreq = 1
         if( ((nt-1)/ifreq)*ifreq .eq. nt-1 )
     &    write(istdout,100)
         if ((nt/niter)*niter .eq. nt .or. nt.eq.1) then
          write(istdout,110) nt,rmax,dmax,tmax,jmxr,kmxr,nmxr,
     &                del_cputime,cd,cl
          if (DEBUG) then
           print*, ' resmax,resmax0,rmax = ',resmax,resmax0,rmax
           print*, ' jmxr,kmxr,nmxr = ',jmxr,kmxr,nmxr
           print*, ' dmax,tmax = ',dmax,tmax
           print*
          endif
         endif
100      format(/'     nt   resmax     divmax     turres      j     k  ' 
     &           '  nz     cputime       cd          cl'/,
     &           '   ----  ---------  ---------  ---------  ----  ---- '
     &           ' ----  ----------  ----------  ----------')
110      format('XX',i5,1p3e11.3,3i6,1p3e12.4)
c*********************************
      endif
c*********************************

c output solution file in plot3d format
c+++++++++++++++++++++++++++++++
      if (ntime.eq.ntmax.and.write_plot3d_q_file) then
c+++++++++++++++++++++++++++++++
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c+++++++++++++++++++++++++++++++
       if (nodeid.eq.0) then
c+++++++++++++++++++++++++++++++
        write(6,*) ' '
        write(6,*) ' >> WRITING 2-D FLOW VARIABLES TO FILE >Q2D.p3d< IN 
     1PLOT3D FORMAT'

c OUTPUT: open the Q2D.p3d output file used for flow variables; in PLOT3D format
        filename='Q2D.p3d'
        filename=trim(filename)
        open(25,file=filename,iostat=ierr_unit25,form='formatted',
     1   status='unknown')
        if (ierr_unit25.ne.0) then
          write(6,*) ' '
          write(6,*) ' ERROR: in main routine smac2d.f'
          write(6,*) '   File Q2D.p3d could not be opened.'
          write(6,*) ' '
          write(6,*) '   PROGRAM STOPPED: Problem with i/o file.'
          write(6,*) ' '
          call MPI_ABORT(MPI_COMM_WORLD,1)
          stop 'smac2d.f: problem with Q2D.p3d file'
        endif
c
c get all jmax,kmax values from all nodes and write to Q2D.p3d file
      write(25,*) numprocs
      write(25,*) jmax,kmax
      jmax_zone(1) = jmax
      kmax_zone(1) = kmax
      do 400 i=2,numprocs
       nodeid_send = i-1
       if (DEBUG) then
        print*
        print*, ' RECV from nodeid',nodeid_send,':'
       endif
       call MPI_RECV(ivar,2,MPI_INTEGER,nodeid_send,nodeid_send,
     1    MPI_COMM_WORLD,stat,ierr)
       write(25,*) ivar(1),ivar(2)
       jmax_zone(i) = ivar(1)
       kmax_zone(i) = ivar(2)
       if (DEBUG) then
        print*, ' jmax_zone(i),kmax_zone(i) =',jmax_zone(i),kmax_zone(i)
       endif
400   continue

c+++++++++++++++++++++++++++++++
      else
c+++++++++++++++++++++++++++++++
c
        ivar(1) = jmax
        ivar(2) = kmax
        call MPI_SEND(ivar,2,MPI_INTEGER,0,nodeid,
     1     MPI_COMM_WORLD,ierr)

c+++++++++++++++++++++++++++++++
       endif
c+++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++
      endif
c+++++++++++++++++++++++++++++++
c
c now get and write all q values (flow field variables) to Q2D.p3d
c+++++++++++++++++++++++++++++++
      if (ntime.eq.ntmax.and.write_plot3d_q_file) then
c+++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++
       if (nodeid.eq.0) then
c+++++++++++++++++++++++++++++++

c for node 0...
        write(25,205) 0.0,alpha,reynum,cputime
        allocate(q_temp(jmax,kmax,4))
        do 405 k=1,kmax
        do 405 j=1,jmax
          pp = q(j,k,1)
          uu = q(j,k,2)
          vv = q(j,k,3)
          q_temp(j,k,1) = 1.0
          q_temp(j,k,2) = uu
          q_temp(j,k,3) = vv
          q_temp(j,k,4) = pp/0.4 + 0.5*(uu*uu + vv*vv)
405     continue
        write(25,206) (((q_temp(j,k,n),j=1,jmax),k=1,kmax),n=1,4)
        deallocate(q_temp)
c
c for all other nodes...
      do 410 i=2,numprocs
       nodeid_send = i-1
       allocate (q_temp(jmax_zone(i),kmax_zone(i),4))

#ifdef D_PRECISION
       call MPI_RECV(q_temp,jmax_zone(i)*kmax_zone(i)*4,MPI_REAL8,
     1   nodeid_send,nodeid_send,MPI_COMM_WORLD,stat,ierr)
#else
       call MPI_RECV(q_temp,jmax_zone(i)*kmax_zone(i)*4,MPI_REAL,
     1   nodeid_send,nodeid_send,MPI_COMM_WORLD,stat,ierr)
#endif

       write(25,205) 0.0,alpha,reynum,cputime
      write(25,206) (((q_temp(j,k,n),
     1  j=1,jmax_zone(i)),k=1,kmax_zone(i)),n=1,4)
      deallocate(q_temp)
410   continue

c+++++++++++++++++++++++++++++++
      else
c+++++++++++++++++++++++++++++++
       allocate (q_temp(jmax,kmax,4))
       do 407 k=1,kmax
       do 407 j=1,jmax
         pp = q(j,k,1)
         uu = q(j,k,2)
         vv = q(j,k,3)
         q_temp(j,k,1) = 1.0
         q_temp(j,k,2) = uu
         q_temp(j,k,3) = vv
         q_temp(j,k,4) = pp/0.4 + 0.5*(uu*uu + vv*vv)
407      continue

#ifdef D_PRECISION
       call MPI_SEND(q_temp,jmax*kmax*4,MPI_REAL8,0,nodeid,
     1   MPI_COMM_WORLD,ierr)
#else
       call MPI_SEND(q_temp,jmax*kmax*4,MPI_REAL,0,nodeid,
     1   MPI_COMM_WORLD,ierr)
#endif

       deallocate (q_temp)
        
c+++++++++++++++++++++++++++++++
       endif
c+++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++
      endif
c+++++++++++++++++++++++++++++++

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if (nodeid.eq.0.and.ntime.eq.ntmax.and.
     &   (.not.write_plot3d_q_file)) then
       write(6,*) ' '
       write(6,*) ' '
       write(6,*) 
     & ' >> 2-D FLOW VARIABLES WERE NOT WRITTEN TO FILE >Q2D.p3d< IN PLO
     &T3D FORMAT'
       write(6,'(a,l2)') ' write_plot3d_q_file =',write_plot3d_q_file
       write(6,*) ' '
      endif 



c+++++++++++++++++++++++++++++++
      if (nodeid.eq.0.and.ntime.eq.ntmax) then
c+++++++++++++++++++++++++++++++
c-----
c  Output summary of results for airfoil computations
c  Useful in parametric studies
c-----
         write(istdout,300) alpha, reynum, cd, cl, cm
         write(istdout,310) rmax, dmax, abs(cd-cdold), abs(cl-clold)
c output forces and moment
         write(istdout,210) cd,cl,cm
c output to yaml file
         write(yaml_unit,*) 'airfoil:'
         write(yaml_unit,*) '   cd: ',cd
         write(yaml_unit,*) '   cl: ',cl
         write(yaml_unit,*) '   cm: ',cm
         write(yaml_unit,*) 'residual: ',rmax
c         write(yaml_unit,*) 'divergence: ',dmax
c         write(yaml_unit,*) 'delta_cl: ',abs(cl-clold)
c         write(yaml_unit,*) 'delta_cd: ',abs(cd-cdold)
c+++++++++++++++++++++++++++++++
      endif
c+++++++++++++++++++++++++++++++
c
c format statements
205   format(1pf5.2,1x,1pe13.5,1x,1pe13.5,1x,1pe13.5)
206   format(6f10.6)
210   format(/'airfoil cd,cl,cm: ',3e12.4)
300   format(/'**************************************************',/,
     &        ' summary_of_results:alpha,re,cd,cl,cm: ',
     &          f10.2,e12.4,3e14.6)
310   format( ' summary_of_convergence:resmax,divmax,delta(cd,cl): ',
     &          4e12.5,
     &       /'**************************************************',/)

c-----
c  End of iterout
c-----
      return
      end
