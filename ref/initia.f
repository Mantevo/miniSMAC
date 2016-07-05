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
c  1. initia
c  2. inizone 
c
c************************************************************
c
c
c**********************************************************
      subroutine initia(geometry,title)
c**********************************************************

c  All ranks are initialized with the following variables
c
c----------------------------------------------------------
#include "precis.h"
#include "mpif.h"
#include "mpi_params.f"
c
      common/impli/impsch,njsp,nksp,nbcimp,ntjsp,ntksp
      common/mass/pin,pout
      common/timm/beta,dtau,dt,time
      common/visc/reynum,vnu
      common/iparmi/nt,niter,ntmax,ntime,iflxo,ivis,iturb
      common/airfr/alpha,clift,cdrag,cmom
      common/implr/underr,dcoef2,dcoef4
      common/unitno/istdout
      common/gridmods/kadd
      common/write_p3d_q/write_plot3d_q_file
      common/bcmain_location/bcmain_directory,bcmain_filename
      common/xygrid_location/xygrid_directory,xygrid_filename

      character*60 title,geometry, namelist_input_file
      character*100 bcmain_directory,bcmain_filename
      character*100 xygrid_directory,xygrid_filename

      logical DEBUG
      logical zz_run_default_input
      logical write_plot3d_q_file
      double precision, parameter :: pi = 4.d0*datan(1.d0)
      integer nodeid,numprocs,ierr,stat(MPI_STATUS_SIZE)
      __REAL rvar(6)
      __INTEGER ivar(4)

      namelist/input/
     &  alpha,beta,dcoef2,dcoef4,dt,dtau,geometry,impsch,iturb,ivis,
     &  kadd,ntmax,reynum,title,underr,zz_run_default_input,
     &  write_plot3d_q_file,bcmain_directory,bcmain_filename,
     &  xygrid_directory,xygrid_filename
c********************************************************

      DEBUG = .false.
c      DEBUG = .true.

      if (DEBUG) then
       print*
       print*, ' >>> Entering sub initia - node = ',nodeid,' <<<'
       print*, '   nodeid = ',nodeid
       print*, '   numprocs = ',numprocs
      endif

c Default values for select variables 
c ... ints
      nbcimp = 1        ! line relaxation sweep implicit bc update frequency
      niter = 1         ! output at every niter iterations
      istdout = 6       ! standard output handle
c ... reals
      pin = 1.5D0
      pout = 1.0D0
c ... number of grid points added to outer boundary; used for increasing the distance
c     of the grid outer boundary and for scaling studies
      kadd=0

c Read input namelist from file
c-------------------------------
      if (nodeid.eq.0) then
c-------------------------------
       namelist_input_file = 'smac2d.in'
       open (unit=900,file=namelist_input_file,err=998)
       goto 999
998    continue
       write(istdout,11) namelist_input_file
       write(*,11) namelist_input_file
11     format(/'ERROR: namelist input file not found:',/,
     & '   file = ',a,/,
     & '    Program terminating.'/)
       call MPI_ABORT(MPI_COMM_WORLD,1)
       stop 'stopping: namelist input file not found'

999    continue

       read(900,input,err=991)
       goto 992
991    continue
       write(istdout,12) namelist_input_file
       write(*,12) namelist_input_file
12     format(/'ERROR: cannot read from namelist input file:',/,
     &  '   file = ',a,/,
     &  '    Program terminating.'/)
       call MPI_ABORT(MPI_COMM_WORLD,1)
       stop 'stopping: cannot read from namelist input file'

992    continue
       print*
       print*, 'From namelist file:',namelist_input_file
       print*, '-----------------------------------------------'
       write(istdout,input)
       close(900)

c DEFAULT VALUES SPECIFIC TO AIRFOIL (airfoil)
c-------------------
       if (zz_run_default_input) then
c-------------------
c ... use default values for airfoil
        alpha = 8.0d0             ! degs, airfoil angle of attack
        beta = 5.D0               ! artificial compressibility factor
        dcoef2 = 0.00d0		  ! 2nd-order smoothing coefficient
        dcoef4 = 0.00d0           ! 4th-order smoothing coefficient
        dt = 1.0d12               ! time step for turbulence model
        dtau=0.1d0                ! pseudo-time step for steady state case 
        geometry='airfoil'        ! geometry under consideration
        impsch = 1                ! implicit scheme, 1=line relaxation; 2=lusgs
        iturb = 1		  ! use Spallart-Allmaras turbulence model
        ivis = 1                  ! variable viscosity and non-orthogonal grid
        kadd=0                    ! add layers (points) to grid outer boundary
        ntmax = 500               ! number of time steps to take
        reynum = 1.50D6           ! reynolds number of 1,500,000
        title='SMAC2D -- Sandia MiniAero Code 2D - DEFAULT VALUES' !title
        underr = 0.01D0           ! under-relaxation parameter, 0 < underr < 1
        write_plot3d_q_file=.false. ! do not write plot3d q file, may be big
        bcmain_directory = '.'
        bcmain_filename = 'bcmain.dat_original'
        xygrid_directory = '.'
        xygrid_filename = 'xy.dat_original'
        print*
        print*, 'DEFAULT AIRFOIL VALUES ARE BEING USED FOR THIS RUN:'
        if (DEBUG) then
         write(*,955) alpha,beta,dt,dtau,trim(geometry),impsch,iturb,
     &    ivis,ntmax,reynum,trim(title),underr,write_plot3d_q_file,
     &    bcmain_directory,bcmain_filename,
     &    xygrid_directory,xygrid_filename
955      format('alpha = ',f6.3,/,
     &     'beta = ',f6.3,/,
     &     'dt = ',1pd14.7,/,
     &     'dtau = ',1pd14.7,/,
     &     'geometry = ',a,/,
     &     'impsch = ',i1,/,
     &     'iturb = ',i1,/,
     &     'ivis = ',i1,/,
     &     'kadd = ',i10,/,
     &     'ntmax = ',i7,/,
     &     'reynum = ',1pd14.7,/,
     &     'title = ',a,/,
     &     'underr = ',1pd11.4,/,
     &     'write_plot3d_q_file = ',l2,/,
     &     'bcmain_directory = ',a,/,
     &     'bcmain_filename = ',a,/,
     &     'xygrid_directory = ',a,/,
     &     'xygrid_filename = ',a,/
     &     )
       endif
      endif
c a few checks on input values
c ... ivis
      if (ivis.ne.0.and.ivis.ne.1) then
       print*
       print*, 'ERROR: invalid value for ivis on nodeid =',nodeid
       print*, '  ivis must be either 0 (inviscid) or 1 (viscous)'
       print*, '  Current value of ivis =',ivis
       print*, 'Program terminating.'
       print*
       call MPI_ABORT(MPI_COMM_WORLD,1)
       stop 'stop: invalid value for ivis'
      endif
c ... iturb
      if (iturb.ne.0.and.iturb.ne.1) then
       print*
       print*, 'ERROR: invalid value for iturb on nodeid =',nodeid
       print*, '  iturb must be either 0 (laminar) or 1 (turbulent)'
       print*, '  Current value of iturb =',iturb
       print*, 'Program terminating.'
       print*
       call MPI_ABORT(MPI_COMM_WORLD,1)
       stop 'stop: invalid value for iturb'
      endif
c ... kadd
      if (kadd.lt.0) then
       print*
       print*, 'ERROR: invalid value for kadd on nodeid =',nodeid
       print*, '  kadd must be equal to or greater than zero.'
       print*, '  Current value of kadd =',kadd
       print*, 'Program terminating.'
       print*
       call MPI_ABORT(MPI_COMM_WORLD,1)
       stop 'stop: invalid value for kadd'
      endif
c ... underr
      if (underr.le.0.0d0.or.underr.gt.1.0d0) then
       print*
       print*, ' ERROR: invalid value for underr on nodeid =',nodeid
       print*, '   underr must be between 0 and 1'
       print*, '   Current value of underr =',underr
       print*, 'Program terminating.'
       print*
       call MPI_ABORT(MPI_COMM_WORLD,1)
       stop 'stop: invalid value for underr'
      endif
c Send values to other nodes
c ... pack the namelist input variables
c ...   integers
      ivar(1) = impsch
      ivar(2) = iturb
      ivar(3) = ivis
      ivar(4) = ntmax
c ...   reals
      rvar(1) = alpha
      rvar(2) = beta
      rvar(3) = dt
      rvar(4) = dtau
      rvar(5) = reynum
      rvar(6) = underr
c ... send integers
      do 888 i=2,numprocs
       node_out = i-1
c ... send integers
       call MPI_SEND(ivar,4,MPI_INTEGER,node_out,node_out,
     &  MPI_COMM_WORLD,ierr)

c ... send reals
#ifdef D_PRECISION
       call MPI_SEND(rvar,6,MPI_REAL8,node_out,node_out+1,
     &  MPI_COMM_WORLD,ierr)
#else
       call MPI_SEND(rvar,6,MPI_REAL,node_out,node_out+1,
     &  MPI_COMM_WORLD,ierr)
#endif

888   continue

c-------------------
      else
c-------------------
c for nodes other than node 0
c ... receive integers
       call MPI_RECV(ivar,4,MPI_INTEGER,0,nodeid,
     &  MPI_COMM_WORLD,stat,ierr)
c ...   unpack integers
       impsch = ivar(1)
       iturb = ivar(2)
       ivis = ivar(3)
       ntmax = ivar(4)

c ... receive reals
#ifdef D_PRECISION
       call MPI_RECV(rvar,6,MPI_REAL8,0,nodeid+1,
     &  MPI_COMM_WORLD,stat,ierr)
#else
       call MPI_RECV(rvar,6,MPI_REAL,0,nodeid+1,
     &  MPI_COMM_WORLD,stat,ierr)
#endif

c ...   unpack reals
       alpha = rvar(1)
       beta = rvar(2)
       dt = rvar(3)
       dtau = rvar(4)
       reynum = rvar(5)
       underr = rvar(6)
c-------------------
      endif
c-------------------

      if (DEBUG) then
       call flush(6)
       call flush(istdout)
c       call MPI_ABORT(MPI_COMM_WORLD,1)
c       stop 'stop: after reading namelist input file'
      endif

c-----
c  Initialize following variables on all nodes
c-----
c ... kinemative viscosity
      vnu = 1.0D0/reynum
c ... initial conditions
      ntime = 0
      nt = 0
      time = 0.0
      clift = 0.0

c-----
c  Print initialized variables
c-----
c********************
      if (nodeid.eq.0) then
c********************
       print*
       print*, '+++++++++++++++++++++++++++++++++++++++++++++++++++++'
       print*, '       Flow Conditions and Initialized Variables     '
       print*, '+++++++++++++++++++++++++++++++++++++++++++++++++++++'
       write(*,550) title
550    format(' title: ',a)
       write(*,555) geometry
555    format(' geometry: ',a)
       print*, 'nbcimp: ',nbcimp
       print*, 'niter: ',niter
       print*, 'ntmax: ',ntmax
       print*, 'impsch: ',impsch
       print*, 'iturb: ',iturb
       print*, 'ivis: ',ivis
       print*
#ifdef D_PRECISION
       print*, 'fp_precision: double'
#else
       print*, 'fp_precision: single'
#endif
       print*
       print*, 'beta: ',beta
       print*, 'dtau: ',dtau
       print*, 'dt: ',dt
       print*
       print*, 'kadd: ',kadd
       print*
       print*, 'pin: ',pin
       print*, 'pout: ',pout
       print*, 'reynum: ',reynum
       print*, 'underr: ',underr
       print*, 'vnu = 1./reynum: ',vnu
       print*, 'write_plot3d_q_file: ',write_plot3d_q_file
       print*, 'bcmain_directory: ',bcmain_directory
       print*, 'bcmain_filename: ',bcmain_filename
       print*, 'xygrid_directory: ',xygrid_directory
       print*, 'xygrid_filename: ',xygrid_filename
       print*, '+++++++++++++++++++++++++++++++++++++++++++++++++++++'
       print*
c********************
      endif
c********************

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if (DEBUG) then
       print*
       print*, ' Finished in initia.f/initia 100 - node =',nodeid
       call flush(6)
       print*, ' This is nodeid =',nodeid
       call MPI_ABORT(MPI_COMM_WORLD,1)
       stop 'stop: at end of initia.f/initia'
      endif

c-----
c  End of initia
c-----
      return
      end
c
c
c***********************************************************************
      subroutine inizone
c***********************************************************************
c  Initializes variables specific to one particular zone.
c
c------------------------------------------------------------------
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      integer nodeid,numprocs,ierr,stat(MPI_STATUS_SIZE)

      logical kpr
      logical DEBUG

      DEBUG=.false.
c      DEBUG = .true.

c-----
c DEFAULT SETTINGS 
c-----
      kpr = .false.
      njswp = 2
      nkswp = 2
      ntjsp = 2
      ntksp = 0
      iflxodr = 3
      itran = 0
      xtr1 = 0.0
      xtr2 = 0.0
      cdiss = 1.0
   
      njsp = njswp
      nksp = nkswp
      kend2 = 2
      kendm = kmax - 1
      iflxo = iflxodr
      itrans = itran
      cdis = cdiss
      xtran1 = xtr1
      xtran2 = xtr2

      if (nodeid.eq.0) then
       print*
       print*, ' Settings from initia.f/inizone:'
       print*, '  kpr = ',kpr
       print*, '  njswp = ',njswp
       print*, '  nkswp = ',nkswp
       print*, '  ntjsp = ',ntjsp
       print*, '  ntksp = ',ntksp
       print*, '  iflxodr = ',iflxodr
       print*, '  itran = ',itran
       print*, '  xtr1 = ',xtr1
       print*, '  xtr2 = ',xtr2
       print*, '  cdiss = ',cdiss
       print*, '  njsp = ',njsp
       print*, '  nksp = ',nksp
       print*, '  kend2 = ',kend2
       print*, '  kendm = ',kendm
       print*, '  iflxo = ',iflxo
       print*, '  itrans = ',itrans
       print*, '  cdis = ',cdis
       print*, '  xtran1 = ',xtran1
       print*, '  xtran2 = ',xtran2
       print*, '    ---  end ---'
       print*
      endif


      if (DEBUG) then
         print*
         print*, ' >> Node ',nodeid,': finished in initia.f/inizone '
      endif

c-----
c  End inizone
c-----
      return
      end
c
c
