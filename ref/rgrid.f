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
c  1. rgridh
c  2. rgrid
c
c************************************************************
c
c
c************************************************************************
      subroutine rgridh(iunit,jmax,kmax,form,mzone,iblank,
     1   nzne,nzone,jkmax,istdout,jkmx,iselect)
c************************************************************************
c  This routine attempts to find the grid file, a plot3d 2D type file,
c  and determines whether it is formatted or unformatted, single or
c  multiple zone, and with or without an iblank array.  It reads in the
c  dimensions of the grid and performs parameter checks. 
c  The grid file is left open after determining grid dimensions. 
c
c------------------------------------------------------------------------
#include "precis.h"
#include "mpif.h"
#include "mpi_params.f"
      common/xygrid_location/xygrid_directory,xygrid_filename
      character*100 xygrid_directory,xygrid_filename
      integer nodeid,numprocs,ierr,stat(MPI_STATUS_SIZE)

      character aline*132,alinedim*1000,
     &          word1*16,word2*8,word3*7
      character*110 filename
      logical yes,form,mzone,iblank,first,threed
      logical DEBUG
      logical lvar(3)
      integer ivar(3)
      integer, pointer :: jdum(:)
      integer, pointer :: kdum(:)

      DEBUG = .false.
c      DEBUG = .true.

c-----
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      allocate (jdum(numprocs),kdum(numprocs))

      if (DEBUG) then
        print*
        print*, '>>> MODULE: rgrid.f/rgridh'
        print*
        if (nodeid.eq.0) then
          print*
          print*, ' >>> Node 0: reached MPI_BARRIER for rgridh 001'
        endif
      endif

c++++++++++++++++++++++++++++++++++++++++++++
      if (nodeid.eq.0) then
c++++++++++++++++++++++++++++++++++++++++++++

      write(istdout,5) 
5     format(/'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',
     &       /'                Grid Information (sub rgrid.f/rgridh)',/)

      print*
c      print*, ' stopping at rgrid.f/rgridh 50'
c      print*
c      call MPI_ABORT(MPI_COMM_WORLD,1)
c      stop ' stoppping at rgrid.f/rgridh 50'

c-----
c  Number of columns to try to process
c-----
      ncols = 132
      ncolsdim = 1000
c-----
c  name of grid file
c-----

c read original grid file
      filename = trim(xygrid_directory)//'/'//trim(xygrid_filename)
      print*
      print*, 'filename = ',filename
      print*, '  string length = ',len_trim(filename),'/110 characters'
      print*
c
      if (DEBUG) then
       print*, ' inquiring about grid file = ',filename
       print*
      endif

      inquire(file=filename,exist=yes)
      if(yes) then
         goto 30
      endif
15    continue
      write(istdout,20)
      write(istdout,25) (filename)
      write(*,*) ' '
      write(*,*) 'ERROR: Could not find grid file'
      write(*,*) '       grid file name: ',filename
      write(*,*) ' Program terminated.'
      call MPI_ABORT(MPI_COMM_WORLD,1)
      stop '-- STOP: no grid file found; stopped in rgrid.f'
20    format(' ERROR: No grid file found.  Please supply a 2D plot3d',/
     &       '        type grid file using the following name:')
25    format(10x,a)

c-----
c  Got file name.
c  Test formatted vs. unformatted by assuming it is formatted, and read
c  the first line as characters and check for invalid characters.
c-----
30    continue

      write(istdout,31) nodeid,trim(filename)
31    format('Node ',i4,': Opening file named "',a,'"')

      if (DEBUG) then
        print*
        print*, ' nodeid 0: opening formatted file = ',filename
        print*, '   ... iunit = ',iunit
        print*
      endif

c open file
      open(iunit,file=filename,form='formatted',status='old')

c check if formatted or unformatted
      form = .false.
      read(iunit,35,err=45) aline
35    format(a80)
      form = .true.
      do 40 ic=1,ncols
         if(aline(ic:ic).ne.' ' .and. aline(ic:ic).ne.',' .and.
     &     (aline(ic:ic).lt.'0' .or.  aline(ic:ic).gt.'9')) then
            form = .false.
            goto 45
         endif
40    continue
45    continue
      if(.not. form) then
        if (DEBUG) then
          print*
          print*, ' Closing iunit = ',iunit
          print*
        endif
        close(iunit)
      endif

      if (DEBUG) then
       print*
       print*, ' nodeid 0: form = ',form
       print*
      endif

c-----
c  Formatted file
c-----
      if(form) then
c
c  single or multiple grid file
c
         ints = 0
         first = .true.
         do 50 ic=1,ncols
            if(aline(ic:ic) .ne. ' ' .and.
     &         aline(ic:ic) .ne. ',' .and. first) then
               ints=ints+1
               first=.false.
            elseif(aline(ic:ic).eq.' '.or.aline(ic:ic).eq.',')then
               first=.true.
            endif
50       continue

         if (DEBUG) then
          print*
          print*, ' ints in sub rgrid.f/rgridh = ',ints
          print*
         endif
c
         if(ints.eq.1)then
            mzone=.true.
         elseif(ints.eq.2)then
            mzone=.false.
            threed=.false.
         elseif(ints.eq.3)then
            mzone=.false.
            threed=.true.
         else
            write(istdout,*)'ERROR: reading grid file.... stopping.'
            call exit(1)
         endif
c
c  2d vs. 3d for multi-zone file
c
         if(mzone)then
            rewind(iunit)
            read(iunit,*) numgrd
            if (DEBUG) then
             print*
             print*, ' nodeid 0: mzone = ',mzone,' numgrd = ',numgrd
             print*
            endif
c
            ints=0
55          continue
               intt=0
               first=.true.
               read(iunit,60) alinedim
60             format(a1000)
               do 65 ic=1,ncolsdim
                  if(alinedim(ic:ic).eq.'.')then
                     goto 70
                  elseif(alinedim(ic:ic).ne.' '.and.alinedim(ic:ic)
     &               .ne.','.and.first)then
                     intt=intt+1
                     first=.false.
                  elseif(alinedim(ic:ic).eq.' '.or.
     &                   alinedim(ic:ic).eq.',') then
                     first=.true.
                  endif
65             continue
               ints=ints+intt
            goto 55
70          continue
c
            threed=.true.
            if(ints/numgrd.eq.2)threed=.false.
         endif
c
c  must have a 2d file
c
         if(threed) then
           write(istdout,*) ' ERROR: input grid file appears to be a'
           write(istdout,*) '        three-dimensional file. Must use'
           write(istdout,*) '        a two-dimensional file. Stopping.'
           write(istdout,*) '        ints, numgrd: ',ints,numgrd
           print*
           print*, ' nodeid 0 has reached rgrid.f/rgridh 100'
           print*
           call MPI_ABORT(MPI_COMM_WORLD,1)
           stop ' stopping at rgrid.f/rgridh 100'
         endif
c
c  iblank vs. no iblank
c
         rewind(iunit)
         if(mzone)then
            read(iunit,*)numgrd
            if (DEBUG) then
             print*
             print*, ' nodeid 0: numgrd = ',numgrd
             print*
            endif
         else
            numgrd=1
         endif
c
         read(iunit,*)idim,jdim,(idum,idum,ng=1,numgrd-1)
         read(iunit,*) ((dum,i=1,idim),j=1,jdim),
     &                 ((dum,i=1,idim),j=1,jdim)
c
         iblank = .false.
         read(iunit,35,end=79)aline
         iblank=.true.
         do 77 ic=1,ncols
            if(aline(ic:ic).eq.'.')iblank=.false.
77       continue
79       continue

         if (DEBUG) then
          print*
          print*, ' nodeid 0: iblank = ',iblank
          print*
         endif

         if (iblank) then
           write(istdout,*) ' Grid file appears to have an iblank'
           write(istdout,*) '      array. MiniSMAC2D can only handle'
           write(istdout,*) '      patched grid problems for now.'
           write(istdout,*) ' '
           print*
           print*, ' stopping at rgrid.f/rgridh 110'
           print*
           call MPI_ABORT(MPI_COMM_WORLD,1)
           stop ' stopping at rgrid.f/rgridh 110'
        endif

c-----
c  Unformatted file
c-----

      else

         open(unit=iunit,file=filename,form='unformatted',status='old')
c
c  single or multiple grid file
c
         read(iunit,err=80)i,j,k
            mzone=.false.
            threed=.true.
            goto 130
80       rewind(iunit)
            read(iunit,err=90)i,j
            mzone=.false.
            threed=.false.
            goto 130
90       rewind(iunit)
            read(iunit,err=100)numgrd
            mzone=.true.
            goto 110
100      write(istdout,*)' ERROR: bad grid...  stopping'
         print*
         print*, ' stopping at rgrid.f/rgridh 120'
         print*
         call MPI_ABORT(MPI_COMM_WORLD,1)
         stop ' stopping at rgrid.f/rgridh 120'
c
c  2d or 3d
c
110      continue
         read(iunit,err=120)(i,j,k,ng=1,numgrd)
            threed=.true.
            goto 130
120      continue
         threed=.false.
130      continue
c
c  must have a 2d file
c
         if(threed) then
           write(istdout,*) ' ERROR: input grid file appears to be a'
           write(istdout,*) '        three-dimensional file. Must use'
           write(istdout,*) '        a two-dimensional file. Stopping.'
           print*
           print*, ' stopping at rgrid.f/rgridh 120'
           print*
           call MPI_ABORT(MPI_COMM_WORLD,1)
           stop ' stopping at rgrid.f/rgridh 120'
         endif
c
c  iblank vs. no iblank
c
         rewind(iunit)

         if(mzone) then
            read(iunit) numgrd
         else
            numgrd = 1
         endif
c
         iblank = .true.

         if(iblank) then
           write(istdout,*) ' ERROR: input grid is iblanked. This'
           write(istdout,*) '    is not allowed for smac2d. Stopping.'
           print*
           print*, ' stopping at rgrid.f/rgridh 130'
           print*
           call MPI_ABORT(MPI_COMM_WORLD,1)
           stop ' stopping at rgrid.f/rgridh 130'
          endif

         iblank = .false.
150      continue

c end of 'formatted or unformatted'
      endif

c-----
c  Read in dimensions
c-----

      close(iunit)

      if(form) then

c       if (DEBUG) then
         print*
         print*, ' nodeid 0: opening iunit = ',iunit,
     &      ' filename = ',filename
         print*, ' mzone = ',mzone
         print*
c       endif

         open(iunit,file=filename,form='formatted',status='old')
         print*, ' file is formatted'
         print*

         rewind(iunit)

         if(mzone) then
            read(iunit,*) nzone
         else
            nzone = 1
         endif
         print*, '   nzone = ',nzone
         print*
         print*, ' nodeid 0: reading jdum,kdum, nz=1,nzone'
         print*, '    ... nzone = ',nzone
         read(iunit,*) (jdum(nz),kdum(nz),nz=1,nzone)
         print*, '  grid file header has been read in sub. rgridh'
         print*

      else
         open(iunit,file=filename,form='unformatted',status='old')
         if(mzone) then
            read(iunit) nzone
         else
            nzone = 1
         endif
         read(iunit) (jdum(nz),kdum(nz),nz=1,nzone)
      endif


c-----
c  Perform parameter checks
c-----

      if(nzone .gt. numprocs) then
         write(istdout,160) nzone,numprocs
160      format(' ERROR: number of zones in grid file = ',i4,/
     &          '        is greater than parameter numprocs = ',i4,/
     &          '        Stopping.')
         call exit(1)
      endif
c
      imaxtot = 0
      jkmx = 0
      jmax= jdum(1)
      kmax= kdum(1)
      imaxtot = jmax*kmax
      jkmx = max(jkmx,jmax,kmax)

      if(jkmx .gt. jkmax) then
         write(istdout,190) jkmx,jkmax
         call exit(1)
      endif
190   format(' Error in code dimensions: jkmax too small:',/
     &       '    Largest dimension in input grid: ',i6,/
     &       '    Dimensioning parameter in code jkmax = ',i6)
c-----
c  Output grid info
c-----
      if(form) then
         word1 = 'a formatted, '
      else
         word1 = 'an unformatted, '
      endif

      if(mzone) then
         word2 = 'multiple'
      else
         word2 = 'single'
      endif

      if(iblank) then
         word3 = 'with'
      else
         word3 = 'without'
      endif

      write(istdout,200) word1,word2,word3
200   format(/'========================================',/
     & ' This is ',a16,a8,' zone grid file ',a7,
     & ' an iblank array.'/
     & '========================================'/)

      write(istdout,215)
      do 210 nz=1,nzone
         write(istdout,220) nodeid,nz,jdum(nz),kdum(nz)

210   continue

215   format('Original dimensions of grid:')
220   format('nodeid ',i4,': nz = ',i4,' jmax = ',i5,' kmax = ',i5)

      write(*,250)

c********************************************************

c++++++++++++++++++++++++++++++++++++++++++++
c end 'if (nodeid.eq.0)'
      endif
c++++++++++++++++++++++++++++++++++++++++++++
c return if don't want to send grid header to other nodes
      if (iselect.eq.0) return

c********************************************************
c send grid info to other nodes

c********************************************************

c+++++++++++++++++++++++++++++++++++++
      if (ng.gt.1) then
      if (nodeid.eq.0) then
c+++++++++++++++++++++++++++++++++++++
      print *, ' '
      print *, '========================================'

      lvar(1) = form
      lvar(2) = mzone
      lvar(3) = iblank
      write(*,222) nodeid, jmax,kmax,mzone
222   format(' Nodeid=',i5,', jmax=',i4,', kmax=',i4,', mzone=',l10)
      do 205 i=2,nzone
        node_out = i-1
        ivar(1) = jdum(i)
        ivar(2) = kdum(i)
        ivar(3) = numprocs
        if (DEBUG) then
          print*
          print*, ' ivar for sending to node = ',node_out
          print *, ivar(1),ivar(2),ivar(3)
          print*
          print*, ' lvar for sending to node = ',node_out
          print *, lvar(1),lvar(2),lvar(3)
        endif
        call MPI_SEND(ivar,3,MPI_INTEGER,node_out,node_out,
     &    MPI_COMM_WORLD,ierr)
        call MPI_SEND(lvar,3,MPI_LOGICAL,node_out,node_out+100,
     &    MPI_COMM_WORLD,ierr)
c        write(*,221) i,node_out,jdum(i),kdum(i),mzone
221   format(i3,'. nz=',i4,' jmax=',i5,' kmax=',i5,' nodeid=',i5)
205   continue

      print *, '========================================'
      print *, ' '
c
      write(istdout,250)
c
c      deallocate (jdum,kdum)

      write(istdout,230) nodeid,imaxtot,jkmx
230   format(' node ',i4,' total number of grid points = ',i7,/,
     &       ' largest single dimension    = ',i5)
      write(istdout,250)
250   format('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

c for all nodes other than node 0
c+++++++++++++++++++++++++++++++++++++
      else
c+++++++++++++++++++++++++++++++++++++

      if (DEBUG) then
        print *
        print *, ' NODEID ',nodeid,': CALL MPI_RECV'
        print *
      endif
      ivar(1) = 0
      ivar(2) = 0
      ivar(3) = 0
      call MPI_RECV(ivar,3,MPI_INTEGER,0,nodeid,
     &  MPI_COMM_WORLD,stat,ierr)
c unpack ivar
      jmax = ivar(1)
      kmax = ivar(2)
      numprocs = ivar(3)
      if (DEBUG) then
        print*
        print*, ' ivar received from node 0 for node = ',nodeid
        print*, '  ivar(1),ivar(2),ivar(3) = '
        print*, ivar(1),ivar(2),ivar(3)
      endif

      call MPI_RECV(lvar,3,MPI_LOGICAL,0,nodeid+100,
     &  MPI_COMM_WORLD,stat,ierr)
c unpack lvar
      form = lvar(1)
      mzone = lvar(2)
      iblank = lvar(3)
      if (DEBUG) then
        print*
        print*, ' lvar received from node 0 for node = ',nodeid
        print *, lvar(1),lvar(2),lvar(3)
      endif

      nzone = 1
      imaxtot = jmax*kmax
      jkmx = max(jmax,kmax)
      if(jkmx.gt.jkmax) then
        write(*,191) jkmx,jkmax
        call MPI_ABORT(MPI_COMM_WORLD,1)
        stop '2. jkmx > jkmax'
      endif
191   format(' Error in code dimensions: jkmax too small:',/
     & '    Largest dimension in input grid: ',i6,/
     & '    Dimensioning parameter in code jkmax = ',i6)

c+++++++++++++++++++++++++++++++++++++++
c end of 'if (nodeid.eq.0)' 
      endif
c end of 'if (nz.gt.1)'
      endif
c+++++++++++++++++++++++++++++++++++++++

      deallocate (jdum,kdum)

      if (DEBUG) then
       print*
       print *, ' -- Calling MPI_BARRIER at end of rgridh, node = ',
     &  nodeid
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      endif

c-----
c  End of rgridh
c-----
      return
      end
