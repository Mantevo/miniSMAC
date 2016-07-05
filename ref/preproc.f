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
c  1. readgo
c  2. addk2d
c  3. cube
c  4. partishn
c  5. average
c  6. break
c  7. patch_2d
c  8. min_surf_vol_ratio
c  9. subgrid_dimensions
c 10. grid_point_comparison
c 11. global_indices
c 12. write_subgrids
c 13. get_global_index
c 14. out1planetoplt3d_break
c
c************************************************************
c
c
c************************************************************************
      subroutine readgo(iunit,x,y,jdim,kdim,kadd,mzone)
c************************************************************************
#include "forttype.h"
c
c Purpose: This subroutine reads the original x,y grid values 
c
c Author: D. W. Barnette, SNL
c
      __REAL x(jdim,kdim),y(jdim,kdim)
      logical DEBUG
      logical mzone
c
c   ^ y,k
c   |
c   |
c   |
c   -----------> x,j
c where x-y is the geometry plane

      DEBUG = .false.
c      DEBUG = .true.

      if (DEBUG) then
       print*
       print*, '>> In subroutine readgo'
       print*
       print*, ' jdim,kdim,kadd = ',jdim,kdim,kadd
c       stop 'stop: at beginning of preproc.f/readgo'
      endif
c
c read 2d array into x and y
       rewind(iunit)
       if (mzone) read(iunit,*) idum
       read(iunit,*) idum1,idum2
       read(iunit,*)
     1  ((x(j,k),j=1,jdim),k=1,kdim-kadd),
     2  ((y(j,k),j=1,jdim),k=1,kdim-kadd)

      if (DEBUG) then
        print*
        print*, '=================================================='
        print*
        print*, ' In sub readgo, grid corner values are:'
        print*, ' jdim,kdim-kadd = ',
     1    jdim,kdim-kadd
        print*, ' 1. xy(1,1) = ',
     1    x(1,1),y(1,1)

        print*, ' 2. xy(jmx,1) = ',
     1    x(jdim,1),y(jdim,1)

        print*, ' 3. xy(1,kmx-kadd) = ',
     1    x(1,kdim-kadd),y(1,kdim-kadd)

        print*, ' 4. xy(1,1) = ',
     1    x(1,1),y(1,1)

        print*, ' 5. xy(jmx,kmx-kadd) = ',
     1    x(jdim,kdim-kadd),y(jdim,kdim-kadd)

        print*, ' 6. xy(jmx,kmx-kadd) = ',
     1    x(jdim,kdim-kadd),y(jdim,kdim-kadd)

        print*
        print*, '=================================================='
      endif
c
      return
      end  
c
c
c======================================================================|
      subroutine extend_grid_k(x,y,jdim,kdim,kadd)
c
c Purpose: This subroutine adds points in the radial or k direction to
c          arbitrarily increase the grid size for scaling studies.
c          The number of grid points in the j direction is unchanged.
c
      __REAL x(jdim,kdim),y(jdim,kdim) 
      logical DEBUG
c
      DEBUG = .false.
c      DEBUG = .true.

      do 10 j=1,jdim
c calculate delta-x at outer grid boundary
        delta_x = x(j,kdim-kadd) - x(j,kdim-kadd-1)
c calculate delta_y at outer grid boundary
        delta_y = y(j,kdim-kadd) - y(j,kdim-kadd-1)
c calculate slope delta_y/delta_x of outermost radial grid line
        slope = delta_y/delta_x
      do 10 k=1,kadd
c add grid points using straight line 2-point formula
        delta_x_total = k*delta_x
        x(j,kdim-kadd+k) = x(j,kdim-kadd) + delta_x_total
        y(j,kdim-kadd+k) = y(j,kdim-kadd) + slope*delta_x_total
10    continue

      return
      end
c
c
c-----------------------------------------------------------------
      subroutine cube(x,y,jdim,kdim)
c
c  Purpose: This subroutine writes out select grid points
c           for verification
c   Author: D. W. Barnette, Org. 1556, SNL
c  Written: 1-28-91
c Comments:
c
#include "forttype.h"
c
      __REAL x(jdim,kdim),y(jdim,kdim)
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.
c
c
      print *
      print *,'Block corner locations with/without surrounding points:'
      print *
      print *,'     Jmax, Kmax, = ',jdim, kdim
      print *
      print *,'          2 __________________ 3'
      print *,'           |                  |'
      print *,'           |                  |'
      print *,'           |                  |'
      print *,'       (K) ^                  |'
      print *,'           |                  |'
      print *,'           |                  |'
      print *,'         1 |------>-----------|4'
      print *,'                   (J)'
      print *
      print *,'  point   j     k         x          y     '
      print *,' ------ ----- -----  ---------- ----------'
c
c Outer block indices
      j0=1
      k0=1
      jdimm=jdim
      kdimm=kdim
c
      print *,'Corner grid points:'
      write(6,500) 1,j0,k0,x(j0,k0),y(j0,k0)
      write(6,500) 2,j0,kdimm,x(j0,kdimm),y(j0,kdimm)
      write(6,500) 3,jdimm,kdimm,x(jdimm,kdimm),y(jdimm,kdimm)
      write(6,500) 4,jdimm,k0,x(jdimm,k0),y(jdimm,k0)
c
c-------------------------
      if (DEBUG) then
c-------------------------
c write out more points only if desired
      print *
      print *,' POINT VALUES ALONG BLOCK EDGES'
      print *,'          j     k       x          y     '
      print *,' ------ ----- ----- ---------- ----------'
c
      print *
      print *,'All points from 1 to 4:'
      do 13 j=1,jdim
      do 13 k=1,1
      write(6,502) j,k,x(j,k),y(j,k)
13    continue
c
      print *
      print *,'All points from 1 to 2:'
      do 15 j=1,1
      do 15 k=1,kdim
      write(6,502) j,k,x(j,k),y(j,k)
15    continue
c
      print *
      print *,'All points from 4 to 3:'
      do 17 j=jdim,jdim
      do 17 k=1,kdim
      write(6,502) j,k,x(j,k),y(j,k)
17    continue
c
      print *
      print *,'All points from 2 to 3:'
      do 19 j=1,jdim
      do 19 k=kdim,kdim
      write(6,502) j,k,x(j,k),y(j,k)
19    continue

c-------------------------
      endif
c-------------------------
c
c
500   format(i3,5x,i5,1x,i5,1x,g10.3,1x,g10.3)
502   format(8x,i5,1x,i5,1x,g10.3,1x,g10.3)
c
      return
      end
c
c
c-----------------------------------------------------------------
      subroutine partishn(form,mzone,xorig,yorig,
     & jdim,kdim,idimj,idimk,jk,nodeid,numprocs,
     & j_subgrid,k_subgrid,jindex_global,kindex_global,ioverlap
     & )
c
c Purpose:
c     partitions original grid into user-specified (numprocs)
c     number of subgrids which are load balanced and optimized for maximum
c     volume and minimum surface area; subrids are overlapped for
c     subgrid-to-subgrid communication
c
c Output:
c     file bcpzn.dat which tells which subgrid communicates with which
c     other subgrid
c
#include "forttype.h"
c     
      __REAL xorig(jk),yorig(jk)
      dimension j_subgrid(numprocs),k_subgrid(numprocs)
      dimension jindex_global(numprocs),kindex_global(numprocs)
      dimension ipnum(numprocs)
      character*10 meshname
c
c
c calculate total grid points, average number of grid points per 
c  grid, etc.
      nzonetot = 1
      nzne = 1
      meshname='2D_mesh'
      call average(nzonetot,nzne,meshname,numprocs,jdim,kdim)
c
c breakup grids
      call break(nzonetot,meshname,nodeid,numprocs,form,
     1 ioverlap,ipnum(1),jdim,kdim,idimj,idimk,
     2 j_subgrid(1),k_subgrid(1),jk,
     3 xorig(1),yorig(1),
     4 jindex_global(1),kindex_global(1)
     5 )
c   
c link patched grids to subgrids
       call link_overlap(ioverlap,ipnum(1),nodeid,numprocs,
     1 jdim,kdim,idimj,idimk,
     2 j_subgrid(1),k_subgrid(1))
c

      write(6,*) '                                                     '
      write(6,*) '------------------------------------------------     '
      write(6,*) '  Subgrids output file: GRIDS2D.p3d (PLOT3D format)  '
      write(6,*) '     and boundary condition link file: 2D_PATCH_TABLE'
      write(6,*) '------------------------------------------------     '
      write(6,*) '                                                     '
      write(6,*) '   >> Grid partitioning finished <<                  '
      write(6,*) '                                                     '
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c NONE                                                                  |
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end

c======================================================================|
      subroutine average(nzone,ndim,meshname,numprocs,jdim,kdim)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c       This subroutine takes in all grids and computes average number
c       of grid points per subgrid.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c       D. W. Barnette
c       Sandia National Laboratories
c       Org. 9221, MS 1111
c       Albuquerque, NM 87112
c       email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c       breakup_for_fun
c       breakup_grids
c       chooser
c
c Called routines, in order of appearance:
c       NONE
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c       02-13-96        Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
#include "forttype.h"
c
      common/averages/isum,gpavg,igptotal
c
      character*10 meshname
c
c
      igptotal=jdim*kdim
      isum = igptotal
c
      gpavg=float(igptotal)/float(numprocs)
      write(6,100)
c
      write(6,110) meshname,jdim,kdim,igptotal
c
      write(6,115) isum,numprocs,gpavg
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
100   format( //,3x,'Grid',5x,'Jmax',3x,'Kmax',6x,'Total')              |
110   format(a10,2x,i4,3x,i4,3x,i8)                                     |
115   format(/,'Total number of grid points = ',i10,//,                 |
     1 'Total number of final zones (or subgrids) = ',i3,//,            |
     2 'Average number of grid points/grid = ',f12.2,//)                |
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c
c=======================================================================
      subroutine break(nzone,meshname,nodeid,numprocs,form,
     1 ioverlap,ipnum,jdim,kdim,idimj,idimk,
     2 j_subgrid,k_subgrid,jkmax,xorig,yorig,
     3 jindex_global,kindex_global
     4 )
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c       This subroutine breaks up grids into a user-specified number
c       of zones.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c       D. W. Barnette
c       Sandia National Laboratories
c       Org. 9221, MS 1111
c       Albuquerque, NM 87112
c       email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c       partishn 
c
c Called routines, in order of appearance:
c       min_surf_vol_ratio
c       subgrid_dimensions
c       grid_point_comparison
c       global_indices
c       get_global_index
c       write_subgrids
c       get_global_index
c       write_subgrids
c       out1planetoplt3d_break
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c       02-13-96        Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c
c   Guess at these; error catching will inform user if not correct
c      parameter(iicntr_max=maxdim*500,
c     1 ibcntr_max=iicntr_max)
c
#include "forttype.h"

      dimension j_subgrid(numprocs),k_subgrid(numprocs)
      dimension ipnum(numprocs)
      dimension jindex_global(numprocs),kindex_global(numprocs)
      __REAL xorig(jkmax),yorig(jkmax)

      common/averages/isum,gpavg,igptotal
      common/parallel/ippzone
c
      character*10 meshname
c number of subgrids (or cubes) formed in original grid 
      integer idimj,idimk
c
      logical form
c
c original grid dimensions
      jmax_compout = jdim
      kmax_compout = kdim

      write(6,5) numprocs
      print*, ' jkmax = ',jkmax
      write(6,*) '                                                     '
c Calculate grid points per processor node
      gppproc=isum/numprocs
c Print parameters
      do 41 i=1,nzone
      write(6,42) i,jmax_compout*kmax_compout
41    continue
      write(6,43) isum,gpavg,gppproc

c
c Calculate processors per zone
      isump=0
      do 15 i=1,nzone
       ippzone=nint(float(igptotal)/real(gppproc))
       if(ippzone.eq.0) ippzone=1
       isump=isump+ippzone
       write(6,25) i,ippzone
15    continue
      write(6,80) numprocs,isump

c
c Determine exact number jxkxl per processor for each zone.
c Optimize surface to volume ratio for efficient 'communications
c  to compute' ratio, given the number of processors for each zone.
c
      write(6,*) '                                                     '
      write(6,*) ' Calculate j*k values for subgrids for minimum     '
      write(6,*) '   area-to-volume ratio.                             '
c
c
c The following logic is applicable to dividing the grid into cubes.
c  This is optimal for load-balancing 3D grids.
c
c
       write(6,301) meshname
c
c Calculate minimum surface/volume ratio, given the number of desired 
c  processors
c   (determines idimj, idimk, etc.)
      call min_perimeter_to_area_ratio(1,ippzone,jmax_compout,
     1 kmax_compout,
     2 jpts_per_subgrid,kpts_per_subgrid,nzone,
     3 nodeid,numprocs,idimj,idimk,jdim,kdim)
c
c Divide the grid and determine dimensions
c  (determines j_subgrid, k_subgrid, etc.
       call subgrid_dimensions(1,jmax_compout,kmax_compout,
     1 jpts_per_subgrid,kpts_per_subgrid,
     2 nodeid,numprocs,idimj,idimk,
     3 j_subgrid(1),k_subgrid(1))
c
c Sum of grid points in the subdivided grids must add up to number in 
c  original grid before subdividing
       call grid_point_comparison(1,jmax_compout,kmax_compout,
     1 nodeid,numprocs,idimj,idimk,
     2 j_subgrid(1),k_subgrid(1))
c
300   continue
c
c assign processor numbers for all zones in the array ipnum()
      itemp=0
      do 302 j=1,idimj
      do 302 k=1,idimk
       index=(j-1)*idimk+(k-1)+1
       ipnum(index)=itemp+1
       itemp=ipnum(index)
302   continue
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Set global indices
      call global_indices(nzone,meshname,nodeid,numprocs,
     & j_subgrid(1),k_subgrid(1),idimj,idimk,
     2 jindex_global(1),kindex_global(1))
c
c establish new grids for each zone, based on new dimensions
c
      write(6,*) '                                                     '
      write(6,*) '                                                     '
      write(6,*) ' Breaking up grids into subgrids WITH NO overlap:    '
      write(6,295)
c
      iproc=0
c
       ii=0
       do 325 j=1,idimj
       do 325 k=1,idimk
        index=(j-1)*idimk+(k-1)+1
c
c set ioverlap
      ioverlap=0
c
c get global indices
        call get_global_index(1,j,k,1,1,
     1  jp3d_min,jp3d_max,kp3d_min,kp3d_max,
     2  ioverlap,nodeid,numprocs,
     3  j_subgrid(1),k_subgrid(1),
     4  idimj,idimk,jdim,kdim,
     5  jindex_global(1),kindex_global(1))
c
        itotal=j_subgrid(index)*k_subgrid(index)
c
c
c number the subgrids
        iproc=iproc+1
        ii=ii+1
c
        write(6,333) 1,ii,iproc,j_subgrid(index),
     1  k_subgrid(index),itotal,jp3d_min,
     2  jp3d_max,kp3d_min,kp3d_max
325    continue
c
       iproc_total=iproc+1
c
       if(1.lt.nzone) then
        write(6,*) '----------------------------------------------------
     1----------------------------'
       endif
c
320   continue
c
      iproc_total=iproc_total-1
c
c
c Write out Zone, Subzone, J-range, K-range, (relative values)
c   for writing grids to file GRIDS
c
      write(6,*) '                                                     '
      write(6,*) '                                                     '
      write(6,*) ' Breaking up grids into subgrids WITH overlap:       '
      write(6,295)
c
      iproc=0
c
      ioverlap=1
c
      ii=0
      do 425 j=1,idimj
      do 425 k=1,idimk
       index=(j-1)*idimk+(k-1)+1
c
c get global indices
       call get_global_index(1,j,k,1,1,jp3d_min,jp3d_max,
     1 kp3d_min,kp3d_max,ioverlap,nodeid,numprocs,
     2  j_subgrid(1),k_subgrid(1),idimj,idimk,
     3  jdim,kdim,
     4  jindex_global(1),kindex_global(1))


c j values
       if (j.lt.idimj) then
        j_subgrid(index)=j_subgrid(index)+ioverlap
       endif
c k values
       if (k.lt.idimk) then
        k_subgrid(index)=k_subgrid(index)+ioverlap
       endif
c
       itotal=(jp3d_max - jp3d_min + 1)*
     1        (kp3d_max - kp3d_min + 1)
c
c
c number the subgrids
       iproc=iproc+1
       ii=ii+1
c
       write(6,433) 1,ii,iproc,j_subgrid(index),
     1 k_subgrid(index),itotal,jp3d_min,
     2 jp3d_max,kp3d_min,kp3d_max
425     continue
c
       iproc_total=iproc+1
c
       if(1.lt.nzone) then
        write(6,*) '----------------------------------------------------
     1---------------------------'
       endif
c
420   continue

c
      iproc_total=iproc_total-1
c 
c+++++++++++++++++++++
c
c Write new subgrid gridpoint values x,y to one file in PLOT3D format
       write(6,304)
c
       call write_subgrids(nzone,form,nodeid,numprocs,
     1  ioverlap,idimj,idimk,jdim,kdim,
     2  j_subgrid(1),k_subgrid(1),
     3  ippzone,jkmax,
     4  xorig(1),yorig(1),
     5  jindex_global(1),kindex_global(1))
c
c
       write(6,*) '                                                    '
       write(6,303)
c+++++++++++++++++++++
c
c
c Now, do some output
c
      write(6,*) '                                                     '
      write(6,*) '       Output ONE plane of data for each subgrid for '
      write(6,*) '       two-dimensional representation. Each subgrid  '
      write(6,*) '       plane written to same file, >GRIDS2D.p3d<.    '
      write(6,*) '       This file will have multiblock 2-D PLOT3D     '
      write(6,*) '       format.                                       '
      write(6,*) '                                                     '
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
5     format(' >> Number of processors to be used: ',i6)                |
25    format(' For zone ',i3,' use ',i3,' processors.')                 |
40    format(//,' Parallel processing parameters:',//,                  |
     1'    Number of processors to be used: ',i9,/)                     |
42    format(/,' Number of grid points in zone ',i3,': ',i9)            |
43    format(/' Total number of grid points over all grids: ',i12,/,    |
     1'     Average number of grid points per grid: ',f12.2,/,          |
     2'Average number of grid points per processor: ',f12.2,/)          |
80    format(/,' User-specified processorss to be used = ',i5,/,        |
     1' Nodes calculated to be used = ',i5,/)                           |
100   format(//,' Is this ok? If not, you will be asked to enter',/     |
     110x,'a different number of processors. Enter (y/n): ',$)          |
200   format(' Zone = ',i3,', Total grid points = ',i9)                 |
205   format(/,' Warning: ',i3,' grids will not fit optimally on ',     |
     1i5,' processor(s)!',/)                                            |
210   format(' Optimal grid points per processor = ',f15.1,/)           |
295   format(/,' Zone',' Subzone',' Prc#','  Jmax','  Kmax',            |
     1'     Points','     J-rel range','     K-rel range')              |
301   format(//,'==========================================='           |
     1,/,' For mesh ',a,':',/)                                          |
303   format(/,' Subgrids written to file GRIDS2D.p3d in PLOT3D format.'|
     1 ,/)                                                              |
304   format(///,' Writing subgrids to file.',//)                       |
333   format(i5,i6,i6,i7,i6,i12,2x,i5,' to ',i5,2x,i5,' to ',i5)        |
433   format(i5,i6,i6,i7,i6,i12,2x,i5,' to ',i5,2x,i5,' to ',i5)        |
660   format('     Enter choice (1-5): ',$)                             |
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c
c======================================================================|
      subroutine link_overlap(ioverlap,ipnum,nodeid,numprocs,
     1 jdim,kdim,idimj,idimk,j_subgrid,k_subgrid)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c       This subroutine links overlapped (patched) subgrids in a zone.
c       This subroutine is not called for overlapped regions for which
c       interpolation coefficients have been generated. The subgrids
c       are assumed to align point-to-point.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c       D. W. Barnette
c       Sandia National Laboratories
c       Org. 9221, MS 1111
c       Albuquerque, NM 87112
c       email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c       chooser
c
c Called routines, in order of appearance:
c       write_unravel
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c       02-13-96        Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      dimension ipnum(numprocs)
      dimension j_subgrid(numprocs),k_subgrid(numprocs)

      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.
c
      if (DEBUG) then
       print*
       print*, ' >>> In preproc.f/link_overlap <<<'
       print*
       print*, ' ipnum(1), ipnum(2) = ',ipnum(1),ipnum(2)
       print*, ' idimj,idimk = ',idimj,idimk
       print*, ' numprocs = ',numprocs
       print*, ' ioverlap = ',ioverlap
      endif
c
c link overlapped portions of single grids
c
c Notes:
c       - Grids assumed to consist of right-handed coord. system
c       - Grids in this section have point-to-point matchup
c
c
      jkmax = jdim*kdim 

      write(6,*) '                                                     '
      write(6,*) ' >> Calculating LINK Table for adjacent sub-grids in e
     1ach zone.'
      write(6,*) '       Output file will be >bcpzn.dat<, which        '
      write(6,*) '        lists which subgrids are linked to which     '
      write(6,*) '        other neighboring subgrids.                  '
c
c
c write to file >3D_PATCH_TABLE<
c
c OUTPUT: open the >bcpzn.dat< output file used to store link info
      open(21,file='bcpzn.dat',iostat=ierr_unit21,form='formatted',
     1 status='unknown')
c     1 status='unknown')
      if(ierr_unit21.ne.0) then
       write(6,*) '                                                    '
       write(6,*) ' ERROR: Subroutine LINK_OVERLAP                     '
       write(6,*) '   File bcpzn.dat could not be opened.              '
       write(6,*) '                                                    '
       write(6,*) '   PROGRAM STOPPED: Problem with i/o file.          '
       write(6,*) '                                                    '
       stop '1: Subroutine LINK_OVERLAP '
      endif
c
c write header for standard output
      write(6,*) '                                                     '
      write(6,*) ' Print first few lines of file >bcpzn.dat<           '
      write(6,*) '                                                     '
c
c write heading
      write(6,540)
      write(21,540)
c
c start counter
      icnt = 0
c
c cycle thru cube indices (1 cube = 1 subgrid)
      do 20 j=1,idimj
      do 20 k=1,idimk
c
      index    =(j-1)*idimk+(k-1)+1
      index_jm1=(j-2)*idimk+(k-1)+1
      index_jp1=(j  )*idimk+(k-1)+1
      index_km1=(j-1)*idimk+(k-2)+1
      index_kp1=(j-1)*idimk+(k  )+1
c
c do j faces
c
c at j=1
      if(j.eq.1) then
       if(j.ne.idimj) then
        jtarget_min=j_subgrid(index)
        jtarget_max=jtarget_min
        ktarget_min=1
        ktarget_max=k_subgrid(index)
        jbase_min=1+ioverlap
        jbase_max=1+ioverlap
        kbase_min=1
        kbase_max=k_subgrid(index_jp1)
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,
     1                               ktarget_min,ktarget_max
        write(21,15) ipnum(index_jp1),jbase_min,jbase_max,
     1                                 kbase_min,kbase_max
        if(icnt.le.10) then
         write(6,15) ipnum(index),jtarget_min,jtarget_max,
     1                               ktarget_min,ktarget_max
         write(6,15) ipnum(index_jp1),jbase_min,jbase_max,
     1                                 kbase_min,kbase_max
        endif
       endif
c at j=jmax
      else if(j.eq.idimj)then
       if(j.ne.1) then
        jtarget_min=1
        jtarget_max=1
        ktarget_min=1
        ktarget_max=k_subgrid(index)
        jbase_min=j_subgrid(index_jm1)-ioverlap
        jbase_max=jbase_min
        kbase_min=1
        kbase_max=k_subgrid(index_jm1)
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,
     1                               ktarget_min,ktarget_max
        write(21,15) ipnum(index_jm1),jbase_min,jbase_max,
     1                                 kbase_min,kbase_max
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,
     1                               ktarget_min,ktarget_max
        write(6,15) ipnum(index_jm1),jbase_min,jbase_max,
     1                                 kbase_min,kbase_max
       endif
       endif
c everything in between j=1 and j=jmax
       else if((j.gt.1).and.(j.lt.idimj)) then
        jtarget_min=j_subgrid(index)
        jtarget_max=jtarget_min
        ktarget_min=1
        ktarget_max=k_subgrid(index)
        jbase_min=1+ioverlap
        jbase_max=1+ioverlap
        kbase_min=1
        kbase_max=k_subgrid(index_jp1)
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,
     1                               ktarget_min,ktarget_max
        write(21,15) ipnum(index_jp1),jbase_min,jbase_max,
     1                                 kbase_min,kbase_max
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,
     1                               ktarget_min,ktarget_max
        write(6,15) ipnum(index_jp1),jbase_min,jbase_max,
     1                                 kbase_min,kbase_max
       endif
        jtarget_min=1
        jtarget_max=1
        ktarget_min=1
        ktarget_max=k_subgrid(index)
        jbase_min=j_subgrid(index_jm1)-ioverlap
        jbase_max=jbase_min
        kbase_min=1
        kbase_max=k_subgrid(index_jm1)
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,
     1                               ktarget_min,ktarget_max
        write(21,15) ipnum(index_jm1),jbase_min,jbase_max,
     1                                 kbase_min,kbase_max
c
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,
     1                               ktarget_min,ktarget_max
        write(6,15) ipnum(index_jm1),jbase_min,jbase_max,
     1                                 kbase_min,kbase_max
       endif
      endif
c
c        
c do k faces
c
c at k=1
      if(k.eq.1) then
       if(k.ne.idimk) then
        jtarget_min=1
        jtarget_max=j_subgrid(index)
        ktarget_min=k_subgrid(index)
        ktarget_max=ktarget_min
        jbase_min=1
        jbase_max=j_subgrid(index_kp1)
        kbase_min=1+ioverlap
        kbase_max=1+ioverlap
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,
     1                               ktarget_min,ktarget_max
        write(21,15) ipnum(index_kp1),jbase_min,jbase_max,
     1                                 kbase_min,kbase_max
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,
     1                               ktarget_min,ktarget_max
        write(6,15) ipnum(index_kp1),jbase_min,jbase_max,
     1                                 kbase_min,kbase_max
       endif
       endif
c at k=kmax
      else if(k.eq.idimk)then
       if(k.ne.1) then
        jtarget_min=1
        jtarget_max=j_subgrid(index)
        ktarget_min=1
        ktarget_max=1
        jbase_min=1
        jbase_max=j_subgrid(index_km1)
        kbase_min=k_subgrid(index_km1)-ioverlap
        kbase_max=kbase_min
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,
     1                               ktarget_min,ktarget_max
        write(21,15) ipnum(index_km1),jbase_min,jbase_max,
     1                                 kbase_min,kbase_max
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,
     1                               ktarget_min,ktarget_max
        write(6,15) ipnum(index_km1),jbase_min,jbase_max,
     1                                 kbase_min,kbase_max
       endif
       endif
c everything in between k=1 and k=kmax
       else if((k.gt.1).and.(k.lt.idimk)) then
        jtarget_min=1
        jtarget_max=j_subgrid(index)
        ktarget_min=k_subgrid(index)
        ktarget_max=ktarget_min
        jbase_min=1
        jbase_max=j_subgrid(index_kp1)
        kbase_min=1+ioverlap
        kbase_max=kbase_min
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,
     1                               ktarget_min,ktarget_max
        write(21,15) ipnum(index_kp1),jbase_min,jbase_max,
     1                                 kbase_min,kbase_max
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,
     1                               ktarget_min,ktarget_max
        write(6,15) ipnum(index_kp1),jbase_min,jbase_max,
     1                                 kbase_min,kbase_max
       endif
        jtarget_min=1
        jtarget_max=j_subgrid(index)
        ktarget_min=1
        ktarget_max=ktarget_min
        jbase_min=1
        jbase_max=j_subgrid(index_km1)
        kbase_min=k_subgrid(index_km1)-ioverlap
        kbase_max=kbase_min
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,
     1                               ktarget_min,ktarget_max
        write(21,15) ipnum(index_km1),jbase_min,jbase_max,
     1                                 kbase_min,kbase_max
c
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,
     1                               ktarget_min,ktarget_max
        write(6,15) ipnum(index_km1),jbase_min,jbase_max,
     1                                 kbase_min,kbase_max
       endif
      endif
c
c        
20    continue
c
c        
      write(6,*) '                                                     '
      write(6,*) '  Link table finished. File >bcpzn.dat< written.'
      write(6,*) '  Number of points written to file: ',icnt

c
c close bcpzn.dat file
       close(21)
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
15    format(2x,i4,2x,i5,2x,i5,1x,i5,1x,i5)                             |
540   format('c This is the link table for overlapped subgrids.',/,     |
     1'c First line of each pair = target points; second, base points.',|
     2 /,                                                               |
     1'c Zone#  J_beg J_end K_beg K_end ')                              |
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end

c======================================================================|
      subroutine min_perimeter_to_area_ratio
     1 (nz,ippzone,jmaxx,kmaxx,jpts_per_subgrid,kpts_per_subgrid,
     2  nzone,nodeid,numprocs,idimj,idimk,jdim,kdim)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c 	This subroutine calculates the 3-factors for a given number of
c   	subgrids, then calculates the surface to volume ratio for each 
c	set of 3-factors given the maximum j,k,l dimensions of each grid. 
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c	D. W. Barnette
c   	Sandia National Laboratories
c   	Org. 9221, MS 1111
c     	Albuquerque, NM 87112
c  	email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c	break
c	breakup_for_fun
c
c Called routines, in order of appearance:
c	NONE
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c	02-13-96	Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c
      parameter (imax=300)
c      
      dimension idim1(imax),idim2(imax)  
c
      jkmax=jdim*kdim
      max_subs = numprocs
c
      write(6,*) ' List all possible factors for number of subgrids    '
      write(6,*) '       for this zone.                                '
      write(6,5)
      fnumprocs=real(numprocs)     
      itag=0
      do 10 i=1,numprocs
       xtry1=fnumprocs/float(i)
       check1=xtry1-int(xtry1)
       if(check1.eq.0) then   
        do 15 j=1,int(xtry1)
         xtry2=xtry1/float(j)
         check2=xtry2-int(xtry2)
         check_answer=i*j
         if ((check2.eq.0.0).and.(check_answer.eq.numprocs)) then
          itag=itag+1
          if (itag.gt.imax) then
             write(6,*) '                                              '
             write(6,*) ' ERROR: Subroutine MIN_SURF_VOL_RATIO         '
             write(6,*) '    ITAG is greater than IMAX!                '
             write(6,*) '      ITAG = ',itag
             write(6,*) '      IMAX = ',imax
             write(6,*) '        Change parameter statement!           '
             write(6,*) '                                              '
             write(6,*) '  PROGRAM STOPPED.                            '
             write(6,*) '                                              '
             stop 'stopping - 1: Subroutine MIN_SURF_VOL_RATIO '
          endif
          idim1(itag) = i
          idim2(itag) = j
          write(6,25) itag,i,j
         endif
15      continue
        endif
10     continue
c
c calculate the min(perim/area)
      write(6,*) '                                                     ' 
      write(6,*) '                                                     '
      write(6,*) ' Finding min(perimeter/area) ratio for:              '
      write(6,*) '      Jmax = ',jmaxx
      write(6,*) '      Kmax = ',kmaxx
c 
c write header for output
      write(6,30)
      do 35 i=1,itag
      perim_over_area = 2*(float(idim1(i))/float(jmaxx)
     1  + float(idim2(i))/float(kmaxx) ) 
c     1  + float(idim2(i))/float(kmaxx) + float(idim3(i))/float(lmaxx)) 
      if(i.gt.1) then
       xmin_save=xmin_perim_over_area
       xmin_perim_over_area = amin1(xmin_perim_over_area, 
     1  perim_over_area) 
       if (xmin_perim_over_area.ne.xmin_save) then 
        imin=i
       endif    
      else
       xmin_perim_over_area = perim_over_area
       imin=1    
      endif  
      x1=perim_over_area
      x2=xmin_perim_over_area
      write(6,40) i,x1,x2
35    continue 
      idimj=idim1(imin)
      idimk=idim2(imin)
c
c error check
      idim_check=idimj*idimk
      if(idim_check.gt.max_subs) then
       imax_zone=0
       do 200 i=1,nzone
        imax_zone=max(ippzone,imax_zone)
200    continue
       write(6,*) '                                                    '
       write(6,*) ' ERROR: Subroutine MIN_SURF_VOL_RATIO               '
       write(6,*) '   More subgrids are needed than allowed for in     '
       write(6,*) '   parameter statement.                             '
       write(6,100) nz,idim_check,max_subs,imax_zone,imax_zone
       write(6,*) '                                                    '
       write(6,*) '   PROGRAM STOPPED.                                 '
       write(6,*) '                                                    '
       stop 'stopping - 2: Subroutine MIN_SURF_VOL_RATIO '
      endif
c
c output grid values
c   jpts_per_subgrid, etc: approx number of grid points per subgrid

      jpts_per_subgrid=(jmaxx-1)/idimj + 1
      kpts_per_subgrid=(kmaxx-1)/idimk + 1
      write(6,*) '                                                     '
      write(6,*) ' Final values for min(perim/area) for this zone are: '
      write(6,*) '    itag = ',imin
      write(6,50) xmin_perim_over_area
      write(6,60) idimj,jmaxx,jpts_per_subgrid
      write(6,70) idimk,kmaxx,kpts_per_subgrid
      write(6,*) '                                                     '
c
c  
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
5     format(/,' ITAG',7x,'J',7x,'K',/                                  |
     1         ' ----',6x,'---',5x,'---')                               |
25    format(i5,6x,i3,5x,i3)                                            |
30    format(/,' itag',7x,'perim/area',5x,'min(perim/area)')            |
40    format(i4,5x,f11.6,5x,f11.6)                                      |
50    format(5x,'min(perim/area) = ',f11.6)                             |
60    format(5x,'N/Jmax = ',i4,' /',i4,' ==> approx. ',i10,' pts/subgrid|
     1 in J direction')                                                 |
70    format(5x,'M/Kmax = ',i4,' /',i4,' ==> approx. ',i10,' pts/subgrid|
     1 in K direction')                                                 |
100   format(/,                                                         |
     1'              Subgrids needed for zone ',i5,': ',i5,/            |
     2'       Number of subgrids currently allowed : ',i5,/             |
     3'              Max number of subgrids needed : ',i5,//            |
     3'   Change max_subs in preproc.f/min_surf_vol_ratio to ',i5/)     |
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c end 'min_surf_vol_ratio'


c======================================================================|
      subroutine subgrid_dimensions(nz,jmaxx,kmaxx,
     1 jpts_per_subgrid,kpts_per_subgrid,
     2 nodeid,numprocs,idimj,idimk,j_subgrid,k_subgrid)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c 	This subroutine calculates the dimensions of the subgrids.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c	D. W. Barnette
c   	Sandia National Laboratories
c   	Org. 9221, MS 1111
c     	Albuquerque, NM 87112
c  	email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c	break
c
c Called routines, in order of appearance:
c	NONE
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c	02-13-96	Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c      
      dimension j_subgrid(numprocs),k_subgrid(numprocs)
c
c check if grids are evenly divisible or if remainder is left over
      j_remainder = mod(jmaxx-1,idimj)
      k_remainder = mod(kmaxx-1,idimk)
      if(j_remainder.ne.0) then
       write(6,*) ' Warning: Jmaxx/N not evenly divisible!'
       write(6,5) j_remainder 
      endif          
      if(k_remainder.ne.0) then
       write(6,*) ' Warning: Kmaxx/M not evenly divisible!'
       write(6,10) k_remainder
      endif          
c
      if((j_remainder.ne.0).or.(k_remainder.ne.0)) then
       write(6,*) '                                                    '
       write(6,*) '  BREAKUP will make adjustments in the size of the  '
       write(6,*) '    subgrids to compensate for the above remainders.'
       write(6,*) '    Adjustments will start at the boundaries and    '
       write(6,*) '    progress inward. Hence, when remainders appear, '
       write(6,*) '    outer subgrids will have their dimensions       '
       write(6,*) '    increased in the corresponding direction.       '
      endif
c 
c fill in approximate values      
      do 20 j=1,idimj
      do 20 k=1,idimk
       index=(j-1)*idimk+(k-1)+1
       j_subgrid(index)=jpts_per_subgrid
       k_subgrid(index)=kpts_per_subgrid
20    continue
c now include remainders
      j_pool=j_remainder
      k_pool=k_remainder
c work on j values 
      if(j_remainder.ne.0) then
       do 30 k=1,idimk
        index_j1=(1-1)*idimk+(k-1)+1
        j_subgrid(index_j1)=j_subgrid(index_j1)+1  
30     continue
       j_pool=j_pool-1
       if(j_pool.ne.0) then
        do 35 k=1,idimk
         index_jmax=(idimj-1)*idimk
     1   +(k-1)+1
         j_subgrid(index_jmax)=j_subgrid(index_jmax)+1
35      continue
        j_pool=j_pool-1
        if(j_pool.ne.0) then
         jmaxcp1=idimj+1
         jmaxcp1_half=jmaxcp1/2
         do 40 j=2,jmaxcp1_half
          do 45 k=1,idimk
           index=(j-1)*idimk+(k-1)+1
           j_subgrid(index)=j_subgrid(index)+1
45        continue
          j_pool=j_pool-1
          if(j_pool.eq.0) goto 55
          jtemp=idimj-j+1
          do 50 k=1,idimk
           index_jtemp=(jtemp-1)*idimk
     1     +(k-1)+1
           j_subgrid(index_jtemp)=j_subgrid(index_jtemp)+1  
50        continue
          j_pool=j_pool-1
          if(j_pool.eq.0) goto 55 
40       continue 
55       continue
        endif
       endif
      endif 
c      
c work on k values
      if(k_remainder.ne.0) then
       do 60 j=1,idimj
        index_k1=(j-1)*idimk+(1-1)+1
        k_subgrid(index_k1)=k_subgrid(index_k1)+1  
60     continue
       k_pool=k_pool-1
       if(k_pool.ne.0) then
        do 65 j=1,idimj
         index_kmax=(j-1)*idimk
     1   +(idimk-1)+1
         k_subgrid(index_kmax)=k_subgrid(index_kmax)+1
65      continue
        k_pool=k_pool-1
        if(k_pool.ne.0) then
         kmaxcp1=idimk+1
         kmaxcp1_half=kmaxcp1/2
         do 70 k=2,kmaxcp1_half
          do 75 j=1,idimj
           index=(j-1)*idimk+(k-1)+1
           k_subgrid(index)=k_subgrid(index)+1
75        continue
          k_pool=k_pool-1
          if(k_pool.eq.0) goto 85
          ktemp=idimk-k+1
          do 80 j=1,idimj
           index_ktemp=(j-1)*idimk
     1     +(ktemp-1)+1
           k_subgrid(index_ktemp)=k_subgrid(index_ktemp)+1  
80        continue
          k_pool=k_pool-1
          if(k_pool.eq.0) goto 85 
70       continue 
85       continue
        endif
       endif
      endif 
c
c write out final dimensions for subgrids 
      write(6,*) '                                                     '
      write(6,*) ' Final subgrid dimensions are:                       '
c
c print out final grid point values 
      write(6,120)
      inum=0
      do 125 j=1,idimj
      do 125 k=1,idimk
       index=(j-1)*idimk+(k-1)+1
       inum=inum+1
       itotl=j_subgrid(index)*k_subgrid(index)
       write(6,130) inum,j,k,j_subgrid(index),k_subgrid(index),
     1 itotl
125   continue
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
5     format(5x,'J-remainder = ',i3)                                    |
10    format(5x,'K-remainder = ',i3)                                    |
15    format(5x,'L-remainder = ',i3)                                    |
120   format(/,' Grid #',4x,'Jcube',4x,'Kcube',4x,'Jmax',5x,            |
     1 'Kmax',7x,'Total',/,                                             |
     2         '-------',4x,'-----',4x,'-----',4x,'----',5x,            |
     3 '----',5x,'---------')                                           |
130   format(i6,5x,i4,5x,i4,5x,i4,5x,i4,2x,i12)                         |
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c end 'subgrid_dimensions'


c======================================================================|
      subroutine grid_point_comparison(nzone,jmaxx,kmaxx,
     1 nodeid,numprocs,idimj,idimk,j_subgrid,k_subgrid)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c 	This subroutine calculates the total number of points for the 
c	original grid, and the total number of grid points for the 
c	subgrids of the original grid, and compares the two. The program 
c	stops if there is a difference. This is to provide some error 
c	checking.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c	D. W. Barnette
c   	Sandia National Laboratories
c   	Org. 9221, MS 1111
c     	Albuquerque, NM 87112
c  	email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c	break
c	breakup_for_fun
c
c Called routines, in order of appearance:
c	NONE
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c	02-13-96	Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c      
c      
      dimension j_subgrid(numprocs),k_subgrid(numprocs)
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.
c 
c calculate total number of grid points and compare with original
      write(6,*) '                                                     '
      write(6,*) ' --- SANITY CHECK ---'
      write(6,*) ' Compare total number of grid points with original   '
      write(6,*) '     in this zone:                                   '
      ipoints_orig=jmaxx*kmaxx 

      if (DEBUG) then
        print*, ' In preproc.f/grid_point_comparison:'
        print*, ' idimj = ',idimj
        print*, ' idimk = ',idimk
        print*, '         i    j_subgrid   k_subgrid'
        print*, '    ------    ---------   ---------'
        do 5 i=1,numprocs
         print*,i,j_subgrid(i),k_subgrid(i)
5       continue
      endif
c      
c First, add all subgrids together
      ipoints_calc=0
      do 10 j=1,idimj
      do 10 k=1,idimk
       index=(j-1)*idimk+(k-1)+1
       ipoints_calc=ipoints_calc+(j_subgrid(index)
     1  *k_subgrid(index))
10    continue
c
c Now, subtract common faces
c   j-direction 
      jflag=0
      do 20 j=1,idimj-1
      jflag=1
      do 20 k=1,idimk
       index=(j-1)*idimk+(k-1)+1
       ipoints_calc=ipoints_calc-k_subgrid(index)
20    continue 
      if(jflag.eq.1) then 
       write(6,*) '   Subtracted common faces in J-direction.          ' 
      endif
c   k-direction
      kflag=0 
      do 30 k=1,idimk-1
      kflag=1
      do 30 j=1,idimj
       index=(j-1)*idimk+(k-1)+1
       if(j.ne.idimj) then 
        ipoints_calc=ipoints_calc
     1  -(j_subgrid(index)-jflag)
       else
        ipoints_calc=ipoints_calc
     1  -j_subgrid(index)
       endif
30    continue
      if(kflag.eq.1) then 
       write(6,*) '   Subtracted common faces in K-direction.          ' 
      endif
c      
      write(6,40) ipoints_orig,ipoints_calc
c     
      if(ipoints_orig.ne.ipoints_calc) then
       write(6,*) '                                                    '
       write(6,*) ' ERROR: Subroutine GRID_POINT_COMPARISON            '
       write(6,*) '   Original number of grid points does not          '
       write(6,*) '   equal the number summed from the subgrids.       '
       write(6,*) '                                                    '
       write(6,*) '  PROGRAM STOPPED.                                  '
       write(6,*) '                                                    '
       stop '1: Subroutine GRID_POINT_COMPARISON '
      endif
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
40    format(5x,'...Original no. of points in this zone = ',i10,/,      |
     15x,'...Calculated no. of points in this zone = ',i10)             |
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c end 'grid_point_comparison'


c======================================================================|
      subroutine global_indices(nzone,meshname,nodeid,numprocs,
     1 j_subgrid,k_subgrid,idimj,idimk,
     2 jindex_global,kindex_global)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c 	This subroutine sets up the global indices needed to break up 
c	the zone. Call this routine after calling subroutine 
c	subgrid_dimensions.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c	D. W. Barnette
c   	Sandia National Laboratories
c   	Org. 9221, MS 1111
c     	Albuquerque, NM 87112
c  	email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c	break
c
c Called routines, in order of appearance:
c	NONE
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c	02-13-96	Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      dimension j_subgrid(numprocs),k_subgrid(numprocs)
      dimension jindex_global(numprocs),kindex_global(numprocs)
c
      integer idimj,idimk
      character*10 meshname
c
      nz = nzone
      max_subs = numprocs
c
      write(6,6)
c
c do checks on dimensions
      idim_max=0
c
      write(6,12)
c
c write original zone dimensions
       write(6,13) nz,idimj,idimk
       idim_max=max(idim_max,idimj*idimk)
c
      if(idim_max.gt.max_subs) then
       write(6,1) idim_max,max_subs
       write(6,14) max_subs,idim_max
       write(6,*) '                                                    '
       write(6,*) '  PROGRAM STOPPED.                                  '
       write(6,*) '                                                    '
       stop '1: Subroutine GLOBAL_INDICES '
      endif
c                          
       do 15 j=1,idimj
       do 15 k=1,idimk
         index=(j-1)*idimk+(k-1)+1
         index_km1=(j-1)*idimk+(k-2)+1
        if(k.eq.1) then
         kindex_global(index)=k_subgrid(index)
        else
         kindex_global(index)=k_subgrid(index)
     1                   +kindex_global(index_km1)-1
        endif
15     continue                                               
c                          
       do 20 k=1,idimk
       do 20 j=1,idimj
       index=(j-1)*idimk+(k-1)+1
       index_jm1=(j-2)*idimk+(k-1)+1
        if(j.eq.1) then
         jindex_global(index)=j_subgrid(index)
        else
         jindex_global(index)=j_subgrid(index)
     1                   +jindex_global(index_jm1)-1
        endif
20     continue 
c
       write(6,*) '                                                    '
       write(6,*) '                                                    '
       write(6,21) nzone,meshname
       write(6,*)  ' Grd# Jcube Kcube Jmax Kmax  Jming Jmaxg K
     1ming Kmaxg'
       write(6,*)  ' ---- ----- ----- ---- ----  ----- ----- -
     1---- -----'
c      
c j,k are cube dimensions in original grid
c 
       icount=0      
       do 25 j=1,idimj
       do 25 k=1,idimk
        index=(j-1)*idimk+(k-1)+1
        icount=icount+1
        jp3d_max=jindex_global(index)
        kp3d_max=kindex_global(index)
        jp3d_min=jp3d_max-j_subgrid(index)+1
        kp3d_min=kp3d_max-k_subgrid(index)+1 
        jlocal=j_subgrid(index)
        klocal=k_subgrid(index)
        write(6,26) icount,j,k,jlocal,klocal,jp3d_min,jp3d_max,
     1  kp3d_min,kp3d_max
25     continue
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
1     format(/,' ERROR: Subroutine GLOBAL_INDICES',/                    |
     1'     Array dimensions exceeded',/                                |
     2'         idim_max = ',i4,/                                       |
     3'         max_subs = ',i4,/                                       |
     4'     Adjust parameter statement for max_subs.')                  |
6     format(//,' -----------------------------------------',/          |
     1' Local and Global Indices for Each Subgrid',/                    |
     2' -----------------------------------------')                     |
12    format(//' Summary of zone dimensions:',//                        |
     1'  Zone   jcube_max  kcube_max',/                                 |
     2' ------  ---------  ---------')                                  |
13    format(i6,i10,2i11)                                               |
14    format(//,' > Parameter statement needs to be updated:',/         |
     1'       Current values: max_subs = ',i5,/                         |
     2' Suggest changing to : max_subs = ',i5,/)                        |
21    format(' Zone =',i4,', Meshname = ',a)                            |
26    format(i5,8i6)                                                    |
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c end 'global_indices'

c======================================================================|
      subroutine write_grid_for_1_processor(jmax,kmax,x,y)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c       This subroutine writes out the grid in PLOT3D format in 2D
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author:
c       D. W. Barnette
c       Sandia National Laboratories
c       Albuquerque, NM 87112
c       email: dwbarne@sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
#include "forttype.h"
c
      __REAL x(jmax*kmax),y(jmax*kmax)
c
      logical form
      logical DEBUG
c
      DEBUG = .false.
c      DEBUG = .true.
c
      if (DEBUG) then
        print*
        print*, ' >>> In preproc.f/write_grid_for_1_processor <<<'
        print*
      endif
c
      write(6,*) ' '
c
       open(75,file='GRIDS2D.p3d',form='formatted',
     1   iostat=ierr_unit75,status='unknown')
      if(ierr_unit75.ne.0) then
       write(6,*) '                                                    '
       write(6,*) ' ERROR: Subroutine WRITE_SUBGRIDS                   '
       write(6,*) '   File GRIDS2D.p3d could not be opened.            '
       write(6,*) '                                                    '
       write(6,*) '                                                    '
       write(6,*) '   PROGRAM STOPPED: Problem with i/o file.          '
       write(6,*) '                                                    '
       stop '1: Subroutine WRITE_SUBGRIDS '
      endif
c
      write(6,*) 'Output file GRIDS2D.p3d, 1 proc, opened as formatted.'
      write(6,*) ' '
c
      write(75,*) 1 
c      stop 'stopped'

c write jmax,kmax 
       write(75,8) jmax,kmax

c formatted
       write(75,1900) 
     1      ((x((kp-1)*jmax + jp),
     2                         jp=1,jmax),
     3                         kp=1,kmax),
     5      ((y((kp-1)*jmax + jp),
     6                         jp=1,jmax),
     7                         kp=1,kmax)
c
c
      close(75)
      num_grid_points = jmax*kmax 

      print*
      print*, ' Total number of grid points in GRIDS2D.p3d: ',
     1 num_grid_points
      print*
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
8     format(2i10)                                                      |
1900  format(3f20.8)                                                    |
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c======================================================================|
      subroutine write_subgrids(nzone,form,nodeid,numprocs,
     1 ioverlap,idimj,idimk,jdim,kdim,
     2 j_subgrid,k_subgrid,
     3 ippzone,jkmax,
     4 x,y,
     5 jindex_global,kindex_global)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c 	This subroutine writes out the subgrids formed in subroutine
c 	SUBGRID_DIMENSIONS and subroutine GLOBAL_INDICES.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c	D. W. Barnette
c   	Sandia National Laboratories
c   	Org. 9221, MS 1111
c     	Albuquerque, NM 87112
c  	email: dwbarne@sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Called routines, in order of appearance:
c	get_global_index
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c	02-13-96	Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c
#include "forttype.h"

      dimension j_subgrid(numprocs),k_subgrid(numprocs)
      dimension jindex_global(numprocs),kindex_global(numprocs)
c
      __REAL x(jkmax),y(jkmax)
c
      logical form
      logical DEBUG
c
      DEBUG = .false.
c      DEBUG = .true.
c
      if (DEBUG) then
        print*
        print*, ' >>> In preproc.f/write_subgrids <<<'
        print*, '  ioverlap = ',ioverlap
        print*, '  jkmax = ',jkmax
        print*
      endif
c
      jmax_compout=jdim
      kmax_compout=kdim

c output is formatted
c
      write(6,*) ' '
c
       open(75,file='GRIDS2D.p3d',form='formatted',
     1   iostat=ierr_unit75,status='unknown')
      if(ierr_unit75.ne.0) then
       write(6,*) '                                                    '
       write(6,*) ' ERROR: Subroutine WRITE_SUBGRIDS                   '
       write(6,*) '   File GRIDS2D.p3d could not be opened.            '
       write(6,*) '                                                    '
       write(6,*) '                                                    '
       write(6,*) '   PROGRAM STOPPED: Problem with i/o file.          '
       write(6,*) '                                                    '
       stop '1: Subroutine WRITE_SUBGRIDS '
      endif
c
      write(6,*) 'Output file GRIDS2D.p3d opened as formatted.         '
c
      write(6,*) '                                                     '
c
c
c write out plot file for subgrids in PLOT3D format
c
c write header

c
      if (DEBUG) then
        print*
        print*, ' number_of_processors = ',numprocs
        print*, ' ioverlap = ',ioverlap
        print*, ' jkmax = ',jkmax
        print*, ' jmax_compout = ',jmax_compout
        print*, ' kmax_compout = ',kmax_compout
        print*, ' idimj,idimk = ',idimj,idimk
        print*, ' nzone = ',nzone
      endif

      write(75,*) numprocs
      print*,'j_subgrid = ',j_subgrid
      print*,'k_subgrid = ',k_subgrid
c      stop 'stopped'
c
       write(75,8) 
     1((j_subgrid((j-1)*idimk+(k-1)+1),
     2   k_subgrid((j-1)*idimk+(k-1)+1),
     4                k=1,idimk),j=1,idimj) 
c
c
c Extract each subgrid from the original PLOT3D file
c
      i_total_count=0
c
      nz = 1
c
      jmaxc=jmax_compout
      kmaxc=kmax_compout

      write(6,*) ' '
      write(6,*) ' jmaxc and jmax_compout = ',jmaxc
      write(6,*) ' kmaxc and kmax_compout = ',kmaxc
      write(6,*) ' '
c
c read zone grid points from original grid file
c
      write(6,800) nz
      write(6,802)
c
1200  continue
c
c construct subgrids for the zone just read in
c
      itotal=idimj*idimk
c
      print*
      print*, ' Writing to GRIDS2D.p3d'
      print*
c
c loop over sub-blocks of grids
      num_grid_points = 0

      do 30 j=1,idimj
      do 30 k=1,idimk
c
       i_total_count=i_total_count+1
c
       index_ref=(j-1)*idimk + (k-1) + 1 
c
       write(6,705) index_ref,itotal,nz,i_total_count,numprocs,ioverlap
c
c
      if (DEBUG) then
       print*
       print*, ' >>> Calling get_global_index'
       print*
      endif

       call get_global_index
     1 (nz,j,k,1,1,
     2 jp3d_min,jp3d_max,kp3d_min,kp3d_max,
     3 ioverlap,nodeid,numprocs,
     4 j_subgrid(1),k_subgrid(1),
     5 idimj,idimk,jdim,kdim,
     6 jindex_global(1),kindex_global(1))
c
       write(6,900) nz,j,k,jp3d_min,jp3d_max,kp3d_min,
     1 kp3d_max
       write(6,901) jp3d_max - jp3d_min + 1, kp3d_max - kp3d_min + 1
c
c formatted
       write(75,1900) 
     1      ((x((kp-1)*jmaxc + jp),
     2                         jp=jp3d_min,jp3d_max),
     3                         kp=kp3d_min,kp3d_max),
     5      ((y((kp-1)*jmaxc + jp),
     6                         jp=jp3d_min,jp3d_max),
     7                         kp=kp3d_min,kp3d_max)
c
c
      if (DEBUG) then
      write(6,1505) j,k,jmaxc,
     1 x((kp3d_min-1)*jmaxc + jp3d_min),
     2 y((kp3d_min-1)*jmaxc + jp3d_min),
     3 x((kp3d_max-1)*jmaxc + jp3d_max),
     4 y((kp3d_max-1)*jmaxc + jp3d_max)
1505  format(/'min-max values in write_subgrids',/,
     1 '  idimj = ',i3,' idimk = ',i3,' jmaxc = ',i5,/,
     2 '   xmin = ',f20.8,' ymin = ',f20.8,/,
     3 '   xmax = ',f20.8,' ymax = ',f20.8
     4 )
      endif

      num_grid_points = num_grid_points + 
     1  (jp3d_max - jp3d_min + 1)*(kp3d_max - kp3d_min + 1)

30    continue  

      close(75)

      print*
      print*, ' Total number of grid points in GRIDS2D.p3d: ',
     1 num_grid_points
      print*

      if (DEBUG) then
       print*
c       print*, ' stop at end of preproc.f/write_subgrids'
c       stop ' stop at end of preproc.f/write_subgrids'
      endif 
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
8     format(2i10)                                                      |
705   format(                                                           |
     1'> Constructing subgrid #',i5,' of',i5,' subgrids in zone ',i4,/  |
     2'   Total count: subgrid #',i5,' of ',i5,' subgrids',/            |
     3'   Grid cells overlap',i2,' cell widths on grid boundaries.')    |
800   format(/,                                                         |
     1         ' ---------------------------------------------------',/ |
     2         ' Reading original grid file, zone #',i5)                |
802   format(  '  This grid does not contain IBLANK data.',/            |
     1         ' ---------------------------------------------------',/)|
900   format(                                                           |
     1'   nz = ',i4,' jcube = ',i4,' kcube = ',i4,/                     |
     2'       jp3d_min = ',i5,'  jp3d_max = ',i5,/                      |
     3'       kp3d_min = ',i5,'  kp3d_max = ',i5 )                      |
901   format(                                                           |
     1'Local block dimensions: ',i5,' x ',i5/                           |
     2)                                                                 |
c1900  format(4e19.11)                                                   |
1900  format(3f20.11)                                                   |
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c end 'write_subgrids'   


c======================================================================|
      subroutine get_global_index
     1 (nz,j,k,jdo,kdo,j1,j2,k1,k2,ioverlap,
     2 nodeid,numprocs,j_subgrid,k_subgrid,
     3 idimj,idimk,jdim,kdim,
     4 jindex_global,kindex_global)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c 	This subroutine is called to give the index range for subgrids
c  	using global indices. It is called from various locations
c  	throughout the code.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c	D. W. Barnette
c   	Sandia National Laboratories
c   	Org. 9221, MS 1111
c     	Albuquerque, NM 87112
c  	email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c	break
c	link_overset
c	out1planetoplt3d_break
c	out1planetoproc_break
c	outallplanestoplt3d_break
c	outallplanestoproc_break
c	write_subgrids
c
c Called routines, in order of appearance:
c	NONE
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c	02-13-96	Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c
c   	Arguments	Description
c 	---------	-----------------------------------------
c Input:
c	nz		zone number
c	j		index for global index array
c	k		index for global index array
c	jdo		if 1, determine j1 and j2 indices
c	kdo		if 1, determine k1 and k2 indices
c       ioverlap        extent of grid cell overlap between subgrids
c Output (global values):
c	j1		min j index value for subgrid j,k,l in zone nz
c	j2		max j index value for subgrid j,k,l in zone nz
c	k1		min k index value for subgrid j,k,l in zone nz
c	k2		max k index value for subgrid j,k,l in zone nz
c
c
      dimension j_subgrid(numprocs),k_subgrid(numprocs)

      dimension jindex_global(numprocs),kindex_global(numprocs)

      integer idimj,idimk,jdim,kdim,ioverlap

      logical DEBUG
c
      DEBUG = .false.
c      DEBUG = .true.
c
      jmax=jdim
      kmax=kdim
      jkmax=jmax*kmax

c check value for ioverlap
      if((ioverlap.ne.0).and.(ioverlap.ne.1).and.(ioverlap.ne.2)) then
       write(3,*) '                                                    '
       write(3,*) ' ERROR: Subroutine GET_GLOBAL_INDEX                 '
       write(3,*) '   ioverlap in sub. get_global_index is not valid   '
       write(3,10) ioverlap
       write(3,*) '                                                    '
       write(3,*) '  PROGRAM STOPPED.                                  '
       write(3,*) '                                                    '
       write(6,*) '                                                    '
       write(6,*) ' ERROR: Subroutine GET_GLOBAL_INDEX                 '
       write(6,*) '   ioverlap in sub. get_global_index is not valid   '
       write(6,10) ioverlap
       write(6,*) '                                                    '
       write(6,*) '  PROGRAM STOPPED.                                  '
       write(6,*) '                                                    '
       stop '1: Subroutine GET_GLOBAL_INDEX '
      endif
c
c set indices
      if (DEBUG) then 
       print*
       print*,' Before set indices'
       print*, ' idimj,idimk,j,k = ',idimj,idimk,j,k
       print*
      endif
      index    =(j-1)*idimk+(k-1)+1
      index_jm1=(j-2)*idimk+(k-1)+1
      index_km1=(j-1)*idimk+(k-2)+1

c
      if (DEBUG) then
        print*
        print*, ' In sub preproc.f/get_global_index:'
        print*, '   numprocs = ',numprocs
        print*, '   index = ',index
        print*, '   index_jm1 = ',index_jm1
        print*, '   index_km1 = ',index_km1
        print*, '   ioverlap = ',ioverlap
        print*, '   j,k = ',j,k
        print*, '   idimj,idimk = ',idimj,idimk
        print*, '   jdim,kdim = ',jdim,kdim
        print*, '   jdo,kdo = ',jdo,kdo
        print*
        print*, 'Max global indices for this subgrid:'
        print*, index,jindex_global(index),kindex_global(index)
      endif
c
      if(jdo.ne.0) then
c j values
       if (idimj.eq.1) then
        j2=jindex_global(index)
        j1=1
       else if(j.eq.1) then
        j2=jindex_global(index)+ioverlap
        j1=1
       else if(j.eq.idimj) then
        j2=jindex_global(index)
c        j1=j2-j_subgrid(index)+1
        j1=jindex_global(index_jm1)
       else
        j2=jindex_global(index)+ioverlap
c        j1=j2-j_subgrid(index)+1
        j1=jindex_global(index_jm1)
       endif
      endif
c
      if(kdo.ne.0) then
c k values
       if (idimk.eq.1) then
        k2=kindex_global(index)
        k1=1
       else if(k.eq.1) then
        k2=kindex_global(index)+ioverlap
        k1=1
       else if(k.eq.idimk) then
        k2=kindex_global(index)
c        k1=k2-k_subgrid(index)+1 
        k1=kindex_global(index_km1)
       else
        k2=kindex_global(index)+ioverlap
c        k1=k2-k_subgrid(index)+1
        k1=kindex_global(index_km1)
       endif
      endif
c
c
      if (DEBUG) then
        print*, ' j1,k1 = ',j1,k1
        print*, ' j2,k2 = ',j2,k2
        print*
      endif

c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
10    format('      ioverlap = ',i10)                                   |
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c end 'get_global_index'
