!cpb   macros for gcyl <-> not gcyl distinction
#ifdef gcyl
#define GCYL_I i,
#define GCYL_J ,j
#else ! not gcyl
#define GCYL_I
#define GCYL_J
#endif ! not gcyl

#ifdef SS_DOUBLE_PREC
#define SS_PREC_CONV
#define SS_PREC dp
#define SS_PREC_BYTE 8
#define SS_MPI_PREC MPI_REAL8
#else ! not SS_DOUBLE_PREC
#define SS_PREC_CONV real
#define SS_PREC sp
#define SS_PREC_BYTE 4
#define SS_MPI_PREC MPI_REAL4
#endif ! not SS_DOUBLE_PREC


!cpb   =================================================================
!cpb   =================================================================
!cpb                     SUBSPACE DATA COLLECTION
!cpb   =================================================================
!cpb   =================================================================
!cpb
!cpb   This file contains subroutines to support the collection and
!cpb   maintenance of data streams during runtime.

!cpb   =================================================================
!cpb   =================================================================
!cpb                Initialization of subspace data collection
!cpb   =================================================================
!cpb   =================================================================
      subroutine subspace_init()

      use SSVAR
      use mpi
      use subroutineso
      use subroutines3d
      use problemcase
      implicit none

!cpb   =================================================================
!cpb   input arguments of this subroutine
!cpb   =================================================================

!cpb   =================================================================
!cpb   output arguments of this subroutine
!cpb   =================================================================

!cpb   =================================================================
!cpb   internal variables
!cpb   =================================================================
      integer l0
      integer ib,ss_index,istart,iend,jstart,jend,key
      integer colour,npoints,subspace_flow_io_skip
      integer blockstart_ijk(4),blockdims_ijk(4),totaldims_ijk(4)
      integer ss_extra
      character*100 dummy,iofile
      logical lex
!cpb   =================================================================
!cpb   main part
!cpb   =================================================================

!cpb   =================================================================
!cpb             read input file for the histogram collection
!cpb   =================================================================
      if(myid.eq.0) write(*,*) 'Subspaces will be written.'
      !check if inputfile exists
      inquire(file='SUBSPACES_INPUT.dat',exist=lex)
      if(lex)then
        open(70,file='SUBSPACES_INPUT.dat',status='old')
        !read in number of subspaces
        read(70,*) ss_num
      else
        ss_num = 0
      endif

      ss_extra = 0

      if (ioflag.gt.0) ss_extra=ss_extra+1

      ss_num = ss_num+ss_extra
      allocate(ss_capt(ss_num))
      allocate(ss_opt(ss_num))
      allocate(ss_tag(ss_num))
      allocate(ss_tag_grid(ss_num))
      allocate(ss_blocks(mbk+1))
      allocate(ss_num_var(ss_num))
      allocate(ss_spec_block(ss_num,mbk+1,9))
      allocate(ss_dim_block(ss_num,4))
      allocate(ss_dim_proc(ss_num,12))
      allocate(ss_dim_spec(ss_num,3))
      allocate(ss_in_subspace(ss_num))
      allocate(ss_communicator(ss_num))
      allocate(ss_rank(ss_num))
      allocate(ss_fh(ss_num))
      allocate(ss_blockid(ss_num))
      allocate(ss_wall_opt(ss_num))
#ifndef SS_BLOCKING_WRITE
      allocate(ss_status(ss_num,MPI_STATUS_SIZE))
      allocate(ss_first_write(ss_num))
#endif
#ifdef SS_MULTIBLOCK
      allocate(ss_blockid_offset(ss_num))
      allocate(ss_blockid_specs((mbk+1)*3+1,ss_num))
#endif ! SS_MULTIBLOCK
      !read subvolume specifications

      do i=1,ss_num-ss_extra
        ! read specs for subspace
        read(70,*) ss_tag(i),ss_opt(i),ss_capt(i),(ss_blocks(j),j=1,mbk+1)
        do j=1,mbk+1
          if(ss_blocks(j)) then
            !k=1,2,3 start,end,skip in xi-direction
            !k=4,5,6 start,end,skip in eta-direction
            !k=7,8,9 start,end,skip in theta-direction
            read(70,*) (ss_spec_block(i,j,k),k=1,9)
#ifndef THREE_D
            ss_spec_block(i,j,7) = 1
            ss_spec_block(i,j,8) = 1
            ss_spec_block(i,j,9) = 1
#endif ! not THREE_D
          else
            ss_spec_block(i,j,1) = 1
            ss_spec_block(i,j,2) = 0
            ss_spec_block(i,j,3) = 1
            ss_spec_block(i,j,4) = 1
            ss_spec_block(i,j,5) = 0
            ss_spec_block(i,j,6) = 1
            ss_spec_block(i,j,7) = 1
            ss_spec_block(i,j,8) = 0
            ss_spec_block(i,j,9) = 1
          endif
        enddo !j
      enddo !i
      close(70)
!cpb   =================================================================
!cpb           END  read input file for the histogram collection
!cpb   =================================================================
      l0 = 1

      if (ioflag.gt.0) then
        if (mod(ioflag,10).eq.1) then
          ss_tag(ss_num-ss_extra+l0)='FLOW_phys' !name of the file
        elseif (mod(ioflag,10).eq.2) then
          ss_tag(ss_num-ss_extra+l0)='FLOW_spec' !name of the file
        else
          ss_tag(ss_num-ss_extra+l0)='FLOW' !name of the file
        endif
        if (mod(ioflag,10).eq.1) then
          ss_opt(ss_num-ss_extra+l0)= 0 ! option giving what is written out
        elseif (mod(ioflag,10).eq.2) then
          ss_opt(ss_num-ss_extra+l0)= 100 ! option giving what is written out
        else
          ss_opt(ss_num-ss_extra+l0)= 0
        endif
        ss_capt(ss_num-ss_extra+l0)= ncapt ! capturing time step
        if(ioflag/10.eq.0) then
          subspace_flow_io_skip=1
        elseif(ioflag/10.eq.1) then
          subspace_flow_io_skip=2
        elseif(ioflag/10.eq.2) then
          subspace_flow_io_skip=4
        else
          call stopnow('subspace_flow I/O wrong input')
        endif
        do j=0,mbk
          ib=j+1
          ss_spec_block(ss_num-ss_extra+l0,ib,1)=1
          ss_spec_block(ss_num-ss_extra+l0,ib,2)=lximb(j)+1
          ss_spec_block(ss_num-ss_extra+l0,ib,3)=subspace_flow_io_skip
          ss_spec_block(ss_num-ss_extra+l0,ib,4)=1
          ss_spec_block(ss_num-ss_extra+l0,ib,5)=letmb(j)+1
          ss_spec_block(ss_num-ss_extra+l0,ib,6)=subspace_flow_io_skip
          if (mod(ioflag,10).eq.1) then
            ss_spec_block(ss_num-ss_extra+l0,ib,7)=1
            ss_spec_block(ss_num-ss_extra+l0,ib,8)=lzemb(j)+1
            ss_spec_block(ss_num-ss_extra+l0,ib,9)=subspace_flow_io_skip
  !        elseif (mod(ioflag,10).eq.2) then
  !          ss_spec_block(ss_num-ss_extra+l0,i,7)=0
  !          ss_spec_block(ss_num-ss_extra+l0,i,8)=kmod/subspace_flow_io_skip
  !          ss_spec_block(ss_num-ss_extra+l0,i,9)=1
          else
            ss_spec_block(ss_num-ss_extra+l0,ib,7)=1
            ss_spec_block(ss_num-ss_extra+l0,ib,8)=lzemb(j)+1
            ss_spec_block(ss_num-ss_extra+l0,ib,9)=subspace_flow_io_skip
          endif
         end do
      endif
!cpb   =================================================================
!cpb                     localization of subspaces
!cpb   =================================================================
#ifdef SS_MULTIBLOCK
      ss_blockid_offset=0
      ss_blockid_specs(1,:)=0
#endif ! SS_MULTIBLOCK
      do ss_index=1,ss_num
      ss_in_subspace(ss_index)=.false.
        if(   (ibegin(sx)  .le. ss_spec_block(ss_index,mb+1,2)) &
       .and. (ibegin(sx+1).gt. ss_spec_block(ss_index,mb+1,1)) &
       .and. (jbegin(sy)  .le. ss_spec_block(ss_index,mb+1,5)) &
       .and. (jbegin(sy+1).gt. ss_spec_block(ss_index,mb+1,4)) &
       )then
          !flag if process is involved
          ss_in_subspace(ss_index)=.true.
!cpb   -----------------------------------------------------------------
!cpb        correction that spectral output is restricted to planes
!cpb   -----------------------------------------------------------------
      if (ss_opt(ss_index).ge.100 &
         .or. ss_opt(ss_index).eq.8 &
         .or. ss_opt(ss_index).eq.10 &
         .or. ss_opt(ss_index).eq.11 &
               ) then
        do ib=1,mbk+1
          if(mb+1.eq.ib)then
#ifdef THREE_D
            ss_dim_spec(ss_index,1) = ss_spec_block(ss_index,ib,7)
            ss_dim_spec(ss_index,2) = ss_spec_block(ss_index,ib,8)
            ss_dim_spec(ss_index,3) = ss_spec_block(ss_index,ib,9)
#else ! not THREE_D
            ss_dim_spec(ss_index,1) = 0
            ss_dim_spec(ss_index,2) = 0
            ss_dim_spec(ss_index,3) = 1
#endif ! not THREE_D
          endif
          ! make sure that excluded blocks are not included again
          if(ss_spec_block(ss_index,ib,8) .lt.ss_spec_block(ss_index,ib,7)) then
            ss_spec_block(ss_index,ib,7) = 1
            ss_spec_block(ss_index,ib,8) = 0
            ss_spec_block(ss_index,ib,9) = 1
          else
            ss_spec_block(ss_index,ib,7) = 1
            ss_spec_block(ss_index,ib,8) = 1
            ss_spec_block(ss_index,ib,9) = 1
          endif

        end do
      endif
!cpb   -----------------------------------------------------------------
!cpb       end correction that spectral output is restricted to planes
!cpb   -----------------------------------------------------------------


#ifdef SS_MULTIBLOCK
        do ib=1,mbk+1
          npoints= ((ss_spec_block(ss_index,ib,2)-ss_spec_block(ss_index,ib,1)) &
                    /ss_spec_block(ss_index,ib,3)+1) &
                 *((ss_spec_block(ss_index,ib,5)-ss_spec_block(ss_index,ib,4)) &
                   /ss_spec_block(ss_index,ib,6)+1) &
                 *((ss_spec_block(ss_index,ib,8)-ss_spec_block(ss_index,ib,7)) &
                   /ss_spec_block(ss_index,ib,9)+1)
          !set offset for multiblock writing
          if(npoints.gt.0) then
            if(ib.le.mb)then
              ss_blockid_offset(ss_index)=ss_blockid_offset(ss_index)+npoints
            endif
            ss_blockid_specs(ss_blockid_specs(1,ss_index)*3+2,ss_index)= &
              ((ss_spec_block(ss_index,ib,2)-ss_spec_block(ss_index,ib,1)) &
                    /ss_spec_block(ss_index,ib,3)+1)
            ss_blockid_specs(ss_blockid_specs(1,ss_index)*3+3,ss_index)= &
              ((ss_spec_block(ss_index,ib,5)-ss_spec_block(ss_index,ib,4)) &
                   /ss_spec_block(ss_index,ib,6)+1)
            ss_blockid_specs(ss_blockid_specs(1,ss_index)*3+4,ss_index)= &
              ((ss_spec_block(ss_index,ib,8)-ss_spec_block(ss_index,ib,7)) &
                   /ss_spec_block(ss_index,ib,9)+1)
            ss_blockid_specs(1,ss_index)=ss_blockid_specs(1,ss_index)+1
          endif
        enddo !ib
#else ! not SS_MULTIBLOCK
        ss_blockid(ss_index) = 0
        do ib=1,mb+1
          npoints= ((ss_spec_block(ss_index,ib,2)-ss_spec_block(ss_index,ib,1)) &
                    /ss_spec_block(ss_index,ib,3)+1) &
                 *((ss_spec_block(ss_index,ib,5)-ss_spec_block(ss_index,ib,4)) &
                   /ss_spec_block(ss_index,ib,6)+1) &
                 *((ss_spec_block(ss_index,ib,8)-ss_spec_block(ss_index,ib,7)) &
                   /ss_spec_block(ss_index,ib,9)+1)
          !set offset for multiblock writing
          if(npoints.gt.0) then
            ss_blockid(ss_index)=ss_blockid(ss_index)+1
          endif
        enddo !ib
#endif ! not SS_MULTIBLOCK
!cpb   -----------------------------------------------------------------
!cpb                set block specifications for sub-volumes
!cpb   -----------------------------------------------------------------
          !total number of points of block in xi-direction
          ss_dim_block(ss_index,1) = (ss_spec_block(ss_index,mb+1,2)-ss_spec_block(ss_index,mb+1,1)) &
          /ss_spec_block(ss_index,mb+1,3)+1
          !total number of points of block in eta-direction
          ss_dim_block(ss_index,2) = (ss_spec_block(ss_index,mb+1,5)-ss_spec_block(ss_index,mb+1,4)) &
          /ss_spec_block(ss_index,mb+1,6)+1
          !total number of points of block in theta-direction
          ss_dim_block(ss_index,3) = (ss_spec_block(ss_index,mb+1,8)-ss_spec_block(ss_index,mb+1,7)) &
          /ss_spec_block(ss_index,mb+1,9)+1
!cpb   -----------------------------------------------------------------
!cpb            set process specifications for sub-volumes
!cpb   -----------------------------------------------------------------
          !total number of points of process in xi-direction

          !start or volume in actual process
          istart = 1+max(0,ss_spec_block(ss_index,mb+1,1)-ibegin(sx))
          iend =lxi+1-max(0,(ibegin(sx+1)-1)-ss_spec_block(ss_index,mb+1,2))
          !if sx is not first process in i direction then shift due to
          !skipping
          if(istart.eq.1) then
            istart=1+ &
           mod( ss_spec_block(ss_index,mb+1,3)- &
               mod( (ibegin(sx)-ss_spec_block(ss_index,mb+1,1)) &
                  ,ss_spec_block(ss_index,mb+1,3)) &
             ,ss_spec_block(ss_index,mb+1,3))
          endif !start=1

          jstart = 1+max(0,ss_spec_block(ss_index,mb+1,4)-jbegin(sy))
          jend =let+1-max(0,(jbegin(sy+1)-1) &
            -ss_spec_block(ss_index,mb+1,5))
          if(jstart.eq.1) then
            jstart=1+ &
           mod( ss_spec_block(ss_index,mb+1,6)- &
                  mod( (jbegin(sy)-ss_spec_block(ss_index,mb+1,4) ) &
                     ,ss_spec_block(ss_index,mb+1,6)) &
             ,ss_spec_block(ss_index,mb+1,6))
          endif !start=1

          !number i-entries local proc
          ss_dim_proc(ss_index,1) = (iend-istart)/ss_spec_block(ss_index,mb+1,3) + 1
          !number j-entries local proc
          ss_dim_proc(ss_index,2) = (jend-jstart)/ss_spec_block(ss_index,mb+1,6) + 1

          !number k-entries local proc
          ss_dim_proc(ss_index,3) = &
         (ss_spec_block(ss_index,mb+1,8) - ss_spec_block(ss_index,mb+1,7)) &
            /ss_spec_block(ss_index,mb+1,9)+1

          !number i-entries before
          if(ibegin(sx).le.ss_spec_block(ss_index,mb+1,1)) then
            ss_dim_proc(ss_index,4) = 0
          else
            ss_dim_proc(ss_index,4) = &
             (ibegin(sx)-1  -ss_spec_block(ss_index,mb+1,1)) &
             /ss_spec_block(ss_index,mb+1,3) &
             +1
            ! check if process does not contain a point due to skipping
            if( ss_dim_block(ss_index,1).le.ss_dim_proc(ss_index,4) ) &
              !if not containing a point => take it out of subspace
             ss_in_subspace(ss_index)=.false.
          endif

           !number j-entries before
          if(jbegin(sy).le.ss_spec_block(ss_index,mb+1,4)) then
            ss_dim_proc(ss_index,5) = 0
          else
            ss_dim_proc(ss_index,5) = &
           (jbegin(sy)-1 -ss_spec_block(ss_index,mb+1,4)) &
           /ss_spec_block(ss_index,mb+1,6) &
           +1
            ! check if process does not contain a point due to skipping
            if( ss_dim_block(ss_index,2).le.ss_dim_proc(ss_index,5) ) &
              !if not containing a point => take it out of subspace
             ss_in_subspace(ss_index)=.false.
          endif

          !number k-entries before
          ss_dim_proc(ss_index,6) = 0

          !start index i in local process
          ss_dim_proc(ss_index,7) = istart

          !start index j in local process
          ss_dim_proc(ss_index,8) = jstart

          !start index k in local process
          ss_dim_proc(ss_index,9) = ss_spec_block(ss_index,mb+1,7)

          !start index i in local process
          ss_dim_proc(ss_index,10) = iend

          !start index j in local process
          ss_dim_proc(ss_index,11) = jend

          !start index k in local process
          ss_dim_proc(ss_index,12) =  ss_spec_block(ss_index,mb+1,8)



!cpb   set respective number of variables
        if(ss_opt(ss_index) .eq. 0) then
          ss_num_var(ss_index) = 5 ! !cpb passive scalars are now option 13  !cds, generalize.
        elseif(ss_opt(ss_index) .le. 6) then
          ss_num_var(ss_index) = 1
        elseif(ss_opt(ss_index) .eq. 7) then
          ss_num_var(ss_index) = 6
        elseif(ss_opt(ss_index) .eq. 8) then
          ss_num_var(ss_index) = 2*(ss_dim_spec(ss_index,2)-ss_dim_spec(ss_index,1)+1) &
         /ss_dim_spec(ss_index,3)
        elseif(ss_opt(ss_index) .eq. 9) then
          ss_num_var(ss_index) = 3
        elseif(ss_opt(ss_index) .eq. 10) then
          ss_num_var(ss_index) = 3*(ss_dim_spec(ss_index,2)-ss_dim_spec(ss_index,1)+1) &
         /ss_dim_spec(ss_index,3)
        elseif(ss_opt(ss_index) .eq. 11) then
          ss_num_var(ss_index) = 1*(ss_dim_spec(ss_index,2)-ss_dim_spec(ss_index,1)+1) &
         /ss_dim_spec(ss_index,3)
        elseif(ss_opt(ss_index) .eq. 12) then
          ss_num_var(ss_index) = 1
          if ((ss_spec_block(ss_index,mb+1,4).eq. &
              ss_spec_block(ss_index,mb+1,5)).and. &
             (ss_spec_block(ss_index,mb+1,4).eq. &
                             letmb(mb+1)+1))  &
                    ss_wall_opt(ss_index) = 1 
          if ((ss_spec_block(ss_index,mb+1,4).eq. &
              ss_spec_block(ss_index,mb+1,5)).and. &
             (ss_spec_block(ss_index,mb+1,4).eq.1))  &
                    ss_wall_opt(ss_index) = 2 
          if ((ss_spec_block(ss_index,mb+1,1).eq. &
              ss_spec_block(ss_index,mb+1,2)).and. &
             (ss_spec_block(ss_index,mb+1,1).eq. &
                             lximb(mb+1)+1))  &
                    ss_wall_opt(ss_index) = 3
          if ((ss_spec_block(ss_index,mb+1,1).eq. &
              ss_spec_block(ss_index,mb+1,2)).and. &
             (ss_spec_block(ss_index,mb+1,1).eq.1))  &
                    ss_wall_opt(ss_index) = 4
        elseif(ss_opt(ss_index) .eq. 13) then
          ss_num_var(ss_index) = 5-5
        elseif(ss_opt(ss_index) .eq. 100) then
          ss_num_var(ss_index) = 5* &
         (ss_dim_spec(ss_index,2) &
          -ss_dim_spec(ss_index,1) &
          +1) &
         /ss_dim_spec(ss_index,3)
        elseif(ss_opt(ss_index) .eq. 101) then

          ss_num_var(ss_index) = 5* &
         (ss_dim_spec(ss_index,2) &
          -ss_dim_spec(ss_index,1) &
          +1) &
         /ss_dim_spec(ss_index,3)
        elseif(ss_opt(ss_index) .eq. 102) then
          ss_num_var(ss_index) = 5* &
         (ss_dim_spec(ss_index,2) &
          -ss_dim_spec(ss_index,1) &
          +1) &
         /ss_dim_spec(ss_index,3)
        else
          call stopnow('Subspace option not defined')
        endif


        else
          ss_in_subspace(ss_index)=.false.
          ss_dim_block(ss_index,1) = 0
          ss_dim_block(ss_index,2) = 0
          ss_dim_block(ss_index,3) = 0
          ss_dim_proc(ss_index,1) = 0
          ss_dim_proc(ss_index,2) = 0
          ss_dim_proc(ss_index,3) = 0
          ss_dim_proc(ss_index,4) = 0
          ss_dim_proc(ss_index,5) = 0
          ss_dim_proc(ss_index,6) = 0
          ss_dim_proc(ss_index,7) = 0
          ss_dim_proc(ss_index,8) = 0
          ss_dim_proc(ss_index,9) = 0
          ss_dim_proc(ss_index,10) = 0
          ss_dim_proc(ss_index,11) = 0
          ss_dim_proc(ss_index,12) = 0
        endif !in subspace
!cpb   =================================================================
!cpb                  set the limits for each block
!cpb   =================================================================

      totaldims_ijk(1)  = ss_dim_block(ss_index,1)
      totaldims_ijk(2)  = ss_dim_block(ss_index,2)
      totaldims_ijk(3)  = ss_dim_block(ss_index,3)
      totaldims_ijk(4)  = ss_num_var(ss_index)
      blockdims_ijk(1)  = ss_dim_proc(ss_index,1)
      blockdims_ijk(2)  = ss_dim_proc(ss_index,2)
      blockdims_ijk(3)  = ss_dim_proc(ss_index,3)
      blockdims_ijk(4)  = ss_num_var(ss_index)
      blockstart_ijk(1) = ss_dim_proc(ss_index,4)
      blockstart_ijk(2) = ss_dim_proc(ss_index,5)
      blockstart_ijk(3) = ss_dim_proc(ss_index,6)
      blockstart_ijk(4) = 0

!cpb   =================================================================
!cpb                  set communicator for each subvolume
!cpb   =================================================================
      key=myid
      if(.not.(ss_in_subspace(ss_index))) then ! those are cores outside
        colour=MPI_UNDEFINED
      else ! the following for cores that are inside volume to be written
        colour=ss_index
         call MPI_TYPE_CREATE_SUBARRAY(4, totaldims_ijk, &
          blockdims_ijk, blockstart_ijk, MPI_ORDER_FORTRAN, &
          SS_MPI_PREC, ss_dim_block(ss_index,4), ierr)
         call MPI_TYPE_COMMIT(ss_dim_block(ss_index,4), ierr)
      endif
#ifndef SS_MULTIBLOCK
      call MPI_COMM_SPLIT(block_comm,colour,key,ss_communicator(ss_index),ierr)
#else ! SS_MULTIBLOCK
      call MPI_COMM_SPLIT(icom,colour,key,ss_communicator(ss_index),ierr)
#endif ! SS_MULTIBLOCK
        !get rank in subvolume communicator
      if( ss_in_subspace(ss_index) ) then ! those are cores outside
        call MPI_COMM_RANK(ss_communicator(ss_index),ss_rank(ss_index),ierr)
!cpb   write parts of the final tag now to reduce wasted runtime
!cpb   and to make more general names possible
        if(ss_index .le. ss_num-ss_extra)  then

#ifndef SS_MULTIBLOCK
!cpb   write GRID tag
          write(ss_tag_grid(ss_index),'(a,a,''_GRID_'',i0,''.xyz'')') &
            TRIM(ss_path),trim(ss_tag(ss_index)),ss_blockid(ss_index)
!cpb   write FIELD tag
          write(iofile,'(a,a,''_'', i0 ,''_var_'', i0)') &
            TRIM(ss_path),trim(ss_tag(ss_index)),ss_blockid(ss_index),ss_num_var(ss_index)
          ss_tag(ss_index)=iofile
#else ! SS_MULTIBLOCK
!cpb   write GRID tag
          write(ss_tag_grid(ss_index),'(a,a,''_GRID.xyz'')') &
            TRIM(ss_path),trim(ss_tag(ss_index))
!cpb   write FIELD tag
          write(iofile,'(a,a,''_var_'', i0)') &
            TRIM(ss_path),trim(ss_tag(ss_index)),ss_num_var(ss_index)
          ss_tag(ss_index)=iofile
#endif ! SS_MULTIBLOCK
        else ! FLOW field
#ifndef SS_MULTIBLOCK
!cpb   write GRID tag
          write(ss_tag_grid(ss_index), '(''FLOW/'',a,''_GRID_'',i0,''.xyz'')') &
            trim(ss_tag(ss_index)),ss_blockid(ss_index)
!cpb   write FIELD tag
          write(iofile,'(''FLOW/'',a,''_'', i0)') &
            trim(ss_tag(ss_index)),ss_blockid(ss_index)
          ss_tag(ss_index)=iofile
#else ! SS_MULTIBLOCK
!cpb   write GRID tag
          write(ss_tag_grid(ss_index),'(''FLOW/'',a,''_GRID.xyz'')') &
            trim(ss_tag(ss_index))
!cpb   write FIELD tag
          write(iofile,'(''FLOW/'',a)') &
            trim(ss_tag(ss_index))
          ss_tag(ss_index)=iofile
#endif ! SS_MULTIBLOCK
        endif
      endif
      enddo!ss_index
!cpb   =================================================================
!cpb                   END localization of subspaces
!cpb   =================================================================

!cpb   =================================================================
!cpb   =================================================================
!cpb           END initialization of subspace data collection
!cpb   =================================================================
!cpb   =================================================================
      end !subroutine ss_init


!cpb   =================================================================
!cpb   =================================================================
!cpb                IO subroutine to write out subspaces
!cpb   =================================================================
!cpb   =================================================================
!cpb
      subroutine subspace_write_data(nit,sample_option)
!cpb   =================================================================
!cpb   This routine writes out
!cpb   subspaces of the computational domain. The specification of these
!cpb   subdomains is done via the input file 'SUBSPACES_INPUT.dat' which
!cpb   has to be in the same directory as the executable.
!cpb
!cpb   The specification is grid dependent!
!cpb
!cpb   Bechlars, November 2012
!cpb   =================================================================

!cpb   =================================================================
!cpb                             load modules
!cpb   =================================================================
      use SSVAR
      use mpi
      use subroutineso
      use subroutines3d
      use problemcase
      implicit none
!cpb   =================================================================
!cpb                         subroutine arguments
!cpb   =================================================================
      integer nit
      integer sample_option
!cpb   =================================================================
!cpb                  additional variables/parameters
!cpb   =================================================================
      integer ivar,iii,jjj,kkk

      !MPI/PLOT3D environment
      integer accmod,ialloc(7)
      integer (kind=MPI_OFFSET_KIND) :: disp,ioff,iread,l2
      integer totaldims_ijk(4),blockdims_ijk(4),blockstart_ijk(4)
!      integer subarray_ijk
      integer status(MPI_STATUS_SIZE), contig4, contig_plot3d
      integer real_mp_size
      real(kind=wp) xxl,xnit,tw
      real(kind=SS_PREC), save ::  riobuf(4)
      integer*4, save :: iiobuf(4)
      real(kind=SS_PREC) ::  riobuf_4(4),temp,normx,normy
#ifdef SS_BLOCKING_WRITE
      real(kind=SS_PREC), allocatable :: ss_qio(:)
#endif
#ifndef SS_MULTIBLOCK
      integer*4 :: iiobuf_plot3d(3)
#else ! SS_MULTIBLOCK
      integer*4 , allocatable :: iiobuf_plot3d(:)
#endif ! SS_MULTIBLOCK
      character*8 chl
      character*70 iofile
      logical lex
      equivalence (riobuf,iiobuf)
      integer ss_start,ss_index
      integer ss_idim,ss_jdim,ss_kdim
      integer ss_istart,ss_iend,ss_iskip
      integer ss_jstart,ss_jend,ss_jskip
      integer ss_kstart,ss_kend,ss_kskip
      real(kind=wp) xdudz,xdvdr,xdudr,xdvdz

!cpb   =================================================================
!cpb                          executing part
!cpb   =================================================================

      ss_start = 1
      do ss_index=1,ss_num
!cpb   ONLY DO THIS FOR CORES THAT ARE WITHIN VOLUME
      if((ss_in_subspace(ss_index)) &
        .and.(mod(n,ss_capt(ss_index)).eq.0) &
        .and.((    n .ne. 0) &
               .or.(n .eq. 0)  &
              .and.(ss_index.eq.ss_num) &
              .and.(ioflag.gt.0) &
             ) &
        .and. &
                ! ss options < 100 are called for sample option 0
              ((    (ss_opt(ss_index).lt.100).and.(0.eq.sample_option)) &
                ! ss option > 100 are called for sample option 1
                .or.((ss_opt(ss_index).ge.100).and.(1.eq.sample_option)) &
              ) &
        ) then
        accmod = ior(MPI_MODE_WRONLY,MPI_MODE_CREATE)
        !name and directory of output files

        write(iofile,'(a,''_'',i0,''.raw'')')trim(ss_tag(ss_index)),n

!cpb   =================================================================
!cpb                  set the limits for each block and proc
!cpb   =================================================================
        ss_idim=ss_dim_proc(ss_index,1)
        ss_jdim=ss_dim_proc(ss_index,2)
        ss_kdim=ss_dim_proc(ss_index,3)

        totaldims_ijk(1)  = ss_dim_block(ss_index,1)
        totaldims_ijk(2)  = ss_dim_block(ss_index,2)
        totaldims_ijk(3)  = ss_dim_block(ss_index,3)
        totaldims_ijk(4)  = ss_num_var(ss_index)
        blockdims_ijk(1)  = ss_idim
        blockdims_ijk(2)  = ss_jdim
        blockdims_ijk(3)  = ss_kdim
        blockdims_ijk(4)  = ss_num_var(ss_index)
        blockstart_ijk(1) = ss_dim_proc(ss_index,4)
        blockstart_ijk(2) = ss_dim_proc(ss_index,5)
        blockstart_ijk(3) = ss_dim_proc(ss_index,6)
        blockstart_ijk(4) = 0

        !limits in i-direction
        ss_istart =ss_dim_proc(ss_index,7)
        ss_iend =ss_dim_proc(ss_index,10)
        ss_iskip =ss_spec_block(ss_index,mb+1,3)

        !limits in j-direction
        ss_jstart =ss_dim_proc(ss_index,8)
        ss_jend =ss_dim_proc(ss_index,11)
        ss_jskip =ss_spec_block(ss_index,mb+1,6)

        !limits in k-direction
        ss_kstart =ss_dim_proc(ss_index,9)
        ss_kend =ss_dim_proc(ss_index,12)
        ss_kskip =ss_spec_block(ss_index,mb+1,9)

!cpb   =================================================================
!cpb                  set subarray for each subvolume
!cpb   =================================================================
!        subarray_ijk=ss_dim_block(ss_index,4)
!cpb   ==================================================================
!cpb                      OUTPUT IN PLOT 3D FORMAT
!cpb   ==================================================================
#ifdef SS_BLOCKING_WRITE
        !allocate writing vector
        allocate (ss_qio( &
                     int(ss_idim,kind=MPI_OFFSET_KIND) &
                    *int(ss_jdim,kind=MPI_OFFSET_KIND) &
                    *int(ss_kdim,kind=MPI_OFFSET_KIND) &
                    *int(ss_num_var(ss_index),kind=MPI_OFFSET_KIND) &
                    ), stat=ialloc(1))
        !check if allocation went fine
        if(ialloc(1) .ne. 0)then
          write(*,*)'allocate error (ss_volume):',ialloc
          call MPI_FINALIZE(ierr)
          stop
        endif
#endif

#ifndef SS_BLOCKING_WRITE
        if (ss_first_write(ss_index)) then
#ifndef SS_NON_COLLECTIVE_WRITE
          call MPI_File_write_all_end(ss_fh(ss_index), ss_qio,ss_status(ss_index,1),ierr)
          call MPI_FILE_CLOSE(ss_fh(ss_index),ierr)
#else
          call MPI_Wait(ss_status(ss_index,1),status,ierr)
          call MPI_FILE_CLOSE(ss_fh(ss_index),ierr)
#endif
!crp     set flag to true so for every subsequent write during this run it
!crp     is waited for the last writing to finish
          ss_first_write(ss_index)=.false.
        end if
#endif !non blocking write

        if(ss_opt(ss_index) .eq. 0) then
          !write data in vector
          !loop over values
          l2=0
          do ivar=1,5   !!cpb: this shall be 1-5 passive scalars are written with option 13 ! cds: general case
          !loop theta direction
          do k=ss_kstart,ss_kend,ss_kskip
          !loop eta direction
          do j=ss_jstart,ss_jend,ss_jskip
            !loop xi direction
            do i=ss_istart,ss_iend,ss_iskip
              l2 = l2 + 1
              ! Q=Q/J
              ss_qio((ss_start-1)+l2) = SS_PREC_CONV(qa(indx3(i-1,j-1,k-1,1),ivar))
            end do
          end do
          end do
          end do
        elseif (ss_opt(ss_index) .eq. 1) then
          !write data in vector
          !loop over values
          l2=0
          !loop theta direction
          do k=ss_kstart,ss_kend,ss_kskip
          !loop eta direction
          do j=ss_jstart,ss_jend,ss_jskip
            !loop xi direction
            do i=ss_istart,ss_iend,ss_iskip
              l2 = l2 + 1
              ! Q=Q/J
              ss_qio((ss_start-1)+l2) = SS_PREC_CONV(qb(indx3(i-1,j-1,k-1,1),2)/qb(indx3(i-1,j-1,k-1,1),1)+umf(1))
            end do
          end do
          end do
        elseif (ss_opt(ss_index) .eq. 2) then
          !write data in vector
          !loop over values
          l2=0
          !loop theta direction
          do k=ss_kstart,ss_kend,ss_kskip
          !loop eta direction
          do j=ss_jstart,ss_jend,ss_jskip
            !loop xi direction
            do i=ss_istart,ss_iend,ss_iskip
              l2 = l2 + 1
              ! Q=Q/J
              ss_qio((ss_start-1)+l2) = SS_PREC_CONV(qa(indx3(i-1,j-1,k-1,1),2))
            end do
          end do
          end do
        elseif (ss_opt(ss_index) .eq. 3) then
          !write data in vector
          !loop over values
          l2=0
          !loop theta direction
          do k=ss_kstart,ss_kend,ss_kskip
          !loop eta direction
          do j=ss_jstart,ss_jend,ss_jskip
            !loop xi direction
            do i=ss_istart,ss_iend,ss_iskip
              l2 = l2 + 1
              ! Q=Q/J
              ss_qio((ss_start-1)+l2) = SS_PREC_CONV(qa(indx3(i-1,j-1,k-1,1),3))
            end do
          end do
          end do
        elseif (ss_opt(ss_index) .eq. 4) then
          !write data in vector
          !loop over values
          l2=0
          !loop theta direction
          do k=ss_kstart,ss_kend,ss_kskip
          !loop eta direction
          do j=ss_jstart,ss_jend,ss_jskip
            !loop xi direction
            do i=ss_istart,ss_iend,ss_iskip
              l2 = l2 + 1
              ! Q=Q/J
              ss_qio((ss_start-1)+l2) = SS_PREC_CONV(qa(indx3(i-1,j-1,k-1,1),4))
            end do
          end do
          end do
        elseif (ss_opt(ss_index) .eq. 5) then
          !write data in vector
          !loop over values
          l2=0
          !loop theta direction
          do k=ss_kstart,ss_kend,ss_kskip
          !loop eta direction
          do j=ss_jstart,ss_jend,ss_jskip
            !loop xi direction
            do i=ss_istart,ss_iend,ss_iskip
              l2 = l2 + 1
              ! Q=Q/J
              ss_qio((ss_start-1)+l2) = SS_PREC_CONV(qa(indx3(i-1,j-1,k-1,1),5))     
            end do
          end do
          end do
        elseif (ss_opt(ss_index) .eq. 6) then
          !write data in vector
          !loop over values
          l2=0
          !loop theta direction
          do k=ss_kstart,ss_kend,ss_kskip
          !loop eta direction
          do j=ss_jstart,ss_jend,ss_jskip
            !loop xi direction
            do i=ss_istart,ss_iend,ss_iskip
              l2 = l2 + 1
              ! Q=Q/J
              ss_qio((ss_start-1)+l2) = SS_PREC_CONV(p(indx3(i-1,j-1,k-1,1)))
            end do
          end do
          end do
        else
          call stopnow('Subspace option not defined')
        endif


        ioff = 0
        disp = 0


        !if rank in subvolume communicator = 0 then overwrite existing file
        if(ss_rank(ss_index).eq.0)then
          inquire(file=iofile,exist=lex)
          if(lex)then
            call MPI_FILE_DELETE(iofile,MPI_INFO_NULL,ierr)
          endif
        endif

        !set iobuffer for plot3d

#ifndef SS_MULTIBLOCK
        iiobuf_plot3d(1)  =  totaldims_ijk(1)
        iiobuf_plot3d(2)  =  totaldims_ijk(2)
        iiobuf_plot3d(3)  =  totaldims_ijk(3)
#else ! SS_MULTIBLOCK
        allocate(iiobuf_plot3d(ss_blockid_specs(1,ss_index)*3+1))
        iiobuf_plot3d(1)  = ss_blockid_specs(1,ss_index)
        do i=1, ss_blockid_specs(1,ss_index)*3
            iiobuf_plot3d(1+i)  =  ss_blockid_specs(1+i,ss_index)
        enddo
#endif ! SS_MULTIBLOCK

#ifndef SS_MULTIBLOCK
        call MPI_TYPE_CONTIGUOUS(3, MPI_INTEGER4,contig_plot3d,ierr)
#else ! SS_MULTIBLOCK
        call MPI_TYPE_CONTIGUOUS(ss_blockid_specs(1,ss_index)*3+1,MPI_INTEGER4, contig_plot3d,ierr)
#endif ! SS_MULTIBLOCK
        call MPI_TYPE_COMMIT(contig_plot3d, ierr)
        call MPI_FILE_OPEN(ss_communicator(ss_index),iofile,accmod,MPI_INFO_NULL,ss_fh(ss_index),ierr)
        call MPI_FILE_SET_VIEW(ss_fh(ss_index), disp, MPI_INTEGER4,contig_plot3d, 'native', MPI_INFO_NULL, ierr)

        ! IF-block added by Jon Gibson, NAG
#ifndef SS_MULTIBLOCK
        if (ss_rank(ss_index) .eq. 0) then
          call MPI_FILE_WRITE(ss_fh(ss_index),iiobuf_plot3d,3, MPI_INTEGER4, status, ierr)
        endif
#else ! SS_MULTIBLOCK
        if (ss_rank(ss_index) .eq. 0) then
          call MPI_FILE_WRITE(ss_fh(ss_index), iiobuf_plot3d,ss_blockid_specs(1,ss_index)*3+1, MPI_INTEGER4, status, ierr)
        endif
#endif ! SS_MULTIBLOCK

#ifndef SS_MULTIBLOCK
        disp = disp + 4*3
#else ! SS_MULTIBLOCK
        disp = disp + 4*(ss_blockid_specs(1,ss_index)*3+1)
        deallocate(iiobuf_plot3d)
#endif ! SS_MULTIBLOCK

        call MPI_TYPE_CONTIGUOUS(4, SS_MPI_PREC, contig4,ierr)
        call MPI_TYPE_COMMIT(contig4, ierr)

        riobuf_4(1)  = SS_PREC_CONV(amach1)
        riobuf_4(2)  = SS_PREC_CONV(0.0d0)
        riobuf_4(3)  = SS_PREC_CONV(reoo)
        riobuf_4(4)  = SS_PREC_CONV(timo)

        call MPI_FILE_SET_VIEW(ss_fh(ss_index), disp,SS_MPI_PREC,contig4,'native', MPI_INFO_NULL, ierr)

        ! IF-block added by Jon Gibson, NAG
        if (ss_rank(ss_index) .eq. 0) then
          call MPI_FILE_WRITE(ss_fh(ss_index), riobuf_4,4, SS_MPI_PREC, status, ierr)
        endif

        disp = disp + SS_PREC_BYTE*4

#ifdef SS_MULTIBLOCK
        !set offset for each block
        disp = disp + ss_num_var(ss_index) * SS_PREC_BYTE * ss_blockid_offset(ss_index)
#endif ! SS_MULTIBLOCK

        call MPI_FILE_SET_VIEW(ss_fh(ss_index), disp, SS_MPI_PREC,ss_dim_block(ss_index,4), 'native', MPI_INFO_NULL, ierr)

        iread =      int(ss_idim,kind=MPI_OFFSET_KIND) &
                    *int(ss_jdim,kind=MPI_OFFSET_KIND) &
                    *int(ss_kdim,kind=MPI_OFFSET_KIND) &
                    *int(ss_num_var(ss_index),kind=MPI_OFFSET_KIND)
#ifdef SS_BLOCKING_WRITE
        call MPI_FILE_WRITE_ALL(ss_fh(ss_index), ss_qio, iread,SS_MPI_PREC, status, ierr)
        call MPI_FILE_CLOSE(ss_fh(ss_index),ierr)
#elif SS_NON_COLLECTIVE_WRITE
        call MPI_FILE_IWRITE(ss_fh(ss_index), ss_qio(ss_start), iread,SS_MPI_PREC, ss_status(ss_index,1), ierr)
        ss_first_write(ss_index)=.true.
#else
        call MPI_File_write_all_begin(ss_fh(ss_index), ss_qio(ss_start), iread   , SS_MPI_PREC,ierr)
        ss_first_write(ss_index)=.true.
#endif

        call MPI_TYPE_FREE(contig4, ierr)
        call MPI_TYPE_FREE(contig_plot3d, ierr)

#ifdef SS_BLOCKING_WRITE
        deallocate (ss_qio)
#else
        !set start index for iobuffer
        ss_start = ss_start + ss_idim*ss_jdim*ss_kdim*ss_num_var(ss_index)
#endif
      end if
      enddo
!cpb   ==================================================================
!cpb                    END OUTPUT IN PLOT 3D FORMAT
!cpb   ==================================================================
      return
      end


      subroutine subspace_write_grid()
!cpb   =================================================================
!cpb   This routine writes out
!cpb   subspaces of the computational domain. The specification of these
!cpb   subdomains is done via the input file 'SUBSPACES_INPUT.dat' which
!cpb   has to be in the same directory as the executable.
!cpb
!cpb   The specification is grid dependent!
!cpb
!cpb   Bechlars, November 2012
!cpb   =================================================================

!cpb   =================================================================
!cpb                             load modules
!cpb   =================================================================
      use SSVAR
      use mpi
      use subroutineso
      use subroutines3d
      use problemcase
      implicit none
!cpb   =================================================================
!cpb                         subroutine arguments
!cpb   =================================================================

!cpb   =================================================================
!cpb                  additional variables/parameters
!cpb   =================================================================
      integer iii,jjj,kkk,ss_index

      !MPI/PLOT3D environment
      integer accmod,ialloc(7)
      integer (kind=MPI_OFFSET_KIND) :: disp,ioff,iread,l2,ss_start
      integer totaldims_ijk(4),blockdims_ijk(4),blockstart_ijk(4)
      integer subarray_ijk
      integer status(MPI_STATUS_SIZE), contig_plot3d
      integer real_mp_size
      real(kind=wp) xxl,xnit,tw
      real(kind=wp), save ::  riobuf(4)
      integer*8, save :: iiobuf(4)
      real(kind=SS_PREC), allocatable :: ss_qio_grid(:)
      real(kind=wp) ::  riobuf_4(4)
#ifndef SS_MULTIBLOCK
      integer*4 :: iiobuf_plot3d(3)
#else ! SS_MULTIBLOCK
      integer*4, allocatable  :: iiobuf_plot3d(:)
#endif ! SS_MULTIBLOCK
      character*8 chl
      character*70 iofile
      character*70 linux_command
      logical lex
      equivalence (riobuf,iiobuf)
      integer ss_idim,ss_jdim,ss_kdim
      integer ss_istart,ss_iend,ss_iskip
      integer ss_jstart,ss_jend,ss_jskip
      integer ss_kstart,ss_kend,ss_kskip


!cpb   =================================================================
!cpb                          executing part
!cpb   =================================================================

      ss_start = 1
      do ss_index=1,ss_num
!cpb   ONLY DO THIS FOR CORES THAT ARE WITHIN VOLUME
      if((ss_in_subspace(ss_index))) then
        !name and directory of output files
        linux_command='mkdir -p FLOW'
        if(ss_rank(ss_index) .eq. 0) CALL system(linux_command)
        !check if ss_path variable is set
        if(ss_path.eq.'') then
         if(ss_rank(ss_index) .eq. 0) write(*,*) 'ss_path is not set => subvolume output in main directory'
        else
          !create new directory if not already existing
          linux_command='mkdir -p '//ss_path
          if(ss_rank(ss_index) .eq. 0) CALL system(linux_command)
          iofile=trim(ss_tag_grid(ss_index ))
        endif
!cpb   =================================================================
!cpb                  set the limits for each block and proc
!cpb   =================================================================
        ss_idim=ss_dim_proc(ss_index,1)
        ss_jdim=ss_dim_proc(ss_index,2)
        ss_kdim=ss_dim_proc(ss_index,3)

        totaldims_ijk(1)  = ss_dim_block(ss_index,1)
        totaldims_ijk(2)  = ss_dim_block(ss_index,2)
        totaldims_ijk(3)  = ss_dim_block(ss_index,3)
        totaldims_ijk(4)  = 3
        blockdims_ijk(1)  = ss_idim
        blockdims_ijk(2)  = ss_jdim
        blockdims_ijk(3)  = ss_kdim
        blockdims_ijk(4)  = 3
        blockstart_ijk(1) = ss_dim_proc(ss_index,4)
        blockstart_ijk(2) = ss_dim_proc(ss_index,5)
        blockstart_ijk(3) = ss_dim_proc(ss_index,6)
        blockstart_ijk(4) = 0

        !limits in i-direction
        ss_istart =ss_dim_proc(ss_index,7)
        ss_iend =ss_dim_proc(ss_index,10)
        ss_iskip =ss_spec_block(ss_index,mb+1,3)

        !limits in j-direction
        ss_jstart =ss_dim_proc(ss_index,8)
        ss_jend =ss_dim_proc(ss_index,11)
        ss_jskip =ss_spec_block(ss_index,mb+1,6)

        !limits in k-direction
        ss_kstart =ss_dim_proc(ss_index,9)
        ss_kend =ss_dim_proc(ss_index,12)
        ss_kskip =ss_spec_block(ss_index,mb+1,9)
!cpb   =================================================================
!cpb                  set subarray for each subvolume
!cpb   =================================================================
         call MPI_TYPE_CREATE_SUBARRAY(4, totaldims_ijk, &
          blockdims_ijk, blockstart_ijk, MPI_ORDER_FORTRAN, &
          SS_MPI_PREC, subarray_ijk, ierr)
         call MPI_TYPE_COMMIT( subarray_ijk, ierr)
  
!cpb   ==================================================================
!cpb                      OUTPUT IN PLOT 3D FORMAT
!cpb   ==================================================================
        !allocate writing vector
        allocate (ss_qio_grid(int(ss_idim,kind=MPI_OFFSET_KIND) &
                             *int(ss_jdim,kind=MPI_OFFSET_KIND) &
                             *int(ss_kdim,kind=MPI_OFFSET_KIND) &
                             *int(3,kind=MPI_OFFSET_KIND)),stat=ialloc(1))
        !check if allocation went fine
        if(ialloc(1) .ne. 0)then
          write(*,*)'allocate error (ss_volume):',ialloc
          call MPI_FINALIZE(ierr)
          stop
        endif

        l2=0
        !write data in vector
        !loop over values
        !loop theta direction
        do k=ss_kstart,ss_kend,ss_kskip
        !loop eta direction
        do j=ss_jstart,ss_jend,ss_jskip
          !loop xi direction
          do i=ss_istart,ss_iend,ss_iskip
            ! Q=Q/J
            l2=l2+1
            ss_qio_grid((ss_start-1)+l2) = SS_PREC_CONV(ss(indx3(i-1,j-1,k-1,1),1))
           end do
        end do
        end do
        !loop theta direction
        do k=ss_kstart,ss_kend,ss_kskip
        !loop eta direction
        do j=ss_jstart,ss_jend,ss_jskip
          !loop xi direction
          do i=ss_istart,ss_iend,ss_iskip
            ! Q=Q/J
            l2=l2+1
            ss_qio_grid((ss_start-1)+l2) = SS_PREC_CONV(ss(indx3(i-1,j-1,k-1,1),2))
           end do
        end do
        end do
        !loop theta direction
        do k=ss_kstart,ss_kend,ss_kskip
        !loop eta direction
        do j=ss_jstart,ss_jend,ss_jskip
          !loop xi direction
          do i=ss_istart,ss_iend,ss_iskip
            ! Q=Q/J
            l2=l2+1
#ifdef THREE_D
            ss_qio_grid((ss_start-1)+l2) = SS_PREC_CONV(ss(indx3(i-1,j-1,k-1,1),3))
#else
            ss_qio_grid((ss_start-1)+l2) = 0.e0
#endif
           end do
        end do
        end do


        ioff = 0
        disp = 0

        accmod = ior(MPI_MODE_WRONLY,MPI_MODE_CREATE)

        !if rank in subvolume communicator = 0 then overwrite existing file
        if(ss_rank(ss_index).eq.0)then
          inquire(file=iofile,exist=lex)
          if(lex)then
            call MPI_FILE_DELETE(iofile,MPI_INFO_NULL,ierr)
          endif
        endif

        !set iobuffer for plot3d

#ifndef SS_MULTIBLOCK
        iiobuf_plot3d(1)  =  totaldims_ijk(1)
        iiobuf_plot3d(2)  =  totaldims_ijk(2)
        iiobuf_plot3d(3)  =  totaldims_ijk(3)
#else ! SS_MULTIBLOCK
        allocate(iiobuf_plot3d(ss_blockid_specs(1,ss_index)*3+1))
        iiobuf_plot3d(1)  = ss_blockid_specs(1,ss_index)
        do i=1, ss_blockid_specs(1,ss_index)*3
            iiobuf_plot3d(1+i)  =  ss_blockid_specs(1+i,ss_index)
        enddo
#endif ! SS_MULTIBLOCK

#ifndef SS_MULTIBLOCK
        call MPI_TYPE_CONTIGUOUS(3, MPI_INTEGER4,contig_plot3d,ierr)
#else ! SS_MULTIBLOCK
        call MPI_TYPE_CONTIGUOUS(iiobuf_plot3d(1)*3+1, MPI_INTEGER4,contig_plot3d,ierr)
#endif ! SS_MULTIBLOCK
        call MPI_TYPE_COMMIT(contig_plot3d, ierr)
        call MPI_COMM_SIZE(block_comm,jjj,ierr)
        call MPI_COMM_SIZE(ss_communicator(ss_index),jjj,ierr)
        call MPI_FILE_OPEN(ss_communicator(ss_index),iofile,accmod,MPI_INFO_NULL,ss_fh(ss_index),ierr)
        call MPI_FILE_SET_VIEW(ss_fh(ss_index), disp, MPI_INTEGER4,contig_plot3d, 'native', MPI_INFO_NULL, ierr)

        ! IF-block added by Jon Gibson, NAG
        if (ss_rank(ss_index) .eq. 0) then
!cpb   changed from write_AT(disp) to write()
#ifndef SS_MULTIBLOCK
          call MPI_FILE_WRITE(ss_fh(ss_index), iiobuf_plot3d,3, MPI_INTEGER4, status, ierr)
#else ! SS_MULTIBLOCK
          call MPI_FILE_WRITE(ss_fh(ss_index),  iiobuf_plot3d,iiobuf_plot3d(1)*3+1, MPI_INTEGER4, status, ierr)
#endif ! SS_MULTIBLOCK
        endif

#ifndef SS_MULTIBLOCK
        disp = disp + 4*3
#else ! SS_MULTIBLOCK
        disp = disp + 4*(ss_blockid_specs(1,ss_index)*3+1)
        deallocate(iiobuf_plot3d)
#endif ! SS_MULTIBLOCK

#ifdef SS_MULTIBLOCK
        !set offset for each block
        disp = disp + 3 * SS_PREC_BYTE * ss_blockid_offset(ss_index)
#endif ! SS_MULTIBLOCK

        call MPI_FILE_SET_VIEW(ss_fh(ss_index), disp, MPI_REAL4,subarray_ijk, 'native', MPI_INFO_NULL, ierr)

        iread =  int(ss_idim,kind=MPI_OFFSET_KIND) &
                *int(ss_jdim,kind=MPI_OFFSET_KIND) &
                *int(ss_kdim,kind=MPI_OFFSET_KIND) &
                *int(3,kind=MPI_OFFSET_KIND)
         call MPI_FILE_WRITE_ALL(ss_fh(ss_index), ss_qio_grid(ss_start) &
         , iread, SS_MPI_PREC, status, ierr)
       call MPI_FILE_CLOSE(ss_fh(ss_index),ierr)


        call MPI_TYPE_FREE(subarray_ijk, ierr)
        call MPI_TYPE_FREE(contig_plot3d, ierr)
        deallocate (ss_qio_grid)

      end if
      enddo
!cpb   ==================================================================
!cpb                    END OUTPUT IN PLOT 3D FORMAT
!cpb   ==================================================================
      return
      end


!cpb   =================================================================
!cpb   =================================================================
!cpb               END IO subroutine to write out subspaces
!cpb   =================================================================
!cpb   =================================================================


