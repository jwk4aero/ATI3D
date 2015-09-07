      subroutine stopnow(msg)
!C
!C  stops program after flushing and closing output files
!C
      use mpi
      use subroutineso
      use subroutines3d
      use problemcase
      character*(*) msg

!crp   write output in red
      write (*,'(a,i5,a,i5,a,i5,a,i5)') 'procid=', myid,   &   
           ' ,blockid=', mb, ' ,sx=',sx, ' ,sy=',sy
      
            print *, achar(27)//'[31m'// msg //achar(27)//'[0m'
!c      write (*,'(/a)') msg
!c  add other files to be closed as appropriate
!      if(nodeout) close(nodeio)

!C******************************
!C  Close down all processors  *
!C******************************
#ifdef MPI
      call MPI_abort(icom,icode,ierr)
#endif

      stop 'Program terminated abnormally'

      end
