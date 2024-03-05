!C
!C***
!C*** INPUT_CNTL
!C***
!C
      subroutine INPUT_CNTL
      use pfem_util

      implicit REAL*8 (A-H,O-Z)

      if (my_rank.eq.0) then
        open (11,file= 'INPUT.DAT', status='unknown')
        read (11,'(a80)') HEADER
        read (11,*) ITER
        read (11,*) COND, QVOL
        read (11,*) RESID
        read (11,*) NFLAG
        close (11)
      endif

      call MPI_BCAST (HEADER, 80, MPI_CHARACTER, 0, MPI_COMM_WORLD,ierr)
      call MPI_BCAST (ITER   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST (COND   , 1, MPI_DOUBLE_PRECISION, 0, 
     &                                            MPI_COMM_WORLD, ierr)
      call MPI_BCAST (QVOL   , 1, MPI_DOUBLE_PRECISION, 0, 
     &                                            MPI_COMM_WORLD, ierr)
      call MPI_BCAST (RESID  , 1, MPI_DOUBLE_PRECISION, 0, 
     &                                            MPI_COMM_WORLD, ierr)

      pfemRarray(1)= RESID
      pfemIarray(1)= ITER

      return
      end
