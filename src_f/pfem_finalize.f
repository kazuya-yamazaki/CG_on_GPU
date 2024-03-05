!C
!C***
!C*** PFEM_FINALIZE
!C***
!C
      subroutine PFEM_FINALIZE
      use pfem_util
      implicit REAL*8 (A-H,O-Z)

      call MPI_FINALIZE (errno)
      if (my_rank.eq.0) stop ' * normal termination'

      return
      end
