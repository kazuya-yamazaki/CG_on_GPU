!C
!C***
!C*** PFEM_INIT
!C***
!C
!C    INIT. PFEM-FEM process's
!C
      subroutine PFEM_INIT
      use pfem_util
#ifdef OPENACC
      use openacc
#endif
      implicit REAL*8 (A-H,O-Z)
      integer ngpu, ierr

      call MPI_INIT      (ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD, PETOT, ierr )
      call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr )
      
#ifdef OPENACC
      ngpu = acc_get_num_devices(acc_device_nvidia)
      call acc_set_device_num(mod(my_rank, ngpu), acc_device_nvidia)
#endif

      pfemRarray= 0.d0
      pfemIarray= 0

      return
      end
