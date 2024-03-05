!C
!C*** 
!C*** module solver_SR
!C***
!C
      module solver_SR
      contains
!C
!C*** SOLVER_SEND_RECV
!C
      subroutine  SOLVER_SEND_RECV                                      &
     &            ( N0, N, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,&
     &                                        EXPORT_INDEX, EXPORT_ITEM,&
     &                  WS, WR, X, my_rank)

      implicit REAL*8 (A-H,O-Z)
      include  'mpif.h'
      include  'precision.inc'

      integer(kind=kint )                , intent(in)   ::  N0, N
      integer(kind=kint )                , intent(in)   ::  NEIBPETOT
      integer(kind=kint ), pointer :: NEIBPE      (:)
      integer(kind=kint ), pointer :: IMPORT_INDEX(:)
      integer(kind=kint ), pointer :: IMPORT_ITEM  (:)
      integer(kind=kint ), pointer :: EXPORT_INDEX(:)
      integer(kind=kint ), pointer :: EXPORT_ITEM  (:)
      real   (kind=kreal), dimension(N), intent(inout):: WS
      real   (kind=kreal), dimension(N), intent(inout):: WR
      real   (kind=kreal), dimension(N), intent(inout):: X
      integer                          , intent(in)   :: my_rank

      integer(kind=kint ), dimension(:,:), save, allocatable :: sta
      integer(kind=kint ), dimension(:  ), save, allocatable :: req

      integer(kind=kint ), save :: NFLAG
      data NFLAG/0/

!C
!C-- INIT.
      if (NFLAG.eq.0) then
        allocate (sta(MPI_STATUS_SIZE,2*NEIBPETOT))
        allocate (req(2*NEIBPETOT))
        NFLAG= 1
      endif

!$acc  data present(IMPORT_ITEM, EXPORT_ITEM, WS, WR, X)
       
!C
!C-- SEND
      do neib= 1, NEIBPETOT
        istart= EXPORT_INDEX(neib-1)
        inum  = EXPORT_INDEX(neib  ) - istart
!$omp parallel do private(k,ii)
!$acc kernels
        do k= istart+1, istart+inum
           ii   = EXPORT_ITEM(k)
           WS(k)= X(ii)
        enddo
!$acc end kernels
      enddo

      do neib= 1, NEIBPETOT
        istart= EXPORT_INDEX(neib-1)
        inum  = EXPORT_INDEX(neib  ) - istart
!$acc host_data use_device(WS)
        call MPI_Isend (WS(istart+1), inum, MPI_DOUBLE_PRECISION,         &
     &                  NEIBPE(neib), 0, MPI_COMM_WORLD, req(neib),       &
     &                  ierr)
!$acc end host_data
      enddo

!C
!C-- RECEIVE
      do neib= 1, NEIBPETOT
        inum  = IMPORT_INDEX(neib  ) - IMPORT_INDEX(neib-1)
        istart= IMPORT_INDEX(neib-1) + N0 + 1
!$acc host_data use_device(X)
        call MPI_Irecv (X(istart), inum, MPI_DOUBLE_PRECISION,            &
     &     NEIBPE(neib), 0, MPI_COMM_WORLD, req(NEIBPETOT+neib),          &
     &     ierr)
!$acc end host_data
      enddo

      call MPI_Waitall (2*NEIBPETOT, req, sta, ierr)

!$acc end data

      end subroutine solver_send_recv
      end module     solver_SR
