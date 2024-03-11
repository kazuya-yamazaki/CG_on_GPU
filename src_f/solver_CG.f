!C
!C*** 
!C*** module solver_CG
!C***
!C
      module solver_CG
      contains
!C
!C*** CG
!C
!C    CG solves the linear system Ax = b using the Conjugate Gradient 
!C    iterative method with the following preconditioners
!C
      subroutine CG                                                     &
     &   (N, NP, NPLU, D, AMAT, index, item, B, X, RESID,               &
     &    ITER, ERROR, my_rank,                                         &
     &    NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,                 &
     &                       EXPORT_INDEX, EXPORT_ITEM)

      use  solver_SR

      implicit REAL*8(A-H,O-Z)
      include  'precision.inc'
      include  'mpif.h'

      integer(kind=kint ), intent(in):: N, NP, NPLU, my_rank
      integer(kind=kint ), intent(in):: NEIBPETOT
      integer(kind=kint ), intent(inout):: ITER, ERROR
      real   (kind=kreal), intent(inout):: RESID

      real(kind=kreal), dimension(NP)   , intent(inout):: B, X, D
      real(kind=kreal), dimension(NPLU), intent(inout):: AMAT

      integer(kind=kint ), dimension(0:NP),intent(in) :: index
      integer(kind=kint ), dimension(NPLU),intent(in) :: item

      integer(kind=kint ), pointer :: NEIBPE(:)
      integer(kind=kint ), pointer :: IMPORT_INDEX(:), IMPORT_ITEM(:)
      integer(kind=kint ), pointer :: EXPORT_INDEX(:), EXPORT_ITEM(:)

      real(kind=kreal), dimension(:),    allocatable:: WS, WR
      real(kind=kreal), dimension(:,:),  allocatable:: WW

      integer(kind=kint), parameter ::  R= 1
      integer(kind=kint), parameter ::  Z= 2
      integer(kind=kint), parameter ::  Q= 2
      integer(kind=kint), parameter ::  P= 3
      integer(kind=kint), parameter :: DD= 4

      integer(kind=kint ) :: MAXIT
      real   (kind=kreal) :: TOL, W, SS
      
      integer(kind=kint) :: i,j,k

!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      COMMtime= 0.d0
      COMPtime= 0.d0

      ERROR= 0

      allocate (WW(NP,4),WR(NP),WS(NP))

      MAXIT  = ITER
       TOL   = RESID           

!$acc  data create(WW, WS, WR) present(B, D, X, AMAT)
!$acc+      copyin(IMPORT_ITEM, EXPORT_ITEM)

!$omp parallel do private(i)
!$acc kernels
      do i= 1, N
        WW(i,1)= 0.d0
        WW(i,2)= 0.d0
        WW(i,3)= 0.d0
        WW(i,4)= 0.d0
        WS(i  )= 0.d0
        WR(i  )= 0.d0
        X (i  )= 0.d0
      enddo

!$omp parallel do private(i)
      do i= N+1, NP
        WW(i,1)= 0.d0
        WW(i,2)= 0.d0
        WW(i,3)= 0.d0
        WW(i,4)= 0.d0
        WS(i  )= 0.d0
        WR(i  )= 0.d0
        X (i  )= 0.d0
      enddo
!$acc end kernels
!C===

!C
!C +-----------------------+
!C | {r0}= {b} - [A]{xini} |
!C +-----------------------+
!C===

!C
!C-- INTERFACE data EXCHANGE
      call SOLVER_SEND_RECV                                             &
     &   ( N, NP, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,         &
     &     EXPORT_INDEX, EXPORT_ITEM, WS, WR, X , my_rank)

!$omp parallel do private(j,k,i,WVAL)
!$acc kernels
      do j= 1, N
        WW(j,DD)= 1.d0/D(j)
        WVAL= B(j) - D(j)*X(j)
        do k= index(j-1)+1, index(j)
          i= item(k)
          WVAL= WVAL - AMAT(k)*X(i)
        enddo
        WW(j,R)= WVAL
      enddo
!$acc end kernels


      BNRM20= 0.d0
!$omp parallel do private(i) reduction (+:BNRM20)
!$acc kernels loop           reduction (+:BNRM20)
      do i= 1, N
        BNRM20= BNRM20 + B(i)**2
      enddo
!$acc end kernels
      call MPI_Allreduce (BNRM20, BNRM2, 1, MPI_DOUBLE_PRECISION,       &
     &                    MPI_SUM, MPI_COMM_WORLD, ierr)

      if (BNRM2.eq.0.d0) BNRM2= 1.d0
      ITER = 0
      
!C===

      do iter= 1, MAXIT
!C
!C************************************************* Conjugate Gradient Iteration

!C
!C +----------------+
!C | {z}= [Minv]{r} |
!C +----------------+
!C===

!$omp parallel do private(i)
!$acc kernels
      do i= 1, N
        WW(i,Z)= WW(i,R) * WW(i,DD)
      enddo
!$acc end kernels
!C===
      
!C
!C +---------------+
!C | {RHO}= {r}{z} |
!C +---------------+
!C===
      RHO0= 0.d0

!$omp parallel do private(i) reduction (+:RHO0)
!$acc kernels loop           reduction (+:RHO0)
      do i= 1, N
        RHO0= RHO0 + WW(i,R)*WW(i,Z)
      enddo
!$acc end kernels

      call MPI_Allreduce (RHO0, RHO, 1, MPI_DOUBLE_PRECISION,       
     &                    MPI_SUM, MPI_COMM_WORLD, ierr)
!C===

!C
!C +-----------------------------+
!C | {p} = {z} if      ITER=1    |
!C | BETA= RHO / RHO1  otherwise |
!C +-----------------------------+
!C===
      if ( ITER.eq.1 ) then

!$omp parallel do private(i)
!$acc kernels
        do i= 1, N
          WW(i,P)= WW(i,Z)
        enddo
!$acc end kernels
       else
         BETA= RHO / RHO1
!$omp parallel do private(i)
!$acc kernels
         do i= 1, N
           WW(i,P)= WW(i,Z) + BETA*WW(i,P)
         enddo
!$acc end kernels
      endif
!C===

!C
!C +-------------+
!C | {q}= [A]{p} |
!C +-------------+
!C===        
!C
!C-- INTERFACE data EXCHANGE
      call SOLVER_SEND_RECV                                             &
     &   ( N, NP, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,         &
     &     EXPORT_INDEX, EXPORT_ITEM, WS, WR, WW(1,P), my_rank)

!$omp parallel do private(j,k,i,WVAL)
!$acc kernels loop private(i,WVAL)
      do j= 1, N
        WVAL= D(j)*WW(j,P)
!$acc loop seq
        do k= index(j-1)+1, index(j)
           i= item(k)
          WVAL= WVAL + AMAT(k)*WW(i,P)
        enddo
        WW(j,Q)= WVAL
      enddo
!$acc end kernels
!C===

!C
!C +---------------------+
!C | ALPHA= RHO / {p}{q} |
!C +---------------------+
!C===
      C10= 0.d0
!$omp parallel do private(i) reduction (+:C10)
!$acc kernels loop           reduction (+:C10)
      do i= 1, N
        C10= C10 + WW(i,P)*WW(i,Q)
      enddo
!$acc end kernels
      call MPI_Allreduce (C10, C1, 1, MPI_DOUBLE_PRECISION,     
     &                    MPI_SUM, MPI_COMM_WORLD, ierr)

      ALPHA= RHO / C1
!C===

!C
!C +----------------------+
!C | {x}= {x} + ALPHA*{p} |
!C | {r}= {r} - ALPHA*{q} |
!C +----------------------+
!C===

!$omp parallel do private(i)
!$acc kernels
      do i= 1, N
         X(i)  = X (i)   + ALPHA * WW(i,P)
        WW(i,R)= WW(i,R) - ALPHA * WW(i,Q)
      enddo
!$acc end kernels

      DNRM20= 0.d0
!$omp parallel do private(i) reduction (+:DNRM20)
!$acc kernels loop           reduction (+:DNRM20)
      do i= 1, N
        DNRM20= DNRM20 + WW(i,R)**2
      enddo
!$acc end kernels
      call MPI_Allreduce (DNRM20, DNRM2, 1, MPI_DOUBLE_PRECISION,     
     &                    MPI_SUM, MPI_COMM_WORLD, ierr)

      RESID= dsqrt(DNRM2/BNRM2)

!C##### ITERATION HISTORY
      if (my_rank.eq.0) then
        write (*, 1000) ITER, RESID
      endif
 1000   format (i5, 1pe16.6)
! 1010   format (1pe16.6)
!C#####

        if ( RESID.le.TOL   ) exit
        if ( ITER .eq.MAXIT ) ERROR= -300
        
        RHO1 = RHO                                                             
      enddo
!C===

   30 continue


!C
!C-- INTERFACE data EXCHANGE
      call SOLVER_SEND_RECV                                             &
     &   ( N, NP, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,         &
     &     EXPORT_INDEX, EXPORT_ITEM, WS, WR, X , my_rank)

!$acc end data

      deallocate (WW,WR,WS)

      end subroutine        CG
      
      end module     solver_CG
