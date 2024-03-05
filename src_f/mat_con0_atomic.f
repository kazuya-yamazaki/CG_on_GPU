!C
!C***
!C*** MAT_CON0
!C***
!C
      subroutine MAT_CON0
      use pfem_util
      use util
      use nvtxwrap
      implicit REAL*8 (A-H,O-Z)

      NLU= 26

      allocate (INLU(NP), IALU(NP,NLU))

      call mynvtxStartRange("NODES")
!$acc data copyout(INLU, IALU) copyin(ICELNOD)
!$acc kernels
      INLU(:) = 0
      IALU(:,:) = 0
!$acc loop independent collapse(3) private(in1, in2)
      do icel= 1, ICELTOT
        do i1 = 1, 8
          do i2 = 1, 8
            if (i1 == i2) cycle
            in1 = ICELNOD(icel,i1)
            in2 = ICELNOD(icel,i2)
            call FIND_TS_NODE (INLU,IALU,NP,NLU,in1,in2)
          enddo
        enddo
      enddo
!$acc end kernels
!$acc end data
      call mynvtxEndRange

      call mynvtxStartRange("sort")
      do in= 1, N
        NN= INLU(in)
        do k= 1, NN
          NCOL1(k)= IALU(in,k)
        enddo
          call mSORT (NCOL1, NCOL2, NN)
        do k= NN, 1, -1
          IALU(in,NN-k+1)= NCOL1(NCOL2(k))
        enddo        
      enddo
      call mynvtxEndRange

      contains
!C
!C***
!C*** FIND_TS_NODE
!C***
!C
      subroutine FIND_TS_NODE (INLU,IALU,NP,NLU,ip1,ip2)
!$acc routine seq
      implicit none
      
      integer, intent(in) :: NP, NLU, ip1, ip2
      integer, intent(inout) :: INLU(NP), IALU(NP, NLU)

      ! local variables
      integer :: inlu_ip1, kk, ks, icou, ind
      logical :: retry

      ks = 1
      do while (.true.)

!$acc atomic read
        inlu_ip1 = INLU(ip1)
!$acc end atomic
        if (inlu_ip1 > NLU) then
          ! If multiple threads increment INLU(ip1) at the same time, it may
          ! exceed the upper limit. In that case, start over.
          cycle
        endif
        
        retry = .false.
        do kk= ks, inlu_ip1
!$acc atomic read
          ind = IALU(ip1,kk)
!$acc end atomic
          if (ind == 0) then
            ! It is possible that I might read the array value before other thread
            ! store the value in IALU(ip1,kk). In that case, start over.
            retry = .true.
            exit
          else if (ip2.eq.ind) then
            return
          end if
          ks = kk + 1
        enddo
        if (retry) cycle

!$acc atomic capture
        INLU(ip1) = INLU(ip1) + 1
        icou = INLU(ip1)
!$acc end atomic

        if (icou .ne. inlu_ip1 + 1) then
          ! Other threads may have incremented INLU(ip1).
          ! In that case, start over.
!$acc atomic update
          INLU(ip1) = INLU(ip1) - 1
!$acc end atomic
          cycle
        end if
        
!$acc atomic write
        IALU(ip1,icou)= ip2
!$acc end atomic
        exit

      end do

      return

      end subroutine FIND_TS_NODE
      end subroutine MAT_CON0
