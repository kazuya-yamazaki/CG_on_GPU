      program heat3Dp

      use solver11
      use pfem_util
      use nvtxwrap

      implicit REAL*8(A-H,O-Z)
!C
!C +-------+
!C | INIT. |
!C +-------+
!C=== 
!$acc kernels
      dummy = 0.0
!$acc end kernels
      call mynvtxStartRange("INIT_INPUT")
      call PFEM_INIT
      call INPUT_CNTL
      call INPUT_GRID
      call mynvtxEndRange()
!C===

!C
!C +---------------------+
!C | matrix connectivity |
!C +---------------------+
!C===
      START_TIME= MPI_WTIME()

      call mynvtxStartRange("MAT_CON0")
      call MAT_CON0
      call mynvtxEndRange()
      call mynvtxStartRange("MAT_CON1")
      call MAT_CON1
      call mynvtxEndRange()

      END_TIME= MPI_WTIME()
      if (my_rank.eq.0) then
        write (*, '("*** matrix conn. ", 1pe16.6, " sec.")')            &
     &         END_TIME-START_TIME
      endif
!C===
      
!C
!C +-----------------+
!C | MATRIX assemble |
!C +-----------------+
!C===
      START_TIME= MPI_WTIME()
      
      call mynvtxStartRange("MAT_ASS_MAIN")
!$acc data copyin(index, item)
      call MAT_ASS_MAIN
      call mynvtxEndRange()
      call mynvtxStartRange("MAT_ASS_BC")
      call MAT_ASS_BC
      call mynvtxEndRange()

      END_TIME= MPI_WTIME()
      if (my_rank.eq.0) then
        write (*, '("*** matrix ass.  ", 1pe16.6, " sec.",/)')          &
     &         END_TIME-START_TIME
      endif
!C===

!C
!C +--------+
!C | SOLVER |
!C +--------+
!C===
      START_TIME= MPI_WTIME()
      call mynvtxStartRange("SOLVE11")
      call SOLVE11
!$acc end data
!$acc exit data delete(B, D, AMAT) copyout(X)
      call mynvtxEndRange()
      END_TIME= MPI_WTIME()

      if (my_rank.eq.0) then
        write (*, '("*** real  COMP.  ", 1pe16.6, " sec.",/)')          &
     &         END_TIME-START_TIME
      endif
!C===

!C
!C +--------+
!C | OUTPUT |
!C +--------+
!C===

!      call OUTPUT_UCD
      do i= 1, N
        if (XYZ(i,1).eq.0.d0.and.XYZ(i,2).eq.0.d0.
     &                       and.XYZ(i,3).eq.0.d0) then
          write (*,'("Value on edge (",2i8,")=",1pe16.6)')              &
     &            my_rank, i, X(i)
        endif
      enddo

!C===

      call PFEM_FINALIZE

      end program heat3Dp

!C
!C***
!C*** ERROR_EXIT
!C***
!C
      subroutine ERROR_EXIT (IFLAG, my_rank)

      if (my_rank.eq.0) then

      if (IFLAG.eq. 101) then
        write (*,'(a)') '#### PFEM-SOL-E0101: PEsmpTOT must be >0'
        call MPI_FINALIZE (errno)
        stop
      endif

      if (IFLAG.eq. 102) then
        write (*,'(a)') '#### PFEM-SOL-E0102: METHOD must be 1 or 2'
        call MPI_FINALIZE (errno)
        stop
      endif

      if (IFLAG.eq. 103) then
        write (*,'(a)') '#### PFEM-SOL-E0103: PRECOND must be 0 or 1'
        call MPI_FINALIZE (errno)
        stop
      endif

      if (IFLAG.eq. 104) then
        write (*,'(a)') '#### PFEM-SOL-E0104: ITER must be >0'
        call MPI_FINALIZE (errno)
        stop
      endif

      if (IFLAG.eq. 111) then
        write (*,'(a)') '#### PFEM-SOL-W0111: iterPREmax must be >=1'
        write (*,'(a)') '                     iterPREmax set to   =1'
        return
      endif

      if (IFLAG.eq. 112) then
        write (*,'(a)') '#### PFEM-SOL-W0112: iterPREmax must be =<4'
        write (*,'(a)') '                     iterPREmax set to   =4'
        return
      endif

      endif

      if (IFLAG.eq. 201) then
        write (*,'(a)') '#### PFEM-SOL-E0201: too many PEs specified'
        write (*,'(a)') '                       in MPIRUN.'
        call MPI_ABORT (errno)
        stop
      endif

      if (IFLAG.eq. 202) then
        write (*,'(a,i8)') '#### PFEM-SOL-E0202: invalid mesh data',    &
     &                                           my_rank
        call MPI_FINALIZE (errno)
        stop
      endif

      return
      end
      

