      module SOLVER11

      contains
      subroutine SOLVE11

        use pfem_util
        use solver_CG

        implicit REAL*8 (A-H,O-Z)

        integer :: ERROR, ICFLAG
        character(len=char_length) :: BUF

        data ICFLAG/0/

!C
!C +------------+
!C | PARAMETERs |
!C +------------+
!C===
        ITER = pfemIarray(1)
        RESID= pfemRarray(1)
!C===

!C
!C +------------------+
!C | ITERATIVE solver |
!C +------------------+
!C===
        call CG                                                         &
     &     ( N, NP, NPLU, D, AMAT, index, item, B, X,  RESID,           &
     &       ITER, ERROR, my_rank,                                      &
     &       NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,              &
     &                          EXPORT_INDEX, EXPORT_ITEM) 
        ITERactual= ITER
!C===

      end subroutine SOLVE11
      end module SOLVER11
