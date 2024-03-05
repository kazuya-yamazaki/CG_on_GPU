!C
!C***
!C*** MAT_CON1
!C***
!C
      subroutine MAT_CON1
      use pfem_util
      implicit REAL*8 (A-H,O-Z)

      allocate (index(0:NP))
      index= 0

      do i= 1, NP
        index(i)= index(i-1) + INLU(i)
      enddo

      NPLU= index(NP)

      allocate (item(NPLU))
      if (NFLAG.eq.0) then
        item= 0
       else
!$omp parallel do private(j,k)
        do j= 1, N
          do k= index(j-1)+1, index(j)
            item(k)= 0
          enddo
        enddo
      endif

      do i= 1, NP
        do k= 1, INLU(i)
               kk = k + index(i-1)
          item(kk)=      IALU(i,k)
        enddo
      enddo

      deallocate (INLU, IALU)

      end subroutine MAT_CON1
