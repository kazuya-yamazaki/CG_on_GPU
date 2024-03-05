!C
!C***
!C*** MAT_ASS_BC
!C***
!C
      subroutine MAT_ASS_BC
      use pfem_util
      implicit REAL*8 (A-H,O-Z)

      allocate (IWKX(NP,2))

!C
!C== Z=Zmax

!$acc data create(IWKX) present(B,D,AMAT, index,item)
!$acc& copyin(NODGRP_INDEX)

!$omp parallel do private(in)
!$acc kernels
      IWKX(:,:)= 0
!$acc end kernels

      ib0= -1
      do ib0= 1, NODGRPtot
        if (NODGRP_NAME(ib0).eq.'Zmax') exit
      enddo

!$omp parallel do private(ib,in)      
!$acc kernels
      do ib= NODGRP_INDEX(ib0-1)+1, NODGRP_INDEX(ib0)
        in= NODGRP_ITEM(ib)
        IWKX(in,1)= 1
      enddo
!$acc end kernels

!$omp parallel do private(ib,in,iS,iE,k)      
!$acc kernels
      do in= 1, NP
        if (IWKX(in,1).eq.1) then
          B(in)= 0.d0
          D(in)= 1.d0
          
          iS= index(in-1) + 1
          iE= index(in  )
          do k= iS, iE
            AMAT(k)= 0.d0
          enddo
        endif
      enddo
!$acc end kernels

!$omp parallel do private(k,in,iS,iE)  
!$acc kernels    
      do in= 1, NP
        iS= index(in-1) + 1
        iE= index(in  )
        do k= iS, iE
          if (IWKX(item(k),1).eq.1) then
            AMAT(k)= 0.d0
          endif
        enddo
      enddo
!$acc end kernels

!$acc end data

!C==
      return
      end
