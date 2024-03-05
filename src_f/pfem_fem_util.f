!C
!C***
!C*** pfem_fem_util
!C***
!C
      module pfem_fem_util
      use pfem_util

      real(kind=kreal), dimension(2,2,8) :: PNQ, PNE, PNT
      real(kind=kreal), dimension(2)     :: WEI, POS
      integer(kind=kint), dimension(100) :: NCOL1, NCOL2

      real(kind=kreal), dimension(2,2,2,8) :: SHAPE
      real(kind=kreal), dimension(2,2,2,8) :: PNX, PNY, PNZ
      real(kind=kreal), dimension(2,2,2  ) :: DETJ

      end module pfem_fem_util
