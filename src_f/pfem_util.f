!C
!C***
!C*** pfem_util
!C***

      module pfem_util
      include 'mpif.h'
      include 'precision.inc'
!C
!C +--------------+
!C | MPI settings |
!C +--------------+
!C===
      integer :: PEsmpTOT
      integer :: PETOT, my_rank, errno

      character(len= 4 ) ::  penum, penum_left
      character(len= 80) ::  fname, HEADER
!C===

!C
!C +-------------------------+
!C | DISTRIBUTED MESH FILE's |
!C +-------------------------+
!C===
!C
!C-- COMM. TABLE
      integer          :: NEIBPETOT
      integer, pointer :: IMPORT_INDEX(:), IMPORT_ITEM(:)
      integer, pointer :: EXPORT_INDEX(:), EXPORT_ITEM(:)
      integer, pointer :: NEIBPE(:)

!C
!C-- CONNECTIVITIES & BOUNDARY nodes
      integer:: ICELTOT, ICELTOT_INT, NODGRPtot, ELMGRPtot, SUFGRPtot

      real(kind=kreal)  , dimension(:,:),allocatable :: XYZ
      integer(kind=kint), dimension(:,:),allocatable :: iCELNOD
      integer(kind=kint), dimension(:,:),allocatable :: NODE_ID
      integer(kind=kint), dimension(:,:),allocatable :: ELEM_ID
      integer(kind=kint), dimension(:),  allocatable :: intELEM_list

      integer(kind=kint), dimension(:),  allocatable ::                 &
     &  NODGRP_INDEX, NODGRP_ITEM, ELMGRP_INDEX, ELMGRP_ITEM,           &
     &  SUFGRP_INDEX, SUFGRP_ITEM

      character(len=80), dimension(:),  allocatable ::                  &
     &  NODGRP_NAME, ELMGRP_NAME, SUFGRP_NAME
!C===

!C
!C +-----------------+
!C | MATRIX & SOLVER |
!C +-----------------+
!C===

!C
!C-- MATRIX SCALARs
      integer(kind=kint) :: N, NP, N2, NL, NU, NPLU, ELMCOLORtot

!C
!C-- MATRIX arrays
      real(kind=kreal), dimension(:), allocatable ::  D, B, X
      real(kind=kreal), dimension(:), allocatable ::  AMAT

      integer(kind=kint), dimension(:), allocatable :: index, item

      integer(kind=kint), dimension(:),  allocatable :: INLU
      integer(kind=kint), dimension(:,:),allocatable :: IALU

      integer(kind=kint), dimension(:,:),allocatable :: IWKX
      integer(kind=kint), dimension(:),allocatable   :: ELMCOLORindex
      integer(kind=kint), dimension(:),allocatable   :: ELMCOLORitem
!C
!C-- PARAMETER's for LINEAR SOLVER
      integer(kind=kint) :: ITER, ITERactual
      real   (kind=kreal) :: RESID, SIGMA_DIAG, SIGMA
!C===

!C
!C +-------------+
!C | PARAMETER's |
!C +-------------+
!C===

!C
!C-- GENERAL PARAMETER's
      integer(kind=kint ), dimension(100) :: pfemIarray
      real   (kind=kreal), dimension(100) :: pfemRarray
      integer(kind=kint )                 :: NFLAG

      real   (kind=kreal), parameter :: O8th= 0.125d0
!C
!C-- PARAMETER's for FEM
      real(kind=kreal), dimension(2,2,8) :: PNQ, PNE, PNT
      real(kind=kreal), dimension(2)     :: WEI, POS
      integer(kind=kint), dimension(100) :: NCOL1, NCOL2

      real(kind=kreal), dimension(2,2,2,8) :: SHAPE
      real(kind=kreal), dimension(2,2,2,8) :: PNX, PNY, PNZ
      real(kind=kreal), dimension(2,2,2  ) :: DETJ
!C
!C-- PROBLEM PARAMETER's
      real(kind=kreal):: COND, QVOL, QVC
!C===
      end module pfem_util
