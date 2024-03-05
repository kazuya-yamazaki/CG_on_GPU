      subroutine OUTPUT_UCD 
      use pfem_util

      implicit REAL*8 (A-H,O-Z)

      real(kind=kreal), dimension(:), allocatable :: XYZ_G
      real(kind=kreal), dimension(:), allocatable :: VAL_G, WSarray
      integer, dimension(:), allocatable :: NODEflag
      integer, dimension(:), allocatable :: ICELNOD_L, ICELNOD_G
      integer, dimension(:), allocatable :: PEnode, PEelem

      integer, dimension(:), allocatable :: rcountsN, displsN      
      integer, dimension(:), allocatable :: rcountsE, displsE      

      character(len=6) :: ETYPE

!C
!C +--------------------+
!C | INTERNAL ELEMENT's |
!C +--------------------+
!C===
      allocate (NODEflag(NP), ICELNOD_L(8*ICELTOT_INT))

      ICELNOD_L= 0
      NODEflag = 0

      icou= 0
      do icel0= 1, ICELTOT_INT
        icel= intELEM_list(icel0)

        in1= ICELNOD(icel,1)
        in2= ICELNOD(icel,2)
        in3= ICELNOD(icel,3)
        in4= ICELNOD(icel,4)
        in5= ICELNOD(icel,5)
        in6= ICELNOD(icel,6)
        in7= ICELNOD(icel,7)
        in8= ICELNOD(icel,8)

        NODEflag(in1)= 1
        NODEflag(in2)= 1
        NODEflag(in3)= 1
        NODEflag(in4)= 1
        NODEflag(in5)= 1
        NODEflag(in6)= 1
        NODEflag(in7)= 1
        NODEflag(in8)= 1
      enddo

      NElocal= ICELTOT_INT

      icou= 0
      do i= 1, NP
        if (NODEflag(i).ne.0) then
          icou= icou + 1
          NODEflag(i)= icou
        endif
      enddo
!C===

!C
!C +---------------+
!C | GLOBAL ARRAYs |
!C +---------------+
!C===
      allocate (rcountsN(PETOT), displsN(0:PETOT))
      allocate (rcountsE(PETOT), displsE(0:PETOT))

      rcountsN= 0
      rcountsE= 0
      displsN = 0
      displsE = 0

      NNlocal= icou
      call MPI_Allgather                                                & 
     &    (NNlocal, 1, MPI_INTEGER, rcountsN, 1, MPI_INTEGER,           &
     &     MPI_COMM_WORLD, ierr)
      call MPI_Allgather                                                & 
     &    (NElocal, 1, MPI_INTEGER, rcountsE, 1, MPI_INTEGER,           &
     &     MPI_COMM_WORLD, ierr)

      do ip= 1, PETOT
        displsN(ip)= displsN(ip-1) + rcountsN(ip)
        displsE(ip)= displsE(ip-1) + rcountsE(ip)
      enddo

      iS0= displsN(my_rank)
      do icel0= 1, ICELTOT_INT
        icel= intELEM_list(icel0)

        in1= NODEflag(ICELNOD(icel,1))
        in2= NODEflag(ICELNOD(icel,2))
        in3= NODEflag(ICELNOD(icel,3))
        in4= NODEflag(ICELNOD(icel,4))
        in5= NODEflag(ICELNOD(icel,5))
        in6= NODEflag(ICELNOD(icel,6))
        in7= NODEflag(ICELNOD(icel,7))
        in8= NODEflag(ICELNOD(icel,8))

        ICELNOD_L(8*icel0-7)= in1 + iS0
        ICELNOD_L(8*icel0-6)= in2 + iS0
        ICELNOD_L(8*icel0-5)= in3 + iS0
        ICELNOD_L(8*icel0-4)= in4 + iS0
        ICELNOD_L(8*icel0-3)= in5 + iS0
        ICELNOD_L(8*icel0-2)= in6 + iS0
        ICELNOD_L(8*icel0-1)= in7 + iS0
        ICELNOD_L(8*icel0  )= in8 + iS0
      enddo

      NNtotG= displsN(PETOT)
      NEtotG= displsE(PETOT)
!C===

!C
!C +--------------+
!C | ALL_GATHER_V |
!C +--------------+
!C===
      allocate (ICELNOD_G(8*NEtotG))
      allocate (VAL_G(NNtotG), XYZ_G(3*NNtotG), WSarray(3*NNtotG))

!C
!C-- ELEMENT CONNECTIVITY
      do ip= 1, PETOT
         displsE(ip)= 8* displsE(ip)
        rcountsE(ip)= 8*rcountsE(ip)
      enddo

      call MPI_Allgatherv                                               &
     &    (ICELNOD_L, NElocal*8,            MPI_INTEGER,                &
     &     ICELNOD_G, rcountsE, displsE(0), MPI_INTEGER,                &
     &     MPI_COMM_WORLD, ierr)

!C
!C-- NODE COORDINATE/VAL
      do ip= 1, PETOT
         displsN(ip)= 3* displsN(ip)
        rcountsN(ip)= 3*rcountsN(ip)
      enddo

      do i= 1, NP
        if (NODEflag(i).ne.0) then
          in= NODEflag(i)
          WSarray(3*in-2)= XYZ(i,1)
          WSarray(3*in-1)= XYZ(i,2)
          WSarray(3*in  )= XYZ(i,3)
        endif
      enddo

      call MPI_Allgatherv                                               &
     &    (WSarray, 3*NNlocal,            MPI_DOUBLE_PRECISION,         &
     &     XYZ_G  , rcountsN, displsN(0), MPI_DOUBLE_PRECISION,         &
     &     MPI_COMM_WORLD, ierr)

      do ip= 1, PETOT
         displsN(ip)=  displsN(ip) / 3
        rcountsN(ip)= rcountsN(ip) / 3
      enddo

      do i= 1, NP
        if (NODEflag(i).ne.0) then
          in= NODEflag(i)
          WSarray(in)= X(i)
        endif
      enddo

      call MPI_Allgatherv                                               &
     &    (WSarray, NNlocal,              MPI_DOUBLE_PRECISION,         &
     &     VAL_G  , rcountsN, displsN(0), MPI_DOUBLE_PRECISION,         &
     &     MPI_COMM_WORLD, ierr)
!C===
      
!C
!C +----------+
!C | AVS file |
!C +----------+
!C===
      if (my_rank.eq.0) then

        open (21 ,file= 'test.inp', status='unknown')

        N0= 0
        N1= 1
        N3= 3
        N4= 4
        ZERO= 0.d0

        write (21,'(5i8)')  NNtotG, NEtotG, N1, N0, N0
        do i= 1, NNtotG
          XX= XYZ_G(3*i-2)
          YY= XYZ_G(3*i-1)
          ZZ= XYZ_G(3*i  )
          write (21,'(i8,3(1pe16.6))') i, XX, YY, ZZ
        enddo
        do ie= 1, NEtotG
          ETYPE= 'hex   '
          in1= ICELNOD_G(8*ie-7)
          in2= ICELNOD_G(8*ie-6)
          in3= ICELNOD_G(8*ie-5)
          in4= ICELNOD_G(8*ie-4)
          in5= ICELNOD_G(8*ie-3)
          in6= ICELNOD_G(8*ie-2)
          in7= ICELNOD_G(8*ie-1)
          in8= ICELNOD_G(8*ie  )

          write (21,'(i8,i3,1x,a6,1x,8i8)')                             &
     &      ie, N1, ETYPE, in1, in2, in3, in4, in5, in6, in7, in8

        enddo

        write (21,'(10i3)')  N1, N1
        write (21,'(a  )') 'temp, temp'

        igmax    = 0
        VAL_G_max= 0.d0
        do i= 1, NNtotG
          write (21,'(i8, 1pe16.6)')  i, VAL_G(i)
        enddo
        close (21)

      endif
!C===
      end subroutine output_ucd

