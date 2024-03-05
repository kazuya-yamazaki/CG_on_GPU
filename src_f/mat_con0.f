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

      INLU= 0
      IALU= 0
      
      call mynvtxStartRange("NODES")
      do icel= 1, ICELTOT
        in1= ICELNOD(icel,1)
        in2= ICELNOD(icel,2)
        in3= ICELNOD(icel,3)
        in4= ICELNOD(icel,4)
        in5= ICELNOD(icel,5)
        in6= ICELNOD(icel,6)
        in7= ICELNOD(icel,7)
        in8= ICELNOD(icel,8)

        call FIND_TS_NODE (in1,in2)
        call FIND_TS_NODE (in1,in3)
        call FIND_TS_NODE (in1,in4)
        call FIND_TS_NODE (in1,in5)
        call FIND_TS_NODE (in1,in6)
        call FIND_TS_NODE (in1,in7)
        call FIND_TS_NODE (in1,in8)

        call FIND_TS_NODE (in2,in1)
        call FIND_TS_NODE (in2,in3)
        call FIND_TS_NODE (in2,in4)
        call FIND_TS_NODE (in2,in5)
        call FIND_TS_NODE (in2,in6)
        call FIND_TS_NODE (in2,in7)
        call FIND_TS_NODE (in2,in8)

        call FIND_TS_NODE (in3,in1)
        call FIND_TS_NODE (in3,in2)
        call FIND_TS_NODE (in3,in4)
        call FIND_TS_NODE (in3,in5)
        call FIND_TS_NODE (in3,in6)
        call FIND_TS_NODE (in3,in7)
        call FIND_TS_NODE (in3,in8)

        call FIND_TS_NODE (in4,in1)
        call FIND_TS_NODE (in4,in2)
        call FIND_TS_NODE (in4,in3)
        call FIND_TS_NODE (in4,in5)
        call FIND_TS_NODE (in4,in6)
        call FIND_TS_NODE (in4,in7)
        call FIND_TS_NODE (in4,in8)

        call FIND_TS_NODE (in5,in1)
        call FIND_TS_NODE (in5,in2)
        call FIND_TS_NODE (in5,in3)
        call FIND_TS_NODE (in5,in4)
        call FIND_TS_NODE (in5,in6)
        call FIND_TS_NODE (in5,in7)
        call FIND_TS_NODE (in5,in8)

        call FIND_TS_NODE (in6,in1)
        call FIND_TS_NODE (in6,in2)
        call FIND_TS_NODE (in6,in3)
        call FIND_TS_NODE (in6,in4)
        call FIND_TS_NODE (in6,in5)
        call FIND_TS_NODE (in6,in7)
        call FIND_TS_NODE (in6,in8)

        call FIND_TS_NODE (in7,in1)
        call FIND_TS_NODE (in7,in2)
        call FIND_TS_NODE (in7,in3)
        call FIND_TS_NODE (in7,in4)
        call FIND_TS_NODE (in7,in5)
        call FIND_TS_NODE (in7,in6)
        call FIND_TS_NODE (in7,in8)

        call FIND_TS_NODE (in8,in1)
        call FIND_TS_NODE (in8,in2)
        call FIND_TS_NODE (in8,in3)
        call FIND_TS_NODE (in8,in4)
        call FIND_TS_NODE (in8,in5)
        call FIND_TS_NODE (in8,in6)
        call FIND_TS_NODE (in8,in7)
      enddo
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
        subroutine FIND_TS_NODE (ip1,ip2)

          do kk= 1, INLU(ip1)
            if (ip2.eq.IALU(ip1,kk)) return
          enddo
          icou= INLU(ip1) + 1
          IALU(ip1,icou)= ip2
          INLU(ip1     )= icou
          return

        end subroutine FIND_TS_NODE
      end subroutine MAT_CON0
