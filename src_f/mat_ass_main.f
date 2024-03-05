!C
!C***
!C*** MAT_ASS_MAIN
!C***
!C
      subroutine MAT_ASS_MAIN
      use pfem_util
      use util
      use nvtxwrap
      implicit REAL*8 (A-H,O-Z)
      integer(kind=kint), dimension(  8) :: nodLOCAL

!C
!C +------------------+
!C | ELEMENT Coloring |
!C +------------------+
!C===
      allocate (ELMCOLORindex(0:NP))
      allocate (ELMCOLORitem (ICELTOT))
      if (allocated (IWKX)) deallocate (IWKX)
      allocate (IWKX(0:NP,3))
      
      call mynvtxStartRange("Coloring")

      IWKX= 0
      icou= 0
      do icol= 1, NP
        do i= 1, NP
          IWKX(i,1)= 0
        enddo
        do icel= 1, ICELTOT
          if (IWKX(icel,2).eq.0) then
            in1= ICELNOD(icel,1)
            in2= ICELNOD(icel,2)
            in3= ICELNOD(icel,3)
            in4= ICELNOD(icel,4)
            in5= ICELNOD(icel,5)
            in6= ICELNOD(icel,6)
            in7= ICELNOD(icel,7)
            in8= ICELNOD(icel,8)

            ip1= IWKX(in1,1)
            ip2= IWKX(in2,1)
            ip3= IWKX(in3,1)
            ip4= IWKX(in4,1)
            ip5= IWKX(in5,1)
            ip6= IWKX(in6,1)
            ip7= IWKX(in7,1)
            ip8= IWKX(in8,1)

            isum= ip1 + ip2 + ip3 + ip4 + ip5 + ip6 + ip7 + ip8
            if (isum.eq.0) then 
              icou= icou + 1
              IWKX(icol,3)= icou
              IWKX(icel,2)= icol
              ELMCOLORitem(icou)= icel

              IWKX(in1,1)= 1
              IWKX(in2,1)= 1
              IWKX(in3,1)= 1
              IWKX(in4,1)= 1
              IWKX(in5,1)= 1
              IWKX(in6,1)= 1
              IWKX(in7,1)= 1
              IWKX(in8,1)= 1
              if (icou.eq.ICELTOT) goto 100            
            endif
          endif
        enddo
      enddo

 100  continue
      ELMCOLORtot= icol
      IWKX(0          ,3)= 0
      IWKX(ELMCOLORtot,3)= ICELTOT

      do icol= 0, ELMCOLORtot
        ELMCOLORindex(icol)= IWKX(icol,3)
      enddo

!      write (*,'(a,2i8)') '### Number of Element Colors', 
!     &                     my_rank, ELMCOLORtot
      deallocate (IWKX)
!C===
   
      allocate (AMAT(NPLU))
      allocate (B(NP), D(NP), X(NP))
!$acc enter data create(B, D, X, AMAT)   

!$acc kernels present(B, D, X, AMAT) 
      AMAT(:) = 0.d0
         B(:) = 0.d0
         X(:) = 0.d0
         D(:) = 0.d0
!$acc end kernels

      WEI(1)= +1.0000000000D+00
      WEI(2)= +1.0000000000D+00

      POS(1)= -0.5773502692D+00
      POS(2)= +0.5773502692D+00
      
      call mynvtxEndRange()

!C
!C-- INIT.
!C     PNQ   - 1st-order derivative of shape function by QSI
!C     PNE   - 1st-order derivative of shape function by ETA
!C     PNT   - 1st-order derivative of shape function by ZET
!C
      call mynvtxStartRange("Init coefs")
      do kp= 1, 2
      do jp= 1, 2
      do ip= 1, 2
        QP1= 1.d0 + POS(ip)
        QM1= 1.d0 - POS(ip)
        EP1= 1.d0 + POS(jp)
        EM1= 1.d0 - POS(jp)
        TP1= 1.d0 + POS(kp)
        TM1= 1.d0 - POS(kp)
        SHAPE(ip,jp,kp,1)= O8th * QM1 * EM1 * TM1
        SHAPE(ip,jp,kp,2)= O8th * QP1 * EM1 * TM1
        SHAPE(ip,jp,kp,3)= O8th * QP1 * EP1 * TM1
        SHAPE(ip,jp,kp,4)= O8th * QM1 * EP1 * TM1
        SHAPE(ip,jp,kp,5)= O8th * QM1 * EM1 * TP1
        SHAPE(ip,jp,kp,6)= O8th * QP1 * EM1 * TP1
        SHAPE(ip,jp,kp,7)= O8th * QP1 * EP1 * TP1
        SHAPE(ip,jp,kp,8)= O8th * QM1 * EP1 * TP1
        PNQ(jp,kp,1)= - O8th * EM1 * TM1
        PNQ(jp,kp,2)= + O8th * EM1 * TM1
        PNQ(jp,kp,3)= + O8th * EP1 * TM1
        PNQ(jp,kp,4)= - O8th * EP1 * TM1
        PNQ(jp,kp,5)= - O8th * EM1 * TP1
        PNQ(jp,kp,6)= + O8th * EM1 * TP1
        PNQ(jp,kp,7)= + O8th * EP1 * TP1
        PNQ(jp,kp,8)= - O8th * EP1 * TP1
        PNE(ip,kp,1)= - O8th * QM1 * TM1
        PNE(ip,kp,2)= - O8th * QP1 * TM1
        PNE(ip,kp,3)= + O8th * QP1 * TM1
        PNE(ip,kp,4)= + O8th * QM1 * TM1
        PNE(ip,kp,5)= - O8th * QM1 * TP1
        PNE(ip,kp,6)= - O8th * QP1 * TP1
        PNE(ip,kp,7)= + O8th * QP1 * TP1
        PNE(ip,kp,8)= + O8th * QM1 * TP1
        PNT(ip,jp,1)= - O8th * QM1 * EM1
        PNT(ip,jp,2)= - O8th * QP1 * EM1
        PNT(ip,jp,3)= - O8th * QP1 * EP1
        PNT(ip,jp,4)= - O8th * QM1 * EP1
        PNT(ip,jp,5)= + O8th * QM1 * EM1
        PNT(ip,jp,6)= + O8th * QP1 * EM1
        PNT(ip,jp,7)= + O8th * QP1 * EP1
        PNT(ip,jp,8)= + O8th * QM1 * EP1
      enddo
      enddo
      enddo
      call mynvtxEndRange()
      
!$acc  data copyin(SHAPE, PNQ, PNE, PNT, ELMCOLORindex, ELMCOLORitem)
!$acc+ copyin(ICELNOD, XYZ, WEI) present(B,D,AMAT)
!$acc+ create(nodLOCAL, PNX, PNY, PNZ, DETJ)

      call mynvtxStartRange("Main loop")
      do icol= 1, ELMCOLORtot
!$acc parallel loop private(PNX,PNY,PNZ,DETJ,nodLOCAL)
!$omp parallel do private (icel0,icel,in1,in2,in3,in4,in5,in6,in7,in8)  &
!$omp&            private (nodLOCAL,ie,je,ip,jp,kp,kk,iiS,iiE,k)        &
!$omp&            private (DETJ,PNX,PNY,PNZ)                            &
!$omp&            private (PNXi,PNYi,PNZi,PNXj,PNYj,PNZj,COEFij,SHi)    &
!$omp&            private (X1,X2,X3,X4,X5,X6,X7,X8)                     &
!$omp&            private (Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)                     &
!$omp&            private (Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8)                     &
!$omp&            private (QVC,QV0,ipn,jpn,kpn,coef)                    &
      do icel0= ELMCOLORindex(icol-1)+1, ELMCOLORindex(icol)
        icel= ELMCOLORitem(icel0)

        in1= ICELNOD(icel,1)
        in2= ICELNOD(icel,2)
        in3= ICELNOD(icel,3)
        in4= ICELNOD(icel,4)
        in5= ICELNOD(icel,5)
        in6= ICELNOD(icel,6)
        in7= ICELNOD(icel,7)
        in8= ICELNOD(icel,8)
!C
!C== JACOBIAN & INVERSE JACOBIAN
        nodLOCAL(1)= in1
        nodLOCAL(2)= in2
        nodLOCAL(3)= in3
        nodLOCAL(4)= in4
        nodLOCAL(5)= in5
        nodLOCAL(6)= in6
        nodLOCAL(7)= in7
        nodLOCAL(8)= in8

        X1= XYZ(in1,1)
        X2= XYZ(in2,1)
        X3= XYZ(in3,1)
        X4= XYZ(in4,1)
        X5= XYZ(in5,1)
        X6= XYZ(in6,1)
        X7= XYZ(in7,1)
        X8= XYZ(in8,1)

        Y1= XYZ(in1,2)
        Y2= XYZ(in2,2)
        Y3= XYZ(in3,2)
        Y4= XYZ(in4,2)
        Y5= XYZ(in5,2)
        Y6= XYZ(in6,2)
        Y7= XYZ(in7,2)
        Y8= XYZ(in8,2)

        QVC= O8th * (X1+X2+X3+X4+X5+X6+X7+X8+
     &               Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8)
        Z1= XYZ(in1,3)
        Z2= XYZ(in2,3)
        Z3= XYZ(in3,3)
        Z4= XYZ(in4,3)
        Z5= XYZ(in5,3)
        Z6= XYZ(in6,3)
        Z7= XYZ(in7,3)
        Z8= XYZ(in8,3)

        call JACOBI (DETJ, PNQ, PNE, PNT, PNX, PNY, PNZ,                &
     &               X1, X2, X3, X4, X5, X6, X7, X8,                    &
     &               Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8,                    &
     &               Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8 )
!C
!C== CONSTRUCT the GLOBAL MATRIX
        do ie= 1, 8
          ip = nodLOCAL(ie)
        do je= 1, 8
          jp = nodLOCAL(je)

          kk= 0
          if (jp.ne.ip) then
            iiS= index(ip-1) + 1
            iiE= index(ip  )
            do k= iiS, iiE
              if ( item(k).eq.jp ) then
                kk= k
                exit
              endif
            enddo
          endif

          QV0   = 0.d0
          COEFij= 0.d0
          do kpn= 1, 2
          do jpn= 1, 2
          do ipn= 1, 2
            coef= WEI(ipn)*WEI(jpn)*WEI(kpn)

            PNXi= PNX(ipn,jpn,kpn,ie)
            PNYi= PNY(ipn,jpn,kpn,ie)
            PNZi= PNZ(ipn,jpn,kpn,ie)

            PNXj= PNX(ipn,jpn,kpn,je)
            PNYj= PNY(ipn,jpn,kpn,je)
            PNZj= PNZ(ipn,jpn,kpn,je)

            COEFij= COEFij + coef * COND * 
     &                      (PNXi*PNXj+PNYi*PNYj+PNZi*PNZj) *
     &                       dabs(DETJ(ipn,jpn,kpn))
            SHi= SHAPE(ipn,jpn,kpn,ie)
            QV0= QV0 + SHi * QVOL * coef * dabs(DETJ(ipn,jpn,kpn))
          enddo
          enddo
          enddo

          if (jp.eq.ip) then
            D(ip)= D(ip) + COEFij
            B(ip)= B(ip) + QV0*QVC
           else
            AMAT(kk)= AMAT(kk) + COEFij
          endif
        enddo
        enddo
      enddo
!$acc end parallel
      enddo   
!$acc end data
      call mynvtxEndRange()   

      return
      end
