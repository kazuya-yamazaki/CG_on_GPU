!C
!C***
!C*** INPUT_GRID
!C***
!C
      subroutine INPUT_GRID
      use pfem_util
      implicit REAL*8 (A-H,O-Z)
      
      call define_file_name (HEADER, fname, my_rank)
      open (11, file= fname, status= 'unknown', form= 'formatted')

!C
!C-- NEIB-PE
      read (11,'(10i10)') kkk
      read (11,'(10i10)') NEIBPETOT
      allocate (NEIBPE(NEIBPETOT))

      read (11,'(10i10)') (NEIBPE(i), i= 1, NEIBPETOT)

      do i= 1, NEIBPETOT
        if (NEIBPE(i).gt.PETOT-1) then
          call ERROR_EXIT (202, my_rank)
        endif
      enddo

!C
!C-- NODE
      read (11,'(10i10)')  NP, N
      allocate (XYZ(NP,3), NODE_ID(NP,2))
      XYZ= 0.d0
        do i= 1, NP
          read (11,*) NODE_ID(i,1), NODE_ID(i,2), (XYZ(i,kk),kk=1,3)
        enddo

!C
!C-- ELEMENT
      read (11,*)  ICELTOT, ICELTOT_INT
      
      allocate (ICELNOD(ICELTOT,8), intELEM_list(ICELTOT))
      allocate (ELEM_ID(ICELTOT,2))
      read (11,'(10i10)') (NTYPE, i= 1, ICELTOT)
      do icel= 1, ICELTOT
        read (11,'(i10,2i5,8i10)') (ELEM_ID(icel,jj),jj=1,2),           &
     &                       IMAT, (ICELNOD(icel,k), k= 1, 8)
      enddo

      read (11,'(10i10)') (intELEM_list(ic0), ic0= 1, ICELTOT_INT)

!C
!C-- COMMUNICATION table
      allocate (IMPORT_INDEX(0:NEIBPETOT))
      allocate (EXPORT_INDEX(0:NEIBPETOT))

      IMPORT_INDEX= 0
      EXPORT_INDEX= 0
      
      if (PETOT.ne.1) then
      read (11,'(10i10)') (IMPORT_INDEX(i), i= 1, NEIBPETOT)
      nn= IMPORT_INDEX(NEIBPETOT)
      allocate (IMPORT_ITEM(nn))
      do i= 1, nn
        read (11,*) IMPORT_ITEM(i)
      enddo

      read (11,'(10i10)') (EXPORT_INDEX(i), i= 1, NEIBPETOT)
      nn= EXPORT_INDEX(NEIBPETOT)
      allocate (EXPORT_ITEM(nn))
      do i= 1, nn
        read (11,*) EXPORT_ITEM(i)
      enddo
      endif

!C
!C-- NODE grp. info.
      read (11,'(10i10)') NODGRPtot
      allocate (NODGRP_INDEX(0:NODGRPtot),NODGRP_NAME(NODGRPtot))
      NODGRP_INDEX= 0

      read (11,'(10i10)') (NODGRP_INDEX(i), i= 1, NODGRPtot)
      nn= NODGRP_INDEX(NODGRPtot)
      allocate (NODGRP_ITEM(nn))
      
      do k= 1, NODGRPtot
        iS= NODGRP_INDEX(k-1) + 1
        iE= NODGRP_INDEX(k  )
        read (11,'(a80)') NODGRP_NAME(k)
        nn= iE - iS + 1
        if (nn.ne.0) then
          read (11,'(10i10)') (NODGRP_ITEM(kk),kk=iS, iE)
        endif
      enddo

      close (11)

      return
      end
