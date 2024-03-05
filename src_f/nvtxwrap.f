      module nvtxwrap
      contains
      
      subroutine mynvtxStartRange(rangename)
#ifdef OPENACC
        use nvtx
#endif
        implicit none
        character(*) :: rangename
        
#ifdef OPENACC
        call nvtxStartRange(rangename)
#endif
        return
      end subroutine mynvtxStartRange
      
      subroutine mynvtxEndRange()
#ifdef OPENACC
        use nvtx
#endif
        implicit none
        
#ifdef OPENACC
        call nvtxEndRange()
#endif
        return
      end subroutine mynvtxEndRange
      
      end module nvtxwrap
