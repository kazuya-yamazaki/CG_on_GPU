!C
!C***
!C*** DEFINE_FILE_NAME
!C***
!C
      subroutine DEFINE_FILE_NAME (HEADERo, filename, my_rank)

      character (len=80) ::  HEADERo, filename
      character (len=80) ::  HEADER
      character (len= 1) ::  SUBindex1
      character (len= 2) ::  SUBindex2
      character (len= 3) ::  SUBindex3
      character (len= 4) ::  SUBindex4
      character (len= 5) ::  SUBindex5
      character (len= 6) ::  SUBindex6
      integer:: LENGTH, ID
   
      HEADER= adjustL (HEADERo)
      LENGTH= len_trim(HEADER)

      if (my_rank.le.9) then
        ID= 1
        write(SUBindex1 ,'(i1.1)') my_rank
       else if (my_rank.le.99) then
        ID= 2
        write(SUBindex2 ,'(i2.2)') my_rank
       else if (my_rank.le.999) then
        ID= 3
        write(SUBindex3 ,'(i3.3)') my_rank
       else if (my_rank.le.9999) then
        ID= 4
        write(SUBindex4 ,'(i4.4)') my_rank
       else if (my_rank.le.99999) then
        ID= 5
        write(SUBindex5 ,'(i5.5)') my_rank
       else if (my_rank.le.999999) then
        ID= 6
        write(SUBindex6 ,'(i6.6)') my_rank
      endif

      if (ID.eq.1) filename= HEADER(1:LENGTH)//'.'//SUBindex1
      if (ID.eq.2) filename= HEADER(1:LENGTH)//'.'//SUBindex2
      if (ID.eq.3) filename= HEADER(1:LENGTH)//'.'//SUBindex3
      if (ID.eq.4) filename= HEADER(1:LENGTH)//'.'//SUBindex4
      if (ID.eq.5) filename= HEADER(1:LENGTH)//'.'//SUBindex5
      if (ID.eq.6) filename= HEADER(1:LENGTH)//'.'//SUBindex6


      end subroutine define_file_name
