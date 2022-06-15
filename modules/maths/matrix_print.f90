module matrix_print

    !==============================================================
    ! This code is part of FCclasses2
    !==============================================================
    !
    ! Despription
    ! -----------
    ! Subroutines to perform matrix manipulations:
    ! *Matrix printing
    !
    !==============================================================

    CONTAINS


   !=========================================
   ! SUBROUTINES FOR PRINTING:
   ! Source: old FCclasses code. Uses 
   !         labels for format and GOTOs 
   !         and does not follow indentation
   !=========================================

      subroutine print_vector(unt,A,N,name)

        !==================================
        ! Description
        ! -----------
        ! Source: adapted from original FCclasses
        !         code for matrix printing (MAT1)
        ! Prints a vector A
        !==================================

        integer,intent(in) :: unt
        real(8),dimension(:),intent(in) :: A
        integer,intent(in) :: N
        character(len=*),intent(in),optional :: name
        !Local
        integer :: i
        integer :: len_line
        character(len=100) :: line 

      !Set the table line accoring to its size
      len_line = 18
      len_line = max(len_line,len_trim(name)+6)
      do i=1,len_line
          line(i:i) = "-"
      enddo

      if (present(name)) then
      write(unt,'(/,X,A)') line(1:len_line)
      write(unt,'(3X,A)') name
      endif 
      write(unt,'(X,A)') line(1:len_line)

      do i=1,N
          write(unt,'(X,I5,F12.5)') i, A(i)
      enddo
      write(unt,'(X,A,/)') line(1:len_line)
! C
      return

      end subroutine print_vector


      subroutine MAT0(unt,AA,NR,NC,name)

        !==================================
        ! Description
        ! -----------
        ! Prints a matrix AA (without 
        ! any associated vector)
        !==================================

        integer,intent(in) :: unt
        real(8),dimension(:,:),intent(in) :: AA
        integer,intent(in) :: NC, NR
        character(len=*),intent(in),optional :: name
        !Local
        integer :: I, J, N, LA, LB, LC, KA, KB, KC
        integer :: len_line
        character(len=100) :: line 
        
      !Set the table line accoring to its size
      len_line = 9+(min(6,NC)*12)
      len_line = max(len_line,len_trim(name)+6)
      do i=1,len_line
          line(i:i) = "-"
      enddo

      if (present(name)) then
      write(unt,'(/,X,A)') line(1:len_line)
      write(unt,'(3X,A)') name
      endif 
      write(unt,'(X,A)') line(1:len_line)

      KA=1
      KC=6
   10 KB=MIN0(KC,NC)
      WRITE (unt,50) (I,I=KA,KB)
      WRITE (unt,70)
      LA=1
      LC=40
   20 LB=MIN0(LC,NR)
      N=0
      DO 30 I=LA,LB
         WRITE (unt,80) I,(AA(I,J),J=KA,KB)
         N=N+1
         IF (N.LT.10) GO TO 30
         WRITE (unt,70)
         N=0
   30 CONTINUE
      IF (LB.EQ.NR) GO TO 40
      LA=LC+1
      LC=LC+40
      WRITE (unt,90)
      GO TO 20
   40 IF (KB.EQ.NC) then
          write(unt,'(X,A,/)') line(1:len_line)
          RETURN
      endif
      KA=KC+1
      KC=KC+6
      IF (NR.GT.25) WRITE (unt,90)
      GO TO 10
! C
   50 FORMAT (/9H         ,I5,9I12)
   70 FORMAT (2H  )
   80 FORMAT (I5,10F12.5)
   90 FORMAT (1H1)
! C
      return

      end subroutine MAT0


      subroutine MAT1(unt,AA,B,NR,NC,name)

        !==================================
        ! Description
        ! -----------
        ! Prints eigenvectors in matrix AA
        ! along with eigenvalues in B
        !==================================

        integer,intent(in) :: unt
        real(8),dimension(:,:),intent(in) :: AA
        real(8),dimension(:),intent(in) :: B
        integer,intent(in) :: NC, NR
        character(len=*),intent(in),optional :: name
        !Local
        integer :: I, J, N, LA, LB, LC, KA, KB, KC
        integer :: len_line
        character(len=81) :: line 

      !Set the table line accoring to its size
      len_line = 9+(min(6,NC)*12)
      len_line = max(len_line,len_trim(name)+6)
      do i=1,len_line
          line(i:i) = "-"
      enddo

      if (present(name)) then
      write(unt,'(/,X,A)') line(1:len_line)
      write(unt,'(3X,A)') name
      endif 
      write(unt,'(X,A)') line(1:len_line)

      KA=1
      KC=6
   10 KB=MIN0(KC,NC)
      WRITE (unt,50) (I,I=KA,KB)
      WRITE (unt,60) (B(I),I=KA,KB)
      WRITE (unt,70)
      LA=1
      LC=40
   20 LB=MIN0(LC,NR)
      N=0
      DO 30 I=LA,LB
         WRITE (unt,80) I,(AA(I,J),J=KA,KB)
         N=N+1
         IF (N.LT.10) GO TO 30
         WRITE (unt,70)
         N=0
   30 CONTINUE
      IF (LB.EQ.NR) GO TO 40
      LA=LC+1
      LC=LC+40
      WRITE (unt,90)
      GO TO 20
   40 IF (KB.EQ.NC) then
          write(unt,'(X,A,/)') line(1:len_line)
          RETURN
      endif
      KA=KC+1
      KC=KC+6
      IF (NR.GT.25) WRITE (unt,90)
      GO TO 10
! C
   50 FORMAT (/9H ROOT NO.,I5,9I12)
   60 FORMAT (/4X,10F12.6)
   70 FORMAT (2H  )
   80 FORMAT (I5,10F12.5)
   90 FORMAT (1H1)
! C
      return

      end subroutine MAT1


end module matrix_print
