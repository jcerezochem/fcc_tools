module general_math

    contains


      function fctorial(nn,mm)
!c
!c  reports the elapsed cpu time
!c
      integer,intent(in) :: nn,mm
      real(8) :: fctorial 
      ! Local
      integer :: i

      fctorial=1.d0
      do i=1,mm
      fctorial=fctorial*dfloat(nn+1-i)/dfloat(i)
!c      write(6,*) nn+1-mm,mm,fctorial
      enddo
      return
      end function fctorial

end module general_math