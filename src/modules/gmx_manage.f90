module gmx_manage

    !==============================================================
    ! This code is part of FCC_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get molecular information
    !   from Psi4 out files
    !    
    ! Notes  
    !  All subroutines rewind the file after using it
    !==============================================================

    !Common declarations:
    !===================
    use constants
    use line_preprocess
    use alerts
    implicit none

    contains

    subroutine read_g96_natoms(unt,Nat)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Get geometry and atom names from xyz. The number of atoms
        ! is also taken
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file 
        ! Nat     (out) int /scalar    Number of atoms
        !
        !==============================================================

        integer,intent(in)  :: unt
        integer,intent(out) :: Nat

        !local
        integer::i, natoms, ii, ios
        character(len=260) :: line
        !The use of get_structure=.false. fails probably due to
        ! memory issues (change without being directly accesed)
        logical :: not_get_structure=.true.

        !The file is organized in sections. Each begining with a given name
        do

            read(unt,'(A)',iostat=ios) line
            if (ios /= 0) exit

            if (adjustl(line) == "POSITION" ) then
                not_get_structure=.false.
                i=0
                do 
                     read(unt,'(A)') line
                     if (adjustl(line) == "END" ) exit
                     i = i + 1
                enddo
                Nat = i
                return
            elseif (adjustl(line) == "POSITIONRED" ) then
            !this section only has info about coordinates (no atom info!)
            write(0,*) "NOTE: No Atom Names in g96. Masses will not be assigned."
                not_get_structure=.false.
                i=0
                do 
                     read(unt,'(A)') line
                     if (adjustl(line) == "END" ) exit
                enddo
                Nat = i
                return
            else
                cycle
            endif

        enddo

        if (not_get_structure) then
            write(0,*) "ERROR: No structure read in g96 file"
            stop
        endif

        return

     end subroutine read_g96_natoms



    subroutine read_g96_geom(unt,Nat,AtName,X,Y,Z)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Get geometry and atom names from xyz. The number of atoms
        ! is also taken
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file 
        ! Nat     (out) int /scalar    Number of atoms
        ! AtName  (out) char/vertor    Atom names
        ! X,Y,Z   (out) real/vectors   Coordinate vectors (ANGSTRONG)
        ! 
        !==============================================================

        integer,intent(in)  :: unt
        integer,intent(out) :: Nat
        character(len=*), dimension(:), intent(out) :: AtName
        real(kind=8), dimension(:), intent(out) :: X,Y,Z
        !local
        integer::i, natoms, ii, ios
        character(len=260) :: line
        character :: dummy_char
        !The use of get_structure=.false. fails probably due to
        ! memory issues (change without being directly accesed)
        logical :: not_get_structure=.true.

        !The file is organized in sections. Each begining with a given name
        do

            read(unt,'(A)',iostat=ios) line
            if (ios /= 0) exit

            if (adjustl(line) == "POSITION" ) then
                not_get_structure=.false.
                i=0
                do 
                     read(unt,'(A)') line
                     if (adjustl(line) == "END" ) exit
                     i=i+1
                     read(line,*) ii,                    & !resseq, &
                                  dummy_char,            & !resname,&
                                  AtName(i),             &
                                  ii,                    & !serial
                                  X(i),                  &
                                  Y(i),                  &
                                  Z(i)
                     ! nm to \AA
                     X(i) = X(i)*10.d0
                     Y(i) = Y(i)*10.d0
                     Z(i) = Z(i)*10.d0
                enddo
                Nat = i
            elseif (adjustl(line) == "POSITIONRED" ) then
            !this section only has info about coordinates (no atom info!)
            write(0,*) "NOTE: No Atom Names in g96. Masses will not be assigned."
                not_get_structure=.false.
                i=0
                do 
                     read(unt,'(A)') line
                     if (adjustl(line) == "END" ) exit
                     i=i+1
                     read(line,*) X(i),                  &
                                  Y(i),                  &
                                  Z(i)
                     ! nm to \AA
                     X(i) = X(i)*10.d0
                     Y(i) = Y(i)*10.d0
                     Z(i) = Z(i)*10.d0
                enddo
                Nat = i
            else
                cycle
            endif

        enddo

        if (not_get_structure) then
            write(0,*) "ERROR: No structure read in g96 file"
            stop
        endif

        return

    end subroutine read_g96_geom


    subroutine read_gmx_hess(unt,Nat,Hlt,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Read Hessian from ascii file generated from an .mtx file (with gmxdump)
        ! 
        ! 
        !Arguments
        ! unt   (inp) scalar   unit for the file
        ! Nat   (inp) scalar   Number of atoms
        ! Hlt   (out) vector   Lower triangular part of Hessian matrix (AU)
        ! error_flag (out) scalar  error_flag :
        !                                 0 : Success
        !                                 1: not square matrix
        !  
        !==============================================================

        integer,intent(in) :: unt
        integer,intent(in) :: Nat
        real(8), dimension(:), intent(out) :: Hlt
        integer,intent(out) :: error_flag

        !Local stuff
        !=============
        character :: cnull
        !Counter
        integer :: N,M
        integer :: i, j
         
        read(unt,'(A)') cnull
        read(unt,*) N, M
        if (N /= M) then
            write(0,*) "Hessian matrix is not square (gmx)"
            error_flag = 1
            stop
        else if (N /= 3*Nat) then
            write(0,*) "Incorrect Hessian dimension (gmx)"
            error_flag = 2
            stop
        endif

        !Read in triangular form
        j=1
        do i=1,3*Nat
            read(unt,*) Hlt(j:j+i-1)
            j=j+i
        enddo

        ! UNIT CONVERSION                                                              ! GROMACS       --> Atomic Units  
        Hlt(1:3*Nat*(3*Nat+1)/2)=Hlt(1:3*Nat*(3*Nat+1)/2)/CALtoJ/HtoKCALM*BOHRtoNM**2  ! KJ/mol * nm-2 --> Hartree * bohr-2

        return

    end subroutine read_gmx_hess


end module gmx_manage
