module molden_manage

    !==============================================================
    ! This code is part of FCC_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get molecular information
    !   from MOLCAS (out) files
    !    
    ! Notes  
    !  All subroutines rewind the file after using it
    !==============================================================

    !Common declarations:
    !===================
    use constants
    use line_preprocess
    implicit none

    contains

    subroutine read_molden_natoms(unt,Nat,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Reads coordinates and atom names from MOLCAS. The coordinates
        ! are retuned as a 3Nat vector

        !Description
        ! Get geometry and atom names from MOLCAS. The number of atoms
        ! is also taken
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file 
        ! Nat     (out) int /scalar    Number of atoms
        ! io_flag  (io ) flag          Error flag:
        !                                   0 : Success
        !                                  -i : Read error on line i
        ! 
        !Note
        ! Need to understand better the MOLCAS output
        !==============================================================

        integer,intent(in)  :: unt
        integer,intent(out) :: Nat
        integer,intent(out) :: error_flag

        !Local variables
        !=============
        character(len=240) :: line=""
        character :: cnull
        !I/O
        integer :: IOstatus
        !Counters
        integer :: i, ii
        
        ! Search section
        ii = 0
        do
            ii = ii + 1
            read(unt,'(A)',IOSTAT=IOstatus) line
            ! Two possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            ! 2) Found what looked for!      
            if ( adjustl(line) == "[N_ATOMS]" ) then
                read(unt,*) Nat 
                exit
            endif

        enddo

        rewind(unt)
        return

    end subroutine read_molden_natoms

    subroutine read_molden_geom(unt,Nat,AtName,X,Y,Z,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Reads coordinates and atom names from MOLCAS. The coordinates
        ! are retuned as a 3Nat vector

        !Description
        ! Get geometry and atom names from MOLCAS. The number of atoms
        ! is also taken
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file 
        ! Nat     (out) int /scalar    Number of atoms
        ! AtName  (out) char/vertor    Atom names
        ! X,Y,Z   (out) real/vectors   Coordinate vectors (ANGSTRONG)
        ! io_flag  (io ) flag          Error flag:
        !                                   0 : Success
        !                                  -i : Read error on line i
        ! 
        !Note
        ! Need to understand better the MOLCAS output
        !==============================================================

        integer,intent(in)  :: unt
        integer,intent(out) :: Nat
        character(len=*), dimension(:), intent(out) :: AtName
        real(kind=8), dimension(:), intent(out) :: X,Y,Z
        integer,intent(out) :: error_flag

        !Local variables
        !=============
        character(len=10)  :: units
        character(len=240) :: line=""
        character :: cnull
        !I/O
        integer :: IOstatus
        !Counters
        integer :: i, ii, AtNum
        
        
        ! Search section
        ii = 0
        do
            ii = ii + 1
            read(unt,'(A)',IOSTAT=IOstatus) line
            ! Two possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            ! 2) Found what looked for!      
            if ( adjustl(line) == "[N_ATOMS]" ) then
                read(unt,*) Nat 
                exit
            endif

        enddo
        rewind(unt)

        ! Search section
        ii = 0
        do
            ii = ii + 1
            read(unt,'(A)',IOSTAT=IOstatus) line
            ! Two possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            ! 2) Found what looked for!      
            if ( index(line,"[ATOMS]") /= 0 ) then
                read(line,*) cnull, units
                exit
            endif

        enddo
        do i=1,Nat
            read(unt,*) AtName(i), cnull, cnull, X(i), Y(i), Z(i)
        enddo
        if (trim(adjustl(units))=='(AU)') then
            !Transform to AA
            X(1:Nat) = X(1:Nat)*BOHRtoANGS
            Y(1:Nat) = Y(1:Nat)*BOHRtoANGS
            Z(1:Nat) = Z(1:Nat)*BOHRtoANGS
        endif

        rewind(unt)
        return

    end subroutine read_molden_geom


end module molden_manage

