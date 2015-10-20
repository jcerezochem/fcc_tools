module molpro_manage

    !==============================================================
    ! This code is part of FCC_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get molecular information
    !   from molpro out files
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

    subroutine read_molpro_natoms(unt,Nat,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Reads coordinates and atom names from molpro. The coordinates
        ! are retuned as a 3Nat vector

        !Description
        ! Get geometry and atom names from molpro. The number of atoms
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
        ! Need to understand better the molpro output
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
        error_flag = 0
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
            if ( adjustl(line) == "Current geometry (xyz format, in Angstrom)" ) then
                !One empty line
                read(unt,'(A)',IOSTAT=IOstatus) line
                exit
            endif
        enddo

        ! Overpass lines until reaching the target table
        !Read Table lines
        read(unt,*,IOSTAT=IOstatus) Nat

        rewind(unt)
        return

    end subroutine read_molpro_natoms


    subroutine read_molpro_geom(unt,Nat,AtName,X,Y,Z,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Reads coordinates and atom names from molpro. The coordinates
        ! are retuned as a 3Nat vector

        !Description
        ! Get geometry and atom names from molpro. The number of atoms
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
        ! Need to understand better the molpro output
        !==============================================================

        integer,intent(in)  :: unt
        integer,intent(out) :: Nat
        character(len=*), dimension(:), intent(out) :: AtName
        real(kind=8), dimension(:), intent(out) :: X,Y,Z
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
        error_flag = 0
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
            if ( adjustl(line) == "Current geometry (xyz format, in Angstrom)" ) then
                !One empty line
                read(unt,'(A)',IOSTAT=IOstatus) line
                exit
            endif
        enddo

        ! Overpass lines until reaching the target table
        !Read Table lines
        read(unt,*,IOSTAT=IOstatus) Nat
        !Title
        read(unt,'(A)',IOSTAT=IOstatus) line
        !Start reading geometry
        do i=1,Nat
            read(unt,*)  AtName(i), &
                         X(i),      &
                         Y(i),      &
                         Z(i)
            if ( IOstatus < 0 ) then
                print*, "Unexpected end of file while reading Geometry"
                error_flag = 1
                rewind(unt)
                return
            endif
        enddo

        rewind(unt)
        return

    end subroutine read_molpro_geom


    subroutine read_molpro_hess(unt,Nat,Hlt,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Read Hessian from molpro output. Returs the triangular part of the
        ! Hessian matrix in AU
        ! 
        !Arguments
        ! unt   (inp) scalar   unit for the file
        ! Nat   (inp) scalar   Number of atoms
        ! Hlt   (out) vector   Lower triangular part of Hessian matrix (AU)
        ! error_flag (out) scalar  error_flag :
        !                                 0 : Success
        !                                -i : Read error on line i
        !                                 2 : Wrong number of elements for Hlt
        ! Notes
        !==============================================================

        integer,intent(in) :: unt
        integer,intent(in) :: Nat
        real(kind=8), dimension(:), intent(out) :: Hlt
        integer,intent(out) :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        character(len=1)   :: cnull
        integer :: N
        !I/O
        integer :: IOstatus
        !Counters
        integer :: i, j, k, ii, jini, jfin, &
                   iblock, nblocks, icols
        !Auxiliar arrays
        real(kind=8),dimension(3*Nat,3*Nat) :: Hpart
        
        
        !Use N to store 3*Nat
        N = 3*Nat

        ! Search section
        ii = 0
        error_flag = 0
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
                if ( INDEX(line,"Force Constants (Second Derivatives of the Energy) in [a.u.]") /= 0 ) then
                    exit
                endif
        enddo

        !Hessian elements arranged in blocks of 5 columns each
        !Only triangular part is shown
        nblocks = N/5
        if (N /= nblocks*5) nblocks=nblocks+1
        do iblock=1,nblocks
            !Rirst line is header
            read(unt,'(A)') line
            jini=(iblock-1)*5+1
            do i=jini,N
                jfin=min(i,iblock*5)
                read(unt,*) cnull, Hpart(i,jini:jfin)
            enddo
        enddo

        !Get Hlt from the half matrix in Hpart
        k = 0
        do i=1,N
            do j=1,i
                k = k + 1
                Hlt(k) = Hpart(i,j)
            enddo
        enddo
        if (k /= (N*(N+1)/2)) then
            error_flag = 2
        endif

        rewind(unt)
        return

    end subroutine read_molpro_hess



end module molpro_manage

