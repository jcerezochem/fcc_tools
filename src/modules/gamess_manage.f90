module gamess_manage

    !==============================================================
    ! This code is part of FCC_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get molecular information
    !   from GAMESS (out) files
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

    subroutine read_gamess_natoms(unt,Nat,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Reads coordinates and atom names from GAMESS. The coordinates
        ! are retuned as a 3Nat vector

        !Description
        ! Get geometry and atom names from GAMESS. The number of atoms
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
        ! Need to understand better the GAMESS output
        !==============================================================

        integer,intent(in)           :: unt
        integer,intent(out)          :: Nat
        integer,intent(out),optional :: error_flag

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
            if ( adjustl(line) == "ATOM      ATOMIC                      COORDINATES (BOHR)" ) then
                read(unt,'(A)') line !part of the header
                exit
            endif
        enddo
        
        i=0
        do
            ii = ii + 1
            read(unt,'(A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            !The geometry ends with a blank line
            if (len_trim(line) == 0 ) exit
            i=i+1
        enddo
        Nat=i

        rewind(unt)
        return

    end subroutine read_gamess_natoms

    subroutine read_gamess_geom(unt,Nat,AtName,X,Y,Z,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Reads coordinates and atom names from GAMESS. The coordinates
        ! are retuned as a 3Nat vector

        !Description
        ! Get geometry and atom names from GAMESS. The number of atoms
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
        ! Need to understand better the GAMESS output
        !==============================================================

        integer,intent(in)                          :: unt
        integer,intent(out)                         :: Nat
        character(len=*), dimension(:), intent(out) :: AtName
        real(kind=8), dimension(:), intent(out)     :: X,Y,Z
        integer,intent(out),optional                :: error_flag

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
            if ( adjustl(line) == "ATOM      ATOMIC                      COORDINATES (BOHR)" ) then
                read(unt,'(A)') line !part of the header
                exit
            endif
        enddo
        
        i=0
        do
            ii = ii + 1
            read(unt,'(A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            !The geometry ends with a blank line
            if (len_trim(line) == 0 ) exit
            i=i+1
            read(line,*) AtName(i), cnull, X(i), Y(i), Z(i)
        enddo
        Nat=i
        !Transform to AA
        X(1:Nat) = X(1:Nat)*BOHRtoANGS
        Y(1:Nat) = Y(1:Nat)*BOHRtoANGS
        Z(1:Nat) = Z(1:Nat)*BOHRtoANGS

        rewind(unt)
        return

    end subroutine read_gamess_geom

    subroutine read_gamess_hess(unt,Nat,Hlt,error_flag)


        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Read Hessian from GAMESS output. Returs the triangular part of the
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
        !  Using variable format as suggested in:
        !  http://stackoverflow.com/questions/9881186/fortran-output-format-dependent-on-a-variable
        !  (Answer by eriktous)
        !==============================================================

        integer,intent(in)                      :: unt
        integer,intent(in)                      :: Nat
        real(kind=8), dimension(:), intent(out) :: Hlt
        integer,intent(out),optional            :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        integer :: N
        !I/O
        integer :: IOstatus
        character(len=100) :: fmt
        !Counters
        integer :: i, j, k, ii, imax, imin, jmin, &
                   ib, nblocks, icols
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
                if ( adjustl(line) == "CARTESIAN FORCE CONSTANT MATRIX" ) then
                    read(unt,*) line
                    exit
                endif
        enddo

        !Organized in blocks, displaying 6 columns each (2atoms)
        nblocks = N/6
        if (6*nblocks /= N) nblocks=nblocks+1
        do ib=1,nblocks
            imin = (ib-1)*6 + 1
            imax = ib    *6
            imax=min(imax,N)
            icols = 1 + (imax-imin)
            write(fmt,'(a,i0,a)') '(20X,',icols,'F9.6)'
            !Pass headers
            read(unt,'(A)') line !blank
            read(unt,'(A)') line !index
            read(unt,'(A)') line !atom names
            read(unt,'(A)') line !axis labels
            !Parse hessian elements
            read(unt,fmt) Hpart(imin:imax,imin:N)
        enddo

        !Hpart is incomplete, but the upper tringular part
        ! is all there (as for the read procedure)
        k = 0
        do i=1,N
            do j=1,i
                k = k + 1
                Hlt(k) = Hpart(j,i)
            enddo
        enddo
        if (k /= (N*(N+1)/2)) then
            error_flag = 2
        endif

        rewind(unt)
        return

    end subroutine read_gamess_hess


end module gamess_manage

