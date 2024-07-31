module cp2k_manage

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

    subroutine read_cp2k_natoms(unt,Nat,error_flag)

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
            if ( index(line,"                             - Atoms:") /= 0 ) then
                exit
            endif
        enddo
        
        read(line,*) cnull, cnull, Nat

        rewind(unt)
        return

    end subroutine read_cp2k_natoms

    subroutine read_cp2k_geom(unt,Nat,AtName,X,Y,Z,error_flag)

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
        ! Nat     (inp) int /scalar    Number of atoms
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
        integer,intent(in)                          :: Nat
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


        print*, "No method to get geometry from CP2K output. You may use MOLDEN file instead."
        error_flag = -99

        rewind(unt)
        return

    end subroutine read_cp2k_geom

    subroutine read_cp2k_grad(unt,Nat,Grad,error_flag)
        !Description
        ! Reads gradient of atoms from CP2K engrad file
        !Arguments
        ! unt     (inp) int /scalar    unit for the file
        ! Nat     (inp) int/scalar     Number of atoms
        ! Grad    (out) real/vector    Gradient
        ! io_flag  (io ) flag          Error flag:
        !                                   0 : Success
        !                                  -i : Read error on line i
        !                          -1000000-i : Error in header
        !==============================================================
        integer,intent(in)               :: unt
        integer,intent(in)               :: Nat
        real(8),dimension(:),intent(out) :: Grad
        integer,intent(out),optional     :: error_flag

        character(len=240) :: line=""
        integer :: IOstatus
        integer :: i, ii

        ii = 0
        do
            ii = ii + 1
            read(unt,'(A)',iostat=IOstatus) line
            if (IOStatus < 0) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            if (adjustl(line) == "Minimum Structure - Energy and Forces:") then
                ii = ii + 1
                read(unt,*,iostat=IOstatus) line ! energy
                read(unt,*,iostat=IOstatus) line ! labels
                exit
                if (IOStatus < 0) then
                    error_flag = -ii
                    rewind(unt)
                    return
                endif
                exit
            endif
        enddo

        read(unt,'(29X,3F17.9)') Grad(1:3*Nat)

        error_flag = 0
        rewind(unt)
        return
    end subroutine read_cp2k_grad

    subroutine read_cp2k_hess(unt,Nat,Hlt,error_flag)


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
        character(len=2) :: AtName
        character :: cnull
        integer :: N
        !I/O
        integer :: IOstatus
        character(len=100) :: fmt
        !Counters
        integer :: i, j, k, ii,jj, imax, imin, jmin, &
                   ib, nblocks, nitems, icols
        !Auxiliar arrays
        real(kind=8),dimension(3*Nat,3*Nat) :: H
        real(kind=8),dimension(Nat) :: Mass
        
        
        !Use N to store 3*Nat
        N = 3*Nat

        ii = 0
        do
            ii = ii + 1
            read(unt,'(A)',iostat=IOstatus) line
            if (IOStatus < 0) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            if (adjustl(line) == "Minimum Structure - Energy and Forces:") then
                ii = ii + 1
                read(unt,*,iostat=IOstatus) line ! energy
                read(unt,*,iostat=IOstatus) line ! labels
                exit
                if (IOStatus < 0) then
                    error_flag = -ii
                    rewind(unt)
                    return
                endif
                exit
            endif
        enddo

        do i=1,Nat
            read(unt,*) cnull, AtName, cnull, cnull, cnull
            ii = atnum_from_atname(AtName)
            Mass(i) = atmass_from_atnum(ii) * AMUtoAU
        enddo
        rewind(unt)

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
                if ( adjustl(line) == "VIB| Hessian in cartesian coordinates" ) then
                    exit
                endif
        enddo

        !Organized in blocks, displaying 5 columns each
        nitems = 5
        nblocks = N/nitems
        if (nitems*nblocks /= N) nblocks=nblocks+1
        do ib=1,nblocks
            imin = (ib-1)*nitems + 1
            imax = ib    *nitems
            imax=min(imax,N)
            icols = 1 + (imax-imin)
            write(fmt,'(a,i0,a)') '(12X,',icols,'F13.6)'
            !Pass headers
            read(unt,'(A)') line !blank
            read(unt,'(A)') line !index
            !Parse hessian elements
            read(unt,fmt) H(imin:imax,1:N)
        enddo

        !Unmassweight the Hessian
        ! By trial and error, we realized that the data in the file
        ! correspond to the Hessian in MWC (not simply Cartesian)
        k=0
        do i=1,3*Nat
        do j=1,i
            k=k+1
            ii = (i-1)/3+1
            jj = (j-1)/3+1
            ! Note we need to divide each element by 1e6 (just by trial and error)
            H(i,j) = H(i,j)*sqrt(Mass(ii)*Mass(jj)) / 1.d6
            H(j,i) = H(i,j)
        enddo
        enddo

        ! is all there (as for the read procedure)
        print*, "Symmetrizing C2PK Hessian"
        k = 0
        do i=1,N
            do j=1,i
                k = k + 1
                Hlt(k) = 0.5d0 * (H(j,i) + H(i,j))
            enddo
        enddo
        if (k /= (N*(N+1)/2)) then
            error_flag = 2
        endif

        rewind(unt)
        return

    end subroutine read_cp2k_hess


end module cp2k_manage

