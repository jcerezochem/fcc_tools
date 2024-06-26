module orca_manage
    use constants
    use line_preprocess
    implicit none

    contains

    subroutine read_orca_natoms(unt,Nat,error_flag)
        !Description
        ! Reads number of atoms from ORCA output
        !Arguments
        ! unt     (inp) int /scalar    unit for the file
        ! Nat     (out) int /scalar    Number of atoms
        ! io_flag  (io ) flag          Error flag:
        !                                   0 : Success
        !                                  -i : Read error on line i
        !                          -1000000-i : Error in header
        !==============================================================
        integer,intent(in)           :: unt
        integer,intent(out)          :: Nat
        integer,intent(out),optional :: error_flag

        character(len=240) :: line=""
        integer :: IOstatus
        integer :: ii

        ii = 0
        do
            ii = ii + 1
            read(unt,'(A)',iostat=IOstatus) line
            if (IOStatus < 0) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            if (adjustl(line) == "$atoms") then
                ii = ii + 1
                read(unt,*,iostat=IOstatus) Nat
                if (IOStatus < 0) then
                    error_flag = -ii
                    rewind(unt)
                    return
                endif
                exit
            endif
        enddo

        error_flag = 0
        rewind(unt)
        return
    end subroutine read_orca_natoms


    subroutine read_orca_geom(unt,Nat,AtName,X,Y,Z,Mass,error_flag)
        !Description
        ! Get geometry and atom names from section 4. The number of atoms
        ! is also taken (since they are needed), using read_gausslog_natoms
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file
        ! Nat     (out) int /scalar    Number of atoms
        ! AtName  (out) char/vector    Atom names
        ! X,Y,Z   (out) real/vectors   Coordinate vectors (ANGSTRONG)
        ! Mass    (out) real/vectors   Coordinate vectors (ANGSTRONG)
        ! io_flag  (io ) flag          Error flag:
        !                                   0 : Success
        !                                   2 : Unkonwn format
        !==============================================================
        integer,intent(in)                        :: unt
        integer,intent(out)                       :: Nat
        real(kind=8),dimension(:),intent(out)     :: X,Y,Z,Mass
        character(len=*),dimension(:),intent(out) :: AtName
        integer,intent(out),optional              :: error_flag

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
            if (adjustl(line) == "$atoms") then
                ii = ii + 1
                read(unt,*,iostat=IOstatus) Nat
                if (IOStatus < 0) then
                    error_flag = -ii
                    rewind(unt)
                    return
                endif
                exit
            endif
        enddo

        if (Nat /= size(AtName)) then
            error_flag = 1
            rewind(unt)
            return
        endif

        do i=1,Nat
            ii = ii + 1
            read(unt,'(A)',iostat=IOstatus) line
            if (IOStatus < 0) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            read(line,*,iostat=IOstatus) AtName(i),Mass(i),X(i),Y(i),Z(i)
            if (IOStatus < 0) then
                error_flag = -ii - 1000000
                rewind(unt)
                return
            endif
        enddo

        X(1:Nat) = X(1:Nat) * BOHRtoANGS
        Y(1:Nat) = Y(1:Nat) * BOHRtoANGS
        Z(1:Nat) = Z(1:Nat) * BOHRtoANGS

        error_flag = 0
        rewind(unt)
        return

    end subroutine read_orca_geom

    subroutine read_orca_hess(unt,Nat,Hlt,error_flag)
        !Description
        ! Read Hessian from ORCA output. Returs the triangular part of the
        ! Hessian matrix in AU
        !
        !Arguments
        ! unt   (inp) scalar   unit for the file
        ! Nat   (inp) scalar   Number of atoms
        ! Hlt   (out) vector   Lower triangular part of Hessian matrix (AU)
        ! error_flag (out) scalar  error_flag :
        !                                 0 : Success
        !                                -i : Read error on line i
        !==============================================================
        integer,intent(in)                      :: unt
        integer,intent(in)                      :: Nat
        real(kind=8), dimension(:), intent(out) :: Hlt
        integer,intent(out),optional            :: error_flag

        character(len=240) :: line=""
        integer :: cols_per_block
        integer,dimension(:),allocatable :: columns
        integer :: IOstatus
        integer :: N
        integer :: i, j, k, ii, imax, imin, &
                   ib, nblocks, icols
        real(kind=8),dimension(3*Nat,3*Nat) :: Hpart
        character :: cnull

        ii = 0
        do
            ii = ii + 1
            read(unt,'(A)',iostat=IOstatus) line
            if (IOStatus < 0) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            if (adjustl(line) == "$hessian") then
                ii = ii + 1
                read(unt,*,iostat=IOstatus) N
                if (IOStatus < 0) then
                    error_flag = -ii
                    rewind(unt)
                    return
                endif
                exit
            endif
        enddo

        if (N /= 3 * Nat) then
            error_flag = 1
            rewind(unt)
            return
        endif

        ! First get cols_per_block from header
        ! as it may vary with the nr of decimals
        read(unt,'(A)',iostat=IOstatus) line ! indices
        allocate(columns(10))
        call string2ivector(line,columns,cols_per_block,' ')
        deallocate(columns)
        allocate(columns(cols_per_block))

        nblocks = N / cols_per_block
        if (cols_per_block * nblocks /= N) nblocks = nblocks + 1
        do ib=1,nblocks
            imin = (ib-1)*cols_per_block + 1
            imax = ib    *cols_per_block
            imax = min(imax,N)
            icols = 1 + (imax-imin)
            !Pass header
            ii = ii + 1
            if (ib > 1) then
                read(unt,'(A)',iostat=IOstatus) line ! indices
                if (IOStatus < 0) then
                    error_flag = -ii
                    rewind(unt)
                    return
                endif
            endif
            read(line,*) columns(1:icols)
            if (columns(1) /= imin - 1) then
                error_flag = -ii - 1000000
                rewind(unt)
                return
            endif
            if (columns(icols) /= imax - 1) then
                error_flag = -ii - 2000000
                rewind(unt)
                return
            endif
            !Parse hessian elements
            do i=1,N
                ii = ii + 1
                read(unt,'(A)',iostat=IOStatus) line
                if (IOStatus < 0) then
                    error_flag = -ii
                    rewind(unt)
                    return
                endif
                read(line,*) j, Hpart(imin:imax,i)
                if (j /= i - 1) then
                    error_flag = -ii - 3000000
                    rewind(unt)
                    return
                endif
            enddo
        enddo

        !===================
        k = 0
        do i=1,N
            do j=1,i
                k = k + 1
                Hlt(k) = Hpart(j,i)
            enddo
        enddo
        if (k /= (N*(N+1)/2)) then
            error_flag = 2
            rewind(unt)
            return
        endif

        error_flag = 0
        rewind(unt)
        return

    end subroutine read_orca_hess
    
    subroutine read_orca4_hess(unt,Nat,Hlt,error_flag)
        !Description
        ! Read Hessian from ORCA4 output. Returs the triangular part of the
        ! Hessian matrix in AU
        !
        !Arguments
        ! unt   (inp) scalar   unit for the file
        ! Nat   (inp) scalar   Number of atoms
        ! Hlt   (out) vector   Lower triangular part of Hessian matrix (AU)
        ! error_flag (out) scalar  error_flag :
        !                                 0 : Success
        !                                -i : Read error on line i
        !==============================================================
        integer,intent(in)                      :: unt
        integer,intent(in)                      :: Nat
        real(kind=8), dimension(:), intent(out) :: Hlt
        integer,intent(out),optional            :: error_flag

        integer,parameter :: cols_per_block = 5

        character(len=240) :: line=""
        integer,dimension(cols_per_block) :: columns
        integer :: IOstatus
        integer :: N
        integer :: i, j, k, ii, imax, imin, &
                   ib, nblocks, icols
        real(kind=8),dimension(3*Nat,3*Nat) :: Hpart
        character :: cnull

        ii = 0
        do
            ii = ii + 1
            read(unt,'(A)',iostat=IOstatus) line
            if (IOStatus < 0) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            if (adjustl(line) == "$hessian") then
                ii = ii + 1
                read(unt,*,iostat=IOstatus) N
                if (IOStatus < 0) then
                    error_flag = -ii
                    rewind(unt)
                    return
                endif
                exit
            endif
        enddo

        if (N /= 3 * Nat) then
            error_flag = 1
            rewind(unt)
            return
        endif

        nblocks = N / cols_per_block
        if (cols_per_block * nblocks /= N) nblocks = nblocks + 1
        do ib=1,nblocks
            imin = (ib-1)*cols_per_block + 1
            imax = ib    *cols_per_block
            imax = min(imax,N)
            icols = 1 + (imax-imin)
            !Pass header
            ii = ii + 1
            read(unt,'(A)',iostat=IOstatus) line ! indices
            if (IOStatus < 0) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            read(line,*) columns(1:icols)
            if (columns(1) /= imin - 1) then
                error_flag = -ii - 1000000
                rewind(unt)
                return
            endif
            if (columns(icols) /= imax - 1) then
                error_flag = -ii - 2000000
                rewind(unt)
                return
            endif
            !Parse hessian elements
            do i=1,N
                ii = ii + 1
                read(unt,'(A)',iostat=IOStatus) line
                if (IOStatus < 0) then
                    error_flag = -ii
                    rewind(unt)
                    return
                endif
                read(line,*) j, Hpart(imin:imax,i)
                if (j /= i - 1) then
                    error_flag = -ii - 3000000
                    rewind(unt)
                    return
                endif
            enddo
        enddo

        !===================
        k = 0
        do i=1,N
            do j=1,i
                k = k + 1
                Hlt(k) = Hpart(j,i)
            enddo
        enddo
        if (k /= (N*(N+1)/2)) then
            error_flag = 2
            rewind(unt)
            return
        endif

        error_flag = 0
        rewind(unt)
        return

    end subroutine read_orca4_hess

    !--------------------------------------
    ! Using .engrad files
    !--------------------------------------

    subroutine read_orcaengrad_natoms(unt,Nat,error_flag)
        !Description
        ! Reads number of atoms from ORCA engrad file
        !Arguments
        ! unt     (inp) int /scalar    unit for the file
        ! Nat     (out) int /scalar    Number of atoms
        ! io_flag  (io ) flag          Error flag:
        !                                   0 : Success
        !                                  -i : Read error on line i
        !                          -1000000-i : Error in header
        !==============================================================
        integer,intent(in)           :: unt
        integer,intent(out)          :: Nat
        integer,intent(out),optional :: error_flag

        character(len=240) :: line=""
        integer :: IOstatus
        integer :: ii

        ii = 0
        do
            ii = ii + 1
            read(unt,'(A)',iostat=IOstatus) line
            if (IOStatus < 0) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            if (adjustl(line) == "# Number of atoms") then
                ii = ii + 1
                read(unt,*,iostat=IOstatus) line ! line with # only
                read(unt,*,iostat=IOstatus) Nat
                if (IOStatus < 0) then
                    error_flag = -ii
                    rewind(unt)
                    return
                endif
                exit
            endif
        enddo

        error_flag = 0
        rewind(unt)
        return
    end subroutine read_orcaengrad_natoms

    subroutine read_orcaengrad_ener(unt,E,error_flag)
        !Description
        ! Reads number of atoms from ORCA engrad file
        !Arguments
        ! unt     (inp) int /scalar    unit for the file
        ! E       (out) real/scalar    Current energy
        ! io_flag  (io ) flag          Error flag:
        !                                   0 : Success
        !                                  -i : Read error on line i
        !                          -1000000-i : Error in header
        !==============================================================
        integer,intent(in)           :: unt
        real(8),intent(out)          :: E
        integer,intent(out),optional :: error_flag

        character(len=240) :: line=""
        integer :: IOstatus
        integer :: ii

        ii = 0
        do
            ii = ii + 1
            read(unt,'(A)',iostat=IOstatus) line
            if (IOStatus < 0) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            if (adjustl(line) == "# The current total energy in Eh") then
                ii = ii + 1
                read(unt,*,iostat=IOstatus) line ! line with # only
                read(unt,*,iostat=IOstatus) E
                if (IOStatus < 0) then
                    error_flag = -ii
                    rewind(unt)
                    return
                endif
                exit
            endif
        enddo

        error_flag = 0
        rewind(unt)
        return
    end subroutine read_orcaengrad_ener

    subroutine read_orcaengrad_grad(unt,Nat,Grad,error_flag)
        !Description
        ! Reads gradient of atoms from ORCA engrad file
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
            if (adjustl(line) == "# The current gradient in Eh/bohr") then
                ii = ii + 1
                read(unt,*,iostat=IOstatus) line ! line with # only
                exit
                if (IOStatus < 0) then
                    error_flag = -ii
                    rewind(unt)
                    return
                endif
                exit
            endif
        enddo

        do i=1,3*Nat
            read(unt,*) Grad(i)
        enddo

        error_flag = 0
        rewind(unt)
        return
    end subroutine read_orcaengrad_grad

    subroutine read_orcaengrad_geom(unt,Nat,AtName,X,Y,Z,Mass,error_flag)
        !Description
        ! Reads geometry from ORCA engrad file
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file
        ! Nat     (out) int /scalar    Number of atoms
        ! AtName  (out) char/vector    Atom names
        ! X,Y,Z   (out) real/vectors   Coordinate vectors (ANGSTRONG)
        ! Mass    (out) real/vectors   Coordinate vectors (ANGSTRONG)
        ! io_flag  (io ) flag          Error flag:
        !                                   0 : Success
        !                                   2 : Unkonwn format
        !==============================================================
        integer,intent(in)                        :: unt
        integer,intent(inout)                     :: Nat
        real(kind=8),dimension(:),intent(out)     :: X,Y,Z,Mass
        character(len=*),dimension(:),intent(out) :: AtName
        integer,intent(out),optional              :: error_flag

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
            if (adjustl(line) == "# The atomic numbers and current coordinates in Bohr") then
                ii = ii + 1
                read(unt,*,iostat=IOstatus) line ! line with # only
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
            read(unt,*) ii, X(i), Y(i), Z(i)
            X(i) = X(i) * BOHRtoANGS
            Y(i) = Y(i) * BOHRtoANGS
            Z(i) = Z(i) * BOHRtoANGS
            AtName(i) = atname_from_atnum(ii)
        enddo

        error_flag = 0
        rewind(unt)
        return
    end subroutine read_orcaengrad_geom

    !--------------------------------------
    ! Using .hess files
    !--------------------------------------
    ! NOTE: mosts parsers are already implemented as read_orca_XX

    subroutine read_orcahess_irdip(unt,DipD,error_flag)
        !Description
        ! Get geometry and atom names from section 4. The number of atoms
        ! is also taken (since they are needed), using read_gausslog_natoms
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file
        ! DipD    (out) real/vector    Dipole derivatives (3x3Nat)
        ! io_flag  (io ) flag          Error flag:
        !                                   0 : Success
        !                                   2 : Unkonwn format
        !==============================================================
        integer,intent(in)                      :: unt
        real(kind=8),dimension(:),intent(out)   :: DipD
        integer,intent(out),optional            :: error_flag

        character(len=240) :: line=""
        integer :: IOstatus
        integer :: i, ii, j, N

        ii = 0
        do
            ii = ii + 1
            read(unt,'(A)',iostat=IOstatus) line
            if (IOStatus < 0) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            if (adjustl(line) == "$dipole_derivatives") then
                ii = ii + 1
                read(unt,*,iostat=IOstatus) N
                if (IOStatus < 0) then
                    error_flag = -ii
                    rewind(unt)
                    return
                endif
                exit
            endif
        enddo

        do i=1,N
            ii = ii + 1
            read(unt,'(A)',iostat=IOstatus) line
            if (IOStatus < 0) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            j = 3*(i-1)
            read(line,*,iostat=IOstatus) DipD(j+1), DipD(j+2), DipD(j+3)
            if (IOStatus < 0) then
                error_flag = -ii - 1000000
                rewind(unt)
                return
            endif
        enddo

        error_flag = 0
        rewind(unt)
        return

    end subroutine read_orcahess_irdip

end module orca_manage

