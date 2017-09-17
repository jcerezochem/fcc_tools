module cfour_manage

    !==============================================================
    ! This code is part of FCC_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get molecular information
    !   from CFOUR (out) files
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

    subroutine read_cfour_natoms(unt,Nat,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Reads coordinates and atom names from CFOUR. The coordinates
        ! are retuned as a 3Nat vector

        !Description
        ! Get geometry and atom names from CFOUR. The number of atoms
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
        ! Need to understand better the CFOUR output
        !==============================================================

        integer,intent(in)           :: unt
        integer,intent(out)          :: Nat
        integer,intent(out),optional :: error_flag

        !Local variables
        !=============
        character(len=240) :: line="", hline
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
            if ( adjustl(line) == "Z-matrix   Atomic            Coordinates (in bohr)" ) then
                read(unt,'(A)') line !part of the header
                ! The table stats with a horizontal line
                read(unt,'(A)') hline !part of the header
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
            !The geometry ends with an horizontal line
            if (adjustl(line) == adjustl(hline)) exit
            i=i+1
        enddo
        Nat=i

        rewind(unt)
        return

    end subroutine read_cfour_natoms

    subroutine read_cfour_geom(unt,Nat,AtName,X,Y,Z,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Reads coordinates and atom names from CFOUR. The coordinates
        ! are retuned as a 3Nat vector

        !Description
        ! Get geometry and atom names from CFOUR. The number of atoms
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
        ! Need to understand better the CFOUR output
        !==============================================================

        integer,intent(in)                          :: unt
        integer,intent(out)                         :: Nat
        character(len=*), dimension(:), intent(out) :: AtName
        real(kind=8), dimension(:), intent(out)     :: X,Y,Z
        integer,intent(out),optional                :: error_flag

        !Local variables
        !=============
        character(len=240) :: line="", hline
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
            if ( adjustl(line) == "Z-matrix   Atomic            Coordinates (in bohr)" ) then
                read(unt,'(A)') line !part of the header
                ! The table stats with a horizontal line
                read(unt,'(A)') hline !part of the header
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
            !The geometry ends with an horizontal line
            if (adjustl(line) == adjustl(hline)) exit
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

    end subroutine read_cfour_geom

    subroutine read_cfour_hess(unt,Nat,Hlt,error_flag)


        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Read Hessian from CFOUR output. Returs the triangular part of the
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

        integer,intent(in)                      :: unt
        integer,intent(in)                      :: Nat
        real(kind=8), dimension(:), intent(out) :: Hlt
        integer,intent(out),optional            :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        character :: cnull
        integer :: N
        !I/O
        integer :: IOstatus
        !Counters
        integer :: i, j, k, ii, imax, imin, jmin, &
                   nblocks, icols, nblank, jini, jfin, &
                   iblock
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
                if ( index(line,"               Molecular hessian") /=0 ) then
                    read(unt,*) line  ! ----------- (hline)
                    exit
                endif
        enddo

        !First, symmetry blocks are shown (always, or only if symm is used?)
        ! [but the values are not the same as in the "real" Hessian...]
        !The the whole hessian is shown after all symmetry section.
        ! It is separated from the last symmetry section by 4 blank lines
        ! while the symmetry sections are separated by 2 blank lines
        nblank=0
        do while (nblank<4)
            read(unt,'(A)') line
            if (len_trim(line)==0) then
                nblank=nblank+1
            else
                ! Restart counter
                nblank=0
            endif
        enddo
        
        !Hessian elements arranged in blocks of 6 columns each
        !Only triangular part is shown
        nblocks = N/6
        if (N /= nblocks*6) nblocks=nblocks+1
        do iblock=1,nblocks
            !Firrst: line is the header (and some blank lines)
            read(unt,'(A)') line
            do while (len_trim(line) == 0)
                read(unt,'(A)') line
            enddo
            jini=(iblock-1)*6+1
            do i=jini,N
                jfin=min(i,iblock*6)
                ! There is a blank line to separate atoms
                read(unt,'(A)') line
                if ( len_trim(line) == 0 ) read(unt,'(A)') line
                read(line,*) cnull,cnull,cnull, Hpart(i,jini:jfin)
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

    end subroutine read_cfour_hess

    subroutine read_FCM_hess(unt,Nat,Hlt,error_flag)


        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Read Hessian from CFOUR output. Returs the triangular part of the
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

        integer,intent(in)                      :: unt
        integer,intent(in)                      :: Nat
        real(kind=8), dimension(:), intent(out) :: Hlt
        integer,intent(out),optional            :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        character :: cnull
        integer :: N
        !I/O
        integer :: IOstatus
        !Counters
        integer :: i, j, k
        !Auxiliar arrays
        real(kind=8),dimension(3*Nat,3*Nat) :: Hpart
        
        
        !Use N to store 3*Nat
        read(unt,*) N
        if (N /= Nat) then
            print*, "Error in read_FCM_hess: Inconsistent nr of atom"
            if (present(error_flag)) then
                error_flag=-1
                return
            else
                stop
            endif
        endif
        N = 3*Nat

        read(unt,'(3(F20.10))') Hpart(1:N,1:N)

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

    end subroutine read_FCM_hess
    
    subroutine read_cfour_grad(unt,Nat,Grad,error_flag)


        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Read Hessian from CFOUR output. Returs gradient as 3N vector
        ! 
        !Arguments
        ! unt   (inp) scalar   unit for the file
        ! Nat   (inp) scalar   Number of atoms
        ! Grad  (out) vector   Gradient vector (AU)
        ! error_flag (out) scalar  error_flag :
        !                                 0 : Success
        !                                -i : Read error on line i
        !                                 2 : Wrong number of elements for Hlt
        ! Notes
        !==============================================================

        integer,intent(in)                      :: unt
        integer,intent(in)                      :: Nat
        real(kind=8), dimension(:), intent(out) :: Grad
        integer,intent(out),optional            :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        integer :: N
        !I/O
        integer :: IOstatus
        character :: cnull
        !Counters
        integer :: i, j, k, ii, imax, imin, jmin, &
                   ib, nblank
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
                if ( index(line,"               Molecular gradient") /=0  ) then
                    read(unt,*) line ! part of the header
                    read(unt,*) line ! part of the header
                    exit
                endif
        enddo

        !First, ther is a column where only non-zero terms seem to be shown
        ! [but the values are not the same as in the "real" gradient...]
        !The the whole gradient is shown after that section.
        ! It is separated from the former section by 2 blank lines
        ! while the separation of the header with former section is 1 blank lines
        nblank=0
        do while (nblank<2)
            read(unt,'(A)') line
            if (len_trim(line)==0) then
                nblank=nblank+1
            else
                ! Restart counter
                nblank=0
            endif
        enddo
        
        ! Read gradient
        read(unt,'(3(9X,F15.10))') Grad(1:N)

        rewind(unt)
        return

    end subroutine read_cfour_grad


end module cfour_manage

