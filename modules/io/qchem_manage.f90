module qchem_manage

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

    subroutine read_qchem_natoms(unt,Nat,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Reads coordinates and atom names from Psi4. The coordinates
        ! are retuned as a 3Nat vector

        !Description
        ! Get geometry and atom names from Psi4. The number of atoms
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
        ! Need to understand better the Psi4 output
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
        integer :: i, ii, k
        !Units flag
        character(len=10) :: DIST_UNITS

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
            if ( INDEX(line,"Standard Nuclear Orientation (") /= 0 ) then
                exit
            endif
        enddo

        ! Overpass lines until reaching the target table
        !Read Headers
        read(unt,'(A)',IOSTAT=IOstatus) line
        !Read table lines
        read(unt,'(A)',IOSTAT=IOstatus) line
        !Start reading geometry
        i=0
        read(unt,'(X,A)',IOSTAT=IOstatus) line
        do while (index(line,'--------------------------')==0)
            i=i+1
            read(line,*) k
            if (k/=i) call alert_msg('fatal','Reading coordinates in read_qchem_natoms')
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) then
                print*, "Unexpected end of file while reading Geometry"
                error_flag = 1
                rewind(unt)
                return
            endif
        enddo
        Nat = i

        rewind(unt)
        return

    end subroutine read_qchem_natoms


    subroutine read_qchem_geom(unt,Nat,AtName,X,Y,Z,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Reads coordinates and atom names from QChem. The coordinates
        ! are retuned as a 3Nat vector
        ! Data are taken from FIRST entry for "Standard Nuclear Orientation"
        ! (uncomment code to get LAST entry)
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
        ! Need to understand better the Psi4 output
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
        integer :: i, ii, k
        !Units flag
        character(len=10) :: DIST_UNITS
        logical :: job_ended

        ! Search section
        error_flag = 0
        ii = 0
        ! Loop to get last entry
!         job_ended = .false.
!         do while (.not.job_ended)
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
            if ( INDEX(line,"Standard Nuclear Orientation (") /= 0 ) then
                exit
            endif
!             ! 3) Job ended
!             if ( index(line,'Total job time') /= 0) then
!                 job_ended = .true.
!                 exit
!             endif
        enddo
!         if (job_ended) exit
            
        !Get units
        call split_line(line,'(',cnull,DIST_UNITS)
        call split_line(line,')',DIST_UNITS,cnull)
        call set_word_upper_case(DIST_UNITS)

        ! Overpass lines until reaching the target table
        !Read Headers
        read(unt,'(A)',IOSTAT=IOstatus) line
        !Read table lines
        read(unt,'(A)',IOSTAT=IOstatus) line
        !Start reading geometry
        i=0
        read(unt,'(X,A)',IOSTAT=IOstatus) line
        do while (index(line,'--------------------------')==0)
            i=i+1
            read(line,*) k,         &
                        AtName(i), &
                        X(i),      &
                        Y(i),      &
                        Z(i)
            if (k/=i) call alert_msg('fatal','Reading coordinates in read_qchem_geom')
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) then
                print*, "Unexpected end of file while reading Geometry"
                error_flag = 1
                rewind(unt)
                return
            endif
        enddo
        Nat = i
        
        !Manage lenght units
        if (adjustl(DIST_UNITS) == "BOHR") then
            X(1:Nat) = X(1:Nat) * BOHRtoANGS
            Y(1:Nat) = Y(1:Nat) * BOHRtoANGS
            Z(1:Nat) = Z(1:Nat) * BOHRtoANGS
        endif
        
!         ! Loop to get last entry (if only Freq, the first is also valid, as we do for Grad!)
!         enddo

        rewind(unt)
        return

    end subroutine read_qchem_geom
    

    subroutine read_qchem_grad(unt,Nat,Grad,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Read Gradient from QChem output (first entry) 
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

        integer,intent(in) :: unt
        integer,intent(in) :: Nat
        real(kind=8), dimension(:), intent(out) :: Grad
        integer,intent(out) :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        integer :: N
        !I/O
        integer :: IOstatus
        character(len=100) :: fmt
        !Counters
        integer :: i, j, k, ii, jini, jfin, &
                   iblock, nblocks, icols
        !Auxiliar arrays
        real(kind=8),dimension(3,6) :: Gpart
        
        
        !Use N to store 3*Nat
        N = Nat

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
                if ( INDEX(line,"Gradient of ") /= 0 ) then
                    exit
                endif
        enddo

        !Hessian elements arranged in blocks of 6 columns each
        nblocks = N/6
        if (N /= nblocks*6) nblocks=nblocks+1
        do iblock=1,nblocks
            !First line contains indices (useless)
            read(unt,'(A)') line
            jini=(iblock-1)*6+1
            jfin=min(N,iblock*6)
            do i=1,3
                read(unt,*) k, Gpart(i,1:jfin-jini+1)
            enddo
            do i=jini,jfin 
                Grad(3*(i-1)+1:3*(i-1)+3) = Gpart(1:3,i-jini+1)
            enddo
        enddo

        rewind(unt)
        return

    end subroutine read_qchem_grad


    subroutine read_qchem_hess(unt,Nat,Hlt,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Read Hessian from QChem output. Returs the triangular part of the
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

        integer,intent(in) :: unt
        integer,intent(in) :: Nat
        real(kind=8), dimension(:), intent(out) :: Hlt
        integer,intent(out) :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        integer :: N
        !I/O
        integer :: IOstatus
        character(len=100) :: fmt
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
                if ( INDEX(line,"Final Hessian.") /= 0 .or. &
                     INDEX(line,"Hessian of the SCF Energy") /= 0) then
                    exit
                endif
        enddo

        !Hessian elements arranged in blocks of 6 columns each
        nblocks = N/6
        if (N /= nblocks*6) nblocks=nblocks+1
        do iblock=1,nblocks
            !First line contains indices (useless)
            read(unt,'(A)') line
            jini=(iblock-1)*6+1
            jfin=min(N,iblock*6)
            do i=1,N
                read(unt,*) k, Hpart(i,jini:jfin)
            enddo
        enddo

        !Get Hlt from the whole matrix in Hpart
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

    end subroutine read_qchem_hess


end module qchem_manage

