module turbomol_manage

    !==============================================================
    ! This code is part of FCC_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get molecular information
    !   from TURBOMOL (out) files
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

    subroutine read_turbomol_natoms(unt,Nat,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Reads coordinates and atom names from TURBOMOL. The coordinates
        ! are retuned as a 3Nat vector

        !Description
        ! Get geometry and atom names from TURBOMOL. The number of atoms
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
        ! Need to understand better the TURBOMOL output
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
            ! 2) Found what looked for!  2 different places to look for geom    
!                 if ( adjustl(line) == "| Atomic coordinate, charge and isotop information |" ) then
!                     read(unt,'(A)') line ! bar
!                     read(unt,'(A)') line ! blank
!                     read(unt,'(A)') line ! header
                if ( adjustl(line) == "actual cartesian coordinates" ) then
                    read(unt,'(A)') line ! bar
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
            if (len_trim(line) == 0) exit
            i=i+1
        enddo
        Nat=i

        rewind(unt)
        return

    end subroutine read_turbomol_natoms

    subroutine read_turbomol_geom(unt,Nat,AtName,X,Y,Z,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Reads coordinates and atom names from TURBOMOL. The coordinates
        ! are retuned as a 3Nat vector

        !Description
        ! Get geometry and atom names from TURBOMOL. The number of atoms
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
        ! Need to understand better the TURBOMOL output
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
!                 if ( adjustl(line) == "| Atomic coordinate, charge and isotop information |" ) then
!                     read(unt,'(A)') line ! bar
!                     read(unt,'(A)') line ! blank
!                     read(unt,'(A)') line ! header
                if ( adjustl(line) == "actual cartesian coordinates" ) then
                    read(unt,'(A)') line ! bar
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
            if (len_trim(line) == 0) exit
            i=i+1
            read(line,*) cnull, AtName(i), X(i), Y(i), Z(i)
            !Get atom name from atomic number
            call set_word_upper_case(AtName(i))
        enddo
        Nat=i
        !Transform to AA
        X(1:Nat) = X(1:Nat)*BOHRtoANGS
        Y(1:Nat) = Y(1:Nat)*BOHRtoANGS
        Z(1:Nat) = Z(1:Nat)*BOHRtoANGS

        rewind(unt)
        return

    end subroutine read_turbomol_geom

    subroutine read_turbomol_hess(unt,Nat,Hlt,error_flag)


        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Read Hessian from TURBOMOL output. Returs the triangular part of the
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
        integer :: N, NN
        !I/O
        integer :: IOstatus
        character(len=100) :: fmt
        !Counters
        integer :: i, j, k, ii, imax, imin, jmin, &
                   ib, nblocks, icols, istart, isum
        
        
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
                if ( adjustl(line) == "CARTESIAN FORCE CONSTANT MATRIX (hartree/bohr**2)" ) then
                    read(unt,'(A)') line ! subtitle
                    read(unt,'(A)') line ! bar
                    exit
                endif
        enddo

        !Organized in blocks of 6 cols
        icols = 6
        nblocks = N/icols
        if (mod(N,icols) /= 0) nblocks = nblocks+1
        do ib = 1,nblocks
            ! Place the tip on the correct line
            read(unt,'(A)') line ! blank
            read(unt,'(A)') line ! header1
            read(unt,'(A)') line ! header2
            ! Initialize auxiliar counters
            istart = (ib-1)*icols
            isum = 0
            do i=1,N
                ! Accumunalte terms (isum) but cycle if the block is not read
                isum = isum + i
                if (i<=istart) cycle
                ! Determine the parts we are goint to read
                ! The j-th row is from 1:sum(1···j)
                ! But at the current round we get for the j-th col 
                !   istart:min(istart+i,istart+6)
                imin = isum-i + istart+1
                ii = min(i-istart,icols)
                imax = isum-i + istart+ii
                read(unt,'(14X,1000F10.7)') Hlt(imin:imax)
            enddo
        enddo

        rewind(unt)
        return

    end subroutine read_turbomol_hess


!     subroutine read_turbomol_grad(unt,Nat,Grad,error_flag,symm)
! 
! 
!         !==============================================================
!         ! This code is part of FCC_TOOLS
!         !==============================================================
!         !Description
!         ! Read Hessian from TURBOMOL output. Returs the triangular part of the
!         ! Hessian matrix in AU
!         ! 
!         !Arguments
!         ! unt   (inp) scalar   unit for the file
!         ! Nat   (inp) scalar   Number of atoms
!         ! Grad  (out) vector   Gradient (AU)
!         ! error_flag (out) scalar  error_flag :
!         !                                 0 : Success
!         !                                -i : Read error on line i
!         !                                 2 : Wrong number of elements for Grad
!         ! symm  (inp,opt) scalar symmetry flag
!         !
!         !==============================================================
! 
!         integer,intent(in) :: unt
!         integer,intent(in) :: Nat
!         real(kind=8), dimension(:), intent(out) :: Grad
!         integer,intent(out) :: error_flag
!         character(len=*),intent(in),optional :: symm
! 
!         !Local stuff
!         !=============
!         character(len=240) :: line=""
!         character(len=2)   :: irrep, cnull
!         integer :: N, NN
!         !I/O
!         integer :: IOstatus
!         !Counters
!         integer :: i, ii
!         
!         
!         !Use N to store 3*Nat
!         N = 3*Nat
! 
!         ! Search section
!         ii = 0
!         error_flag = 0
!         do
!                 ii = ii + 1
!                 read(unt,'(A)',IOSTAT=IOstatus) line
!                 ! Two possible scenarios while reading:
!                 ! 1) End of file
!                 if ( IOstatus < 0 ) then
!                     error_flag = -ii
!                     rewind(unt)
!                     return
!                 endif
!                 ! 2) Found what looked for!      
!                 if ( adjustl(line) == "*              Molecular gradients               *" ) then
!                     read(unt,'(A)') line ! header box
!                     read(unt,'(A)') line ! header box
!                     read(unt,'(A)') line ! blank line
!                     read(unt,'(A)') line ! irrep info
!                     ! Get irrep
!                     call split_line(line,":",line,irrep)
!                     print*, "irrep: ", irrep
!                     read(unt,'(A)') line ! blank line
!                     exit
!                 endif
!         enddo
! 
!         if (present(symm).and.adjustl(symm)=="CI") then
!             if (mod(Nat,2) /= 0) then
!                 print*, "ERROR: cannot read Ci with odd number of atoms"
!                 stop 
!             endif
!             ii=0
!             do i=1,Nat/2
!                 ii=ii+1
!                 read(unt,*) cnull, cnull, Grad(3*ii-2)
!                 read(unt,*) cnull, cnull, Grad(3*ii-1)
!                 read(unt,*) cnull, cnull, Grad(3*ii  )
!                 ii=ii+1
!                 Grad(3*ii-2) = -Grad(3*ii-5)
!                 Grad(3*ii-1) = -Grad(3*ii-4)
!                 Grad(3*ii  ) = -Grad(3*ii-3)
!             enddo
! 
!         else 
!             print*, "ERROR: molcas gradient reader under developent. Only molcas-ci available"
!             stop
!         endif
! 
!         return
! 
!     end subroutine read_turbomol_grad

    subroutine read_turbomol_nm(unt,Nvib,Nat,Freq,L,error_flag,withTrRot)


        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Read Hessian from TURBOMOL output. Returs the triangular part of the
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

        integer,intent(in)  :: unt
        integer,intent(in)  :: Nat
        integer,intent(out) :: Nvib
        real(kind=8), dimension(:),   intent(out) :: Freq
        real(kind=8), dimension(:,:), intent(out) :: L
        integer,intent(out) :: error_flag
        logical,intent(in),optional :: withTrRot

        !Local stuff
        !=============
        character(len=240) :: line=""
        character(len=20)  :: check
        integer :: N, NN
        !I/O
        integer :: IOstatus
        character(len=100) :: fmt
        !Counters
        integer :: i, j, k, ii, imax, imin, jmin, &
                   ib, nblocks, icols, istart, isum
        !Auxiliar arrays
        integer,dimension(6) :: modenums
        integer              :: nfreq
        
        
        !Use N to store 3*Nat-6
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
                if ( adjustl(line) == "NORMAL MODES and VIBRATIONAL FREQUENCIES (cm**(-1))" ) then
                    do i=1,12
                        read(unt,'(A)') line
                    enddo
                    exit
                endif
        enddo

        !Organized in blocks 
        do
            ! Place the tip on the correct line
            read(unt,'(A)') line ! blank
            read(unt,'(A20,A)') check, line ! mode number
            if (trim(adjustl(check)) /= "mode") exit
            call string2vector_int(line,modenums,nfreq," ")
            imin = modenums(1)
            imax = modenums(nfreq)
            read(unt,'(A)') line ! blank

            ! Frequencies
            read(unt,'(20X,1000(X,F8.2))') Freq(imin:imax)
            ! --

            read(unt,'(A)') line ! blank
            read(unt,'(A)') line ! Symmetry
            read(unt,'(A)') line ! blank
            read(unt,'(A)') line ! IR
            read(unt,'(A)') line ! dDIP/dQ
            read(unt,'(A)') line ! intensity(km/mol)
            read(unt,'(A)') line ! intensity(%)
            read(unt,'(A)') line ! blank
            read(unt,'(A)') line ! RAMAN
            read(unt,'(A)') line ! blank
            do i=1,N
                read(unt,'(20X,1000(X,F8.5))') L(i,imin:imax)
            enddo
            read(unt,'(A)') line ! blank
            read(unt,'(A)') line ! Reduced mass
            read(unt,'(A)') line ! blank
        enddo
        ! Take Nvib from last imax
        Nvib = imax

        ! Turbomol reports 3N modes. Tr+Rot with freq=0.00, so we remove Tr+Rot
        ! unless requested
        if (present(withTrRot) .and. withTrRot) then
            rewind(unt)
            return
        else
            k=0
            NN=Nvib
            do i=1,Nvib
                k=k+1
                if (Freq(k) == 0.d0) then
                    do j=k,NN-1
                        Freq(j)      = Freq(j+1)
                        L(1:3*Nat,j) = L(1:3*Nat,j+1)
                    enddo
                    k  = k-1
                    NN = NN-1
                endif
            enddo
            Nvib = NN
        endif

        rewind(unt)
        return

    end subroutine read_turbomol_nm


end module turbomol_manage

