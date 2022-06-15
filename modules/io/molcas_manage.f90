module molcas_manage

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

    subroutine read_molcasUnSym_natoms(unt,Nat,error_flag)

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
                if ( adjustl(line) == "*BEGIN COORDINATES" ) then
                    read(unt,'(A16,I20)') line !Header (does it depend on the JobType?)
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
            if (adjustl(line) == "*END COORDINATES") exit
            i=i+1
        enddo
        Nat=i

        rewind(unt)
        return

    end subroutine read_molcasUnSym_natoms

    subroutine read_molcasUnSym_geom(unt,Nat,AtName,X,Y,Z,error_flag)

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
                if ( adjustl(line) == "*BEGIN COORDINATES" ) then
                    read(unt,'(A16,I20)') line !Header (does it depend on the JobType?)
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
            if (adjustl(line) == "*END COORDINATES") exit
            i=i+1
            read(line,*) cnull, X(i), Y(i), Z(i), AtNum
            !Get atom name from atomic number
            AtName(i) = atname_from_atnum(AtNum)
        enddo
        Nat=i
        !Transform to AA
        X(1:Nat) = X(1:Nat)*BOHRtoANGS
        Y(1:Nat) = Y(1:Nat)*BOHRtoANGS
        Z(1:Nat) = Z(1:Nat)*BOHRtoANGS

        rewind(unt)
        return

    end subroutine read_molcasUnSym_geom

    subroutine read_molcasUnSym_hess(unt,Nat,Hlt,error_flag)


        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Read Hessian from MOLCAS output. Returs the triangular part of the
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
                if ( adjustl(line) == "*BEGIN HESSIAN" ) then
                    read(unt,'(A16,I20)') line, NN
                    if ( NN /= N) then
                        print*, "Error reading Hessian (MOLCAS)"
                        error_flag = 1
                        rewind(unt)
                        return
                    endif
                    exit
                endif
        enddo

        !Organized in blocks by columns (Header + Colum elements)
        do i=1,N
            read(unt,*) line
            read(unt,*) Hpart(1:N,i)
        enddo

        !Get triangular part (we have it whole)
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

    end subroutine read_molcasUnSym_hess


    subroutine read_molcas_grad(unt,Nat,Grad,error_flag,symm)


        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Read Hessian from MOLCAS output. Returs the triangular part of the
        ! Hessian matrix in AU
        ! 
        !Arguments
        ! unt   (inp) scalar   unit for the file
        ! Nat   (inp) scalar   Number of atoms
        ! Grad  (out) vector   Gradient (AU)
        ! error_flag (out) scalar  error_flag :
        !                                 0 : Success
        !                                -i : Read error on line i
        !                                 2 : Wrong number of elements for Grad
        ! symm  (inp,opt) scalar symmetry flag
        !
        !==============================================================

        integer,intent(in) :: unt
        integer,intent(in) :: Nat
        real(kind=8), dimension(:), intent(out) :: Grad
        integer,intent(out) :: error_flag
        character(len=*),intent(in),optional :: symm

        !Local stuff
        !=============
        character(len=240) :: line=""
        character(len=2)   :: irrep, cnull
        integer :: N, NN
        !I/O
        integer :: IOstatus
        !Counters
        integer :: i, ii
        
        
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
                if ( adjustl(line) == "*              Molecular gradients               *" ) then
                    read(unt,'(A)') line ! header box
                    read(unt,'(A)') line ! header box
                    read(unt,'(A)') line ! blank line
                    read(unt,'(A)') line ! irrep info
                    ! Get irrep
                    call split_line(line,":",line,irrep)
                    print*, "irrep: ", irrep
                    read(unt,'(A)') line ! blank line
                    exit
                endif
        enddo

        if (present(symm).and.adjustl(symm)=="CI") then
            if (mod(Nat,2) /= 0) then
                print*, "ERROR: cannot read Ci with odd number of atoms"
                stop 
            endif
            ii=0
            do i=1,Nat/2
                ii=ii+1
                read(unt,*) cnull, cnull, Grad(3*ii-2)
                read(unt,*) cnull, cnull, Grad(3*ii-1)
                read(unt,*) cnull, cnull, Grad(3*ii  )
                ii=ii+1
                Grad(3*ii-2) = -Grad(3*ii-5)
                Grad(3*ii-1) = -Grad(3*ii-4)
                Grad(3*ii  ) = -Grad(3*ii-3)
            enddo

!         else 
!             print*, "ERROR: molcas gradient reader under developent. Only molcas-ci available"
!             stop
        endif

        return

    end subroutine read_molcas_grad

end module molcas_manage

